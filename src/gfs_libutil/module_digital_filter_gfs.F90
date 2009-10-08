      module module_digital_filter_gfs
!
! a generic digital filter for any model under ESMF 
!
! March 2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
!
      use esmf_mod

      implicit none

! ---------
! dynamics
! ---------
      real, allocatable, save :: dyn_array_save(:,:)
      real             , save :: totalsum
      character(20), allocatable, save :: dyn_name(:)
      integer,       allocatable, save :: dyn_dim(:)
      integer,                    save :: dyn_items, kstep, nstep
      character(20)	state_name
! ---------
! physics
! ---------
      type(esmf_state) , save :: phy_state_save
      character(20), allocatable, save :: phy_name(:)
      integer		phy_items


      contains

! ---------------------------------------------------------------
! subroutine for dynamics
! ---------------------------------------------------------------
      subroutine digital_filter_dyn_init_gfs(dyn_state,ndfistep)
!
      implicit none
      type(esmf_state), intent(in) ::	dyn_state
      integer,          intent(in) ::   ndfistep
!
      TYPE(ESMF_VM)                :: vm
      type(esmf_array)             :: tmp_array
      integer                      :: tmp_rank, dim, dim_max
      integer                      :: m,n,rc
      integer                      :: me, nodes

      INTEGER, DIMENSION(:, :), POINTER :: tmp_counts
!
      CALL ESMF_VMGetCurrent(vm, rc = rc)
      CALL ESMF_VMGet(vm, localpet = me, petcount = nodes, rc = rc)
      ALLOCATE(tmp_counts(10, nodes))

      nstep = ndfistep 
      kstep = - nstep -1
      call esmf_stateget(state     = dyn_state 				&
                        ,name      = state_name 			&
                        ,itemcount = dyn_items 				&
                        ,rc=rc)
      print *,' in digital_filter_dyn_init itemcount = ',dyn_items
      allocate(dyn_name(dyn_items))
      allocate(dyn_dim (dyn_items))

      call esmf_stateget(state=dyn_state 				&
                        ,itemnamelist = dyn_name 			&
                        ,rc=rc)
      print *,' in digital_filter_dyn_init itemname = ',dyn_name
      dim_max=0
      do n=1,dyn_items
        CALL ESMF_StateGet(dyn_state, dyn_name(n), tmp_array, rc = rc)

        call esmf_arrayget(tmp_array	   				&
                          ,rank              = tmp_rank                 &
                          ,indexCountPDimPDe = tmp_counts               &
                          ,rc                = rc)
        dim=tmp_counts(1, me + 1)
        if( dim.gt.1 ) then
          do m=2,tmp_rank 
            dim=dim*tmp_counts(m, me + 1)
          enddo
          dyn_dim(n)=dim
        else
          dyn_dim(n)=1
        endif
        dim_max=max(dyn_dim(n),dim_max)
      enddo

      print *,' in digital_filter_dyn_init max dimension = ',dyn_dim

      allocate(dyn_array_save(dim_max,dyn_items))

      totalsum=0.0
      dyn_array_save=0.0
      
     DEALLOCATE(tmp_counts)

      end subroutine digital_filter_dyn_init_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_sum_gfs(dyn_state)
!
      implicit none
      type(esmf_state), intent(in)  :: dyn_state
!
      TYPE(ESMF_Array)		    :: tmp_array
      real(ESMF_KIND_R8), dimension(:,:), pointer :: tmp_ptr
      real                          :: sx, wx, digfil
      integer                       :: n, i, rc

        kstep = kstep + 1
        sx     = acos(-1.)*kstep/nstep
        wx     = acos(-1.)*kstep/(nstep+1)
        if( kstep.ne.0 ) then
            digfil = sin(wx)/wx*sin(sx)/sx
        else
            digfil=1
        endif 

!        print *,' in digital_filter_sum digfil = ',digfil

        totalsum = totalsum + digfil
        
        do n=1,dyn_items
          if( dyn_dim(n).gt.1 ) then
          CALL ESMF_StateGet(dyn_state, dyn_name(n), tmp_array, rc = rc)

          nullify(tmp_ptr)
          CALL ESMF_ArrayGet(tmp_array, 0, tmp_ptr, rc = rc)

          do i=1,dyn_dim(n)
          dyn_array_save(i,n)=dyn_array_save(i,n)+digfil*tmp_ptr(i,1)
          enddo
          endif
        enddo

      end subroutine digital_filter_dyn_sum_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_average_gfs(dyn_state)
!
      implicit none
      type(esmf_state), intent(inout) :: dyn_state
!
      TYPE(ESMF_Array)                :: tmp_array
      TYPE(ESMF_DistGrid)             :: tmp_distgrid
      CHARACTER(ESMF_Maxstr)          :: name
      real(ESMF_KIND_R8), dimension(:,:), pointer   :: tmp_ptr
      real                            :: totalsumi
      integer                         :: n, i, rc
!
      totalsumi = 1.0 / totalsum

      do n=1,dyn_items
        if( dyn_dim(n).gt.1 ) then
        do i=1,dyn_dim(n)
        dyn_array_save(i,n)=dyn_array_save(i,n)*totalsumi
        enddo

        CALL ESMF_StateGet(dyn_state, dyn_name(n), tmp_array, rc = rc)

        nullify(tmp_ptr)
        CALL ESMF_ArrayGet(tmp_array, 0, tmp_ptr, rc = rc) 

        do i=1,dyn_dim(n)
        tmp_ptr(i,1) = dyn_array_save(i,n)
        enddo

        CALL ESMF_StateAdd(dyn_state, tmp_array, rc = rc)
        endif
      enddo

      deallocate(dyn_name)
      deallocate(dyn_dim)
      deallocate(dyn_array_save)

      end subroutine digital_filter_dyn_average_gfs

! ---------------------------------------------------------------
! subroutines for physics
! ---------------------------------------------------------------
      subroutine digital_filter_phy_init_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(in) :: phy_state
!
      integer 			rc

      call esmf_stateget(state     = phy_state,				&
                         itemcount = phy_items,				&
                         rc=rc)
      allocate(phy_name(phy_items))
      call esmf_stateget(state=phy_state,				&
                         itemnamelist = phy_name,			&
                         rc=rc)
      phy_state_save=esmf_statecreate(statename="digital filter phy"	&
                                     ,statetype=esmf_state_unspecified	&
                                     ,rc       =rc)
!
      end subroutine digital_filter_phy_init_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_phy_save_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(in) :: phy_state
!
      TYPE(ESMF_Array)             :: tmp_array
      integer                      :: n, rc
!
      do n=1,phy_items
        CALL ESMF_StateGet(phy_state, phy_name(n), tmp_array, rc = rc)

        CALL ESMF_StateAdd(phy_state_save, tmp_array, rc = rc)
      enddo
      end subroutine digital_filter_phy_save_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_phy_restore_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(inout) :: phy_state
!
      TYPE(ESMF_Array)                :: tmp_array
      integer                         :: n, rc
!
      do n=1,phy_items
        CALL ESMF_StateGet(phy_state_save, phy_name(n), tmp_array, rc = rc)

        CALL ESMF_StateAdd(phy_state, tmp_array, rc = rc)
      enddo
      call esmf_statedestroy(phy_state_save,rc=rc)

      end subroutine digital_filter_phy_restore_gfs


      end module module_digital_filter_gfs
