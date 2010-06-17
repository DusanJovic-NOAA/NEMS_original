!-----------------------------------------------------------------------
!
      MODULE module_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!
! a generic digital filter for any model under ESMF 
!
!-----------------------------------------------------------------------
! March 2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
!-----------------------------------------------------------------------
!
      use esmf_mod
      use module_include
      implicit none

! ---------
! dynamics
! ---------
      real, allocatable, save :: array_save_2d(:,:,:)
      real, allocatable, save :: array_save_3d(:,:,:,:)
      real, allocatable, save :: array_save_4d(:,:,:,:,:)
      real             , save :: totalsum
      character(20), allocatable, save :: name_save_2d(:)
      character(20), allocatable, save :: name_save_3d(:)
      character(20), allocatable, save :: name_save_4d(:)
      integer,                    save :: tot_rank_2d=0 
      integer,                    save :: tot_rank_3d=0
      integer,                    save :: tot_rank_4d=0
      integer,                    save :: kstep, nstep
      
! ---------
! physics
! ---------
      type(esmf_state) , save :: phy_state_save
      character(20), allocatable, save :: phy_name(:)
      integer		phy_items


      contains

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_init_nmm(dyn_state                  &
                                            ,ndfistep                   &
                                            ,num_water                  &
                                            ,num_tracers)
!-----------------------------------------------------------------------
!
      use module_dm_parallel
!
      type(esmf_state), intent(in) :: dyn_state   
      integer, intent(in)          :: ndfistep
      integer, intent(in)          :: num_water,num_tracers
!
      type(esmf_field)             :: tmp_field 
      integer                      :: tmp_rank,dyn_items
      integer                      :: spec_max,rc,N
      character(20), allocatable   :: dyn_name(:)
      character(20)                :: state_name
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      nstep = ndfistep 
      kstep = - nstep -1

      call esmf_stateget(state     = dyn_state                          &
                        ,name      = state_name                         &
                        ,itemcount = dyn_items                          &
                        ,rc=rc)
!
      allocate(dyn_name(dyn_items))
      allocate(name_save_2d(dyn_items))
      allocate(name_save_3d(dyn_items))
      allocate(name_save_4d(dyn_items))
      call esmf_stateget(state     = dyn_state                          &
                        ,itemnamelist = dyn_name                        &
                        ,rc=rc)
      DO N=1,dyn_items
!
        CALL ESMF_StateGet(dyn_state, dyn_name(N), tmp_field, rc=rc)

        call esmf_fieldget(tmp_field                                    &
                          ,dimCount          = tmp_rank                 &
                          ,rc                = rc)
!
        IF (tmp_rank == 2) THEN
          tot_rank_2d=tot_rank_2d+1
          name_save_2d(tot_rank_2d)=dyn_name(N)
        ENDIF
!
        IF (tmp_rank == 3) THEN
          tot_rank_3d=tot_rank_3d+1
          name_save_3d(tot_rank_3d)=dyn_name(N)
        ENDIF
!
        IF (tmp_rank ==4 ) THEN
          tot_rank_4d=tot_rank_4d+1  
          name_save_4d(tot_rank_4d)=dyn_name(N)
        ENDIF
!
      ENDDO 
!
      SPEC_MAX=MAX(NUM_TRACERS,NUM_WATER)
     
      IF (tot_rank_2d > 0) THEN
      	allocate(array_save_2d(ITS:ITE,JTS:JTE,tot_rank_2d))
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
      	allocate(array_save_3d(ITS:ITE,JTS:JTE,LM,tot_rank_3d))
      ENDIF
!
      IF (tot_rank_4d > 0) THEN
      	allocate(array_save_4d(ITS:ITE,JTS:JTE,LM,SPEC_MAX,tot_rank_4d))
      ENDIF
!
      deallocate(dyn_name)
      array_save_2d=0.
      array_save_3d=0.
      array_save_4d=0.
      totalsum=0.
!
!-----------------------------------------------------------------------

      end subroutine digital_filter_dyn_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_sum_nmm(dyn_state                   &
                                           ,mean_on                     &
                                           ,num_water                   &
                                           ,num_tracers)
!-----------------------------------------------------------------------
      use module_dm_parallel
!
      type(esmf_state),intent(in)  :: dyn_state 
      integer(kind=kint),intent(in) :: mean_on,num_water,num_tracers
!
      integer(kind=kint) :: i,ii,j,jj,l,n,num_spec,p,rc,rc_upd
      real(kind=kfpt) :: digfil,prod,sx,wx
      real(kind=kfpt),dimension(:,:)    ,pointer :: hold_2d
      real(kind=kfpt),dimension(:,:,:)  ,pointer :: hold_3d
      real(kind=kfpt),dimension(:,:,:,:),pointer :: hold_4d
!
      character(20) :: field_name
!
      type(ESMF_Field) :: hold_field
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc    =esmf_success
      rc_upd=esmf_success
!
      kstep = kstep + 1
      sx     = acos(-1.)*kstep/nstep
      wx     = acos(-1.)*kstep/(nstep+1)
  
      if( kstep/=0)then
        digfil = sin(wx)/wx*sin(sx)/sx
      else
        digfil=1.
      endif 
        
      if(mean_on>0)then
        digfil=1.
      endif
!
      write(0,*)' In digital_filter_sum digfil = ',digfil

      totalsum = totalsum + digfil
!
      if(tot_rank_2d>0) then
!
        do n=1,tot_rank_2d
          field_name=name_save_2d(N)
          nullify(hold_2d) 

          call ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Field
                            ,itemName=FIELD_NAME                        &  !<-- Name of the
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
          call ESMF_FieldGet(field  =HOLD_FIELD                         &  !<-- Field that holds the data pointer
                            ,localDe=0                                  &
                            ,farray =HOLD_2D                            &  !<-- Put the pointer here
                            ,rc     =RC)
          do j=jts,jte
          do i=its,ite
            array_save_2d(i,j,n)=array_save_2d(i,j,n)+digfil*hold_2d(i,j)
          enddo
          enddo
        enddo
      endif
!
      if(tot_rank_3d>0)then
        do n=1,tot_rank_3d
          field_name=name_save_3d(N)
          nullify(hold_3d)
!
          call ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Array
                            ,itemName=FIELD_NAME                        &  !<-- Name of the Field to extract
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
          call ESMF_FieldGet(field  =HOLD_FIELD                         &  !<-- Field that holds the data pointer
                            ,localDe=0                                  &
                            ,farray =HOLD_3D                            &  !<-- Put the pointer here
                            ,rc     =RC)
!
          do l=1,lm  
            do j=jts,jte
            do i=its,ite
              array_save_3d(i,j,l,n)=array_save_3d(i,j,l,n)+digfil*hold_3d(i,j,l)
            enddo
            enddo
         enddo
        enddo
!
      endif
!
      if(tot_rank_4d>0)then
        do n=1,tot_rank_4d
          field_name=name_save_4d(N)
          nullify(hold_4d)
          if (field_name == 'TRACERS') then
            num_spec=num_tracers
          else if (field_name == 'WATER') then
            num_spec=num_water
          endif
!
          call ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Field
                            ,itemName=FIELD_NAME                        &  !<-- Name of the Field to extract
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
          call ESMF_FieldGet(field  =HOLD_FIELD                         &  !<-- Field that holds the data pointer
                            ,localDe=0                                  &
                            ,farray =HOLD_4D                            &  !<-- Put the pointer here
                            ,rc     =RC)
!
          do p=1,num_spec
            do l=1,lm
              do j=jts,jte
              do i=its,ite
                array_save_4d(i,j,l,p,n)=array_save_4d(i,j,l,p,n)+digfil*hold_4d(i,j,l,p)
              enddo
              enddo
            enddo
          enddo
        enddo
!
      endif


      if(rc_upd==esmf_success)then
!       write(0,*)'DYNAMICS UPDATE SUCCEEDED'
      else
        write(0,*)'DYNAMICS UPDATE FAILED RC_UPD=',rc_upd
      endif
!-----------------------------------------------------------------------

      end subroutine digital_filter_dyn_sum_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_average_nmm(dyn_state               &
                                               ,num_water               &
                                               ,num_tracers)
!-----------------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL
!
      type(esmf_state), intent(inout) :: dyn_state
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS
!
      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,P,RC,RC_UPD
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
      REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D
!
      CHARACTER(20)    :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD

      TYPE(ESMF_Field)                :: tmp_field
!
      CHARACTER(ESMF_Maxstr)          :: name
      real, dimension(:,:), pointer   :: tmp_ptr
      real                            :: totalsumi
      integer                         :: NUM_SPEC

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS
!
      totalsumi = 1.0 / totalsum   

      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_2d(I,J,N)=totalsumi*array_save_2d(I,J,N)
          ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          DO L=1,LM
            DO J=JTS,JTE
            DO I=ITS,ITE
              array_save_3d(I,J,L,N)=totalsumi*array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_4d > 0) THEN
        DO N=1,tot_rank_4d
          FIELD_NAME=name_save_4d(N)
          IF (FIELD_NAME == 'TRACERS') THEN
            NUM_SPEC=NUM_TRACERS
          ELSE IF (FIELD_NAME == 'WATER') THEN
            NUM_SPEC=NUM_WATER
          ENDIF
          DO P=1,NUM_SPEC
            DO L=1,LM
              DO J=JTS,JTE
              DO I=ITS,ITE
                array_save_4d(I,J,L,P,N)=totalsumi*array_save_4d(I,J,L,P,N)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          FIELD_NAME=name_save_2d(N)
          NULLIFY(HOLD_2D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Field
                            ,itemName=FIELD_NAME                        &  !<-- Name of the
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field that holds the data pointer
                            ,localDe  =0                                &
                            ,farray   =HOLD_2D                          &  !<-- Put the pointer here
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d(I,J,N)
          ENDDO
          ENDDO
          CALL ESMF_StateAdd(dyn_state,HOLD_FIELD,rc) 
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          FIELD_NAME=name_save_3d(N)
          NULLIFY(HOLD_3D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Field
                            ,itemName=FIELD_NAME                        &  !<-- Name of the
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field that holds the data pointer
                            ,localDe  =0                                &
                            ,farray   =HOLD_3D                          &  !<-- Put the pointer here
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO L=1,LM
            DO J=JTS,JTE
            DO I=ITS,ITE
              HOLD_3D(I,J,L)=array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO
          CALL ESMF_StateAdd(dyn_state,HOLD_FIELD,rc)
        ENDDO
      ENDIF
!
      IF (tot_rank_4d > 0) THEN 
        DO N=1,tot_rank_4d 
          FIELD_NAME=name_save_4d(N)
          IF (FIELD_NAME == 'TRACERS') THEN
            NUM_SPEC=NUM_TRACERS
          ELSE IF (FIELD_NAME == 'WATER') THEN
            NUM_SPEC=NUM_WATER
          ENDIF
          NULLIFY(HOLD_4D)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =DYN_STATE                         &  !<-- State that holds the Field
                            ,itemName=FIELD_NAME                        &  !<-- Name of the
                            ,field   =HOLD_FIELD                        &  !<-- Put extracted Field here
                            ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field that holds the data pointer
                            ,localDe  =0                                &
                            ,farray   =HOLD_4D                          &  !<-- Put the pointer here
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO P=1,NUM_SPEC
            DO L=1,LM
              DO J=JTS,JTE
              DO I=ITS,ITE
                HOLD_4D(I,J,L,P)=array_save_4d(I,J,L,P,N)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL ESMF_StateAdd(dyn_state,HOLD_FIELD,rc)
        ENDDO
      ENDIF


      deallocate(name_save_2d)
      deallocate(name_save_3d)
      deallocate(name_save_4d)
      deallocate(array_save_2d)
      deallocate(array_save_3d)
      deallocate(array_save_4d)
      tot_rank_2d=0
      tot_rank_3d=0
      tot_rank_4d=0
      kstep=0
      nstep=0
!-----------------------------------------------------------------------
      end subroutine digital_filter_dyn_average_nmm
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_init_nmm(phy_state)
!-----------------------------------------------------------------------
!
      type(esmf_state), intent(in) :: phy_state
!
      integer :: rc
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

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
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_save_nmm(phy_state)
!-----------------------------------------------------------------------
!
      type(esmf_state), intent(in) :: phy_state
!
      TYPE(ESMF_Field)             :: tmp_field
      integer                      :: n, rc
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      do n=1,phy_items
        CALL ESMF_StateGet(state   =phy_state                           &
                          ,itemName=phy_name(n)                         &
                          ,field   =tmp_field                           &
                          ,rc      =rc)

        CALL ESMF_StateAdd(state=phy_state_save                         &
                          ,field=tmp_field                              &
                          ,rc   =rc)
      enddo
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_save_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_restore_nmm(phy_state)
!-----------------------------------------------------------------------
!
      type(esmf_state), intent(inout) :: phy_state
!
      TYPE(ESMF_Field)                :: tmp_field
      integer                         :: n, rc
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      do n=1,phy_items
        CALL ESMF_StateGet(state=phy_state_save                         &
                          ,itemName=phy_name(n)                         &
                          ,field=tmp_field                              &
                          ,rc = rc)

        CALL ESMF_StateAdd(state=phy_state                              &
                          ,field=tmp_field                              &
                          ,rc   =rc)
      enddo
      deallocate(phy_name)

!-----------------------------------------------------------------------

      end subroutine digital_filter_phy_restore_nmm

!-----------------------------------------------------------------------
!
      end module module_digital_filter_nmm
!
!-----------------------------------------------------------------------
