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
      use module_exchange,only: halo_exch

!      type(esmf_config),save :: cf_1                                !<-- The config object

      implicit none

! ---------
! dynamics
! ---------
      real, allocatable, save :: array_save_2d(:,:,:)
      real, allocatable, save :: array_save_3d(:,:,:,:)
      real, allocatable, save :: array_save_4d(:,:,:,:,:)
      real, allocatable, save :: dolph_wgts(:)
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
                                            ,dt_int,dt_num,dt_den       &
                                            ,num_water                  &
                                            ,num_tracers)
!-----------------------------------------------------------------------
!
      use module_dm_parallel
!
      type(esmf_state), intent(in) :: dyn_state   
      integer, intent(in)          :: ndfistep
      integer, intent(in)          :: num_water,num_tracers
      integer, intent(in)          :: dt_int,dt_num,dt_den
!
      type(esmf_field)             :: tmp_field 
      integer                      :: tmp_rank,dyn_items                &
                                     ,dfihr
      integer                      :: spec_max,rc,n,m
      character(20), allocatable   :: dyn_name(:)
      character(20)                :: state_name
      real                         :: taus, dt
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      nstep = ndfistep 
      kstep = - nstep -1

      if (allocated(dolph_wgts)) deallocate(dolph_wgts)
      allocate(dolph_wgts(-nstep:nstep))

      dt=float(dt_int)+float(dt_num)/float(dt_den)

! hardwiring cutoff frequency based on length of filtering window

      taus=float(2*ndfistep)*dt
!     write(0,*) 'using pulled dt, taus (200% of tdfi) value is: ', dt, taus
!
      call dolph(dt,taus,nstep,dolph_wgts)
!

      call esmf_stateget(state     = dyn_state                          &
                        ,name      = state_name                         &
                        ,itemcount = dyn_items                          &
                        ,rc=rc)

      tot_rank_2d=0
      tot_rank_3d=0
      tot_rank_4d=0

      if (.not. allocated(dyn_name))                                    &
      allocate(dyn_name(dyn_items))

      if (.not. allocated(name_save_2d))                                &
      allocate(name_save_2d(dyn_items))

      if (.not. allocated(name_save_3d))                                &
      allocate(name_save_3d(dyn_items))

      if (.not. allocated(name_save_4d))                                &
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
!     
      IF (tot_rank_2d > 0 .and. .not. allocated(array_save_2d)) THEN
      	allocate(array_save_2d(ITS:ITE,JTS:JTE,tot_rank_2d))
      ENDIF
!
      IF (tot_rank_3d > 0 .and. .not. allocated(array_save_3d)) THEN
      	allocate(array_save_3d(ITS:ITE,JTS:JTE,LM,tot_rank_3d))
      ENDIF
!
      IF (tot_rank_4d > 0 .and. .not. allocated(array_save_4d)) THEN
      	allocate(array_save_4d(ITS:ITE,JTS:JTE,LM,SPEC_MAX,tot_rank_4d))
      ENDIF
!
      deallocate(dyn_name)
      array_save_2d=0.
      array_save_3d=0.
      array_save_4d=0.
      totalsum=0.
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
      logical :: dolph
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

!-----------------------------------------------------------------------
!     Future task: make this dolph switch logical a configure file item
!-----------------------------------------------------------------------

      dolph=.true.

      IF (dolph) THEN
        digfil=dolph_wgts(kstep)

      ELSE 
        sx = acos(-1.)*kstep/nstep
        wx = acos(-1.)*kstep/(nstep+1)
        if( kstep/=0)then
          digfil = sin(wx)/wx*sin(sx)/sx
        else
          digfil=1.
        endif 

        if(mean_on>0)then
          digfil=1.
        endif

      ENDIF
!
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

          CALL HALO_EXCH(hold_2d,1,2,2)


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

          CALL HALO_EXCH(hold_3d,LM,2,2)

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


      IF (tot_rank_2d > 0) THEN 
        deallocate(name_save_2d)
        deallocate(array_save_2d)
      ENDIF

      IF (tot_rank_3d > 0) THEN 
        deallocate(name_save_3d)
        deallocate(array_save_3d)
      ENDIF

      IF (tot_rank_4d > 0) THEN 
        deallocate(name_save_4d)
        deallocate(array_save_4d)
      ENDIF

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

      if (.not. allocated(phy_name)) then
        allocate(phy_name(phy_items))
      endif

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
        CALL ESMF_StateGet(state   =phy_state_save                      &
                          ,itemName=phy_name(n)                         &
                          ,field   =tmp_field                           &
                          ,rc      =rc)
        CALL ESMF_StateAdd(phy_state                                    &
                          ,tmp_field                                    &
                          ,rc = rc)
      enddo
      deallocate(phy_name)

!-----------------------------------------------------------------------

      end subroutine digital_filter_phy_restore_nmm

!-----------------------------------------------------------------------
!
      end module module_digital_filter_nmm
!
!-----------------------------------------------------------------------

   SUBROUTINE dolph(deltat, taus, m, window)

!     calculation of dolph-chebyshev window or, for short,
!     dolph window, using the expression in the reference:
!
!     antoniou, andreas, 1993: digital filters: analysis,
!     design and applications. mcgraw-hill, inc., 689pp.
!
!     the dolph window is optimal in the following sense:
!     for a given main-lobe width, the stop-band attenuation
!     is minimal; for a given stop-band level, the main-lobe
!     width is minimal.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)                  ::  m
      REAL, DIMENSION(0:2*M), INTENT(OUT)    ::  window
      REAL, INTENT(IN)                     :: deltat, taus

      ! local data
      integer, PARAMETER        :: NMAX = 5000
      REAL, dimension(0:NMAX)   :: t, w, time
      real, dimension(0:2*nmax) :: w2
      INTEGER                   :: NPRPE=0        ! no of pe
      CHARACTER*80              :: MES

      real    :: pi, thetas, x0, term1, term2, rr, r,db, sum, arg, sumw
      integer :: n, nm1, i, nt

      PI = 4*ATAN(1.D0)

!      print *, 'in dfcoef, deltat = ', deltat, 'taus=',taus

      N = 2*M+1
      NM1 = N-1

      THETAS = 2*PI*ABS(DELTAT/TAUS)
      X0 = 1/COS(THETAS/2)
      TERM1 = (X0 + SQRT(X0**2-1))**(FLOAT(N-1))
      TERM2 = (X0 - SQRT(X0**2-1))**(FLOAT(N-1))
      RR = 0.5*(TERM1+TERM2)
      R = 1/RR
      DB = 20*LOG10(R)


!      WRITE(0,'(1X,''DOLPH: M,N='',2I8)')M,N
!      WRITE(0,'(1X,''DOLPH: THETAS (STOP-BAND EDGE)='',F10.3)')THETAS
!      WRITE(0,'(1X,''DOLPH: R,DB='',2F10.3)')R, DB

      DO NT=0,M
         SUM = 1
         DO I=1,M
            ARG = X0*COS(I*PI/N)
            CALL CHEBY(T,NM1,ARG)
            TERM1 = T(NM1)
            TERM2 = COS(2*NT*PI*I/N)
            SUM = SUM + R*2*TERM1*TERM2
         ENDDO
         W(NT) = SUM/N
         TIME(NT) = NT
      ENDDO
!     fill in the negative-time values by symmetry.
      DO NT=0,M
         W2(M+NT) = W(NT)
         W2(M-NT) = W(NT)
      ENDDO

!     fill up the array for return
      SUMW = 0.
      DO NT=0,2*M
         SUMW = SUMW + W2(NT)
      ENDDO
!      WRITE(0,'(1X,''DOLPH: SUM OF WEIGHTS W2='',F10.4)')SUMW

      DO NT=0,2*M
         WINDOW(NT) = W2(NT)
      ENDDO

      RETURN

   END SUBROUTINE dolph


   SUBROUTINE cheby(t, n, x)

!     calculate all chebyshev polynomials up to order n
!     for the argument value x.

!     reference: numerical recipes, page 184, recurrence
!         t_n(x) = 2xt_{n-1}(x) - t_{n-2}(x) ,  n>=2.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: n
      REAL, INTENT(IN)     :: x
      REAL, DIMENSION(0:N) :: t

      integer  :: nn

      T(0) = 1
      T(1) = X
      IF(N.LT.2) RETURN
      DO NN=2,N
         T(NN) = 2*X*T(NN-1) - T(NN-2)
      ENDDO

      RETURN

   END SUBROUTINE cheby
