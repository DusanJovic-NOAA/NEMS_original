      module module_digital_filter_nmm
!
! a generic digital filter for any model under ESMF 
!
! March 2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
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

! ---------------------------------------------------------------
! subroutine for dynamics
! ---------------------------------------------------------------
      subroutine digital_filter_dyn_init_nmm(dyn_state,ndfistep,NUM_WATER,NUM_TRACERS)
      USE MODULE_DM_PARALLEL
      type(esmf_state), intent(in) :: dyn_state   
      INTEGER, intent(in)          :: ndfistep
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS
      type(esmf_array)             :: tmp_array
      integer                      :: tmp_rank,dyn_items
      integer                      :: SPEC_MAX,rc,N
      character(20), allocatable   :: dyn_name(:)
      character(20)                :: state_name
      nstep = ndfistep 
      kstep = - nstep -1

      call esmf_stateget(state     = dyn_state                          &
                        ,name      = state_name                         &
                        ,itemcount = dyn_items                          &
                        ,rc=rc)
      allocate(dyn_name(dyn_items))
      allocate(name_save_2d(dyn_items))
      allocate(name_save_3d(dyn_items))
      allocate(name_save_4d(dyn_items))
      call esmf_stateget(state     = dyn_state                          &
                        ,itemnamelist = dyn_name                        &
                        ,rc=rc)
      DO N=1,dyn_items
        CALL ESMF_StateGet(dyn_state, dyn_name(N), tmp_array, rc=rc)

        call esmf_arrayget(tmp_array                                    &
                          ,rank              = tmp_rank                 &
                          ,rc                = rc)
      IF (tmp_rank == 2) THEN
         tot_rank_2d=tot_rank_2d+1
         name_save_2d(tot_rank_2d)=dyn_name(N)
      ENDIF
      IF (tmp_rank == 3) THEN
         tot_rank_3d=tot_rank_3d+1
         name_save_3d(tot_rank_3d)=dyn_name(N)
      ENDIF
      IF (tmp_rank ==4 ) THEN
         tot_rank_4d=tot_rank_4d+1  
         name_save_4d(tot_rank_4d)=dyn_name(N)
      ENDIF
      ENDDO 
      SPEC_MAX=MAX(NUM_TRACERS,NUM_WATER)
     
      IF (tot_rank_2d .GT. 0) THEN
      	allocate(array_save_2d(ITS:ITE,JTS:JTE,tot_rank_2d))
      ENDIF
      IF (tot_rank_3d .GT. 0) THEN
      	allocate(array_save_3d(ITS:ITE,JTS:JTE,LM,tot_rank_3d))
      ENDIF
      IF (tot_rank_4d .GT. 0) THEN
      	allocate(array_save_4d(ITS:ITE,JTS:JTE,LM,SPEC_MAX,tot_rank_4d))
      ENDIF
      deallocate(dyn_name)
      array_save_2d=0.
      array_save_3d=0.
      array_save_4d=0.
      totalsum=0.

      end subroutine digital_filter_dyn_init_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_sum_nmm(dyn_state,MEAN_FLAG,NUM_WATER,NUM_TRACERS)
      USE MODULE_DM_PARALLEL
!
      implicit none
          
      type(esmf_state), intent(in)  :: dyn_state 
      INTEGER(KIND=KINT), intent(in)          :: NUM_WATER,NUM_TRACERS
      LOGICAL(KIND=KLOG), intent(in)          :: MEAN_FLAG
      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,P,RC,RC_UPD
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
      REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D
!
      CHARACTER(20)    :: ARRAY_NAME
!
      TYPE(ESMF_Array) :: HOLD_ARRAY
!
      real                          :: sx, wx, digfil,prod
      integer			    :: NUM_SPEC
      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS
        kstep = kstep + 1
        sx     = acos(-1.)*kstep/nstep
        wx     = acos(-1.)*kstep/(nstep+1)
  
        if( kstep.ne.0 ) then
            digfil = sin(wx)/wx*sin(sx)/sx
        else
            digfil=1.
        endif 

        
        IF (MEAN_FLAG) THEN
            digfil=1.
        ENDIF
   
        print *,' in digital_filter_sum digfil = ',digfil

        totalsum = totalsum + digfil
  

        
        

      IF (tot_rank_2d  .GT. 0) THEN
      DO N=1,tot_rank_2d
      ARRAY_NAME=name_save_2d(N)
      NULLIFY(HOLD_2D) 

      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
      
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_2D                              &  !<-- Put the pointer here
                        ,rc       =RC)

        DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_2d(I,J,N)=array_save_2d(I,J,N)+digfil*HOLD_2D(I,J)
          ENDDO
        ENDDO
      ENDDO
      ENDIF

      IF (tot_rank_3d  .GT. 0 ) THEN
      DO N=1,tot_rank_3d
      ARRAY_NAME=name_save_3d(N)
      NULLIFY(HOLD_3D)

      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
      
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
       
      DO L=1,LM  
        DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_3d(I,J,L,N)=array_save_3d(I,J,L,N)+digfil*HOLD_3D(I,J,L)
          ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDIF

      IF (tot_rank_4d  .GT. 0) THEN
      DO N=1,tot_rank_4d
      ARRAY_NAME=name_save_4d(N)
      NULLIFY(HOLD_4D)
      IF (ARRAY_NAME == 'TRACERS') THEN
      NUM_SPEC=NUM_TRACERS
      ELSE IF (ARRAY_NAME == 'WATER') THEN
      NUM_SPEC=NUM_WATER
      ENDIF

      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
      
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_4D                              &  !<-- Put the pointer here
                        ,rc       =RC)
       
      DO P=1,NUM_SPEC
       DO L=1,LM
        DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_4d(I,J,L,P,N)=array_save_4d(I,J,L,P,N)+digfil*HOLD_4D(I,J,L,P)
          ENDDO
        ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDIF


      IF(RC_UPD==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYNAMICS UPDATE SUCCEEDED'
      ELSE
        WRITE(0,*)'DYNAMICS UPDATE FAILED RC_UPD=',RC_UPD
      ENDIF

      end subroutine digital_filter_dyn_sum_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_average_nmm(dyn_state,NUM_WATER,NUM_TRACERS)
!
      USE MODULE_DM_PARALLEL
      implicit none
      type(esmf_state), intent(inout) :: dyn_state
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS
      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,P,RC,RC_UPD
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
      REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D
!
      CHARACTER(20)    :: ARRAY_NAME
!
      TYPE(ESMF_Array) :: HOLD_ARRAY

      TYPE(ESMF_Array)                :: tmp_array
      TYPE(ESMF_DistGrid)             :: tmp_distgrid
      CHARACTER(ESMF_Maxstr)          :: name
      real, dimension(:,:), pointer   :: tmp_ptr
      real                            :: totalsumi
      integer                         :: NUM_SPEC

      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS


!
      totalsumi = 1.0 / totalsum   

      IF (tot_rank_2d  .GT. 0) THEN
      DO N=1,tot_rank_2d
         DO J=JTS,JTE
          DO I=ITS,ITE
             array_save_2d(I,J,N)=totalsumi*array_save_2d(I,J,N)
          ENDDO
         ENDDO
      ENDDO
      ENDIF
      IF (tot_rank_3d  .GT. 0) THEN
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
      IF (tot_rank_4d  .GT. 0) THEN
      DO N=1,tot_rank_4d
      ARRAY_NAME=name_save_4d(N)
      IF (ARRAY_NAME == 'TRACERS') THEN
      NUM_SPEC=NUM_TRACERS
      ELSE IF (ARRAY_NAME == 'WATER') THEN
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


      IF (tot_rank_2d  .GT. 0) THEN
      DO N=1,tot_rank_2d
      ARRAY_NAME=name_save_2d(N)
      NULLIFY(HOLD_2D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_2D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d(I,J,N)
          ENDDO
        ENDDO
        CALL ESMF_StateAdd(dyn_state,HOLD_ARRAY,rc) 
      ENDDO
      ENDIF
      IF (tot_rank_3d  .GT. 0) THEN
      DO N=1,tot_rank_3d
      ARRAY_NAME=name_save_3d(N)
      NULLIFY(HOLD_3D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO L=1,LM
        DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_3D(I,J,L)=array_save_3d(I,J,L,N)
          ENDDO
        ENDDO
      ENDDO
        CALL ESMF_StateAdd(dyn_state,HOLD_ARRAY,rc)
      ENDDO
      ENDIF
      IF (tot_rank_4d  .GT. 0) THEN 
      DO N=1,tot_rank_4d 
      ARRAY_NAME=name_save_4d(N)
      IF (ARRAY_NAME == 'TRACERS') THEN
      NUM_SPEC=NUM_TRACERS
      ELSE IF (ARRAY_NAME == 'WATER') THEN
      NUM_SPEC=NUM_WATER
      ENDIF
      NULLIFY(HOLD_4D)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =DYN_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_4D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
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
        CALL ESMF_StateAdd(dyn_state,HOLD_ARRAY,rc)
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


      end subroutine digital_filter_dyn_average_nmm
!----------------------------------------------------------------------------
      subroutine digital_filter_phy_init_nmm(phy_state)
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
      end subroutine digital_filter_phy_init_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_phy_save_nmm(phy_state)
!
      implicit none
      type(esmf_state), intent(in) :: phy_state
!
      TYPE(ESMF_Array)             :: tmp_array
      integer                      :: n, rc
!
      do n=1,phy_items
        CALL ESMF_StateGet(state=phy_state, itemName=phy_name(n), array=tmp_array, rc = rc)

        CALL ESMF_StateAdd(phy_state_save, tmp_array, rc = rc)
      enddo
      end subroutine digital_filter_phy_save_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_phy_restore_nmm(phy_state)
!
      implicit none
      type(esmf_state), intent(inout) :: phy_state
!
      TYPE(ESMF_Array)                :: tmp_array
      integer                         :: n, rc
!
      do n=1,phy_items
        CALL ESMF_StateGet(state=phy_state_save, itemName=phy_name(n), array=tmp_array, rc = rc)

        CALL ESMF_StateAdd(phy_state, tmp_array, rc = rc)
      enddo
      deallocate(phy_name)

      end subroutine digital_filter_phy_restore_nmm

      end module module_digital_filter_nmm
