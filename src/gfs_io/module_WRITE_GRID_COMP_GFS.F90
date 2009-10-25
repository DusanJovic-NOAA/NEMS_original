!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_GRID_COMP_GFS
!
!-----------------------------------------------------------------------
!***  THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!***  DATA WAS PUT INTO THIS COMPONENT'S IMPORT STATE DESTINED FOR
!***  HISTORY OUTPUT.  THIS COMPONENT EXTRACTS THAT INFORMATION
!***  FROM THE IMPORT STATE WHOSE CONTENTS ARE SEEN ONLY BY THE
!***  FORECAST TASKS AND TRANSFERS 2D DATA TO GROUPS OF WRITE TASKS
!***  WHERE IT IS PARTIALLY REASSEMBLED.  THE WRITE TASKS THEN
!***  TRANSFER THEIR SUBSECTIONS TO THE LEAD WRITE TASK WHICH
!***  ASSEBMLES THE 2D DATA ONTO THE FULL DOMAIN AND WRITES OUT
!***  ALL SCALAR/1D/2D DATA TO A HISTORY FILE.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions in CPL_REGISTER  
!                                and CPL_INITIALIZE
!       14 Aug 2007:  T. Black - Revised CPL_RUN for general output
!                                selection and added documentation
!                                for users.
!       12 Sep 2007:  T. Black - Replaced the write component and the
!                                write gridded component with only
!                                a gridded component that contains
!                                quilting.
!          Mar 2008:  R. Vasic - Convert from ESMF 3.0.1 to 3.1.0
!       15 Aug 2008:  J. Wang  - Revised for addition of NEMS-IO
!       16 Sep 2008:  J. Wang  - Output array reverts from 3-D to 2-D
!       14 Oct 2008:  R. Vasic - Add restart capability
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_WRITE_INTERNAL_STATE_GFS
      USE MODULE_WRITE_ROUTINES_GFS,ONLY : FIRST_PASS_GFS               &
                                      ,WRITE_NEMSIO_OPEN 
!
      USE MODULE_GFS_MPI_DEF,ONLY    : NUM_PES_FCST,NUM_PES_WRT         &
                                      ,N_GROUP,PETLIST_WRITE            &
                                      ,MPI_COMM_COMP                    &
                                      ,MPI_COMM_INTER_ARRAY
      USE MODULE_GET_CONFIG_WRITE_GFS
      USE MODULE_ERR_MSG       ,ONLY : ERR_MSG,MESSAGE_CHECK
      USE MODULE_INCLUDE_GFS
      USE MODULE_NEMSIO
!
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_REGISTER_GFS
!
      PUBLIC :: WRITE_SETUP_GFS,WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_LENGTH_I1D=5000                            !<-- Max words in all 1-D integer history variables
      INTEGER,PARAMETER :: MAX_LENGTH_I2D=50000                           !<-- Max I,J points in each integer 2D subdomain
      INTEGER,PARAMETER :: MAX_LENGTH_R1D=25000                           !<-- Max words in all 1-D real history variables
      INTEGER,PARAMETER :: MAX_LENGTH_R2D=700000                          !<-- Max I,J points in each real 2D subdomain
      INTEGER,PARAMETER :: MAX_LENGTH_LOG=MAX_DATA_LOG                    !<-- Max logical variables
!
      INTEGER,SAVE      :: LAST_FCST_TASK                                 !<-- Rank of the last forecast task
      INTEGER,SAVE      :: LEAD_WRITE_TASK                                !<-- Rank of the lead (first) write task in this write group
      INTEGER,SAVE      :: LEAD_WRITE_TASK_INTER                         !<-- Rank of the lead (first) write task in this write group
      INTEGER,SAVE      :: LAST_WRITE_TASK                                !<-- Rank of the last write task the write group
      INTEGER,SAVE      :: NTASKS                                         !<-- # of write tasks in the current group + all forecast tasks
      INTEGER,SAVE      :: NWTPG                                          !<-- # of write tasks (servers) per group 
      LOGICAL,SAVE      :: QUILTING                                       !<-- # of write tasks (servers) per group 
!
      INTEGER,SAVE      :: MAXSIZE_I2D=MAX_LENGTH_I2D*MAX_DATA_I2D        !<-- Max size of the 2D integer history array on each task 
!
      INTEGER,SAVE      :: MAXSIZE_R2D=MAX_LENGTH_R2D*MAX_DATA_R2D        !<-- Max size of the 2D real history array on each task 
!
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: NCURRENT_GROUP             !<-- The currently active write group
!
!-----------------------------------------------------------------------
!
!gfs
      INTEGER         :: NFILE2WRT
      INTEGER,DIMENSION(:),allocatable       :: LFILE2WRT
      TYPE(NEMSIO_GFILE),save      :: NEMSIOFILE
!
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE_GFS),POINTER :: WRT_INT_STATE             ! The internal state pointer.
!
!-----------------------------------------------------------------------
      REAL(KIND=kind_evod)         :: btim,btim0
      REAL(KIND=kind_evod),PUBLIC,SAVE :: write_init_tim                     &
                                    ,write_run_tim                      &
                                    ,write_first_tim                    &
                                    ,write_send_data_tim                &
                                    ,write_get_fields_tim
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_REGISTER_GFS(WRT_COMP,RC_WRT)
! 
!-----------------------------------------------------------------------
!***  REGISTER THE WRITE COMPONENT'S 
!***  INITIALIZE, RUN, AND FINALIZE SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang  - Originator
!       30 Jun 2007:  T. Black - Modified to share same traits as
!                                rest of code.
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: WRT_COMP                     ! The write component
      INTEGER,INTENT(OUT)               :: RC_WRT                       ! Final return code
!     
!----------------------------------------------------------------------
!  
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_WRT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Initialize Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(WRT_COMP                          &  !<-- The write component
                                     ,ESMF_SETINIT                      &  !<-- Predefined subroutine type (INIT)
                                     ,WRT_INITIALIZE_GFS                &  !<-- User's subroutineName
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Run Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(WRT_COMP                          &  !<-- The write component
                                     ,ESMF_SETRUN                       &  !<-- Predefined subroutine type (RUN)
                                     ,WRT_RUN_GFS                       &  !<-- User's subroutineName
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Finalize Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
     CALL ESMF_GridCompSetEntryPoint(WRT_COMP                           &  !<-- The write component
                                    ,ESMF_SETFINAL                      &  !<-- Predefined subroutine type (FINALIZE)
                                    ,WRT_FINALIZE_GFS                   &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!
      IF(RC_WRT==ESMF_SUCCESS)THEN
        WRITE(6,*)"PASS: Write_Register."
      ELSE
        WRITE(6,*)"FAIL: Write_Register."
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_REGISTER_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_INITIALIZE_GFS(WRT_COMP                            &
                               ,IMP_STATE_WRITE                         &
                               ,EXP_STATE_WRITE                         &
                               ,CLOCK                                   &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!       29 Jun 2007:  T. Black - Generalize output types; add comment
!                                descriptions.
!       xx May 2009:  J. WANG  - Changed for GFS
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_GridComp),INTENT(INOUT) :: WRT_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE_WRITE  
!
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
!
      INTEGER,INTENT(OUT)               :: RC_INIT
!
!----------------------------------------------------------------------- 
!***  LOCAL VARIABLES
!----------------------------------------------------------------------- 
!
      INTEGER                                :: RC,ISTAT
      REAL(kind=8)                           :: timef
!
      TYPE(ESMF_VM)                          :: VM
      TYPE(WRITE_WRAP_GFS)                   :: WRAP
      TYPE(WRITE_INTERNAL_STATE_GFS),POINTER :: WRT_INT_STATE
!
!----------------------------------------------------------------------- 
!*********************************************************************** 
!----------------------------------------------------------------------- 
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!----------------------------------------------------------------------- 
!***  INITIALIZE THE WRITE COMPONENT TIMERS.
!----------------------------------------------------------------------- 
!
      write_init_tim=0.
      write_run_tim=0.
      write_first_tim=0.
      write_send_data_tim=0.
      write_get_fields_tim=0.
!
!----------------------------------------------------------------------- 
!***  ALLOCATE THE WRITE COMPONENT'S INTERNAL STATE.
!----------------------------------------------------------------------- 
!
      ALLOCATE(WRT_INT_STATE,stat=RC)
!
!----------------------------------------------------------------------- 
!***  ATTACH THE INTERNAL STATE TO THE WRITE COMPONENT.
!----------------------------------------------------------------------- 
!
      wrap%WRITE_INT_STATE=>WRT_INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the Write Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(WRT_COMP                       &  !<-- The write component
                                        ,WRAP                           &  !<-- Pointer to the internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------------------------- 
!***  RETRIEVE THE LOCAL VM.
!----------------------------------------------------------------------- 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine for this group of tasks
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------------------------- 
!***  INITIALIZE THE VALUE OF THE CURRENTLY ACTIVE WRITE GROUP.
!***  THIS WILL BE USED IN ESMF_Send/Recv AND THEREFORE MUST
!***  BE A CONTIGUOUS DATA ARRAY.
!----------------------------------------------------------------------- 
!
      IF(.NOT.ALLOCATED(NCURRENT_GROUP))THEN
        ALLOCATE(NCURRENT_GROUP(1),stat=ISTAT)
        NCURRENT_GROUP(1)=N_GROUP
      ENDIF
!
!----------------------------------------------------------------------- 
!***  EXTRACT THE TASK IDs AND THE NUMBER OF TASKS PRESENT.
!----------------------------------------------------------------------- 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs and Count from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The local VM
                     ,localPet=wrt_int_state%MYPE                       &  !<-- My task ID
                     ,petCount=wrt_int_state%NTASKS                     &  !<-- Number of MPI tasks present in current group (fcst+write)
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTASKS=wrt_int_state%NTASKS
!
!-----------------------------------------------------------------------
!***  RETRIEVE INFORMATION REGARDING OUTPUT FROM THE CONFIGURATION FILE.
!***  THIS INFORMATION CONTAINS MAINLY THE OUTPUT CHANNELS AND NAMES OF
!***  THE DISK FILES.
!-----------------------------------------------------------------------
!
      CALL GET_CONFIG_WRITE_GFS(WRT_COMP,WRT_INT_STATE,RC)              !<-- User's routine to extract configfile data
!
      NWTPG          =wrt_int_state%WRITE_TASKS_PER_GROUP
      IF(QUILTING) THEN
        LAST_FCST_TASK =NTASKS-NWTPG-1
        LEAD_WRITE_TASK=LAST_FCST_TASK+1
        LAST_WRITE_TASK=NTASKS-1
        LEAD_WRITE_TASK_INTER=NUM_PES_FCST
      ELSE
        LAST_FCST_TASK=NTASKS-1
        LEAD_WRITE_TASK=LAST_FCST_TASK
        LAST_WRITE_TASK=LAST_FCST_TASK
        LEAD_WRITE_TASK_INTER=NUM_PES_FCST-1
      ENDIF
!
!-----------------------------------------------------------------------
!***  ALL TASKS ALLOCATE BUFFER DATA ARRAYS THAT WILL HOLD ALL OF EACH
!***  TYPE OF HISTORY DATA AND WILL BE USED TO Send/Recv THAT DATA
!***  BETWEEN THE FORECAST TASKS THAT KNOW IT INITIALLY TO THE
!***  WRITE TASKS THAT OBTAIN IT FROM THE FORECAST TASKS FOR WRITING.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_I1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_I1D(MAX_LENGTH_I1D,wrt_int_state%num_file),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_I2D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_I2D(MAXSIZE_I2D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_R1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_R1D(MAX_LENGTH_R1D,wrt_int_state%num_file),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_R2D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_R2D(MAXSIZE_R2D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_LOG))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_LOG(MAX_LENGTH_LOG,wrt_int_state%num_file),stat=ISTAT)
      ENDIF
!
!-----------------------------------------------------------------------
!***  ALLOCATE DIMENSIONS AS 1-WORD ARRAYS SINCE ESMF NEEDS
!***  CONTIGUOUS DATA ARRAYS FOR ESMF_Sends/ESMF_Recvs WHEN
!***  THOSE DIMENSIONS ARE TRANSMITTED TO THE WRITE TASKS.
!-----------------------------------------------------------------------
!
        IF(.NOT.ALLOCATED(wrt_int_state%IM))THEN
          ALLOCATE(wrt_int_state%IM(1)                                  &
                  ,wrt_int_state%JM(1)                                  &
                  ,wrt_int_state%LM(1))
        ENDIF
!
!-----------------------------------------------------------------------
!***  THE NUMBER OF Attributes (FOR SCALARS AND 1D ARRAYS) AND
!***  Fields (FOR GRIDDED 2D ARRAYS) IN THE WRITE COMPONENT'S
!***  IMPORT STATE ARE NOT KNOWN A PRIORI.
!
!***  EVEN THOUGH THESE COUNTS ARE JUST SCALAR INTEGERS WE MUST
!***  ALLOCATE THEIR POINTERS TO LENGTH 1 SINCE THEY WILL BE
!***  USED IN ESMF_Send/Recv WHICH REQUIRE THEM TO BE CONTIGUOUS
!***  DATA ARRAYS.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%NCOUNT_FIELDS))THEN
!
        allocate(wrt_int_state%field_name(5000,wrt_int_state%num_file))
!
        allocate(wrt_int_state%NAMES_I1D_STRING(wrt_int_state%num_file))
        allocate(wrt_int_state%NAMES_I2D_STRING(wrt_int_state%num_file))
        allocate(wrt_int_state%NAMES_R1D_STRING(wrt_int_state%num_file))
        allocate(wrt_int_state%NAMES_R2D_STRING(wrt_int_state%num_file))
        allocate(wrt_int_state%NAMES_LOG_STRING(wrt_int_state%num_file))
!
        ALLOCATE(wrt_int_state%NCOUNT_FIELDS(wrt_int_state%num_file),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_I1D(wrt_int_state%num_file),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_I2D(wrt_int_state%num_file),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_R1D(wrt_int_state%num_file),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_R2D(wrt_int_state%num_file),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_LOG(wrt_int_state%num_file),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  ALL INTEGER QUANTITIES (AS 1D ARRAYS) AND 1D AND 2D REAL
!***  QUANTITIES WILL BE STRUNG TOGETHER IN SINGLE ARRAYS OF
!***  EACH PARTICULAR TYPE.  WE NEED TO ALLOCATE THE ARRAYS THAT WILL
!***  HOLD THE LENGTH OF EACH OF THE QUANTITIES IN THESE 'STRINGS'
!***  AS THE 'STRINGS' THEMSELVES.
!-----------------------------------------------------------------------
!
        ALLOCATE(wrt_int_state%LENGTH_DATA_I1D(100,wrt_int_state%num_file),stat=ISTAT)          !<-- Lengths of each individual 1-D integer array
        ALLOCATE(wrt_int_state%LENGTH_DATA_R1D(100,wrt_int_state%num_file),stat=ISTAT)          !<-- Lengths of each individual 1-D real array
        ALLOCATE(wrt_int_state%LENGTH_SUM_I1D(wrt_int_state%num_file),stat=ISTAT)             !<-- Length of string of data of ALL 1-D integer arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_R1D(wrt_int_state%num_file),stat=ISTAT)             !<-- Length of string of data of ALL 1-D real arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_LOG(wrt_int_state%num_file),stat=ISTAT)             !<-- Length of string of data of ALL logical variables
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  SET THE IO_BaseTime TO THE INITIAL CLOCK TIME.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the Output Base Time to the Initial Clock Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK                                &  !<-- The ESMF Clock
                        ,startTime=wrt_int_state%IO_BASETIME            &  !<-- The Clock's starting time
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SET THE FIRST HISTORY FILE'S TIME INDEX.
!-----------------------------------------------------------------------
!
      wrt_int_state%NFHOUR=0
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!!!     WRITE(0,*)"PASS: Write_Initialize."
      ELSE
        WRITE(0,*)"FAIL: Write_Initialize."
      ENDIF
!
!      write_init_tim=timef()-btim0
!
!----------------------------------------------------------------------- 
!
      END SUBROUTINE WRT_INITIALIZE_GFS
!
!----------------------------------------------------------------------- 
!####################################################################### 
!----------------------------------------------------------------------- 
!
      SUBROUTINE WRT_RUN_GFS(WRT_COMP                                   &
                        ,IMP_STATE_WRITE                                &
                        ,EXP_STATE_WRITE                                &
                        ,CLOCK                                          &
                        ,RC_RUN)
!
!----------------------------------------------------------------------- 
!***  THE RUN STEP FOR THE WRITE GRIDDED COMPONENT.  
!***  MOVE DATA INTENDED FOR HISTORY OUTPUT FROM THE IMPORT STATE
!***  TO THE WRITE TASKS.
!----------------------------------------------------------------------- 
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!       14 Aug 2007:  T. Black - Major revisions for generalized
!                                selectable history output with
!                                quilting.  Add descriptive comments.
!
!-----------------------------------------------------------------------
!
      USE ESMF_FieldGetMOD
      use module_gfs_mpi_def, only: WRITE_GROUPS, WRITE_TASKS_PER_GROUP
!
      TYPE(ESMF_GridComp),INTENT(IN) :: WRT_COMP
      TYPE(ESMF_Clock)   ,INTENT(IN) :: CLOCK
! 
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_WRITE                  !<-- The write component export state.
                                                                         !    Although it is loaded up only as output from
                                                                         !    this subroutine, its INTENT needs to be INOUT
                                                                         !    to function properly.
      INTEGER,INTENT(OUT)            :: RC_RUN 
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle) :: FILE_BUNDLE                       !<-- The history output data Bundle
!
      INTEGER,SAVE                          :: IH_INT =MPI_REQUEST_NULL &
                                              ,IH_REAL=MPI_REQUEST_NULL
!
      INTEGER,SAVE                          :: NPOSN_1,NPOSN_2
!
      INTEGER                               :: I,I1,IJ,IJG,IM,IMJM      &
                                              ,J,JM,K,L,JmaP            &
                                              ,MYPE,MYPE_LOCAL          &
                                              ,N,N1,N2,NF,NN            &
                                              ,NN_INTEGER,NN_REAL       &
                                              ,N_POSITION               &
                                              ,NUM_ATTRIB
!
      INTEGER,ALLOCATABLE                   :: IRCNT(:),IDISP(:)
!
      INTEGER                               :: DIM1,DIM2,FIELDSIZE,NBDR
!
      INTEGER                               :: IYEAR_FCST               &
                                              ,IMONTH_FCST              &
                                              ,IDAY_FCST                &
                                              ,IHOUR_FCST               &
                                              ,IMINUTE_FCST             &
                                              ,ISECOND_FCST             &
                                              ,ISECOND_NUM              &
                                              ,ISECOND_DEN
!
      INTEGER(KIND=ESMF_KIND_I8)            :: NTIMESTEP_ESMF
      INTEGER(KIND=kind_io4)                :: NTIMESTEP
!
      INTEGER                               :: NF_HOURS                 &
                                              ,NF_MINUTES               &
                                              ,NSECONDS                 &
                                              ,NSECONDS_NUM             &
                                              ,NSECONDS_DEN
!
      INTEGER                               :: ID_DUMMY                 &
                                              ,ID_RECV                  &
                                              ,NFCST_TASKS              &
                                              ,NFIELD                   &
                                              ,NUM_FIELD_NAMES    
!
      INTEGER                               :: ID_START,ID_END          &
                                              ,ISTART,IEND              &
                                              ,JSTART,JEND
!
      INTEGER                               :: JROW_FIRST,JROW_LAST     &
                                              ,JSTA_WRITE,JEND_WRITE
!
      INTEGER                               :: KOUNT_I2D                &
                                              ,KOUNT_I2D_DATA           &
                                              ,KOUNT_R2D                &
                                              ,KOUNT_R2D_DATA
!
      INTEGER                               :: RST_KOUNT_I2D            &
                                              ,RST_KOUNT_I2D_DATA       &
                                              ,RST_KOUNT_R2D            &
                                              ,RST_KOUNT_R2D_DATA
!
      INTEGER                               :: LENGTH                   &
                                              ,MPI_COMM,MPI_COMM2       &
                                              ,target_wrttask,mpi_commun
!
      INTEGER                               :: IO_HST_UNIT,IO_RST_UNIT
!
      INTEGER                               :: IERR,ISTAT,RC
!jw
      INTEGER                               :: nnext,nstart,nlat,jlat    &
                                              ,nfcst_tasks_send
      INTEGER                               :: NSTART_I2D,NSTART_R2D
      INTEGER                               :: JFCST_I2D,JFCST_R2D
      INTEGER                               :: NLAT_I2D,NLAT_R2D
      INTEGER                               :: NBDL,NFILE         
      INTEGER                               :: itr
      INTEGER                               :: ipt_lats_node_a,lats_node_a
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE)    :: JSTAT
!
      INTEGER,DIMENSION(:)  ,POINTER        :: WORK_ARRAY_I1D
      INTEGER,DIMENSION(:,:),POINTER        :: WORK_ARRAY_I2D
!
      REAL                                  :: NF_SECONDS               &
                                              ,SECOND_FCST
      REAL                                  :: DEGRAD
      REAL(esmf_kind_r8)                               :: zhour(1),pdryini(1)
      REAL(kind_evod)                       :: timef
!
      REAL(KIND=kind_evod),DIMENSION(:)  ,POINTER  :: WORK_ARRAY_R1D
!      REAL(KIND=kind_evod),DIMENSION(:,:),POINTER  :: WORK_ARRAY_R2D
      REAL(KIND=kind_io4),DIMENSION(:,:),POINTER  :: WORK_ARRAY_R2D
      REAL(KIND=kind_io4),DIMENSION(:)  ,POINTER  :: GLAT1D,GLON1D,TMP
!
      LOGICAL                               :: WRITE_LOGICAL
      LOGICAL,SAVE                          :: FIRST=.TRUE.
      LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: FILE_FIRST
!
      CHARACTER(ESMF_MAXSTR)                :: NAME,GFNAME
!
      TYPE(ESMF_Logical)                    :: WORK_LOGICAL
!
      TYPE(ESMF_VM)                         :: VM
      TYPE(ESMF_Grid),SAVE                  :: GRID1
      TYPE(ESMF_DELayout)                   :: MY_DE_LAYOUT
      TYPE(ESMF_Field)                      :: FIELD_WORK1
!
      TYPE(ESMF_InternArray)                  :: ARRAY_WORK
      TYPE(WRITE_WRAP_GFS)                        :: WRAP
      TYPE(WRITE_INTERNAL_STATE_GFS),POINTER      :: WRT_INT_STATE
      TYPE(ESMF_LOGICAL),DIMENSION(:),POINTER :: FIRST_IO_PE
      TYPE(ESMF_Time)                         :: CURRTIME
!
      TYPE(ESMF_TypeKind)      :: DATATYPE
!
!-----------------------------------------------------------------------
!
      real(kind=kind_evod) :: wait_time
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!
      RC    =ESMF_SUCCESS
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IT IS IMPORTANT TO NOTE THAT WHILE THE TASKS EXECUTING THIS
!***  STEP INCLUDE ALL FORECAST TASKS PLUS THOSE WRITE TASKS IN
!***  THIS WRITE GROUP, THE IMPORT STATE WAS FILLED IN THE DYNAMICS
!***  AND PHYSICS STEPS (DYN_INITIALIZE AND PHY_INITIALIZE) ONLY BY
!***  THE FORECAST TASKS.  THEREFORE ANY INFORMATION EXTRACTED FROM
!***  THE IMPORT STATE THAT IS NEEDED BY THE WRITE TASKS MUST BE
!***  SENT TO THEM VIA ESMF_Send/Recv.
!
!***  ALSO NOTE THAT HISTORY DATA CONSISTING OF SCALARS OR 1D ARRAYS
!***  ARE PRESENT IN THE WRITE COMPONENT'S IMPORT STATE AS Attributes.
!***  ALL 2D (GRIDDED) HISTORY DATA ARE PRESENT AS Fields.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE WRITE COMPONENT'S ESMF INTERNAL STATE.
!-----------------------------------------------------------------------
!
      btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Write Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(WRT_COMP                       &  !<-- The write component
                                        ,WRAP                           &  !<-- Pointer to internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      WRT_INT_STATE=>wrap%WRITE_INT_STATE                                  !<-- Local working pointer to internal state
!
!-----------------------------------------------------------------------
!***  GET THE CURRENT LOCAL VM.
!***  THIS COMES FROM THE PetList USED TO CREATE
!***  THE WRITE COMPONENTS IN ATM_INITIALIZE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Current VM for WRT_RUN"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(VM,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      MYPE=wrt_int_state%MYPE
!-----------------------------------------------------------------------
!***  INCREMENT THE WRITE GROUP SO WE KNOW WHICH ONE IS ACTIVE.
!***  ONLY THE FORECAST TASKS ENTER THIS RUN STEP OF THE WRITE
!***  COMPONENT EVERY OUTPUT TIME THEREFORE ONLY THEY CAN PROPERLY
!***  INCREMENT THE NUMBER OF THE CURRENT WRITE GROUP. 
!***  LET FORECAST TASK 0 BROADCAST THE CURRENT GROUP NUMBER
!***  TO ALL TASKS IN THE GROUP INCLUDING THE CURRENT WRITE TASKS.
!-----------------------------------------------------------------------
!
      NCURRENT_GROUP(1)=N_GROUP
!jw      IF(QUILTING) THEN
!jw        LAST_FCST_TASK =NTASKS-NWTPG-1
!jw       LEAD_WRITE_TASK=LAST_FCST_TASK+1
!jw       LAST_WRITE_TASK=NTASKS-1
!jw        LEAD_WRITE_TASK_INTER=NUM_PES_FCST
!jw      ELSE
!jw        LAST_FCST_TASK=NTASKS-1
!jw        LEAD_WRITE_TASK=LAST_FCST_TASK
!jw        LAST_WRITE_TASK=LAST_FCST_TASK
!jw        LEAD_WRITE_TASK_INTER=NUM_PES_FCST-1
!jw      ENDIF
!
      MYPE_LOCAL=MYPE
      if(wrt_int_state%quilting .and. MYPE>=LAST_FCST_TASK.and.NWTPG>0) &
         MYPE_LOCAL=MOD(MYPE-NUM_PES_FCST,NWTPG)
      if(.not.wrt_int_state%quilting .and. MYPE==LEAD_WRITE_TASK) &
         MYPE_LOCAL=0
!
!jw      write(0,*)'in write_run, NWTPG=',NWTPG,'LAST_FCST_TASK=',         &
!jw       LAST_FCST_TASK,'LEAD_WRITE_TASK=',LEAD_WRITE_TASK,               &
!jw       'LAST_WRITE_TASK=',LAST_WRITE_TASK,'mype_local=',MYPE_LOCAL,     &
!jw       'NCURRENT_GROUP(1)=',NCURRENT_GROUP(1),'LEAD_WRITE_TASK_INTER=', &
!jw       LEAD_WRITE_TASK_INTER
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Broadcast Current Write Group Number"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMBroadcast(VM                                          &  !<-- The local VM
                           ,NCURRENT_GROUP                              &  !<-- The current active write group
                           ,1                                           &  !<-- # of elements to broadcast
                           ,0                                           &  !<-- Root sender is fcst task 0
                           ,rc=RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!jw      write(0,*)'after broadcast,NCURRENT_GROUP=',NCURRENT_GROUP(1)
!-----------------------------------------------------------------------
!***  THE WRITE COMPONENT IS EXECUTED AT EACH HISTORY OUTPUT TIME
!***  BY ALL FORECAST TASKS PLUS BY ALL WRITE TASKS IN THE ACTIVE
!***  WRITE GROUP.  A 'FIRST' SWITCH IS EMPLOYED BELOW SO THAT THE
!***  FORECAST TASKS EXTRACT FUNDAMENTAL GRID INFORMATION AND SEND IT
!***  TO THE WRITE TASKS ONLY THE FIRST TIME THAT EACH OF THE WRITE
!***  GROUPS EXECUTES THIS ROUTINE.  FOR EXAMPLE, ASSUME THERE ARE
!***  TWO WRITE GROUPS.  AT THE 1ST OUTPUT TIME, THE FORECAST TASKS
!***  EXTRACT CERTAIN INFORMATION FROM THE IMPORT STATE THAT DOES NOT
!***  CHANGE WITH TIME AND TASK 0 SENDS THAT INFORMATION TO EACH OF 
!***  THE WRITE TASKS IN THE 1ST WRITE GROUP.  AT THE 2ND OUTPUT TIME
!***  THE FORECAST TASKS EXTRACT THAT SAME INFORMATION AND TASK 0
!***  SENDS IT TO THE WRITE TASKS IN THE 2ND WRITE GROUP SINCE THE
!***  WRITE GROUPS ARE CYCLING.  AT THE 3RD OUTPUT TIME AND ALL
!***  SUBSEQUENT OUTPUT TIMES THE EXTRACTION AND SENDING/RECEIVING
!***  OF THIS INFORMATION IS SKIPPED SINCE ALL OF THE WRITE TASKS
!***  (OR AT LEAST THE FIRST WRITE TASK IN EACH WRITE GROUP) ALREADY
!***  HAVE THE INFORMATION.
!-----------------------------------------------------------------------
!
      if(.not.allocated(FILE_FIRST)) then
        allocate(FILE_FIRST(wrt_int_state%num_file))
        file_first=.true.
      endif
!
!-----------------------------------------------------------------------
!*** loop on the files that need to write out
!-----------------------------------------------------------------------
      file_loop: Do NBDL=1, wrt_int_state%num_file
!
!
       if(mype<num_pes_fcst) then
         CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE             &  !<-- The Write component's import state
                             ,name      ='ipt_lats_node_a'             &  !<-- Name of the Attribute to extract
                             ,value =wrt_int_state%ipt_lats_node_a     &  !<-- Extract local subdomain starting I's
                             ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE               &  !<-- The Write component's import state
                             ,name     ='lats_node_a'                  &  !<-- Name of the Attribute to extract
                             ,value = wrt_int_state%lats_node_a        &  !<-- Extract local subdomain ending I's
                             ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE               &  !<-- The Write component's import state
                             ,name     ='im_'//trim(wrt_int_state%filename_base(nbdl)) &  !Name of the Attribute to extract
                             ,value = wrt_int_state%im(1)        &  !<-- Extract local subdomain ending I's
                             ,rc=RC)

        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE               &  !<-- The Write component's import state
                             ,name     ='jm_'//trim(wrt_int_state%filename_base(nbdl)) &  !Name of the Attribute to extract
                             ,value = wrt_int_state%jm(1)              &  !<-- Extract local subdomain ending I's
                             ,rc=RC)

        IF(.NOT.ALLOCATED(wrt_int_state%global_lats_a))THEN
         ALLOCATE(wrt_int_state%global_lats_a(1:wrt_int_state%jm(1)),stat=ISTAT)  !<-- Local starting I for each fcst task's subdomain
        ENDIF
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE              &  !<-- The Write component's import state
                              ,name      ='global_lats_a'              &  !<-- Name of the Attribute to extract
                              ,count     =wrt_int_state%jm(1)          &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%global_lats_a  &  !<-- Extract local subdomain starting J's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE              &  !<-- The Write component's import state
                              ,name      ='zhour'                      &  !<-- Name of the Attribute to extract
                              ,value     =zhour(1)                     &  !<-- Extract local subdomain starting J's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE              &  !<-- The Write component's import state
                              ,name      ='pdryini'                    &  !<-- Name of the Attribute to extract
                              ,value     =pdryini(1)                   &  !<-- Extract local subdomain starting J's
                              ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Broadcast pdryini,zhour"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMBroadcast(VM                                        &  !<-- The local VM
                           ,zhour                                       &  !<-- The current active write group
                           ,1                                           &  !<-- # of elements to broadcast
                           ,0                                           &  !<-- Root sender is fcst task 0
                           ,rc=RC)
        wrt_int_state%zhour=zhour(1)
!
        CALL ESMF_VMBroadcast(VM                                        &  !<-- The local VM
                           ,pdryini                                     &  !<-- The current active write group
                           ,1                                           &  !<-- # of elements to broadcast
                           ,0                                           &  !<-- Root sender is fcst task 0 
                           ,rc=RC)
        wrt_int_state%pdryini=pdryini(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(FIRST_IO_PE(1))                                            !<-- A flag indicating that this is or is not
                                                                          !    this write component
!
      IF( FILE_FIRST(NBDL) )THEN
        FIRST_IO_PE(1)=ESMF_TRUE
      ELSE
        FIRST_IO_PE(1)=ESMF_FALSE
      ENDIF
!
!-----------------------------------------------------------------------
!***  HERE THE LEAD WRITE TASK IS TELLING ALL TASKS INCLUDING THE
!***  FORECAST TASKS WHO EXTRACT INFORMATION FROM THE IMPORT STATE
!***  THAT THIS IS OR IS NOT THIS SET OF WRITE TASKS' FIRST PASS
!***  THROUGH THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Broadcast FIRST_PASS Status"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMBroadcast(VM                                          &  !<-- The local VM
                           ,FIRST_IO_PE                                 &  !<-- The 1st write task in this group tells everyone
                                                                           !    if this is or is not the write tasks' first pass
                                                                           !    through this routine.
                           ,1                                           &  !<-- # of elements to broadcast
                           ,LEAD_WRITE_TASK_INTER                       &  !<-- Root sender is the first write task in this group
                           ,rc=RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(FIRST_IO_PE(1)==ESMF_TRUE)THEN
        FIRST=.TRUE.
      ELSE
        FIRST=.FALSE.
      ENDIF
!
      DEALLOCATE(FIRST_IO_PE)
!
!-----------------------------------------------------------------------
!***  CERTAIN WORK NEEDS TO BE DONE ONLY THE FIRST TIME THAT EACH
!***  GROUP OF WRITE TASKS IS INVOKED.  THIS MOSTLY CONSISTS OF
!***  THE FORECAST TASKS' SENDING THE WRITE TASKS HISTORY DATA
!***  FROM THE IMPORT STATE THAT DOES NOT CHANGE WITH FORECAST TIME.
!-----------------------------------------------------------------------
!
      first_block: IF(FIRST)THEN                                           !<-- Execute this routine only if this is the 
                                                                           !    1st pass by the current set of write tasks.
!
!-----------------------------------------------------------------------
!
          CALL FIRST_PASS_GFS(IMP_STATE_WRITE                           &
                       ,WRT_INT_STATE,NTASKS,MYPE,NBDL                  &
                       ,NCURRENT_GROUP(1))
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS PLUS THE WRITE TASKS IN THIS WRITE GROUP
!***  NOW TURN OFF THEIR 'FIRST' SWITCH.  HOWEVER THIS ONLY MATTERS
!***  FOR THE WRITE TASKS SINCE IT IS THEIR VALUE OF 'FIRST' THAT
!***  IS TRANSMITTED TO THE FORECAST TASKS JUST PRIOR TO THIS BLOCK
!***  BY THE 1ST WRITE TASK IN THE WRITE GROUP THAT IS PRESENT.
!-----------------------------------------------------------------------
!
         FILE_FIRST(NBDL)=.FALSE.
!
!-----------------------------------------------------------------------
!
      ENDIF first_block
!
      write_first_tim=write_first_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  THE ELAPSED FORECAST TIME (HOURS) WILL BE APPENDED TO THE NAME
!***  OF EACH HISTORY OUTPUT FILE.  EXTRACT THAT VALUE NOW.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Current Time for Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK                                 &  !<-- The ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current time (ESMF) on the clock
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ESMF TIME DIFFERENCE BETWEEN START TIME AND CURRENT TIME
!-----------------------------------------------------------------------
!
      wrt_int_state%IO_CURRTIMEDIFF=CURRTIME-wrt_int_state%IO_BASETIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Elapsed Forecast Time for Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                               ,h           =wrt_int_state%NFHOUR          &  !<-- The elapsed time in hours (REAL)
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
!
      KOUNT_I2D=wrt_int_state%KOUNT_I2D(NBDL)
      KOUNT_R2D=wrt_int_state%KOUNT_R2D(NBDL)
      NFCST_TASKS=NUM_PES_FCST                                             !<-- Number of fcst tasks sending to this write task
!
!-----------------------------------------------------------------------
      fcst_tasks: IF( MYPE<NFCST_TASKS)THEN      !<-- Only the forecast tasks can see this data so far
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  write out NFILE by getting info from BUNDLE
!***  Bundle name is file name
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract History Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   =wrt_int_state%filename_base(NBDL)     &  !<-- The name of the history data Bundle
                          ,fieldbundle=FILE_BUNDLE                      &  !<-- The history data Bundle inside the import state
                          ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW PULL THE 2D DATA FROM THE IMPORT STATE.
!***  THIS INCLUDES ALL INDIVIDUAL 2D HISTORY QUANTITIES AS WELL AS
!***  ALL MODEL LEVELS OF THE 3D REAL HISTORY ARRAYS.
!-----------------------------------------------------------------------
        NN_INTEGER=0
        NN_REAL   =0
!
        ipt_lats_node_a=wrt_int_state%ipt_lats_node_a
        lats_node_a=wrt_int_state%lats_node_a
!
!-----------------------------------------------------------------------
!***  BE SURE THE INTEGER AND REAL BUFFERS ARE AVAILABLE FOR ISENDs
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL MPI_WAIT(IH_INT,JSTAT,IERR) 
        wait_time=timef()-btim
        if(wait_time>1.e3)write(0,*)' Long integer buffer WAIT =',wait_time*1.e-3
!
        btim=timef()
        CALL MPI_WAIT(IH_REAL,JSTAT,IERR) 
        wait_time=timef()-btim
        if(wait_time>1.e3)write(0,*)' Long real buffer WAIT =',wait_time*1.e-3
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        field_block: DO N=1,wrt_int_state%NCOUNT_FIELDS(NBDL)                 !<-- Loop through all Fields in the import state
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Fields from History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(bundle=FILE_BUNDLE                      &  !<-- The write component's history data Bundle
                                  ,name  =wrt_int_state%FIELD_NAME(N,NBDL)   &  !<-- The ESMF Field's name
                                  ,field =FIELD_WORK1                   &  !<-- The ESMF Field data pointer
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DOES THIS EXTRACTED Field HOLD INTEGER OR REAL DATA?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Check Datatype of Field from History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =FIELD_WORK1                       &  !<-- The ESMF Field
                            ,typekind=DATATYPE                          &  !<-- ESMF specifier of variable type and kind
                            ,rc      =RC)
 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------------------
!                      -- INTEGER FIELDS --
!--------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Integer Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field  =FIELD_WORK1                      &  !<-- The ESMF Field
                              ,localDe=0                                &  !<-- # of DEs in this grid
                              ,farray =WORK_ARRAY_I2D                   &  !<-- Put the 2D integer data from the Field here
                              ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ISTART=LBOUND(WORK_ARRAY_I2D,1)
            IEND  =UBOUND(WORK_ARRAY_I2D,1)
            JSTART=LBOUND(WORK_ARRAY_I2D,2)
            JEND  =UBOUND(WORK_ARRAY_I2D,2)
!*** allow buff_mult has larger subdomain
            if(JEND>wrt_int_state%lats_node_a)JEND=wrt_int_state%lats_node_a
!
            IF(NN_INTEGER+IEND*JEND>MAXSIZE_I2D)THEN
              WRITE(0,*)' WARNING:  YOU MUST INCREASE THE SIZE OF'      &
                       ,' MAX_LENGTH_I2D OR MAX_DATA_I2D BECAUSE'       &
                       ,' YOU HAVE EXCEEDED THE DECLARED SIZE OF'       &
                       ,' ALL_DATA_I2D'
            ENDIF
!
!
!-----------------------------------------------------------------------
!* rearrange data to im*nrecord*lats(WriteTask)) in order to send chunk of 
!* data to write task data(im*nrecord*lats(mytask)
!-----------------------------------------------------------------------
!
            if(NTASKS>1) then
              NNEXT=(IEND-ISTART+1)*KOUNT_I2D
              DO J=JSTART,JEND
                JMAP=wrt_int_state%nwrttask_on_fcst(J)
                N1=(J-JSTART)*NNEXT
                wrt_int_state%ALL_DATA_I2D(NN_INTEGER+1+N1:              &
                  NN_INTEGER+IEND-ISTART+1+N1)=                          &
                  WORK_ARRAY_I2D(ISTART:IEND,JMAP)     !<-- String together this task's 2D real data
              ENDDO
              NN_INTEGER=NN_INTEGER+IEND-ISTART+1   !<-- for next record
!
            else
!
              wrt_int_state%ALL_DATA_I2D(NN_INTEGER+1:NN_INTEGER+IM*JM)= &
                reshape(WORK_ARRAY_I2D,(/im*jm/) )
              NN_INTEGER=NN_INTEGER+IM*JM
            endif
!
!--------------------------------------------------------------------
!                        -- REAL FIELDS --
!--------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           !<-- Extract real gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Real Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field  =FIELD_WORK1                      &  !<-- The ESMF Field
                              ,localDe=0                                &  !<-- # of DEs in this grid
                              ,farray =WORK_ARRAY_R2D                   &  !<-- Put the 2D real data from the Field here
                              ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ISTART=LBOUND(WORK_ARRAY_R2D,1)
            IEND  =UBOUND(WORK_ARRAY_R2D,1)
            JSTART=LBOUND(WORK_ARRAY_R2D,2)
            JEND  =UBOUND(WORK_ARRAY_R2D,2)
!jw allow buff_mult has larger subdomain
            if(JEND>wrt_int_state%lats_node_a)JEND=wrt_int_state%lats_node_a

            IF(NN_REAL+IEND*JEND>MAXSIZE_R2D)THEN
              WRITE(0,*)' WARNING:  YOU MUST INCREASE THE SIZE OF'      &
                       ,' MAX_LENGTH_R2D OR MAX_DATA_R2D BECAUSE'       &
                       ,' YOU HAVE EXCEEDED THE DECLARED SIZE OF'       &
                       ,' ALL_DATA_R2D'
            ENDIF
!
!*** regroup to send to write processes
!
            if(NTASKS>1) then
            NNEXT=(IEND-ISTART+1)*KOUNT_R2D
            DO J=JSTART,JEND
              JMAP=wrt_int_state%nwrttask_on_fcst(J)
              N1=(J-JSTART)*NNEXT
              wrt_int_state%ALL_DATA_R2D(NN_REAL+1+N1:NN_REAL+IEND-ISTART+1+N1)=  &
                WORK_ARRAY_R2D(ISTART:IEND,JMAP)                                  !<-- String together this task's 2D real data
            ENDDO
            NN_REAL=NN_REAL+IEND-ISTART+1

            else
             wrt_int_state%ALL_DATA_R2D(NN_REAL+1:NN_REAL+IM*JM)=          &
               reshape(WORK_ARRAY_R2D,(/im*jm/) )
             NN_REAL=NN_REAL+(IEND-ISTART+1)*(JEND-JSTART+1)
            endif
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO field_block
!
        write_get_fields_tim=write_get_fields_tim+timef()-btim
!
        btim=timef()
!
!-----------------------------------------------------------------------
!
      ifntasks1: IF(NTASKS > 1 .and. mype/= lead_write_task) THEN
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS NOW SEND THEIR STRINGS OF 2D HISTORY DATA
!***  TO THE APPROPRIATE WRITE TASKS.
!*** send data to write tasks
!-----------------------------------------------------------------------
!
         DO N=lead_write_task,last_write_task
!
!-----------------------------------------------------------------------
!***  FIRST THE 2-D INTEGER DATA.
!-----------------------------------------------------------------------
!
            N1=mod(n-nfcst_tasks,nwtpg)+1
            if(wrt_int_state%quilting) then
              target_wrttask=N1-1
              mpi_commun=MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))
            else
              target_wrttask=lead_write_task
              mpi_commun=MPI_COMM_COMP
            endif
!
            NLAT_I2D=wrt_int_state%NLAT_TO_WRITE_TASK(N1)*IM*KOUNT_I2D
            NSTART_I2D=(wrt_int_state%NSTART_TO_WRITE_TASK(N1)-1)*IM*KOUNT_I2D+1
            IF(NLAT_I2D>0)THEN
!
              CALL MPI_ISEND(wrt_int_state%ALL_DATA_I2D(NSTART_I2D)     &  !<-- Fcst tasks' string of 2D integer history data
                            ,NLAT_I2D                                   &  !<-- #of words in the data string
                            ,MPI_INTEGER                                &  !<-- The datatype
                            ,target_wrttask                             &  !<-- The target write task
                            ,wrt_int_state%NFHOUR*1000+mype             &  !<-- An MPI tag
                            ,mpi_commun                                 &  !<-- The MPI intercommunicator between fcst and quilt tasks
                            ,IH_INT                                     &  !<-- MPI communication request handle
                            ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of integer data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
!-----------------------------------------------------------------------
!***  THEN THE 2-D REAL DATA.
!-----------------------------------------------------------------------
!
            NLAT_R2D=wrt_int_state%NLAT_TO_WRITE_TASK(N1)*IM*KOUNT_R2D
            NSTART_R2D=(wrt_int_state%NSTART_TO_WRITE_TASK(N1)-1)*IM*KOUNT_R2D+1
            IF(NLAT_R2D>0)THEN
              CALL MPI_ISEND(wrt_int_state%ALL_DATA_R2D(NSTART_R2D)       &  !<-- Fcst tasks' string of 2D real history data
                            ,NLAT_R2D                                     &  !<-- #of words in the data string
                            ,MPI_REAL                                     &  !<-- The datatype
                            ,target_wrttask                               &  !<-- The target write task
                            ,wrt_int_state%NFHOUR*1000+mype+N1*100        &  !<-- An MPI tag
                            ,mpi_commun                                   &  !<-- The MPI intercommunicator between fcst and quilt tasks
                            ,IH_REAL                                      &  !<-- MPI communication request handle
                            ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of real data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
        ENDDO
        write_send_data_tim=write_send_data_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF  ifntasks1
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks
!
      write_run_tim=write_run_tim+timef()-btim0
!rst-rst ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!-----------------------------------------------------------------------
!***  THE FORECAST TASKS FINISHED SENDING DATA TO WRITE TASKS, THEY
!***  WILL GO TO THE NEXT BUNDLE LOOP CYCLE. BUT FOR 1 TASK, THE FORECAST 
!***  TASK WILL NEED TO WRITE OUT THE OUTPUT FILE.
!-----------------------------------------------------------------------
!
      IF(MYPE<LEAD_WRITE_TASK ) then
         cycle file_loop
      ENDIF
!
!-----------------------------------------------------------------------
!
      ifntasks3:  IF (NTASKS>1) THEN
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK IN THE ACTIVE WRITE GROUP RECEIVES THE
!***  STRINGS OF 2D HISTORY DATA FROM THE APPROPRIATE FCST TASKS.
!-----------------------------------------------------------------------
!
       JSTA_WRITE=wrt_int_state%JSTART_WRITE(MYPE_LOCAL+1)                 !<-- JTS of this write task
       JEND_WRITE=wrt_int_state%JEND_WRITE(MYPE_LOCAL+1)                   !<-- JTE to this write task
!     
       ALLOCATE(wrt_int_state%WRITE_SUBSET_I(1:IM,JSTA_WRITE:JEND_WRITE  &
                                          ,wrt_int_state%KOUNT_I2D(nbdl)))
       ALLOCATE(wrt_int_state%WRITE_SUBSET_R(1:IM,JSTA_WRITE:JEND_WRITE  &
                                          ,wrt_int_state%KOUNT_R2D(nbdl)))
!
!-----------------------------------------------------------------------
       NLAT_I2D=0
       NLAT_R2D=0
       JFCST_I2D=0
       JFCST_R2D=0
       if(wrt_int_state%quilting) then
         mpi_commun=MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))
         NFCST_TASKS_SEND=NFCST_TASKS
       else
         mpi_commun=MPI_COMM_COMP
         NFCST_TASKS_SEND=NFCST_TASKS-1
       endif
!
!-----------------------------------------------------------------------
!*** IF no Quilting, set all_data_I/R2D to write_subset_i/r 
!-----------------------------------------------------------------------
!
       if(.not.wrt_int_state%quilting) then
!
          JSTART=wrt_int_state%nstart_from_fcst_task(lead_write_task+1)
!
          IF(KOUNT_I2D>0) THEN
            NN=0
            DO J=1,wrt_int_state%nlat_from_fcst_task(lead_write_task+1)
              JLAT=wrt_int_state%fcst_lat_on_write_task(JSTART+J-1)
              DO NF=1,KOUNT_I2D                                                 !<-- Loop through all the 2D integer fields
              DO I=1,IM
                NN=NN+1
                wrt_int_state%WRITE_SUBSET_I(I,JLAT,NF)=wrt_int_state%ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(KOUNT_R2D>0) THEN
            NN=0
            DO J=1,wrt_int_state%nlat_from_fcst_task(lead_write_task+1)
             JLAT=wrt_int_state%fcst_lat_on_write_task(JSTART+J-1)
             DO NF=1,KOUNT_R2D                                                 !<-- Loop through all the 2D real fields
             DO I=1,IM
              NN=NN+1
              wrt_int_state%WRITE_SUBSET_R(I,JLAT,NF)=wrt_int_state%ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
             ENDDO
             ENDDO
            ENDDO
          ENDIF
        endif
!-----------------------------------------------------------------------
!*** collect all_data_i/r2d to write pe(s)
!-----------------------------------------------------------------------
!
   from_fcst_tasks: DO N=1,NFCST_TASKS_SEND                               !<-- Loop through fcst tasks sending to this write task
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D INTEGER DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        IF(N==1)THEN
          NSTART_I2D=1
        ELSE
          NSTART_I2D=NSTART_I2D+NLAT_I2D
        ENDIF
        NLAT_I2D=wrt_int_state%nlat_from_fcst_task(N)*IM*KOUNT_I2D
        IF(NLAT_I2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_I2D(NSTART_I2D)          &  !<-- Fcst tasks' string of 2D integer history data
                       ,NLAT_I2D                                        &  !<-- Max #of words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,N-1                                             &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR*1000+N-1+(mype-nfcst_tasks_send+1)*100            &  !<-- An MPI tag
                       ,mpi_commun                                      &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
!
          NN=0
!
          JSTART=wrt_int_state%nstart_from_fcst_task(N)
          DO J=1,wrt_int_state%nlat_from_fcst_task(N)
!
            JFCST_I2D=JFCST_I2D+1
            JLAT=wrt_int_state%fcst_lat_on_write_task(JSTART+J-1)
            DO NF=1,KOUNT_I2D                                                 !<-- Loop through all the 2D integer fields
            DO I=1,IM
              NN=NN+1
              wrt_int_state%WRITE_SUBSET_I(I,JLAT,NF)=wrt_int_state%ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
            ENDDO
            ENDDO
           ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D REAL DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        NLAT_R2D=wrt_int_state%nlat_from_fcst_task(N)*IM*KOUNT_R2D
        IF(NLAT_R2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_R2D(1)                     &  !<-- Fcst tasks' string of 2D real history data
                       ,NLAT_R2D                                          &  !<-- Max #of words in the data string
                       ,MPI_REAL                                          &  !<-- The datatype
                       ,N-1                                               &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR*1000+N-1+(mype-nfcst_tasks_send+1)*100 & !<-- An MPI tag
                       ,mpi_commun                                        &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                             &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK NEEDS TO INSERT THE PIECES OF THE VARIOUS
!***  2D HISTORY ARRAYS RECEIVED FROM THE INDIVIDUAL FCST TASKS
!***  INTO ARRAYS THAT SPAN THE WRITE TASKS' OWN SUBSECTION OF
!***  THE FULL 2D DOMAIN.  THAT SUBSECTION ALWAYS SPANS THE 
!***  ENTIRE EAST-WEST DIMENSION OF THE FULL DOMAIN (SINCE FULL
!***  ROWS OF FCST TASKS ALWAYS SEND TO WRITE TASKS, NEVER 
!***  PARTIAL ROWS) AND AS MUCH OF THE NORTH-SOUTH DIMENSION OF
!***  THE FULL DOMAIN AS COVERED BY THOSE FCST TASKS SENDING TO
!***  A GIVEN WRITE TASK.
!-----------------------------------------------------------------------
!
          NN=0
!
          JSTART=wrt_int_state%nstart_from_fcst_task(N)
          DO J=1,wrt_int_state%nlat_from_fcst_task(N)
!
           JFCST_R2D=JFCST_R2D+1
           JLAT=wrt_int_state%fcst_lat_on_write_task(JSTART+J-1)
           DO NF=1,KOUNT_R2D                                                 !<-- Loop through all the 2D real fields
           DO I=1,IM
            NN=NN+1
            wrt_int_state%WRITE_SUBSET_R(I,JLAT,NF)=wrt_int_state%ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
           ENDDO
           ENDDO
          ENDDO
!
        ENDIF  !end NLAT_R2D>0
!
!-----------------------------------------------------------------------
!
       ENDDO from_fcst_tasks
!
!-----------------------------------------------------------------------
!
      ENDIF ifntasks3
!
!-----------------------------------------------------------------------
!***  AT THIS POINT, ALL WRITE TASKS HAVE RECEIVED ALL OF THE HISTORY
!***  DATA FROM THEIR ASSOCIATED FCST TASKS AND ASSEMBLED IT ONTO 
!***  THEIR OWN SUBSECTIONS OF THE FULL 2D DOMAIN.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  IT IS TIME FOR THE LEAD WRITE TASK TO BEGIN WRITING TO THE
!***  HISTORY FILES.  THE LEAD WRITE TASK ALREADY HOLDS ALL OF THE
!***  SCALAR/1D HISTORY DATA AND CAN GO AHEAD AND WRITE THEM.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_NEMSIOFLAG) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       MESSAGE_CHECK="Lead Write Task Gets Current ESMF Time from Clock"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK                             &  !<-- The ESMF Clock
                            ,currTime=CURRTIME                          &  !<-- The current time (ESMF) on the clock
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE CURRENT FORECAST TIME.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       MESSAGE_CHECK="Lead Write Task Gets Actual Current Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeGet(time=CURRTIME                               &  !<-- The cuurent forecast time (ESMF)
                           ,yy  =IYEAR_FCST                             &  !<-- The current forecast year (integer)
                           ,mm  =IMONTH_FCST                            &  !<-- The current forecast month (integer)
                           ,dd  =IDAY_FCST                              &  !<-- The current forecast day (integer)
                           ,h   =IHOUR_FCST                             &  !<-- The current forecast hour (integer)
                           ,m   =IMINUTE_FCST                           &  !<-- The current forecast minute (integer)
                           ,s   =ISECOND_FCST                           &  !<-- The current forecast second (integer)
                           ,sN  =ISECOND_NUM                            &  !<-- Numerator of current fractional second (integer)
                           ,sD  =ISECOND_DEN                            &  !<-- Denominator of current fractional second (integer)
                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          SECOND_FCST=ISECOND_FCST+REAL(ISECOND_NUM)/REAL(ISECOND_DEN)     !<-- Current forecast seconds (real)
!
!-----------------------------------------------------------------------
!***  ELAPSED FORECAST TIME.
!-----------------------------------------------------------------------
!
          wrt_int_state%IO_CURRTIMEDIFF=CURRTIME-wrt_int_state%IO_BASETIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Actual Elapsed Fcst Time"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                                   ,h           =NF_HOURS               &  !<-- Hours of elapsed time
                                   ,m           =NF_MINUTES             &  !<-- Minutes of elapsed time
                                   ,s           =NSECONDS               &  !<-- Seconds of elapsed time
                                   ,sN          =NSECONDS_NUM           &  !<-- Numerator of fractional elapsed seconds
                                   ,sD          =NSECONDS_DEN           &  !<-- denominator of fractional elapsed seconds
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NF_SECONDS=NSECONDS+REAL(NSECONDS_NUM)/REAL(NSECONDS_DEN)
          wrt_int_state%NFHOUR=NF_HOURS
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  WE WILL NOW ASSEMBLE THE FULL DOMAIN 2-D HISTORY DATA ONTO
!***  THE LEAD WRITE TASK FROM THE SUBSECTIONS ON ALL WRITE TASKS
!***  THEN THE LEAD TASK WILL WRITE EACH 2D FIELD TO THE HISTORY
!***  FILE.
!
!***  NOTE:  THE LEAD WRITE TASK ASSEMBLES AND WRITES TO HISTORY ONLY
!***         ONE 2-D FIELD AT A TIME.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WE WILL NOW ASSEMBLE THE FULL DOMAIN 2-D HISTORY DATA ONTO
!***  THE LEAD WRITE TASK FROM THE SUBSECTIONS ON ALL WRITE TASKS
!***  THEN THE LEAD TASK WILL WRITE EACH 2D FIELD TO THE HISTORY
!***  FILE. 
!
!***  NOTE:  THE LEAD WRITE TASK ASSEMBLES AND WRITES TO HISTORY ONLY
!***         ONE 2-D FIELD AT A TIME.  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST LOOP THROUGH ALL OF THE INTEGER Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WRITE OUT THE HISTORY FILE AS SELECTED.
!-----------------------------------------------------------------------
!
      IF(MYPE_LOCAL==0)THEN
!
        IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
          DEGRAD=90./ASIN(1.)
          CALL WRITE_NEMSIO_OPEN(WRT_INT_STATE                          &
                                      ,NBDL                             &
                                      ,NEMSIOFILE                       &
                                      ,IYEAR_FCST                       &
                                      ,IMONTH_FCST                      &
                                      ,IDAY_FCST                        &
                                      ,IHOUR_FCST                       &
                                      ,IMINUTE_FCST                     &
                                      ,SECOND_FCST                      &
                                      ,NF_HOURS                         &
                                      ,NF_MINUTES                       &
                                      ,NF_SECONDS                       &
                                      ,DIM1,DIM2,NBDR                   &
                                      ,LEAD_WRITE_TASK)
          FIELDSIZE=(DIM1+2*NBDR)*(DIM2+2*NBDR)
          ALLOCATE(TMP(FIELDSIZE))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
    field_loop_int: DO NFIELD=1,KOUNT_I2D                                !<-- Loop through all 2D integer gridded history data
!-----------------------------------------------------------------------
!
        IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%BUFF_INT(NN)=wrt_int_state%WRITE_SUBSET_I(I,J,NFIELD)
          ENDDO
          ENDDO
!
          allocate(idisp(nwtpg),ircnt(nwtpg))
          IDISP(1)=0
          Do I=1,NWTPG
            IRCNT(I)=IM*(wrt_int_state%JEND_WRITE(I)-wrt_int_state%JSTART_WRITE(I)+1)
            if(I>1) IDISP(I)=IDISP(I-1)+IRCNT(I-1)
          ENDDO
!
          CALL MPI_GATHERV(wrt_int_state%BUFF_INT                          &
                          ,NN                                              &
                          ,MPI_INTEGER                                     &
                          ,wrt_int_state%BUFF_INT_TMP                      &
                          ,IRCNT                                           &
                          ,IDISP                                           &
                          ,MPI_INTEGER                                     &
                          ,0                                               &
                          ,MPI_COMM_COMP                                   &
                          ,IERR)

            deallocate(idisp,ircnt)
!
            IF(MYPE_LOCAL==0)THEN                                  !<-- The lead write task
!
              NN=0 
              DO J=1,JM
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%BUFF_INT_TMP(NN) !<-- Lead write task fills its part of full domain
              ENDDO
              ENDDO
            ENDIF
!
         ELSE                                                    !<-- for 1pe
              NN=(NFIELD-1)*(IM*JM)
              DO J=1,JM
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%ALL_DATA_I2D(NN) !<-- Lead write task fills its part of full domain
              ENDDO
              ENDDO

!
          ENDIF
!
          IF(MYPE_LOCAL==0) THEN                                            !<-- WRite task write out data
!
            NPOSN_1=(NFIELD-1)*NAME_MAXSTR+1
            NPOSN_2=NFIELD*NAME_MAXSTR
            NAME=wrt_int_state%NAMES_I2D_STRING(NBDL)(NPOSN_1:NPOSN_2)      !<-- The name of this 2D integer history quantity
! 
!-----------------------------------------------------------------------
!***  FOR NEMSIO FILE
!-----------------------------------------------------------------------
!
           IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: input data dimension ',IM*JM,           &
               ' does not match data size in NEMSIO file ',FIELDSIZE
            ENDIF
!
            TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM),(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,NFIELD,TMP,IRET=IERR)           !<-- Lead write task writes out the 2D int data!
             if(ierr/=0) print *,'rec num=',NFIELD,' write ',trim(NAME), &
               'into file,ierr=',ierr,'NF_hours=',NF_HOURS,'nf_minutes=',&
               NF_MINUTES
!
          ENDIF
!-----------------------------------------------------------------------
!
         ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_int
!
!-----------------------------------------------------------------------
!***  NOW LOOP THROUGH ALL THE REAL Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_real: DO NFIELD=1,wrt_int_state%KOUNT_R2D(NBDL)              !<-- Loop through all 2D real gridded history data
!-----------------------------------------------------------------------
!
        IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%BUFF_REAL(NN)=wrt_int_state%WRITE_SUBSET_R(I,J,NFIELD)
          ENDDO
          ENDDO
!
          allocate(idisp(nwtpg),ircnt(nwtpg))
          IDISP(1)=0
          Do I=1,NWTPG
            IRCNT(I)=IM*(wrt_int_state%JEND_WRITE(I)-wrt_int_state%JSTART_WRITE(I)+1)
            if(I>1) IDISP(I)=IDISP(I-1)+IRCNT(I-1)
          ENDDO
!
          CALL MPI_GATHERV(wrt_int_state%BUFF_REAL                         &
                          ,NN                                              &
                          ,MPI_REAL                                        &
                          ,wrt_int_state%BUFF_REAL_TMP                     &
                          ,IRCNT                                           &
                          ,IDISP                                           &
                          ,MPI_REAL                                        &
                          ,0                                               &
                          ,MPI_COMM_COMP                                   &
                          ,IERR)

            deallocate(idisp,ircnt)
!
            IF(MYPE_LOCAL==0)THEN                                  !<-- The lead write task
!
              NN=0
              DO J=1,JM
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%BUFF_REAL_TMP(NN) !<-- Lead write task fills its part of full domain
              ENDDO
              ENDDO
            ENDIF
!
        ELSE                                                    !<-- for 1pe
          IF(NTASKS>1) then
            DO J=1,JM
            DO I=1,IM
              wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=                         &
                wrt_int_state%WRITE_SUBSET_R(I,J,NFIELD)                        !<-- Lead write task fills its part of full domain
            ENDDO
            ENDDO
          ELSE
            NN=(NFIELD-1)*(IM*JM)
            DO J=1,JM
            DO I=1,IM
              NN=NN+1
              wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%ALL_DATA_R2D(NN) !<-- Lead write task fills its part of full domain
            ENDDO
            ENDDO
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(MYPE_LOCAL==0)THEN                                  !<-- The lead write task
!
          NPOSN_1=(NFIELD-1)*NAME_MAXSTR+1
          NPOSN_2=NFIELD*NAME_MAXSTR
          NAME=wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D real history quantity
!
!-----------------------------------------------------------------------
!***  FOR NEMSIO FILE
!-----------------------------------------------------------------------
!
          IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: data dimension ',IM*JM,                &
               ' does not match data size in NEMSIO file,',FIELDSIZE
            ENDIF
!
!check time average
           itr=-99
           if(INDEX(NAME,"_ave") >0) then
               itr=3
            elseif(INDEX(NAME,"_acc") >0) then
               itr=4
            elseif(INDEX(NAME,"_win") >0) then
               itr=2
            endif
!
            N=NFIELD+wrt_int_state%KOUNT_I2D(1)
            TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
!
!             write(0,*)'after write N=',N,' var=',trim(NAME),'value=',&
!              maxval(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)), &
!              minval(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)), 'itr=',itr
!
            if(itr==-99) then
              CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            else
              CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR,itr=itr,      &
                zhour=wrt_int_state%zhour)
            endif
             if(ierr/=0) print *,'rec num=',N,' write ',trim(NAME), &
               'into file,ierr=',ierr,'NF_hours=',NF_HOURS,'nf_minutes=',&
               NF_MINUTES
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_real
!
      if(NTASKS>1) THEN
        DEALLOCATE(wrt_int_state%WRITE_SUBSET_I)
        DEALLOCATE(wrt_int_state%WRITE_SUBSET_R)
      endif
!------------
!***  NEMSIO
!------------
!
      IF(wrt_int_state%WRITE_NEMSIOFLAG.AND.MYPE_LOCAL==0)THEN
!
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,IERR,gfname=GFNAME)
!
        DEALLOCATE(TMP)
!
        CALL NEMSIO_CLOSE(NEMSIOFILE)
!
        CALL NEMSIO_FINALIZE()
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO file_loop   !nfile
!
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,ierr)
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
        WRITE(0,*)"PASS: WRITE_RUN"
      ELSE
        WRITE(0,*)"FAIL: WRITE_RUN"
      ENDIF
!
!      write_run_tim=write_run_tim+timef()-btim0
!
      IF(MYPE_LOCAL==0)THEN
        WRITE(0,*)' Write Time is ',write_run_tim*1.e-3 &
                 ,' at Fcst ',NF_HOURS,':',NF_MINUTES,':',NF_SECONDS
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRT_RUN_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_FINALIZE_GFS(WRT_COMP                              &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK                                     &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!***  HISTORY
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT):: WRT_COMP
      TYPE(ESMF_State)  ,INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE_WRITE  
      TYPE(ESMF_Clock)  ,INTENT(INOUT) :: CLOCK
!
      INTEGER,INTENT(OUT)              :: RCFINAL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RCFINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      IF(RCFINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'PASS: Write_Finalize.'
      ELSE
        WRITE(0,*)'FAIL: Write_Finalize.'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRT_FINALIZE_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_SETUP_GFS(ATM_GRID_COMP,WRT_COMPS           &
        ,exp_state_dyn,exp_state_phy                               &
        ,imp_state_write,exp_state_write)
! 
      use module_gfs_mpi_def, only : last_fcst_pe
!-----------------------------------------------------------------------
!***  SET UP THE WRITE COMPONENTS WITH THE FORECAST TASKS AND
!***  THE GROUPS OF WRITE TASKS NEEDED FOR QUILTING THE OUTPUT
!***  AND WRITING IT TO HISTORY FILES.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),INTENT(inOUT)      :: WRT_COMPS(:)              !<-- The ATM gridded component
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_DYN
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_PHY
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_state_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_WRITE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
      TYPE(ESMF_VM)          :: VM                                        !<-- The ESMF virtual machine.
!
      INTEGER                :: INPES,JNPES                             & !<-- Number of fcst tasks in I and J directions
                               ,MYPE                                    & !<-- PE ID
                               ,WRITE_GROUPS                            & !<-- Number of groups of write tasks
                               ,WRITE_TASKS_PER_GROUP                     !<-- #of tasks in each write group
      INTEGER,ALLOCATABLE    :: ITMP(:)
      LOGICAL                :: STANDALONE_POST
!
      CHARACTER( 2)          :: MY_WRITE_GROUP
      CHARACTER(6)           :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR) :: WRITE_NAME
!
      INTEGER :: I,J,K,RC,RC_SETUP
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE TASK AND GROUP COUNTS FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Tasks and Groups from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,QUILTING                             &  !<-- Number of write groups from config file
                                  ,label ='quilting:'                   &
                                  ,rc    =RC)

      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_GROUPS                         &  !<-- Number of write groups from config file
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_TASKS_PER_GROUP                &  !<-- Number of write tasks per group from config file
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  HOW MANY FORECAST TASKS DO WE HAVE?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="INPES/JNPES from Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &  !<-- Give the value of this label to INPES
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &  !<-- Give the value of this label to JNPES
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Total number of forecast tasks
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CURRENT VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  WHAT IS MY MPI TASK ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ASSOCIATE ALL OF THE FORECAST TASKS WITH THE WRITE TASKS
!***  IN EACH WRITE GROUP.
!-----------------------------------------------------------------------
!
      IF(QUILTING) THEN
        NUM_PES_WRT=NUM_PES_FCST+WRITE_TASKS_PER_GROUP
      ELSE
        NUM_PES_WRT=NUM_PES_FCST
      ENDIF
      ALLOCATE(PETLIST_WRITE(NUM_PES_WRT,WRITE_GROUPS))                        !<-- Task IDs of all wrt tasks
                                                                               !    plus the write tasks
                                                                               !    by write group
!
!-----------------------------------------------------------------------
!***  COLLECT THE TASK IDs FOR THE WRITE TASKS AND THE ASSOCIATED
!***  FORECAST-WRITE TASKS.
!-----------------------------------------------------------------------
!
      DO I=0,NUM_PES_FCST-1
        DO J=1,WRITE_GROUPS
          PETLIST_WRITE(I+1,J)=I+last_fcst_pe+1-num_pes_fcst            !<-- Collect forecast task IDs to be associated with
                                                                        !    write tasks by write group
        ENDDO
!
      ENDDO
!
      IF(NUM_PES_WRT>NUM_PES_FCST) THEN 
        K=NUM_PES_FCST
!
        DO J=1,WRITE_GROUPS
          DO I=1,WRITE_TASKS_PER_GROUP
            PETLIST_WRITE(NUM_PES_FCST+I,J)=K                           !<-- Append write task IDs to associated forecast task IDs by group
            K=K+1
          ENDDO
        ENDDO
      ENDIF
!      write(0,*)'in write setup, last_fcst_pe=',last_fcst_pe,'num_pes_fcst=', &
!        num_pes_fcst,'WRITE_GROUPS=',WRITE_GROUPS,'petlist_write=',  &
!         PETLIST_WRITE(:,1)
!
!-----------------------------------------------------------------------
!***  CREATE THE WRITE GRIDDED COMPONENT(S).
!***  THERE ARE AS MANY WRITE COMPONENTS AS THERE ARE GROUPS OF
!***  WRITE TASKS SPECIFIED IN THE CONFIGURE FILE.
!***  REGISTER THEIR INIT, RUN, AND FINALIZE STEPS.
!-----------------------------------------------------------------------
!
!---------------------------------
!***  Create the Write components
!---------------------------------
!
      allocate(ITMP(NUM_PES_WRT))
      DO I=1,WRITE_GROUPS
        ITMP(1:NUM_PES_WRT)=PETLIST_WRITE(1:NUM_PES_WRT,I)
        WRITE(MY_WRITE_GROUP,FMT)I
        WRITE_NAME='write_GridComp_'//MY_WRITE_GROUP
!
        WRT_COMPS(I)=ESMF_GridCompCreate(                         &
                                name          =WRITE_NAME         &  !<-- Name of this group's Write gridded component
                               ,configFile    ='atm_namelist.rc'  &  !<-- The configure file for writes
                               ,petList       =ITMP               &  !<-- The task IDs of the write tasks in this group
                                                                     !    provide the local VM information per component.
                               ,rc            =RC)
!
      ENDDO
      deallocate(itmp)
!
!-----------------------------------
!***  Register the Write components
!-----------------------------------
!
      DO I=1,WRITE_GROUPS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Write Components"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(WRT_COMPS(I)         &  !<-- The Write gridded components
                                     ,WRITE_REGISTER_GFS   &  !<-- The user's subroutine name
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Write subcomponent(s)
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_STATE_WRITE=ESMF_StateCreate(statename='Write Import State' &  !<-- Import state name for writes
                                      ,statetype= ESMF_STATE_IMPORT   &
                                      ,rc       = RC)
!
      EXP_STATE_WRITE=ESMF_StateCreate(statename='Write Export State' &  !<-- Export state names for writes
                                      ,statetype= ESMF_STATE_EXPORT   &
                                      ,rc       = RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE WRITE COMPONENTS' IMPORT STATE INTO THE
!***  DYNAMICS' AND PHYSICS' EXPORT STATES SINCE HISTORY
!***  DATA ITSELF MUST COME FROM THE DYNAMICS AND PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Write Import State into Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =EXP_STATE_DYN        & !<-- Dynamics export state receives a state
                        ,nestedState=IMP_STATE_WRITE      & !<-- Add the write components' import state
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Write Import State into Physics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =EXP_STATE_PHY        & !<-- Physics export state receives a state
                        ,nestedState=IMP_STATE_WRITE      & !<-- Add the write components' import state
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_SETUP_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_DESTROY_GFS(ATM_GRID_COMP,WRT_COMPS,             &
        IMP_STATE_WRITE,EXP_STATE_WRITE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  DESTROY ALL OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      USE module_gfs_mpi_def, only: WRITE_GROUPS, WRITE_TASKS_PER_GROUP
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),DIMENSION(:),INTENT(INOUT)      ::WRT_COMPS
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_STATE_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_STATE_WRITE
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,MYPE,N,RC,RC_DES
!
      TYPE(ESMF_VM)                          :: VM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_DES=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CURRENT VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM in Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  WHAT IS MY MPI TASK ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs for Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FINALIZE THE WRITE GRIDDED COMPONENTS IN EACH WRITE GROUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,WRITE_GROUPS              
        IF(MYPE>=PETLIST_WRITE(1,N).AND.                          &
           MYPE<=PETLIST_WRITE(NUM_PES_WRT,N))THEN
!
           CALL ESMF_GridCompFinalize(gridcomp   =WRT_COMPS(N)    &
                                     ,importstate=IMP_STATE_WRITE &
                                     ,exportstate=EXP_STATE_WRITE &
                                     ,clock      =CLOCK_ATM       &
                                     ,phase      =ESMF_SINGLEPHASE &
                                     ,rc         =RC)
        ENDIF
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY THE WRITE COMPONENTS' IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Component Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(IMP_STATE_WRITE,rc=RC)
      CALL ESMF_StateDestroy(EXP_STATE_WRITE,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO J=1,WRITE_GROUPS
!
        CALL ESMF_VMBarrier(vm=VM,rc=RC)
!
        DO I=1,NUM_PES_WRT
          IF(MYPE==PETLIST_WRITE(I,J))THEN
            CALL ESMF_GridCompDestroy(gridcomp=WRT_COMPS(J) &
                                     ,rc      =RC)
          ENDIF
        ENDDO
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_DES==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM FINALIZE STEP FAILED  RC_DES=',RC_DES
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_GRID_COMP_GFS
!
!-----------------------------------------------------------------------
