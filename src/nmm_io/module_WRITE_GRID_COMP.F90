!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_GRID_COMP
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
!       05 Jan 2009:  J. Wang  - Add 10-m wind factor into NMMB
!                                runhistory and restart files
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_WRITE_INTERNAL_STATE
      USE MODULE_WRITE_ROUTINES,ONLY : FIRST_PASS_HST                   &
                                      ,FIRST_PASS_RST                   &
                                      ,OPEN_HST_FILE                    &
                                      ,OPEN_RST_FILE                    &
                                      ,WRITE_RUNHISTORY_OPEN            &
                                      ,WRITE_NEMSIO_RUNHISTORY_OPEN     &
                                      ,WRITE_RUNRESTART_OPEN            &
                                      ,WRITE_NEMSIO_RUNRESTART_OPEN     &
                                      ,TIME_FOR_HISTORY                 &
                                      ,TIME_FOR_RESTART
!
      USE MODULE_DM_PARALLEL   ,ONLY : PARA_RANGE                       &
                                      ,MPI_COMM_COMP                    &
                                      ,MPI_COMM_INTER_ARRAY
      USE MODULE_CONTROL       ,ONLY : TIMEF
      USE MODULE_GET_CONFIG_WRITE
      USE MODULE_ERR_MSG       ,ONLY : ERR_MSG,MESSAGE_CHECK
      USE MODULE_INCLUDE
      USE MODULE_CONSTANTS,ONLY : G
      USE NEMSIO_MODULE
      USE MODULE_BGRID_INTERP, ONLY: V_TO_H_BGRID
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_REGISTER
!
      PUBLIC :: WRITE_SETUP,WRITE_DESTROY
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_LENGTH_I1D=5000                            !<-- Max words in all 1-D integer history variables
      INTEGER,PARAMETER :: MAX_LENGTH_R1D=25000                           !<-- Max words in all 1-D real history variables
      INTEGER,PARAMETER :: MAX_LENGTH_LOG=MAX_DATA_LOG                    !<-- Max logical variables
!
      INTEGER,SAVE      :: LAST_FCST_TASK                                 !<-- Rank of the last forecast task
      INTEGER,SAVE      :: LEAD_WRITE_TASK                                !<-- Rank of the lead (first) write task in this write group
      INTEGER,SAVE      :: LAST_WRITE_TASK                                !<-- Rank of the last write task the write group
      INTEGER,SAVE      :: NTASKS                                         !<-- # of write tasks in the current group + all forecast tasks
      INTEGER,SAVE      :: NWTPG                                          !<-- # of write tasks (servers) per group 
!
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: NCURRENT_GROUP             !<-- The currently active write group
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE                       !<-- The history output data Bundle
      TYPE(ESMF_FieldBundle),SAVE :: RESTART_BUNDLE                       !<-- The restart output data Bundle
!
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE                 ! The internal state pointer.
!
!-----------------------------------------------------------------------
      REAL(KIND=KFPT)             :: btim,btim0
      REAL(KIND=KFPT),PUBLIC,SAVE :: write_init_tim                     &
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
      SUBROUTINE WRITE_REGISTER(WRT_COMP,RC_WRT)
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
                                     ,WRT_INITIALIZE                    &  !<-- User's subroutineName
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
                                     ,WRT_RUN                           &  !<-- User's subroutineName
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
                                    ,WRT_FINALIZE                       &  !<-- User's subroutineName
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
      END SUBROUTINE WRITE_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_INITIALIZE(WRT_COMP                                &
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
!
      TYPE(ESMF_VM)                          :: VM
      TYPE(WRITE_WRAP)                       :: WRAP
      TYPE(WRITE_INTERNAL_STATE),POINTER     :: WRT_INT_STATE
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
        NCURRENT_GROUP(1)=0
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
!***  ALL TASKS ALLOCATE BUFFER DATA ARRAYS THAT WILL HOLD SCALAR/1D
!***  HISTORY/RESTART DATA AND WILL BE USED TO Send/Recv THAT DATA
!***  BETWEEN THE FORECAST TASKS THAT KNOW IT INITIALLY TO THE
!***  WRITE TASKS THAT OBTAIN IT FROM THE FORECAST TASKS FOR WRITING.
!***  LOGICAL DATA BUFFERS ARE ALSO HANDLED HERE.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_I1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_I1D(MAX_LENGTH_I1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_R1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_R1D(MAX_LENGTH_R1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%ALL_DATA_LOG))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_LOG(MAX_LENGTH_LOG),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%RST_ALL_DATA_I1D))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_I1D(MAX_LENGTH_I1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%RST_ALL_DATA_R1D))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_R1D(MAX_LENGTH_R1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ALLOCATED(wrt_int_state%RST_ALL_DATA_LOG))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_LOG(MAX_LENGTH_LOG),stat=ISTAT)
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
        ALLOCATE(wrt_int_state%NCOUNT_FIELDS(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_I1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_I2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_R1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_R2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_LOG(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_NCOUNT_FIELDS(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_I1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%RST_KOUNT_I2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_R1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%RST_KOUNT_R2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_LOG(1),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  ALL INTEGER QUANTITIES (AS 1D ARRAYS) AND 1D AND 2D REAL
!***  QUANTITIES WILL BE STRUNG TOGETHER IN SINGLE ARRAYS OF
!***  EACH PARTICULAR TYPE.  WE NEED TO ALLOCATE THE ARRAYS THAT WILL
!***  HOLD THE LENGTH OF EACH OF THE QUANTITIES IN THESE 'STRINGS'
!***  AS THE 'STRINGS' THEMSELVES.
!-----------------------------------------------------------------------
!
        ALLOCATE(wrt_int_state%LENGTH_DATA_I1D(100),stat=ISTAT)            !<-- Lengths of each individual 1-D integer array
        ALLOCATE(wrt_int_state%LENGTH_DATA_R1D(100),stat=ISTAT)            !<-- Lengths of each individual 1-D real array
        ALLOCATE(wrt_int_state%LENGTH_SUM_I1D(1),stat=ISTAT)               !<-- Length of string of data of ALL 1-D integer arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_R1D(1),stat=ISTAT)               !<-- Length of string of data of ALL 1-D real arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_LOG(1),stat=ISTAT)               !<-- Length of string of data of ALL logical variables
!
        ALLOCATE(wrt_int_state%RST_LENGTH_DATA_I1D(100),stat=ISTAT)        !<-- Lengths of each restart individual 1-D integer array
        ALLOCATE(wrt_int_state%RST_LENGTH_DATA_R1D(100),stat=ISTAT)        !<-- Lengths of each restart individual 1-D real array
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_I1D(1),stat=ISTAT)           !<-- Length of string of restart data of ALL 1-D integer arrays
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_R1D(1),stat=ISTAT)           !<-- Length of string of restart data of ALL 1-D real arrays
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_LOG(1),stat=ISTAT)           !<-- Length of string of restart data of ALL logical variables
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  RETRIEVE INFORMATION REGARDING OUTPUT FROM THE CONFIGURATION FILE.
!***  THIS INFORMATION CONTAINS MAINLY THE OUTPUT CHANNELS AND NAMES OF
!***  THE DISK FILES.
!-----------------------------------------------------------------------
!
      CALL GET_CONFIG_WRITE(WRT_COMP,WRT_INT_STATE,RC)                   !<-- User's routine to extract configfile data
!
      NWTPG          =wrt_int_state%WRITE_TASKS_PER_GROUP
      LAST_FCST_TASK =NTASKS-NWTPG-1
      LEAD_WRITE_TASK=LAST_FCST_TASK+1
      LAST_WRITE_TASK=NTASKS-1
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS THAT HOLD THE LOCAL LIMITS 
!***  OF ALL THE FORECAST TASKS' SUBDOMAINS.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%LOCAL_ISTART))THEN
        ALLOCATE(wrt_int_state%LOCAL_ISTART(0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local starting I for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_IEND  (0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local ending I for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_JSTART(0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local starting J for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_JEND  (0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local ending J for each fcst task's subdomain
      ENDIF
!
!-----------------------------------------------------------------------
!***  EXTRACT THE HISTORY DATA Bundle FROM THE IMPORT STATE
!***  SO THAT IT IS AVAILABLE FOR RETRIEVING DATA FROM IT
!***  DURING THE RUN STEP.
!***  THE Bundle WAS CREATED DURING THE INIT STEP OF THE Dynamics
!***  SINCE SUBROUTINE POINT_DYNAMICS_OUTPUT MUST HAVE IT AVAILABLE
!***  FOR INSERTING DATA POINTERS INTO IT.  
!***  ONLY THE FORECAST TASKS CAN EXTRACT IT PROPERLY SINCE IT WAS
!***  THEY WHO INSERTED IT.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%MYPE<=LAST_FCST_TASK)THEN                           !<-- The forecast tasks
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract History Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   ='Bundle_Output_Data'             &  !<-- The name of the history data Bundle
                          ,fieldbundle=HISTORY_BUNDLE                   &  !<-- The history data Bundle inside the import state
                          ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Restart Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   ='Bundle_Restart_Data'            &  !<-- The name of the restart data Bundle
                          ,fieldbundle=RESTART_BUNDLE                   &  !<-- The restart data Bundle inside the import state
                          ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
      write_init_tim=timef()-btim0
!
!----------------------------------------------------------------------- 
!
      END SUBROUTINE WRT_INITIALIZE
!
!----------------------------------------------------------------------- 
!####################################################################### 
!----------------------------------------------------------------------- 
!
      SUBROUTINE WRT_RUN(WRT_COMP                                       &
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
      INTEGER,SAVE :: IH_INT =MPI_REQUEST_NULL                          &
                     ,IH_REAL=MPI_REQUEST_NULL
!
      INTEGER,SAVE :: RST_IH_INT =MPI_REQUEST_NULL                      &
                     ,RST_IH_REAL=MPI_REQUEST_NULL
!
      INTEGER,SAVE :: ITS,ITE,JTS,JTE
      INTEGER,SAVE :: IMS,IME,JMS,JME
      INTEGER,SAVE :: IHALO,JHALO
!
      INTEGER,SAVE :: N_START,N_END
      INTEGER,SAVE :: NPOSN_1,NPOSN_2
!
      INTEGER      :: I,I1,IJ,IJG,IM,IMJM                               &
                     ,J,JM,K,L                                          &
                     ,MY_LOCAL_ID                                       &
                     ,MYPE,MYPE_ROW                                     &
                     ,N,N1,N2,NF,NN                                     &
                     ,NN_INTEGER,NN_REAL                                &
                     ,N_POSITION                                        &
                     ,NUM_ATTRIB
!
      INTEGER      :: DIM1,DIM2,FIELDSIZE,NBDR
!
      INTEGER      :: IYEAR_FCST                                        &
                     ,IMONTH_FCST                                       &
                     ,IDAY_FCST                                         &
                     ,IHOUR_FCST                                        &
                     ,IMINUTE_FCST                                      &
                     ,ISECOND_FCST                                      &
                     ,ISECOND_NUM                                       &
                     ,ISECOND_DEN
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
      INTEGER(KIND=KINT)         :: NTIMESTEP
!
      INTEGER      :: NF_HOURS                                          &
                     ,NF_MINUTES                                        &
                     ,NSECONDS                                          &
                     ,NSECONDS_NUM                                      &
                     ,NSECONDS_DEN
!
      INTEGER      :: ID_DUMMY                                          &
                     ,ID_RECV                                           &
                     ,NFCST_TASKS                                       &
                     ,NFIELD                                            &
                     ,NPE_WRITE                                         &
                     ,NUM_FIELD_NAMES                                   &
                     ,NUM_PES_FCST
!
      INTEGER     :: ID_START,ID_END                                    &
                    ,ISTART,IEND                                        &
                    ,JSTART,JEND
!
      INTEGER     :: JROW_FIRST,JROW_LAST                               &
                    ,JSTA_WRITE,JEND_WRITE
!
      INTEGER     :: KOUNT_I2D                                          &
                    ,KOUNT_I2D_DATA                                     &
                    ,KOUNT_R2D                                          &
                    ,KOUNT_R2D_DATA
!
      INTEGER     :: RST_KOUNT_I2D                                      &
                    ,RST_KOUNT_I2D_DATA                                 &
                    ,RST_KOUNT_R2D                                      &
                    ,RST_KOUNT_R2D_DATA
!
      INTEGER     :: LENGTH                                             &
                    ,MPI_COMM,MPI_COMM2
!
      INTEGER     :: IO_HST_UNIT,IO_RST_UNIT
!
      INTEGER     :: IERR,ISTAT,RC
      INTEGER,DIMENSION(MPI_STATUS_SIZE)    :: JSTAT
!
      INTEGER,DIMENSION(:)  ,POINTER        :: WORK_ARRAY_I1D
      INTEGER,DIMENSION(:,:),POINTER        :: WORK_ARRAY_I2D
!
      REAL :: DEGRAD                                                    &
             ,NF_SECONDS                                                &
             ,SECOND_FCST
!
      REAL(KIND=KFPT),DIMENSION(:)  ,POINTER  :: WORK_ARRAY_R1D
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER  :: WORK_ARRAY_R2D
      REAL(KIND=KFPT),DIMENSION(:)  ,POINTER  :: GLAT1D,GLON1D,TMP
!
      REAL(KIND=KFPT),DIMENSION(:,:),ALLOCATABLE:: FACT10               &
                                                  ,FACT10TMPU           &
                                                  ,FACT10TMPV           &
                                                  ,HGT
!
      LOGICAL                               :: GLOBAL                   &
                                              ,WRITE_LOGICAL            &
                                              ,OPENED
!
      LOGICAL,SAVE :: FIRST=.TRUE.                                      &
                     ,HST_FIRST=.TRUE.                                  &
                     ,RST_FIRST=.TRUE.
!
      CHARACTER(ESMF_MAXSTR) :: GFNAME,NAME,FILENAME
      CHARACTER(2)           :: MODEL_LEVEL
!
      TYPE(WRITE_WRAP)                   :: WRAP
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE
!
      TYPE(ESMF_Logical)   :: WORK_LOGICAL
!
      TYPE(ESMF_VM)        :: VM
      TYPE(ESMF_Grid),SAVE :: GRID1
      TYPE(ESMF_DELayout)  :: MY_DE_LAYOUT
      TYPE(ESMF_Field)     :: FIELD_WORK1
!
      TYPE(ESMF_InternArray) :: ARRAY_WORK
!
      TYPE(ESMF_Time)     :: CURRTIME
!
      TYPE(ESMF_TypeKind) :: DATATYPE
!
      TYPE(NEMSIO_GFILE)  :: NEMSIOFILE
!
      TYPE(ESMF_LOGICAL),DIMENSION(:),POINTER :: FIRST_IO_PE
!
!-----------------------------------------------------------------------
!
      real(kind=kfpt) :: wait_time
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
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
!
!-----------------------------------------------------------------------
!***  INCREMENT THE WRITE GROUP SO WE KNOW WHICH ONE IS ACTIVE.
!***  ONLY THE FORECAST TASKS ENTER THIS RUN STEP OF THE WRITE
!***  COMPONENT EVERY OUTPUT TIME THEREFORE ONLY THEY CAN PROPERLY
!***  INCREMENT THE NUMBER OF THE CURRENT WRITE GROUP. 
!***  LET FORECAST TASK 0 BROADCAST THE CURRENT GROUP NUMBER
!***  TO ALL TASKS IN THE GROUP INCLUDING THE CURRENT WRITE TASKS.
!-----------------------------------------------------------------------
!
      NCURRENT_GROUP(1)=NCURRENT_GROUP(1)+1
      IF(NCURRENT_GROUP(1)>wrt_int_state%WRITE_GROUPS)THEN
        NCURRENT_GROUP(1)=1
      ENDIF
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
      ALLOCATE(FIRST_IO_PE(1))                                            !<-- A flag indicating that this is or is not
                                                                          !    this write group's first pass through
                                                                          !    this write component
!
      IF(FIRST)THEN
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
                           ,LEAD_WRITE_TASK                             &  !<-- Root sender is the first write task in this group
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
        CALL FIRST_PASS_HST(IMP_STATE_WRITE,HISTORY_BUNDLE              &
                       ,WRT_INT_STATE,NTASKS,MYPE                       &
                       ,NCURRENT_GROUP(1))
!
        CALL FIRST_PASS_RST(IMP_STATE_WRITE,RESTART_BUNDLE              &
                       ,WRT_INT_STATE,NTASKS,MYPE                       &
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
        FIRST=.FALSE.
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
      KOUNT_I2D=wrt_int_state%KOUNT_I2D(1)
      KOUNT_R2D=wrt_int_state%KOUNT_R2D(1)
!
!-----------------------------------------------------------------------
!***  NOW PULL THE 2D DATA FROM THE IMPORT STATE.
!***  THIS INCLUDES ALL INDIVIDUAL 2D HISTORY QUANTITIES AS WELL AS
!***  ALL MODEL LEVELS OF THE 3D REAL HISTORY ARRAYS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      hst_fcst_tasks: IF(TIME_FOR_HISTORY.AND.MYPE<=LAST_FCST_TASK)THEN    !<-- Only the forecast tasks can see this data so far
!-----------------------------------------------------------------------
!
        NN_INTEGER=0
        NN_REAL   =0
!
        ITS=wrt_int_state%LOCAL_ISTART(MYPE)                               !<-- Starting I of this task's integration region
        ITE=wrt_int_state%LOCAL_IEND(MYPE)                                 !<-- Ending I of this task's integration region
        JTS=wrt_int_state%LOCAL_JSTART(MYPE)                               !<-- Starting J of this task's integration region
        JTE=wrt_int_state%LOCAL_JEND(MYPE)                                 !<-- Ending J of this task's integration region
!
        IHALO=wrt_int_state%IHALO                                          !<-- Halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Halo depth in J
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
        field_block: DO N=1,wrt_int_state%NCOUNT_FIELDS(1)                 !<-- Loop through all Fields in the import state
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Fields from History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(bundle=HISTORY_BUNDLE                &  !<-- The write component's history data Bundle
                                  ,name  =wrt_int_state%FIELD_NAME(N)   &  !<-- The ESMF Field's name
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
!
            IF(NN_INTEGER+IEND*JEND>wrt_int_state%NUM_WORDS_SEND_I2D_HST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF INTEGER WORDS YOU'    &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_INTEGER=NN_INTEGER+1
              wrt_int_state%ALL_DATA_I2D(NN_INTEGER)=WORK_ARRAY_I2D(I,J)   !<-- String together this task's 2D integer data
            ENDDO
            ENDDO
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
!
            IF(NN_REAL+IEND*JEND>wrt_int_state%NUM_WORDS_SEND_R2D_HST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF REAL WORDS YOU'       &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_REAL=NN_REAL+1
              wrt_int_state%ALL_DATA_R2D(NN_REAL)=WORK_ARRAY_R2D(I,J)      !<-- String together this task's 2D real data
            ENDDO
            ENDDO
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
!***  ALL FORECAST TASKS NOW SEND THEIR STRINGS OF 2D HISTORY DATA
!***  TO THE APPROPRIATE WRITE TASKS.
!-----------------------------------------------------------------------
!
        KOUNT_I2D_DATA=wrt_int_state%NUM_WORDS_SEND_I2D_HST                !<-- Total #of words in 2D integer history data on this fcst task
        KOUNT_R2D_DATA=wrt_int_state%NUM_WORDS_SEND_R2D_HST                !<-- Total #of words in 2D real history data on this fcst task
!
        MYPE_ROW=MYPE/wrt_int_state%INPES+1                                !<-- Each fcst task's row among all rows of fcst tasks
!
        DO N=1,NWTPG                                                       !<-- Loop through the write tasks in this group
          CALL PARA_RANGE(wrt_int_state%JNPES,NWTPG,N                   &  !<-- Find each write task's first and last rows of
                         ,JROW_FIRST,JROW_LAST)                            !<--   fcst tasks from which it will recv
!
          NPE_WRITE=N-1                                                    !<-- Consider the write task with this local ID
                                                                           !    beginning with 0
!
          IF(MYPE_ROW>=JROW_FIRST.AND.MYPE_ROW<=JROW_LAST)THEN             !<-- This fcst task associated with this write task
!
!-----------------------------------------------------------------------
!***  FIRST THE 2-D INTEGER DATA.
!-----------------------------------------------------------------------
!
            IF(KOUNT_I2D>0)THEN
              CALL MPI_ISEND(wrt_int_state%ALL_DATA_I2D                 &  !<-- Fcst tasks' string of 2D integer history data
                            ,KOUNT_I2D_DATA                             &  !<-- #of words in the data string
                            ,MPI_INTEGER                                &  !<-- The datatype
                            ,NPE_WRITE                                  &  !<-- The target write task
                            ,wrt_int_state%NFHOUR                       &  !<-- An MPI tag
                            ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))    &  !<-- The MPI intercommunicator between fcst and quilt tasks
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
            IF(KOUNT_R2D>0)THEN
              CALL MPI_ISEND(wrt_int_state%ALL_DATA_R2D                   &  !<-- Fcst tasks' string of 2D real history data
                            ,KOUNT_R2D_DATA                               &  !<-- #of words in the data string
                            ,MPI_REAL                                     &  !<-- The datatype
                            ,NPE_WRITE                                    &  !<-- The target write task
                            ,wrt_int_state%NFHOUR                         &  !<-- An MPI tag
                            ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))      &  !<-- The MPI intercommunicator between fcst and quilt tasks
                            ,IH_REAL                                      &  !<-- MPI communication request handle
                            ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of real data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
          ENDIF
!
        ENDDO
        write_send_data_tim=write_send_data_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF hst_fcst_tasks
!
!-----------------------------------------------------------------------
!***  NOW PULL THE 2D RESTART DATA FROM THE IMPORT STATE.
!***  THIS INCLUDES ALL INDIVIDUAL 2D RESTART QUANTITIES AS WELL AS
!***  ALL MODEL LEVELS OF THE 3D REAL RESTART ARRAYS.
!-----------------------------------------------------------------------
!
      RST_KOUNT_I2D=wrt_int_state%RST_KOUNT_I2D(1)
      RST_KOUNT_R2D=wrt_int_state%RST_KOUNT_R2D(1)
!
!-----------------------------------------------------------------------
      rst_fcst_tasks: IF(TIME_FOR_RESTART.AND.MYPE<=LAST_FCST_TASK)THEN    !<-- Only the forecast tasks can see this data so far
!-----------------------------------------------------------------------
!
        NN_INTEGER=0
        NN_REAL   =0
!
        ITS=wrt_int_state%LOCAL_ISTART(MYPE)                               !<-- Starting I of this task's integration region
        ITE=wrt_int_state%LOCAL_IEND(MYPE)                                 !<-- Ending I of this task's integration region
        JTS=wrt_int_state%LOCAL_JSTART(MYPE)                               !<-- Starting J of this task's integration region
        JTE=wrt_int_state%LOCAL_JEND(MYPE)                                 !<-- Ending J of this task's integration region
!
        IHALO=wrt_int_state%IHALO                                          !<-- Halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Halo depth in J
!
!-----------------------------------------------------------------------
!***  BE SURE THE INTEGER AND REAL BUFFERS ARE AVAILABLE FOR ISENDs
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL MPI_WAIT(RST_IH_INT,JSTAT,IERR) 
        wait_time=timef()-btim
        if(wait_time>1.e3)write(0,*)' Long integer buffer WAIT =',wait_time*1.e-3
!
        btim=timef()
        CALL MPI_WAIT(RST_IH_REAL,JSTAT,IERR) 
        wait_time=timef()-btim
        if(wait_time>1.e3)write(0,*)' Long real buffer WAIT =',wait_time*1.e-3
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        rst_field_block: DO N=1,wrt_int_state%RST_NCOUNT_FIELDS(1)         !<-- Loop through all Fields in the import state
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Fields from Restart Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(bundle=RESTART_BUNDLE                  &  !<-- The write component's restart data Bundle
                                  ,name  =wrt_int_state%RST_FIELD_NAME(N) &  !<-- The ESMF Field's name
                                  ,field =FIELD_WORK1                     &  !<-- The ESMF Field data pointer
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
          MESSAGE_CHECK="Check Datatype of Field from Restart Bundle"
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
!
            IF(NN_INTEGER+IEND*JEND>wrt_int_state%NUM_WORDS_SEND_I2D_RST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF INTEGER WORDS YOU'    &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_INTEGER=NN_INTEGER+1
              wrt_int_state%RST_ALL_DATA_I2D(NN_INTEGER)=WORK_ARRAY_I2D(I,J)   !<-- String together this task's 2D integer data
            ENDDO
            ENDDO
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
            RST_KOUNT_R2D=RST_KOUNT_R2D+1                                  !<-- Count # of 2D real Fields (as opposed to 2D integer Fields)
!
            ISTART=LBOUND(WORK_ARRAY_R2D,1)
            IEND  =UBOUND(WORK_ARRAY_R2D,1)
            JSTART=LBOUND(WORK_ARRAY_R2D,2)
            JEND  =UBOUND(WORK_ARRAY_R2D,2)
!
            IF(NN_REAL+IEND*JEND>wrt_int_state%NUM_WORDS_SEND_R2D_RST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF REAL WORDS YOU'       &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_REAL=NN_REAL+1
              wrt_int_state%RST_ALL_DATA_R2D(NN_REAL)=WORK_ARRAY_R2D(I,J)  !<-- String together this task's 2D real data
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO rst_field_block
!
        write_get_fields_tim=write_get_fields_tim+timef()-btim
!
        btim=timef()
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS NOW SEND THEIR STRINGS OF 2D RESTART DATA
!***  TO THE APPROPRIATE WRITE TASKS.
!-----------------------------------------------------------------------
!
        RST_KOUNT_I2D_DATA=wrt_int_state%NUM_WORDS_SEND_I2D_RST            !<-- # of words in 2D integer restart data on this fcst task
        RST_KOUNT_R2D_DATA=wrt_int_state%NUM_WORDS_SEND_R2D_RST            !<-- # of words in 2D real restart data on this fcst task
!
        MYPE_ROW=MYPE/wrt_int_state%INPES+1                                !<-- Each fcst task's row among all rows of fcst tasks
!
        DO N=1,NWTPG                                                       !<-- Loop through the write tasks in this group
          CALL PARA_RANGE(wrt_int_state%JNPES,NWTPG,N                   &  !<-- Find each write task's first and last rows of
                         ,JROW_FIRST,JROW_LAST)                            !<--   fcst tasks from which it will recv
!
          NPE_WRITE=N-1                                                    !<-- Consider the write task with this local ID
                                                                           !    beginning with 0
!
          IF(MYPE_ROW>=JROW_FIRST.AND.MYPE_ROW<=JROW_LAST)THEN             !<-- This fcst task associated with this write task
!
!-----------------------------------------------------------------------
!***  FIRST THE 2-D INTEGER DATA.
!-----------------------------------------------------------------------
!
            IF(RST_KOUNT_I2D>0)THEN
              CALL MPI_ISEND(wrt_int_state%RST_ALL_DATA_I2D             &  !<-- Fcst tasks' string of 2D integer restart data
                            ,RST_KOUNT_I2D_DATA                         &  !<-- #of words in the data string
                            ,MPI_INTEGER                                &  !<-- The datatype
                            ,NPE_WRITE                                  &  !<-- The target write task
                            ,wrt_int_state%NFHOUR                       &  !<-- An MPI tag
                            ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))    &  !<-- The MPI intercommunicator between fcst and quilt tasks
                            ,RST_IH_INT                                 &  !<-- MPI communication request handle
                            ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of integer data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
!-----------------------------------------------------------------------
!***  THEN THE 2-D REAL DATA.
!-----------------------------------------------------------------------
!
            IF(RST_KOUNT_R2D>0)THEN
              CALL MPI_ISEND(wrt_int_state%RST_ALL_DATA_R2D               &  !<-- Fcst tasks' string of 2D real restart data
                            ,RST_KOUNT_R2D_DATA                           &  !<-- #of words in the data string
                            ,MPI_REAL                                     &  !<-- The datatype
                            ,NPE_WRITE                                    &  !<-- The target write task
                            ,wrt_int_state%NFHOUR                         &  !<-- An MPI tag
                            ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))      &  !<-- The MPI intercommunicator between fcst and quilt tasks
                            ,RST_IH_REAL                                  &  !<-- MPI communication request handle
                            ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of real data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
          ENDIF
!
        ENDDO
        write_send_data_tim=write_send_data_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF rst_fcst_tasks
!
      write_run_tim=write_run_tim+timef()-btim0
!rst-rst ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!-----------------------------------------------------------------------
!***  THE FORECAST TASKS ARE COMPLETELY FINISHED WITH HISTORY AND
!***  RESTART OUTPUT NOW SO THEY WILL EXIT THE ROUTINE AND RESUME
!***  THE INTEGRATION.
!-----------------------------------------------------------------------
!
      IF(MYPE<=LAST_FCST_TASK)RETURN
!
!-----------------------------------------------------------------------
!
      history_time: IF(TIME_FOR_HISTORY) THEN
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK IN THE ACTIVE WRITE GROUP RECEIVES THE
!***  STRINGS OF 2D HISTORY DATA FROM THE APPROPRIATE FCST TASKS.
!-----------------------------------------------------------------------
!
      ID_START=wrt_int_state%ID_FTASK_RECV_STA(MYPE)                       !<-- First fcst task that sends to this write task
      ID_END  =wrt_int_state%ID_FTASK_RECV_END(MYPE)                       !<-- Last fcst task that sends to this write task
      NFCST_TASKS=ID_END-ID_START+1                                        !<-- Number of fcst tasks sending to this write task
!
!-----------------------------------------------------------------------
      hst_from_fcst_tasks: DO N=1,NFCST_TASKS                              !<-- Loop through fcst tasks sending to this write task
!-----------------------------------------------------------------------
!
        ID_RECV=ID_START+N-1
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D INTEGER DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        IF(KOUNT_I2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_I2D                      &  !<-- Fcst tasks' string of 2D integer history data
                       ,wrt_int_state%NUM_WORDS_RECV_I2D_HST(ID_RECV)   &  !<-- # of integer words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR                            &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))         &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D REAL DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        IF(KOUNT_R2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_R2D                      &  !<-- Fcst tasks' string of 2D real history data
                       ,wrt_int_state%NUM_WORDS_RECV_R2D_HST(ID_RECV)   &  !<-- # of real words in the data string
                       ,MPI_REAL                                        &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR                            &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))         &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
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
        ITS=wrt_int_state%LOCAL_ISTART(ID_RECV)                            !<-- Local domain integration limits of sending fcst task
        ITE=wrt_int_state%LOCAL_IEND(ID_RECV)                              !<--
        JTS=wrt_int_state%LOCAL_JSTART(ID_RECV)                            !<--
        JTE=wrt_int_state%LOCAL_JEND(ID_RECV)                              !<--
!
        IHALO=wrt_int_state%IHALO                                          !<-- Subdomain halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Subdomain halo depth in J
!
        IMS=ITS-IHALO
        IME=ITE+IHALO
        JMS=JTS-JHALO
        JME=JTE+JHALO
!
        NN=0
!
        DO NF=1,KOUNT_I2D                                                 !<-- Loop through all the 2D integer fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                           !<-- Exclude halo points
            wrt_int_state%WRITE_SUBSET_I(I,J,NF)=wrt_int_state%ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
        NN=0
!
        DO NF=1,KOUNT_R2D                                                 !<-- Loop through all the 2D real fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                           !<-- Exclude halo points
            wrt_int_state%WRITE_SUBSET_R(I,J,NF)=wrt_int_state%ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO hst_from_fcst_tasks
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
!-----------------------------------------------------------------------
      hst_write_begin: IF(MYPE==LEAD_WRITE_TASK)THEN                       !<-- The lead write task
!-----------------------------------------------------------------------
!
        IF(wrt_int_state%WRITE_HST_FLAG.OR.                             &
           wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Current ESMF Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
      ENDIF hst_write_begin
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
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        IF(wrt_int_state%WRITE_HST_FLAG)THEN
!
          CALL WRITE_RUNHISTORY_OPEN(WRT_INT_STATE                      &
                               ,IYEAR_FCST                              &
                               ,IMONTH_FCST                             &
                               ,IDAY_FCST                               &
                               ,IHOUR_FCST                              &
                               ,IMINUTE_FCST                            &
                               ,SECOND_FCST                             &
                               ,NF_HOURS                                &
                               ,NF_MINUTES                              &
                               ,NF_SECONDS                              &
                               ,HST_FIRST                               &
                               ,LEAD_WRITE_TASK)
        ENDIF
!
        IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
          DEGRAD=90./ASIN(1.)
          CALL WRITE_NEMSIO_RUNHISTORY_OPEN(WRT_INT_STATE               &
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
                                      ,DIM1,DIM2,NBDR,GLOBAL            &
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
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
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
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%BUFF_INT                          &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
!
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
          write(0,*)'lead_write_task, send and recv,LAST_WRITE_TASK=',  &
           'LEAD_WRITE_TASK=',LEAD_WRITE_TASK
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%WRITE_SUBSET_I(I,J,NFIELD) !<-- Lead write task fills its part of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%BUFF_INT                      &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%BUFF_INT(NN)  !<-- Insert other write tasks' subsections into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_I2D_STRING(NPOSN_1:NPOSN_2)                        !<-- The name of this 2D integer history quantity
! 
          IF(wrt_int_state%WRITE_HST_FLAG)THEN
!
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)wrt_int_state%OUTPUT_ARRAY_I2D  !<-- Lead write task writes out the 2D real data
!
            IF(HST_FIRST)THEN
              WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_HST_UNIT
            ENDIF
          ENDIF
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
      WRITE(MODEL_LEVEL,'(I2.2)')wrt_int_state%LM(1)
!
!-----------------------------------------------------------------------
      field_loop_real: DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)              !<-- Loop through all 2D real gridded history data
!-----------------------------------------------------------------------
!
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
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
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%BUFF_REAL                         &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
!
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%WRITE_SUBSET_R(I,J,NFIELD) !<-- Lead write task fills its part of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%BUFF_REAL                     &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%BUFF_REAL(NN)   !<-- Insert other write tasks' subsections into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_R2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D real history quantity
!
!-----------------------------------------------------------------------
!***  BEGIN COMPUTATION OF THE 10-M WIND FACTOR FOR GSI.
!-----------------------------------------------------------------------
!
          IF(TRIM(NAME)=='U10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
!
              DO J=1,JM
              DO I=1,IM
                FACT10(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                  &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)*          &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='V10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
              FACT10=0.
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                 &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)*         &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='U_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPU(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPU)
            write(0,*)'fact10tmpu=',maxval(fact10tmpu(1:im,1:jm)),minval(fact10tmpu(1:im,1:jm))
          ENDIF
!
          IF(TRIM(NAME)=='V_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPV(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPV)
            write(0,*)'fact10tmpv=',maxval(fact10tmpv(1:im,1:jm)),minval(fact10tmpv(1:im,1:jm))
          ENDIF
!
          IF(wrt_int_state%WRITE_HST_FLAG)THEN
!
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)wrt_int_state%OUTPUT_ARRAY_R2D   !<-- Lead write task writes out the 2D real data
!
            IF(HST_FIRST)THEN
              WRITE(0,*)'Wrote ',TRIM(NAME)                                &
                       ,' to history file unit ',wrt_int_state%IO_HST_UNIT &
                       ,MAXVAL(wrt_int_state%OUTPUT_ARRAY_R2D)             &
                       ,MINVAL(wrt_int_state%OUTPUT_ARRAY_R2D)
            ENDIF
          ENDIF
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
            IF(TRIM(NAME)=='FIS')wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=          &
                  wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)/G
!
            IF(TRIM(NAME)=='GLAT')THEN
               wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=           &
                 wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
               ALLOCATE(GLAT1D(FIELDSIZE))
               GLAT1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='GLON')THEN
               wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=           &
                 wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
               ALLOCATE(GLON1D(FIELDSIZE))
               GLON1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='VLAT')wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=   &
               wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
!
            IF(TRIM(NAME)=='VLON')wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=   &
               wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
!
            N=NFIELD+wrt_int_state%KOUNT_I2D(1)
            TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
!
            IF(HST_FIRST)THEN
              WRITE(0,*)'Wrote ',TRIM(NAME),' to nemsio history file iret=',ierr
            ENDIF

          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_real
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  COMPLETE COMPUTATION OF 10-M WIND FACTOR AND WRITE IT OUT.
!-----------------------------------------------------------------------
!
      write(0,*)'allocated fact10=',allocated(FACT10),allocated(FACT10TMPU),allocated(FACT10TMPV)
      IF( MYPE==LEAD_WRITE_TASK )THEN
        IF(ALLOCATED(FACT10).AND.ALLOCATED(FACT10TMPU).AND.ALLOCATED(FACT10TMPV)) THEN
          write(0,*)'allocated fact10,fact10tmpu,fact10tmpv'
          DO J=1,JM
          DO I=1,IM
            FACT10TMPV(I,J)=SQRT(FACT10TMPU(I,J)*FACT10TMPU(I,J)+       &
                                 FACT10TMPV(I,J)*FACT10TMPV(I,J))
          ENDDO
          ENDDO
!
          write(0,*)'wind mgn=',maxval(FACT10TMPV(1:IM,1:JM)),minval(FACT10TMPV(1:IM,1:JM))
          write(0,*)'wind10 mgn=',maxval(sqrt(FACT10(1:IM,1:JM))),minval(sqrt(FACT10(1:IM,1:JM)))
!
          DO J=1,JM
          DO I=1,IM
            IF(FACT10TMPV(I,J)/=0) THEN
              FACT10(I,J)=SQRT(FACT10(I,J))/FACT10TMPV(I,J)
            ELSE
              FACT10(I,J)=1.
            ENDIF
          ENDDO
          ENDDO
!
          DEALLOCATE(FACT10TMPU)
          DEALLOCATE(FACT10TMPV)
!
          write(0,*)'WRITE_HST_FLAG=',wrt_int_state%WRITE_HST_FLAG,'NEMSIOFLAG=',wrt_int_state%WRITE_NEMSIOFLAG
          IF(wrt_int_state%WRITE_HST_FLAG)THEN
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)FACT10                    !<-- Lead write task writes out the 2D real data
            write(0,*)'WRITE_HST_FLAG=',wrt_int_state%WRITE_HST_FLAG,'rc=',rc
            IF(HST_FIRST)THEN
              WRITE(0,*)'Wrote FACT10 to history file unit ',wrt_int_state%IO_HST_UNIT &
                       ,maxval(fact10),minval(fact10)
            ENDIF
          ENDIF
!
          IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
            N=N+1
            TMP=RESHAPE(FACT10(1:IM,1:JM),(/FIELDSIZE/))
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            write(0,*)'after nemsio_writerec,n=',n,'fact10=',maxval(tmp),minval(tmp),'iret=',ierr
          ENDIF
!
          DEALLOCATE(FACT10)
!
        ENDIF
!
      ENDIF
!
      HST_FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  CLOSE THE DISK FILE IF NEEDED.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_HST_FLAG.and.MYPE==LEAD_WRITE_TASK)THEN
        CLOSE(wrt_int_state%IO_HST_UNIT)
        write(0,*)' Closed history file with unit=',wrt_int_state%IO_HST_UNIT
      ENDIF
!
!------------
!***  NEMSIO
!------------
!
      IF(wrt_int_state%WRITE_NEMSIOFLAG.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
!
        IF(ASSOCIATED(GLAT1D).AND.ASSOCIATED(GLON1D)) THEN
          CALL NEMSIO_SETFILEHEAD(NEMSIOFILE,IERR,GLAT1D,GLON1D)
          DEALLOCATE(GLAT1D,GLON1D)
        ENDIF
!
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,IERR,gfname=GFNAME)
!
        DEALLOCATE(TMP)
!
        CALL NEMSIO_CLOSE(NEMSIOFILE)
        WRITE(0,*)' Closed nemsio_history file, ', gfname
!
        CALL NEMSIO_FINALIZE()
!
      ENDIF
!
!jw
!-----------------
!*** fcstdone file
!-----------------
      IF(wrt_int_state%WRITE_DONEFILEFLAG .and. wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        write(FILENAME,'(A8,I3.3)' ) 'fcstdone',NF_HOURS
        DO N=51,99
          INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_HST_UNIT=N
            EXIT
          ENDIF
        ENDDO
        OPEN(unit  =IO_HST_UNIT                                         &
          ,file  =trim(FILENAME)                                        &
          ,form='formatted'                                             &
          ,status='new')
        WRITE(IO_HST_UNIT,'(A4)')'DONE'
        CLOSE(IO_HST_UNIT)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF history_time
!
!-----------------------------------------------------------------------
!
      restart_time: IF(TIME_FOR_RESTART) THEN
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK IN THE ACTIVE WRITE GROUP RECEIVES THE
!***  STRINGS OF 2D RESTART DATA FROM THE APPROPRIATE FCST TASKS.
!-----------------------------------------------------------------------
!
      ID_START=wrt_int_state%ID_FTASK_RECV_STA(MYPE)                       !<-- First fcst task that sends to this write task
      ID_END  =wrt_int_state%ID_FTASK_RECV_END(MYPE)                       !<-- Last fcst task that sends to this write task
      NFCST_TASKS=ID_END-ID_START+1                                        !<-- Number of fcst tasks sending to this write task
!
!-----------------------------------------------------------------------
      rst_from_fcst_tasks: DO N=1,NFCST_TASKS                              !<-- Loop through fcst tasks sending to this write task
!-----------------------------------------------------------------------
!
        ID_RECV=ID_START+N-1
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D INTEGER DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        IF(RST_KOUNT_I2D>0)THEN
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_I2D                  &  !<-- Fcst tasks' string of 2D integer restart data
                       ,wrt_int_state%NUM_WORDS_RECV_I2D_RST(ID_RECV)   &  !<-- # of words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR                            &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))         &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  RECEIVE 2-D REAL DATA IF THERE IS ANY.
!-----------------------------------------------------------------------
!
        IF(RST_KOUNT_R2D>0)THEN
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_R2D                  &  !<-- Fcst tasks' string of 2D real restart data
                       ,wrt_int_state%NUM_WORDS_RECV_R2D_RST(ID_RECV)   &  !<-- # of words in the data string
                       ,MPI_REAL                                        &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOUR                            &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP(1))         &  !<-- The MPI intercommunicator between quilt and fcst tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK NEEDS TO INSERT THE PIECES OF THE VARIOUS
!***  2D RESTART ARRAYS RECEIVED FROM THE INDIVIDUAL FCST TASKS
!***  INTO ARRAYS THAT SPAN THE WRITE TASKS' OWN SUBSECTION OF
!***  THE FULL 2D DOMAIN.  THAT SUBSECTION ALWAYS SPANS THE 
!***  ENTIRE EAST-WEST DIMENSION OF THE FULL DOMAIN (SINCE FULL
!***  ROWS OF FCST TASKS ALWAYS SEND TO WRITE TASKS, NEVER 
!***  PARTIAL ROWS) AND AS MUCH OF THE NORTH-SOUTH DIMENSION OF
!***  THE FULL DOMAIN AS COVERED BY THOSE FCST TASKS SENDING TO
!***  A GIVEN WRITE TASK.
!-----------------------------------------------------------------------
!
        ITS=wrt_int_state%LOCAL_ISTART(ID_RECV)                            !<-- Local domain integration limits of sending fcst task
        ITE=wrt_int_state%LOCAL_IEND(ID_RECV)                              !<--
        JTS=wrt_int_state%LOCAL_JSTART(ID_RECV)                            !<--
        JTE=wrt_int_state%LOCAL_JEND(ID_RECV)                              !<--
!
        IHALO=wrt_int_state%IHALO                                          !<-- Subdomain halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Subdomain halo depth in J
!
          IMS=ITS-IHALO
          IME=ITE+IHALO
          JMS=JTS-JHALO
          JME=JTE+JHALO
!
        NN=0
!
        DO NF=1,RST_KOUNT_I2D                                              !<-- Loop through all the 2D integer fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                      !<-- Exclude halo points
            wrt_int_state%RST_WRITE_SUBSET_I(I,J,NF)=wrt_int_state%RST_ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
        NN=0
!
        DO NF=1,RST_KOUNT_R2D                                              !<-- Loop through all the 2D real fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                      !<-- Exclude halo points
            wrt_int_state%RST_WRITE_SUBSET_R(I,J,NF)=wrt_int_state%RST_ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO rst_from_fcst_tasks
!
!-----------------------------------------------------------------------
!***  AT THIS POINT, ALL WRITE TASKS HAVE RECEIVED ALL OF THE RESTART
!***  DATA FROM THEIR ASSOCIATED FCST TASKS AND ASSEMBLED IT ONTO 
!***  THEIR OWN SUBSECTIONS OF THE FULL 2D DOMAIN.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  IT IS TIME FOR THE LEAD WRITE TASK TO BEGIN WRITING TO THE
!***  RESTART FILES.  THE LEAD WRITE TASK ALREADY HOLDS ALL OF THE
!***  SCALAR/1D RESTART DATA AND CAN GO AHEAD AND WRITE THEM.
!-----------------------------------------------------------------------
!
      rst_write_begin: IF(MYPE==LEAD_WRITE_TASK)THEN                       !<-- The lead write task
!
!-----------------------------------------------------------------------
!
        IF(wrt_int_state%WRITE_RST_FLAG.OR.                             &
           wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Current ESMF Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock       =CLOCK                         &  !<-- The ESMF Clock
                            ,currTime    =CURRTIME                      &  !<-- The current time (ESMF) on the clock
                            ,advanceCount=NTIMESTEP_ESMF                &  !<-- # of times the clock has advanced
                            ,rc          =RC)
!
          NTIMESTEP=NTIMESTEP_ESMF
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
!
      ENDIF rst_write_begin
!
!-----------------------------------------------------------------------
!***  WE WILL NOW ASSEMBLE THE FULL DOMAIN 2-D RESTART DATA ONTO
!***  THE LEAD WRITE TASK FROM THE SUBSECTIONS ON ALL WRITE TASKS
!***  THEN THE LEAD TASK WILL WRITE EACH 2D FIELD TO THE RESTART
!***  FILE. 
!
!***  NOTE:  THE LEAD WRITE TASK ASSEMBLES AND WRITES TO RESTART ONLY
!***         ONE 2-D FIELD AT A TIME.  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST LOOP THROUGH ALL OF THE INTEGER Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WRITE OUT THE RESTART FILE AS SELECTED.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        IF(wrt_int_state%WRITE_RST_FLAG)THEN
!
          CALL WRITE_RUNRESTART_OPEN(WRT_INT_STATE                      &
                               ,IYEAR_FCST                              &
                               ,IMONTH_FCST                             &
                               ,IDAY_FCST                               &
                               ,IHOUR_FCST                              &
                               ,IMINUTE_FCST                            &
                               ,SECOND_FCST                             &
                               ,NTIMESTEP                               &
                               ,NF_HOURS                                &
                               ,NF_MINUTES                              &
                               ,NF_SECONDS                              &
                               ,RST_FIRST                               &
                               ,LEAD_WRITE_TASK)
        ENDIF
!
        IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
          DEGRAD=90./ASIN(1.)
          CALL WRITE_NEMSIO_RUNRESTART_OPEN(WRT_INT_STATE               &
                                      ,NEMSIOFILE                       &
                                      ,IYEAR_FCST                       &
                                      ,IMONTH_FCST                      &
                                      ,IDAY_FCST                        &
                                      ,IHOUR_FCST                       &
                                      ,IMINUTE_FCST                     &
                                      ,SECOND_FCST                      &
                                      ,NTIMESTEP                        &
                                      ,NF_HOURS                         &
                                      ,NF_MINUTES                       &
                                      ,NF_SECONDS                       &
                                      ,DIM1,DIM2,NBDR,GLOBAL            &
                                      ,LEAD_WRITE_TASK)
          FIELDSIZE=(DIM1+2*NBDR)*(DIM2+2*NBDR)
          ALLOCATE(TMP(FIELDSIZE))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
      rst_field_loop_int: DO NFIELD=1,RST_KOUNT_I2D                        !<-- Loop through all 2D integer gridded restart data
!-----------------------------------------------------------------------
!
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%RST_BUFF_INT(NN)=wrt_int_state%RST_WRITE_SUBSET_I(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%RST_BUFF_INT                      &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
!
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
          write(0,*)'RST lead_write_task, send and recv,LAST_WRITE_TASK=',  &
           'LEAD_WRITE_TASK=',LEAD_WRITE_TASK
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%RST_OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%RST_WRITE_SUBSET_I(I,J,NFIELD) !<-- Lead write task fills its part
                                                                                                 !    of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%RST_BUFF_INT                  &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%RST_OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%RST_BUFF_INT(NN)   !<-- Insert other write tasks' subsections
                                                                                         !    into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
         NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
         NPOSN_2=NFIELD*ESMF_MAXSTR
         NAME=wrt_int_state%RST_NAMES_I2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D integer restart quantity
!
         IF(wrt_int_state%WRITE_RST_FLAG)THEN
!
          WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)wrt_int_state%RST_OUTPUT_ARRAY_I2D   !<-- Lead write task writes out the 2D integer data
!
          IF(RST_FIRST)THEN
            WRITE(0,*)'Wrote ',TRIM(NAME),' to restart file unit ',wrt_int_state%IO_RST_UNIT, &
             maxval(wrt_int_state%RST_OUTPUT_ARRAY_I2D),minval(wrt_int_state%RST_OUTPUT_ARRAY_I2D)
          ENDIF
         ENDIF
!
!-----------------------------------------------------------------------
!***  FOR NEMSIO FILE
!-----------------------------------------------------------------------
!
         IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
           IF(FIELDSIZE/=IM*JM)THEN
             WRITE(0,*)'WRONG: input data dimension ',IM*JM,                &
              ' does not match data size in NEMSIO file ',FIELDSIZE
           ENDIF
           TMP=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_I2D(1:IM,1:JM),(/FIELDSIZE/))
!
           CALL NEMSIO_WRITEREC(NEMSIOFILE,NFIELD,TMP,IRET=IERR)           !<-- Lead write task writes out the 2D int data!
!
         ENDIF
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO rst_field_loop_int
!
!-----------------------------------------------------------------------
!***  NOW LOOP THROUGH ALL THE REAL Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      rst_field_loop_real: DO NFIELD=1,wrt_int_state%RST_KOUNT_R2D(1)      !<-- Loop through all 2D real gridded restart data
!-----------------------------------------------------------------------
!
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%RST_BUFF_REAL(NN)=wrt_int_state%RST_WRITE_SUBSET_R(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%RST_BUFF_REAL                     &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
!
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%RST_WRITE_SUBSET_R(I,J,NFIELD) !<-- Lead write task fills its part
                                                                                                 !    of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%RST_BUFF_REAL                     &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%RST_BUFF_REAL(NN)   !<-- Insert other write tasks' subsections
                                                                                          !    into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%RST_NAMES_R2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D real restart quantity
!
!-----------------------------------------------------------------------
!***  BEGIN COMPUTATION OF THE 10-M WIND FACTOR FOR GSI.
!-----------------------------------------------------------------------
!
          IF(TRIM(NAME)=='U10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
!
              DO J=1,JM
              DO I=1,IM
                FACT10(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                 &
                          wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)*     &
                          wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='V10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
              FACT10=0.
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
             FACT10(I,J)=FACT10(I,J)+                                  &
                         wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)*      &
                         wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='U_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPU(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPU)
          ENDIF
!
          IF(TRIM(NAME)=='V_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPV(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPV)
          ENDIF

!
          IF(wrt_int_state%WRITE_RST_FLAG)THEN
!
            WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)wrt_int_state%RST_OUTPUT_ARRAY_R2D   !<-- Lead write task writes out the 2D real data
!
            IF(RST_FIRST)THEN
              WRITE(0,*)'Wrote ',TRIM(NAME)                                &
                       ,' to restart file unit ',wrt_int_state%IO_RST_UNIT &
                       ,MAXVAL(wrt_int_state%RST_OUTPUT_ARRAY_R2D)         &
                       ,MINVAL(wrt_int_state%RST_OUTPUT_ARRAY_R2D)
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  FOR NEMSIO FILE
!-----------------------------------------------------------------------
!
          IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: data dimension ',IM*JM,                 &
               ' does not match data size in NEMSIO file,',FIELDSIZE
            ENDIF
            IF(TRIM(NAME)=='FIS') THEN
              IF(.NOT.ALLOCATED(HGT)) ALLOCATE(HGT(1:IM,1:JM))
              HGT(1:IM,1:JM)=wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)/G
            ENDIF
!
            IF(TRIM(NAME)=='GLAT')THEN
               wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)=           &
                 wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
               ALLOCATE(GLAT1D(FIELDSIZE))
               GLAT1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='GLON')THEN
               wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)=           &
                 wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
               ALLOCATE(GLON1D(FIELDSIZE))
               GLON1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='VLAT')wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)=   &
               wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
!
            IF(TRIM(NAME)=='VLON')wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)=   &
               wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)*DEGRAD
!
            N=NFIELD+wrt_int_state%RST_KOUNT_I2D(1)
            TMP=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
!
            IF(RST_FIRST)THEN
              WRITE(0,*)'Wrote ',TRIM(NAME),' to nemsio restart file iret=',ierr
            ENDIF

          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO rst_field_loop_real
!
!-----------------------------------------------------------------------
!***  COMPLETE COMPUTATION OF 10-M WIND FACTOR, AND WRITE FACT10 and HGT OUT.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK) THEN
        IF(ALLOCATED(FACT10).AND.ALLOCATED(FACT10TMPU).AND.ALLOCATED(FACT10TMPV)) THEN
!
          DO J=1,JM
          DO I=1,IM
            FACT10TMPV(I,J)=SQRT(FACT10TMPU(I,J)*FACT10TMPU(I,J)+       &
                                 FACT10TMPV(I,J)*FACT10TMPV(I,J) )
!
            IF(FACT10TMPV(I,J)/=0) THEN
              FACT10(I,J)=SQRT(FACT10(I,J))/FACT10TMPV(I,J)
            ELSE
              FACT10(I,J)=1.
            ENDIF
          ENDDO
          ENDDO
!
          DEALLOCATE(FACT10TMPU)
          DEALLOCATE(FACT10TMPV)
!
          IF(wrt_int_state%WRITE_RST_FLAG)THEN
            WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)FACT10                !<-- Lead write task writes out the 2D real data
            IF(RST_FIRST)THEN
              WRITE(0,*)'Wrote FACT10 to restart file unit ',wrt_int_state%IO_RST_UNIT &
                       ,MAXVAL(fact10),MINVAL(fact10)
            ENDIF
          ENDIF

          IF(wrt_int_state%WRITE_NEMSIOFLAG)THEN
            N=N+1
            TMP=RESHAPE(FACT10(1:IM,1:JM),(/FIELDSIZE/))
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
             write(0,*)'after nemsio_writerec,n=',n,'fact10=',maxval(tmp),minval(tmp),'iret=',ierr
          ENDIF
!
          DEALLOCATE(FACT10)

        ENDIF
!
        IF(wrt_int_state%WRITE_RST_FLAG .and. wrt_int_state%WRITE_NEMSIOFLAG)THEN
          N=N+1
          TMP=RESHAPE(HGT(1:IM,1:JM),(/FIELDSIZE/))
          CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
          DEALLOCATE(HGT)
        ENDIF
!
      ENDIF
!
      RST_FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  CLOSE THE DISK FILE IF NEEDED.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_RST_FLAG.and.MYPE==LEAD_WRITE_TASK)THEN
        CLOSE(wrt_int_state%IO_RST_UNIT)
        write(0,*)' Closed restart file with unit=',wrt_int_state%IO_RST_UNIT
      ENDIF
!
!------------
!***  NEMSIO
!------------
!
      IF(wrt_int_state%WRITE_NEMSIOFLAG.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
!
        IF(ASSOCIATED(GLAT1D).AND.ASSOCIATED(GLON1D)) THEN
          CALL NEMSIO_SETFILEHEAD(NEMSIOFILE,IERR,GLAT1D,GLON1D)
          DEALLOCATE(GLAT1D,GLON1D)
        ENDIF
!  
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,IERR,gfname=GFNAME)
!
        DEALLOCATE(TMP)
!
        CALL NEMSIO_CLOSE(NEMSIOFILE)
        WRITE(0,*)' Closed nemsio_restart file, ', gfname
!
        CALL NEMSIO_FINALIZE()
!
      ENDIF
!
!jw
!-----------------
!*** restartdone file
!-----------------
      IF(wrt_int_state%WRITE_DONEFILEFLAG .and. wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        write(FILENAME,'(A11,I3.3)' ) 'restartdone',NF_HOURS
        DO N=51,99
          INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_RST_UNIT=N
            EXIT
          ENDIF
        ENDDO
        OPEN(unit  =IO_RST_UNIT                                         &
          ,file  =trim(FILENAME)                                        &
          ,form='formatted'                                              &
          ,status='new')
        WRITE(IO_RST_UNIT,'(A4)')'DONE'
        CLOSE(IO_RST_UNIT)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF restart_time
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
        WRITE(0,*)"PASS: WRITE_RUN"
      ELSE
        WRITE(0,*)"FAIL: WRITE_RUN"
      ENDIF
!
      write_run_tim=write_run_tim+timef()-btim0
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
        WRITE(0,*)' Write Time is ',write_run_tim*1.e-3 &
                 ,' at Fcst ',NF_HOURS,':',NF_MINUTES,':',NF_SECONDS
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRT_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_FINALIZE(WRT_COMP                                  &
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
      END SUBROUTINE WRT_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_SETUP(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  SET UP THE WRITE COMPONENTS WITH THE FORECAST TASKS AND
!***  THE GROUPS OF WRITE TASKS NEEDED FOR QUILTING THE OUTPUT
!***  AND WRITING IT TO HISTORY FILES.
!-----------------------------------------------------------------------
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE      
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
      TYPE(ESMF_VM)          :: VM                                        !<-- The ESMF virtual machine.
!
      INTEGER                :: INPES,JNPES                             & !<-- Number of fcst tasks in I and J directions
                               ,MYPE                                    & !<-- My task ID
                               ,NUM_PES_FCST                            & !<-- Number of forecast tasks
                               ,WRITE_GROUPS                            & !<-- Number of groups of write tasks
                               ,WRITE_TASKS_PER_GROUP                     !<-- #of tasks in each write group
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
      atm_int_state%WRITE_GROUPS=WRITE_GROUPS                              !<-- Save for the ATM's Finalize step
      atm_int_state%WRITE_TASKS_PER_GROUP=WRITE_TASKS_PER_GROUP            !<-- Save for the ATM's Finalize step
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
      ALLOCATE(atm_int_state%PETLIST_WRITE(NUM_PES_FCST+WRITE_TASKS_PER_GROUP,WRITE_GROUPS)) !<-- Task IDs of all fcst tasks
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
          atm_int_state%PETLIST_WRITE(I+1,J)=I                          !<-- Collect forecast task IDs to be associated with
                                                                        !    write tasks by write group
        ENDDO
!
      ENDDO
!
      K=NUM_PES_FCST
!
      DO J=1,WRITE_GROUPS
        DO I=1,WRITE_TASKS_PER_GROUP
          atm_int_state%PETLIST_WRITE(NUM_PES_FCST+I,J)=K               !<-- Append write task IDs to associated forecast task IDs by group
          K=K+1
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  CREATE THE WRITE GRIDDED COMPONENT(S).
!***  THERE ARE AS MANY WRITE COMPONENTS AS THERE ARE GROUPS OF
!***  WRITE TASKS SPECIFIED IN THE CONFIGURE FILE.
!***  REGISTER THEIR INIT, RUN, AND FINALIZE STEPS.
!-----------------------------------------------------------------------
!
      IF(WRITE_GROUPS>0)THEN
        ALLOCATE(atm_int_state%WRT_COMPS(WRITE_GROUPS))                 !<-- The Write gridded components
      ENDIF
!
!---------------------------------
!***  Create the Write components
!---------------------------------
!
      DO I=1,WRITE_GROUPS
        WRITE(MY_WRITE_GROUP,FMT)I
        WRITE_NAME='write_GridComp_'//MY_WRITE_GROUP
!
        atm_int_state%WRT_COMPS(I)=ESMF_GridCompCreate(                         &
                                name          =WRITE_NAME                       &  !<-- Name of this group's Write gridded component
                               ,configFile    ='configure_file'                 &  !<-- The configure file for writes
                               ,petList       =atm_int_state%PETLIST_WRITE(:,I) &  !<-- The task IDs of the write tasks in this group
                                                                                   !    provide the local VM information per component.
                               ,rc            =RC)
!
      ENDDO
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
        CALL ESMF_GridCompSetServices(atm_int_state%WRT_COMPS(I)         &  !<-- The Write gridded components
                                     ,WRITE_REGISTER                     &  !<-- The user's subroutine name
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
      atm_int_state%IMP_STATE_WRITE=ESMF_StateCreate(statename='Write Import State' &  !<-- Import state name for writes
                                                    ,statetype= ESMF_STATE_IMPORT   &
                                                    ,rc       = RC)
!
      atm_int_state%EXP_STATE_WRITE=ESMF_StateCreate(statename='Write Export State' &  !<-- Export state names for writes
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
      CALL ESMF_StateAdd(state      =atm_int_state%EXP_STATE_DYN        & !<-- Dynamics export state receives a state
                        ,nestedState=atm_int_state%IMP_STATE_WRITE      & !<-- Add the write components' import state
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
      CALL ESMF_StateAdd(state      =atm_int_state%EXP_STATE_PHY        & !<-- Physics export state receives a state
                        ,nestedState=atm_int_state%IMP_STATE_WRITE      & !<-- Add the write components' import state
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE SUBDOMAIN HORIZONTAL LIMITS FOR ALL TASKS INTO THE
!***  DYNAMICS AND PHYSICS IMPORT STATES SINCE THEY WILL BE NEEDED
!***  FOR SPECIFYING OUTPUT FIELDS WITHIN THE DYNAMICS AND PHYSICS
!***  COMPONENTS.  THIS IS ONLY RELEVANT TO THE FORECAST TASKS
!***  SINCE ONLY THEY ENTER THE DYNAMICS AND PHYSICS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  NOTE!!!
!-----------------------------------------------------------------------
!***  THE LOCAL DOMAIN LIMITS NEED TO HAVE BEEN INSERTED INTO THE
!***  ATM COMPONENT'S INTERNAL STATE IN THE CORE-SPECIFIC SETUP
!***  ROUTINE CALLED BY ATM_INITIALIZE.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- This excludes I/O tasks.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Local Domain Limits in Dynamics Imp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_DYN    &  !<-- Dynamics import state receives an attribute
                              ,name     ='LOCAL_ISTART'                 &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_ISTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_DYN    &  !<-- Dynamics import state receives an attribute
                              ,name     ='LOCAL_IEND'                   &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_IEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_DYN    &  !<-- Dynamics import state receives an attribute
                              ,name     ='LOCAL_JSTART'                 &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_JSTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_DYN    &  !<-- Dynamics import state receives an attribute
                              ,name     ='LOCAL_JEND'                   &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_JEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_PHY    &  !<-- Physics import state receives an attribute
                              ,name     ='LOCAL_ISTART'                 &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_ISTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_PHY    &  !<-- Physics import state receives an attribute
                              ,name     ='LOCAL_IEND'                   &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_IEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_PHY    &  !<-- Physics import state receives an attribute
                              ,name     ='LOCAL_JSTART'                 &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_JSTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =atm_int_state%IMP_STATE_PHY    &  !<-- Physics import state receives an attribute
                              ,name     ='LOCAL_JEND'                   &  !<-- The attribute's name
                              ,count    =NUM_PES_FCST                   &  !<-- The attribute's length
                              ,valueList=atm_int_state%LOCAL_JEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_DESTROY(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  DESTROY ALL OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE      
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
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
      DO N=1,atm_int_state%WRITE_GROUPS              
        IF(MYPE>=atm_int_state%PETLIST_WRITE(1,N).AND.                          &
           MYPE<=atm_int_state%PETLIST_WRITE(atm_int_state%WRITE_TASKS_PER_GROUP,N))THEN
!
           CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%WRT_COMPS(N)    &
                                     ,importstate=atm_int_state%EXP_STATE_WRITE &
                                     ,exportstate=atm_int_state%IMP_STATE_WRITE &
                                     ,clock      =CLOCK_ATM                     &
                                     ,phase      =ESMF_SINGLEPHASE              &
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
      CALL ESMF_StateDestroy(atm_int_state%IMP_STATE_WRITE,rc=RC)
      CALL ESMF_StateDestroy(atm_int_state%EXP_STATE_WRITE,rc=RC)
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
      DO J=1,atm_int_state%WRITE_GROUPS
!
        CALL ESMF_VMBarrier(vm=VM,rc=RC)
!
        DO I=1,atm_int_state%NUM_PES_FCST+atm_int_state%WRITE_TASKS_PER_GROUP
          IF(MYPE==atm_int_state%PETLIST_WRITE(I,J))THEN
            CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%WRT_COMPS(J) &
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
      END SUBROUTINE WRITE_DESTROY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_GRID_COMP
!
!-----------------------------------------------------------------------
