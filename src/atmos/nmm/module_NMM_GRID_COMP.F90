!-----------------------------------------------------------------------
!
      MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This is the NMM-B module.  It will set up one or more Domain
!***  subcomponents then execute their Initialize, Run, and Finalize
!***  steps.
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE module_INCLUDE
!
      USE module_NMM_INTERNAL_STATE,ONLY: NMM_INTERNAL_STATE            &
                                         ,WRAP_NMM_INTERNAL_STATE
!
      USE module_DOMAIN_GRID_COMP,ONLY: DOMAIN_REGISTER                    !<-- The Register routine for DOMAIN_GRID_COMP
!
      USE module_NMM_INTEGRATE,ONLY: NMM_INTEGRATE
!
      USE module_NESTING,ONLY: PARENT_CHILD_COMMS
!
      USE module_PARENT_CHILD_CPL_COMP,ONLY: PARENT_CHILD_CPL_REGISTER  &  !<-- The Register routine for PARENT_CHILD Coupler
                                            ,PARENT_CHILD_COUPLER_SETUP
!
      USE module_CONTROL,ONLY: TIMEF
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NMM_REGISTER
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: MAX_DOMAINS=99                       !<-- There cannot be more domains than this
!
      INTEGER(kind=KINT) :: MYPE                                        &  !<-- Each MPI task ID
                           ,NHOURS_CLOCKTIME                            &  !<-- Hours between clocktime prints
                           ,NPE_PRINT                                   &  !<-- Clocktime diagnostics from this MPI task
                           ,NTIMESTEP                                   &  !<-- The integration timestep
                           ,TIMESTEP_SEC_WHOLE                          &
                           ,TIMESTEP_SEC_NUMERATOR                      &
                           ,TIMESTEP_SEC_DENOMINATOR
!
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: DT                       !<-- The fundamental timestep (s) of the domains
!
      CHARACTER(ESMF_MAXSTR) :: CLOCK_NMM_NAME                             !<-- Name of the NMM's ESMF Clock
!
!
      LOGICAL(kind=KLOG) :: RESTARTED_RUN                               &  !<-- Flag indicating if this is a restarted run
                           ,RST_OUT_00                                     !<-- Shall we write 00h history in restarted run?
!
      TYPE(ESMF_VM),SAVE :: VM                                             !<-- The ESMF virtual machine.
!
      TYPE(ESMF_Time),SAVE :: STARTTIME                                    !<-- The ESMF start time.
!
      TYPE(ESMF_TimeInterval),SAVE :: INTERVAL_CLOCKTIME                &  !<-- ESMF time interval between clocktime prints (h)
                                     ,RUNDURATION                          !<-- The ESMF simulation length (sec)
!
      TYPE(ESMF_TimeInterval),DIMENSION(:),ALLOCATABLE :: INTERVAL_HISTORY &  !<-- ESMF time interval between history output (h)
                                                         ,INTERVAL_RESTART &  !<-- ESMF time interval between restart output (h)
                                                         ,TIMESTEP            !<-- The ESMF timestep (s)
!
      TYPE(ESMF_Clock),DIMENSION(:),ALLOCATABLE :: CLOCK_NMM               !<-- The NMM ESMF Clocks
!
      TYPE(NMM_INTERNAL_STATE),POINTER,SAVE :: NMM_INT_STATE               !<-- The NMM component internal state pointer
!
      TYPE(WRAP_NMM_INTERNAL_STATE),SAVE :: WRAP                           !<-- The F90 wrap of the NMM internal state
!
!---------------------
!***  For NMM Nesting
!---------------------
!
      INTEGER(kind=KINT),SAVE :: NUM_DOMAINS=0
!
      INTEGER(kind=KINT) :: KOUNT_STEPS=0
!
      INTEGER(kind=KINT) :: COMM_FULL_DOMAIN                            &  !<-- Communicator for ALL tasks on domain to be split
                           ,COMM_MY_DOMAIN                              &  !<-- Each domain's local intracommunicator
                           ,COMM_TO_MY_PARENT                           &  !<-- Intercommunicator between a domain and its parent
                           ,MY_DOMAIN_ID                                &  !<-- The ID of each domain
                           ,PARENT_CHILD_TIME_RATIO                        !<-- Ratio of parent timestep to child's
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: COMM_TO_MY_CHILDREN    &  !<-- Array of intercommunicators between domain and its children
                                                ,ID_DOMAINS             &  !<-- IDs of all domains
                                                ,ID_PARENTS             &  !<-- IDs of all domains' parents
                                                ,FTASKS_DOMAIN          &  !<-- # of forecast tasks on each domain excluding descendents
                                                ,NTASKS_DOMAIN          &  !<-- # of tasks on each domain excluding descendents
                                                ,NUM_CHILDREN              !<-- # of children on each domain
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: ID_CHILDREN          &  !<-- IDs of all children of all domains
                                                  ,PETLIST_DOMAIN          !<-- List of task IDs for each domain (DOMAIN Component)
!
      LOGICAL(kind=KLOG) :: NESTING_NMM
!
      TYPE(ESMF_Config),DIMENSION(MAX_DOMAINS),SAVE :: CF                  !<-- The config objects (one per domain)
!
      TYPE(ESMF_Logical),SAVE :: I_AM_A_FCST_TASK                       &  !<-- Am I a forecast task?
                                ,I_AM_A_NEST                               !<-- Am I in a nested domain?
!
      TYPE(ESMF_State),SAVE :: IMP_STATE_CPL_NEST                       &
                              ,EXP_STATE_CPL_NEST
!
      TYPE(ESMF_CplComp),SAVE :: PARENT_CHILD_COUPLER_COMP                 !<-- Coupler component for parent-child/nest exchange
!
      TYPE(ESMF_Alarm),SAVE :: ALARM_CLOCKTIME                          &  !<-- The ESMF Alarm for clocktime prints
                              ,ALARM_RECV_FROM_PARENT                      !<-- The ESMF Alarm for child to recv data from parent
!
!-----------------------------------------------------------------------
!
!---------------------------
!***  For Digital Filtering
!---------------------------
!
      INTEGER(kind=KINT) :: DFIHR
!
!-----------------------------------------------------------------------
!
!-----------
!*** Timing
!-----------
!
      REAL(kind=KDBL) :: btim,btim0                                     &
                        ,atm_drv_init                                   &
                        ,atm_drv_run_1                                  &
                        ,atm_drv_run_2                                  &
                        ,atm_drv_run_cpl1                               &
                        ,atm_drv_run_cpl2
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_REGISTER(NMM_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  Register the NMM component's Initialize, Run, and Finalize steps.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP                   !<-- The NMM component
      INTEGER            ,INTENT(OUT)   :: RC_REG
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      write(0,*) "    NMM_REGISTER"
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Initialize routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETINIT                      &
                                     ,NMM_INITIALIZE                    &
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Run routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETRUN                       &
                                     ,NMM_RUN                           &
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Finalize routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETFINAL                     &
                                     ,NMM_FINALIZE                      &
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_REGISTER succeeded'
      ELSE
        WRITE(0,*)' NMM_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
      write(0,*) "    END OF NMM_REGISTER"
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_INITIALIZE(NMM_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_NEMS                              &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  This routine creates the individual DOMAIN gridded components 
!***  and executes their Initialize step.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP                   !<-- The NMM component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The NMM import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                       !<-- The NMM export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NEMS                      !<-- The NEMS ESMF Clock
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_DOM,ID_X,ISTAT,N
!
      INTEGER(kind=KINT) :: MINUTES_HISTORY                             &  !<-- Hours between history output
                           ,MINUTES_RESTART                             &  !<-- Hours between restart output
                           ,NHOURS_FCST                                 &  !<-- Length of forecast in hours
                           ,NSECONDS_FCST                                  !<-- Length of forecast in seconds
!
      INTEGER(kind=KINT) :: INPES,JNPES,LENGTH,N_TASKS                  &
                           ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
      INTEGER(kind=KINT),DIMENSION(MAX_DOMAINS) :: DOMAIN_ID_TO_RANK=0  &  !<-- The configure file associated with each domain ID
                                                  ,RANK_TO_DOMAIN_ID=0     !<-- The domain ID associated with each configure file
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: CHILD_ID               &
                                                ,PETLIST
!
      LOGICAL(kind=KLOG) :: CFILE_EXIST
!
      CHARACTER(2) :: INT_TO_CHAR
      CHARACTER(6) :: FMT='(I2.2)'
      CHARACTER(7) :: MODE
      CHARACTER(MAX_DOMAINS) :: CONFIG_FILE_NAME
!
      CHARACTER(ESMF_MAXSTR) :: DOMAIN_COMP_BASE='DOMAIN Gridded Component ' &
                               ,DOMAIN_GRID_COMP_NAME
!
      TYPE(ESMF_TimeInterval) :: TIMEINTERVAL_RECV_FROM_PARENT             !<-- ESMF time interval between Recv times from parent
!
      TYPE(ESMF_Logical) :: PHYSICS_ON                                     !<-- Does the integration include physics?
!
      TYPE(ESMF_Config) :: CF_X                                            !<-- Working config object
!
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      write(0,*)' enter NMM_INITIALIZE'
      RC      =ESMF_SUCCESS
      RC_INIT = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Allocate the NMM component's internal state.
!-----------------------------------------------------------------------

      ALLOCATE(NMM_INT_STATE,stat=RC)
      wrap%NMM_INT_STATE=>NMM_INT_STATE
!
!-----------------------------------------------------------------------
!***  Attach the NMM internal state to the NMM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach NMM Internal State to the NMM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(NMM_GRID_COMP                  &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the NMM component.
!***  We need VM now to obtain the MPI task IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from NMM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=NMM_GRID_COMP                      &  !<-- The NMM component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the NMM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain MPI Task IDs from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Each MPI task ID
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load all of the configure objects.  All domains
!***  are functionally equivalent thus each has its own configure
!***  file.  We are counting the configure files as we create the
!***  ESMF configure objects so we will know how many different
!***  domains there are.
!-----------------------------------------------------------------------
!
      DO N=1,MAX_DOMAINS                                                   !<-- Number of config files must not exceed 99
!
        WRITE(INT_TO_CHAR,FMT)N
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Each configure file has a unique number.
!
        CFILE_EXIST=.FALSE.
        INQUIRE(FILE=CONFIG_FILE_NAME,EXIST=CFILE_EXIST)
!
        IF(CFILE_EXIST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Create Temporary Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF_X=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Load the Temp Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF_X                        &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Extract Domain ID From Temp Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_X                      &  !<-- The config object
                                      ,value =ID_X                      &  !<-- The domain's ID
                                      ,label ='my_domain_id:'           &  !<-- Take value from this config labelious variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF(ID_X)=ESMF_ConfigCreate(rc=RC)                                !<-- Domain's ID is its element in the CF array
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Destroy Temporary Config Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigDestroy(config=CF_X                           &  !<-- The temporary config object
                                 ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Load the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(ID_X)                    &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DOMAIN_ID_TO_RANK(ID_X)=N                                        !<-- The configure file rank for a given domain ID
          RANK_TO_DOMAIN_ID(N)=ID_X                                        !<-- The domain ID for a given configure file rank
!
          NUM_DOMAINS=NUM_DOMAINS+1
!
        ELSE
!
          EXIT
!
        ENDIF
!
      ENDDO
!
      NESTING_NMM=.FALSE.
      IF(NUM_DOMAINS>1)NESTING_NMM=.TRUE.                                  !<-- We have nests if more than one domain is present
!
!-----------------------------------------------------------------------
!***  Obtain the global communicator for all tasks in this run.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGET(vm             =VM                                &  !<-- The virtual machine
                     ,mpiCommunicator=COMM_FULL_DOMAIN                  &  !<-- Global intracommunicator for all tasks
                     ,rc             =RC)
!
!-----------------------------------------------------------------------
!***  IF NESTED DOMAINS ARE BEING USED THEN:
!***    (1) Split the MPI Communicator between all domains;
!***    (2) Create a DOMAIN subcomponent for all domains;
!***    (3) Call DOMAIN_INIT recursively for all domains.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      nesting_block_1: IF(NESTING_NMM)THEN                                 !<-- Intercommunicators are needed only for nesting
!
!-----------------------------------------------------------------------
!***  Split the global communicator among all NMM domains and create
!***  Parent-Child intercommunicators.
!-----------------------------------------------------------------------
!
        CALL PARENT_CHILD_COMMS(MYPE                                    &  !<-- This task's global rank (in)
                               ,NUM_DOMAINS                             &  !<-- Total number of domains, all generations (in)
                               ,CF                                      &  !<-- Configure objects for all domains (in)
                               ,COMM_FULL_DOMAIN                        &  !<-- Intracommunicator for ALL tasks (in)
                               ,MY_DOMAIN_ID                            &  !<-- ID of domain on which this task resides (out)
                               ,ID_DOMAINS                              &  !<-- IDs of all domains (out)
                               ,ID_PARENTS                              &  !<-- ID of all domains' parents (out)
                               ,NUM_CHILDREN                            &  !<-- # of children on each domain (out)
                               ,ID_CHILDREN                             &  !<-- IDs of all children of all domains (out)
                               ,COMM_MY_DOMAIN                          &  !<-- Communicators for each individual domain (out)
                               ,COMM_TO_MY_PARENT                       &  !<-- Children's intercommunicators to their parents (out)
                               ,COMM_TO_MY_CHILDREN                     &  !<-- Parents' intercommunicators to their children (out)
                               ,FTASKS_DOMAIN                           &  !<-- # of fcst tasks on each domain excluding descendents (out)
                               ,NTASKS_DOMAIN                           &  !<-- # of tasks on each domain excluding descendents (out)
                               ,PETLIST_DOMAIN                          &  !<-- List of task IDs for each domain (DOMAIN Component) (out)
                               ,RANK_TO_DOMAIN_ID )                        !<-- Domain IDs for each configure file
!
!-----------------------------------------------------------------------
!
      ELSE nesting_block_1                                                 !<-- There is only a single domain
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  How many forecast/write tasks will be active on the domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract INPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =INPES                       &  !<-- The domain's fcst tasks in I
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract JNPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =JNPES                       &  !<-- The domain's fcst tasks in J
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Write_Groups From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_GROUPS                &  !<-- The number of Write groups on this domain
                                    ,label ='write_groups:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Write Tasks Per Group From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_TASKS_PER_GROUP       &  !<-- The number of tasks per Write group
                                    ,label ='write_tasks_per_group:'    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ALLOCATE(NUM_CHILDREN(1))
        NUM_CHILDREN(1)=0
!
        N_TASKS=INPES*JNPES+WRITE_GROUPS*WRITE_TASKS_PER_GROUP             !<-- Total # of tasks on the domain
        ALLOCATE(NTASKS_DOMAIN(1))
        NTASKS_DOMAIN(1)=N_TASKS
        ALLOCATE(PETLIST_DOMAIN(1:N_TASKS,1))
!
        DO N=1,N_TASKS
          PETLIST_DOMAIN(N,1)=N-1                                          !<-- The list of task IDs for the DOMAIN Component
        ENDDO
!
        ALLOCATE(ID_DOMAINS(1))
        ID_DOMAINS(1)=1                                                    !<-- There is a single domain; its ID is 1
        MY_DOMAIN_ID=1
!
        ALLOCATE(ID_CHILDREN(1,1))
        ID_CHILDREN(1,1)=0                                                 !<-- A single domain thus no children
!
        ALLOCATE(ID_PARENTS(1))
        ID_PARENTS(1)=-999                                                 !<-- There is a single domain; it has no valid parent
!
        COMM_TO_MY_PARENT=-999                                             !<-- There is a single domain; it has no parent
!
!-----------------------------------------------------------------------
!
      ENDIF nesting_block_1
!
!-----------------------------------------------------------------------
!***  Allocate the DOMAIN import/export states.
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%IMP_STATE_DOMAIN(1:NUM_DOMAINS),stat=ISTAT)
      ALLOCATE(nmm_int_state%EXP_STATE_DOMAIN(1:NUM_DOMAINS),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  Create the DOMAIN import/export states.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%IMP_STATE_DOMAIN(ID_DOM)=ESMF_StateCreate(        &  !<-- DOMAIN import state
                                       statename='Domain Import State'  &  !<-- DOMAIN import state name
                                      ,statetype= ESMF_STATE_IMPORT     &
                                      ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create the DOMAIN Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
        nmm_int_state%EXP_STATE_DOMAIN(ID_DOM)=ESMF_StateCreate(        &  !<-- DOMAIN export state
                                       statename='Domain Export State'  &  !<-- DOMAIN export state name
                                      ,statetype= ESMF_STATE_EXPORT     &
                                      ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Extract Restart Flag from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                  ,value =RESTARTED_RUN               &  !<-- Logical flag indicating if this is a restarted run
                                  ,label ='restart:'                  &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Rst_out_00 Flag from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                  ,value =RST_OUT_00                  &  !<-- Logical flag indicating if this is a restarted run
                                  ,label ='rst_out_00:'               &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Each task will create Clocks for all domains for simplicity
!***  in executing the major DO loops over the DOMAIN components.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create the domains' clocks with their timesteps, start times,
!***  and run durations.
!-----------------------------------------------------------------------
!
      ALLOCATE(CLOCK_NMM(1:NUM_DOMAINS))
      ALLOCATE(TIMESTEP (1:NUM_DOMAINS))
      ALLOCATE(DT       (1:NUM_DOMAINS))
!
      ALLOCATE(INTERVAL_HISTORY(1:NUM_DOMAINS))
      ALLOCATE(INTERVAL_RESTART(1:NUM_DOMAINS))
!
!-----------------------------------------------------------------------
!***  Extract timestep information and history/restart output frequency
!***  from the config files of all domains.
!-----------------------------------------------------------------------
!
      timeinfo_loop: DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Timestep from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain
                                    ,value =TIMESTEP_SEC_WHOLE          &  !<-- The variable filled (integer part of timestep (sec))
                                    ,label ='dt_int:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain 
                                    ,value =TIMESTEP_SEC_NUMERATOR      &  !<-- The variable filled (numerator of timestep fraction)
                                    ,label ='dt_num:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain 
                                    ,value =TIMESTEP_SEC_DENOMINATOR    &  !<-- The variable filled (denominator of timestep fraction)
                                    ,label ='dt_den:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Establish the timesteps for all of the domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Set Timestep Interval"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP(ID_DOM)         &  !<-- The fundamental timestep on this domain (sec) (ESMF)
                                 ,s           =TIMESTEP_SEC_WHOLE       &
                                 ,sn          =TIMESTEP_SEC_NUMERATOR   &
                                 ,sd          =TIMESTEP_SEC_DENOMINATOR &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DT(ID_DOM)=TIMESTEP_SEC_WHOLE+                                  &  !<-- The domain's fundamental timestep (sec) (REAL)
                   REAL(TIMESTEP_SEC_NUMERATOR)                         &
                  /REAL(TIMESTEP_SEC_DENOMINATOR)
!
!-----------------------------------------------------------------------
!***  Get the NMM history output interval (hours) from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Obtain History Interval from the Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The configure object of this domain
                                    ,value =MINUTES_HISTORY             &  !<-- Fill this variable
                                    ,label ='minutes_history:'          &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ESMF history file output time interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the History Output Time Interval."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_HISTORY(ID_DOM) &  !<-- Time interval between
                                 ,m           =MINUTES_HISTORY          &  !<-- Minutes between history output
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Get the NMM restart output interval (hours) from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Obtain Restart Interval from the Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The configure object of this domain
                                    ,value =MINUTES_RESTART             &  !<-- Fill this variable
                                    ,label ='minutes_restart:'          &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ESMF restart file output time interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Restart Output Time Interval."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_RESTART(ID_DOM) &  !<-- Time interval between restart output (ESMF)
                                 ,m           =MINUTES_RESTART          &  !<-- Minutes between restart output (integer)
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDDO timeinfo_loop
!
!-----------------------------------------------------------------------
!***  Obtain the forecast start time from the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Start Time from NMM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_NEMS                         &  !<-- The NEMS ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time (ESMF)
!!!                     ,runDuration=RUNDURATION                        &  !<-- The simulation run duration (ESMF)
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      clock_loop: DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!***  Obtain the forecast length time from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Forecast Length from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                    &
                                    ,value =NHOURS_FCST                   &
                                    ,label ='nhours_fcst:'                &
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NSECONDS_FCST=NHOURS_FCST*3600                                     !<-- The forecast length (sec) (REAL)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Set the Forecast Length"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION              &  !<-- The forecast length (sec) (ESMF)
                                 ,s           =NSECONDS_FCST            &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  With data from above, create the ESMF Clocks to control
!***  the timestepping within the DOMAIN subcomponent(s).
!***  Each domain will set its own clock in the initialize
!***  step of DOMAIN_GRID_COMP.
!-----------------------------------------------------------------------
!
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        CLOCK_NMM_NAME='CLOCK_NMM_'//INT_TO_CHAR
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clocks for the NMM Domains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CLOCK_NMM(N)=ESMF_ClockCreate(name       =CLOCK_NMM_NAME        &  !<-- The NMM Domain's Clock's name
                                     ,timeStep   =TIMESTEP(ID_DOM)      &  !<-- The fundamental timestep in this component
                                     ,startTime  =STARTTIME             &  !<-- Start time of simulation
                                     ,runDuration=RUNDURATION           &  !<-- Duration of simulation
                                     ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDDO clock_loop
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the DOMAIN gridded component(s).
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%DOMAIN_GRID_COMP(1:NUM_DOMAINS),stat=ISTAT)
!
      IF(ISTAT/=0)THEN
        WRITE(0,*)' ERROR: Failed to allocate DOMAIN_GRID_COMP'
        WRITE(6,*)' ERROR: Failed to allocate DOMAIN_GRID_COMP'
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the DOMAIN gridded component(s).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      domain_comp_create: DO N=1,NUM_DOMAINS
!
!-----------------------------------------------------------------------
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        DOMAIN_GRID_COMP_NAME=DOMAIN_COMP_BASE//INT_TO_CHAR                !<-- Append domain ID to DOMAIN Comp name
        N_TASKS=NTASKS_DOMAIN(ID_DOM)                                      !<-- # of tasks on this domain
        PETLIST=>PETLIST_DOMAIN(1:N_TASKS,ID_DOM)                          !<-- The PETlist for this domain
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create DOMAIN_GRID_COMP"//INT_TO_CHAR
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%DOMAIN_GRID_COMP(ID_DOM)=ESMF_GridCompCreate(     &  !<-- The DOMAIN Component for this domain
                                         name   =DOMAIN_GRID_COMP_NAME  &  !<-- Name of the new DOMAIN gridded component
                                        ,config =CF(ID_DOM)             &  !<-- This domain's configure file
                                        ,petList=PETLIST                &  !<-- The IDs of tasks that will run on this domain
                                        ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Register the DOMAIN components' Init, Run, Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register DOMAIN Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(nmm_int_state%DOMAIN_GRID_COMP(ID_DOM) &  !<-- The DOMAIN component
                                     ,DOMAIN_REGISTER                        &  !<-- User's subroutineName
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        IF(ID_DOM/=MY_DOMAIN_ID)CYCLE                                      !<-- Only need to load Import State properly for my domain
!
!-----------------------------------------------------------------------
!***  Check the configure flag indicating whether or not to run
!***  adiabatically (i.e., with no physics).  Insert the flag
!***  into the DOMAIN import state.
!-----------------------------------------------------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =MODE                        &
                                    ,label ='adiabatic:'                &
                                    ,rc    =rc)
!
        IF(TRIM(MODE)=='true')THEN
          PHYSICS_ON=ESMF_False
          IF(MYPE==0) WRITE(0,*)' NMM will run without physics.'
        ELSE
          PHYSICS_ON=ESMF_True
          IF(MYPE==0) WRITE(0,*)' NMM will run with physics.'
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Physics flag to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component's import state
                              ,name ='PHYSICS_ON'                           &  !<-- The flag indicating if physics is active
                              ,value=PHYSICS_ON                             &  !<-- The value being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the maximum number of domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add MAX_DOMAINS to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component's import state
                              ,name ='MAX_DOMAINS'                          &  !<-- Maximum # of domains
                              ,value=MAX_DOMAINS                            &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the domain IDs into the DOMAIN import state(s) along with
!***  the number of children and the children's domain IDs.
!***  Also insert a flag as to whether the DOMAIN component is a nest.
!
!***  Note that all tasks are aware of all domains' IDs,
!***  number of children, and those children's domain IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Domain IDs to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component's import state
                              ,name ='DOMAIN_ID'                            &  !<-- This DOMAIN Component's domain ID
                              ,value=ID_DOM                                 &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the association of configure file IDs with domain IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Configure File ID Associated With Each Domain ID"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component's import state
                              ,name     ='DOMAIN_ID_TO_RANK'                    &  !<-- Adding Attribute with this name
                              ,count    =MAX_DOMAINS                            &  !<-- Total # of domains
                              ,valueList=DOMAIN_ID_TO_RANK                      &  !<-- Configure file IDs linked to each domain
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Number of Children to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component's import state
                              ,name ='NUM_CHILDREN'                         &  !<-- This DOMAIN Component's # of children
                              ,value=NUM_CHILDREN(ID_DOM)                   &  !<-- Insert this into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(COMM_TO_MY_PARENT<0)THEN
          I_AM_A_NEST=ESMF_FALSE
        ELSE
          I_AM_A_NEST=ESMF_TRUE
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Nest/Not-a-Nest Flag to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component import state
                              ,name ='I-Am-A-Nest Flag'                     &  !<-- Name of Attribute
                              ,value=I_AM_A_NEST                            &  !<-- Logical nest flag
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nesting_block_2: IF(NESTING_NMM)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Domain IDs of Children to the DOMAIN Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          LENGTH=MAX(1,NUM_CHILDREN(ID_DOM))
          CHILD_ID=>ID_CHILDREN(1:LENGTH,ID_DOM)                           !<-- Select only the IDs of this domain's children
!
          CALL ESMF_AttributeSet(state    =nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component import state
                                ,name     ='CHILD_IDs'                            &  !<-- The children's IDs of this DOMAIN Component
                                ,count    =LENGTH                                 &  !<-- Length of inserted array
                                ,valueList=CHILD_ID                               &  !<-- Insert this into the import state
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I_AM_A_NEST==ESMF_TRUE)THEN
!
            PARENT_CHILD_TIME_RATIO=NINT(DT(ID_PARENTS(ID_DOM))/DT(ID_DOM))   !<-- Ratio of parent's timestep to this nest's
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Add Parent-Child Time Ratio to DOMAIN Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=nmm_int_state%IMP_STATE_DOMAIN(ID_DOM) &  !<-- This DOMAIN component import state
                                  ,name ='Parent-Child Time Ratio'              &  !<-- Name of Attribute
                                  ,value=PARENT_CHILD_TIME_RATIO                &  !<-- # of child timesteps per parent timestep
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
        ENDIF nesting_block_2
!
!-----------------------------------------------------------------------
!
      ENDDO domain_comp_create
!
!-----------------------------------------------------------------------
!***  At this point, DOMAIN components for each domain have been 
!***  created and registered.  Now they need to be initialized.
!
!***  The following call will initialize DOMAIN_GRID_COMP for domain #1.
!***  If more than one domain exists, domain #1 is the uppermost and
!***  the remaining domains will be initialized recursively through
!***  the generations of children.  Recursion is necessary because
!***  children must not be initialized before their parents since
!***  a parent might be directed by the user to generate input data
!***  for its children and that must be complete before the parent's
!***  children are initialized and try to read their input data
!***  before it exists.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      CALL CALL_DOMAIN_INITIALIZE(1,CLOCK_NMM)                             !<-- Initiate cascade of DOMAIN Initialize calls for all domains
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Identify the forecast vs. quilt/write tasks since Parent-Child
!***  interaction does not involve any Write tasks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Extract Fcst-or-Write Flag from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) &  !<-- The DOMAIN component export state
                            ,name ='Fcst-or-Write Flag'                         &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_FCST_TASK                             &  !<-- Am I a forecast task?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If there are nests then create a Parent-Child coupler through
!***  which parents will send boundary data to their children.
!***  Load that coupler's import state with the data the parents
!***  need to generate boundary data for their children.
!-----------------------------------------------------------------------
!
      nesting_block_3: IF(NESTING_NMM)THEN                                 !<-- All parents and children create the Coupler.
!
        CALL PARENT_CHILD_COUPLER_SETUP(NUM_DOMAINS                                  &  !
                                       ,MY_DOMAIN_ID                                 &  !
                                       ,NUM_CHILDREN(MY_DOMAIN_ID)                   &  !
                                       ,COMM_TO_MY_CHILDREN                          &  !
                                       ,COMM_TO_MY_PARENT                            &  !
                                       ,COMM_MY_DOMAIN                               &  !
                                       ,DT                                           &  !
                                       ,CHILD_ID                                     &  !     ^
                                       ,nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) &  !     |
                                       ,FTASKS_DOMAIN                                &  !     |
                                       ,ID_PARENTS                                   &  !     |
                                       ,DOMAIN_ID_TO_RANK                            &  !     |
                                       ,MAX_DOMAINS                                  &  !   Input
!                                                                                         ----------
                                       ,IMP_STATE_CPL_NEST                           &  !   Output
                                       ,EXP_STATE_CPL_NEST                           &  !     |
                                       ,PARENT_CHILD_COUPLER_COMP )                     !     v
!
!-----------------------------------------------------------------------
!***  Now we can initialize the Parent_Child coupler subcomponent.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_FCST_TASK==ESMF_TRUE)THEN                                !<-- Only forecast tasks are relevant
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Initialize Parent-Child Coupler"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompInitialize(cplcomp    =PARENT_CHILD_COUPLER_COMP   &  !<-- The dyn_phy coupler component
                                     ,importState=IMP_STATE_CPL_NEST          &  !<-- The dyn-phy coupler import state
                                     ,exportState=EXP_STATE_CPL_NEST          &  !<-- The dyn-phy coupler export state
                                     ,clock      =CLOCK_NMM(MY_DOMAIN_ID)     &  !<-- The DOMAIN Clock
                                     ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Create the Alarm for children to recv from their parents.
!
!***  The uppermost domain has no parent therefore it must never
!***  try to recv from one.
!
!***  The data exchanged between parents and children consists
!***  entirely of prognostic forecast values therefore the
!***  Quilt/Write tasks must never try to recv from the parent.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Alarm for Child to Recv from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(MY_DOMAIN_ID==1.OR.I_AM_A_FCST_TASK==ESMF_FALSE)THEN                !<-- Uppermost domain and quilt/write tasks do not
          TIMEINTERVAL_RECV_FROM_PARENT=TIMESTEP(MY_DOMAIN_ID)*1000000         !     recv from parents
        ELSE
          TIMEINTERVAL_RECV_FROM_PARENT=TIMESTEP(MY_DOMAIN_ID)              &  !<-- Children recv at the end of each parent timestep
                                        *NINT(DT(ID_PARENTS(MY_DOMAIN_ID))  &
                                            /DT(MY_DOMAIN_ID))
        ENDIF
!
        ALARM_RECV_FROM_PARENT=ESMF_AlarmCreate(                         &
                        name             ='ALARM Recv from Parent'       &  !<-- Name of Alarm
                       ,clock            =CLOCK_NMM(MY_DOMAIN_ID)        &  !<-- Each domain's ESMF Clock
                       ,ringTime         =STARTTIME                      &  !<-- First time the Alarm rings (ESMF)
                       ,ringInterval     =TIMEINTERVAL_RECV_FROM_PARENT  &  !<-- Recv from my parent at this frequency (ESMF)
                       ,ringTimeStepCount=1                              &  !<-- The Alarm rings for this many timesteps
                       ,sticky           =.false.                        &  !<-- Alarm does not ring until turned off
                       ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF nesting_block_3
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract ID of the task that will print clocktimes on this domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read MPI Task ID That Provides Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NPE_PRINT                     &  !<-- Fill this variable (this task prints its clocktimes)
                                  ,label ='npe_print:'                  &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the alarm for printing clocktimes used by model routines.
!***  Read in forecast time interval for clocktime output as well as
!***  the selected task ID that will provide the clocktimes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NHOURS_CLOCKTIME              &  !<-- Fill this variable (fcst hrs between clocktime prints)
                                  ,label ='nhours_clocktime:'           &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create ESMF Clocktime Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_CLOCKTIME         &  !<-- Time interval between clocktime writes (h) (ESMF)
                               ,h           =NHOURS_CLOCKTIME           &  !<-- Hours between clocktime writes (INTEGER)
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Clocktime Output Alarm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALARM_CLOCKTIME=ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'           &
                                      ,clock            =CLOCK_NMM(MY_DOMAIN_ID)     &  !<-- DOMAIN Clock
                                      ,ringInterval     =INTERVAL_CLOCKTIME          &  !<-- Time interval between clocktime prints (ESMF)
                                      ,ringTimeStepCount=1                           &  !<-- The Alarm rings for this many timesteps
                                      ,sticky           =.false.                     &  !<-- Alarm does not ring until turned off
                                      ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' NMM_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
      write(0,*) '    END OF NMM_INITIALIZE'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_RUN(NMM_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_NEMS                                     &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  This routine executes the integration timeloop for the NMM
!***  through a call to subroutine NMM_INTEGRATE.
!***  That is preceded by digital filtering if it is requested.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP                   !<-- The NMM component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The NMM import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                       !<-- The NMM export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NEMS                      !<-- The NEMS ESMF Clock
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: FILTER_METHOD,HDIFF_ON,MYPE_LOCAL,NTIMESTEP
!
      INTEGER(kind=KINT) :: RC
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      TYPE(ESMF_Time) :: CURRTIME
!
      TYPE(ESMF_State) :: IMP_STATE_DOMAIN                              &
                         ,EXP_STATE_DOMAIN
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP
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
!
      DOMAIN_GRID_COMP=nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID)
      IMP_STATE_DOMAIN=nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
      EXP_STATE_DOMAIN=nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!***  Obtain current information from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Run: Get current time info from the Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_NMM(MY_DOMAIN_ID)           &
                        ,starttime   =STARTTIME                         &
                        ,currtime    =CURRTIME                          &
                        ,advanceCount=NTIMESTEP_ESMF                    &
                        ,runduration =RUNDURATION                       &
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***  We need the local MPI task ID on the given NMM domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN gridded component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the DOMAIN component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine for this DOMAIN component
                     ,localpet       =MYPE_LOCAL                        &  !<-- Each MPI task ID
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Do we want to run a digital filter?  Extract the flag.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Driver_Run: Get Filter Method from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &
                                  ,value =FILTER_METHOD                 &
                                  ,label ='filter_method:'              &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Run: Put Filter Method into DOMAIN import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                     &  !<-- This DOMAIN component's import state for filter
                            ,name ='Filter_Method'                      &  !<-- Flag for type of digital filter
                            ,value=FILTER_METHOD                        &  !<-- Value of digital filter flag
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set default value for horizontal diffusion flag (1-->ON).
!-----------------------------------------------------------------------
!
      HDIFF_ON=1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Run: Put Horizontal Diffusion Flag into DOMAIN import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                     &  !<-- This DOMAIN component's import state for horiz diff
                            ,name ='HDIFF'                              &  !<-- Flag for diffusion on/off
                            ,value=HDIFF_ON                             &  !<-- Value of horizontal diffusion flag
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If the user has requested digital filtering then proceed with the
!***  selected method before performing the normal forecast integration.
!-----------------------------------------------------------------------
!
      IF(FILTER_METHOD>0)THEN
!
        CALL RUN_DIGITAL_FILTER_NMM                                        !<-- See internal subroutine below.
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Execute the normal forecast integration.
!-----------------------------------------------------------------------
!
      CALL NMM_INTEGRATE(domain_grid_comp  =DOMAIN_GRID_COMP               &
                        ,imp_state_domain  =IMP_STATE_DOMAIN               &
                        ,exp_state_domain  =EXP_STATE_DOMAIN               &
                        ,clock_integrate   =CLOCK_NMM(MY_DOMAIN_ID)        &
                        ,clock_direction   ='Forward'                      &
                        ,currtime          =CURRTIME                       &
                        ,starttime         =STARTTIME                      &
                        ,timestep          =TIMESTEP(MY_DOMAIN_ID)         &
                        ,ntimestep         =NTIMESTEP                      &
                        ,dt                =DT(MY_DOMAIN_ID)               &
                        ,interval_clocktime=INTERVAL_CLOCKTIME             &
                        ,interval_history  =INTERVAL_HISTORY(MY_DOMAIN_ID) &
                        ,interval_restart  =INTERVAL_RESTART(MY_DOMAIN_ID) &
                        ,filter_method     =FILTER_METHOD                  &
                        ,npe_print         =NPE_PRINT                      &
                        ,restarted_run     =RESTARTED_RUN                  &
                        ,rst_out_00        =RST_OUT_00                     &
                        ,i_am_a_fcst_task  =I_AM_A_FCST_TASK               &
                        ,nesting           =NESTING_NMM                    &
                        ,i_am_a_nest       =I_AM_A_NEST                    &
                        ,my_domain_id      =MY_DOMAIN_ID                   &
                        ,comm_to_my_parent =COMM_TO_MY_PARENT              &
                        ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)     &
                        ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP      &
                        ,imp_state_cpl_nest=IMP_STATE_CPL_NEST             &
                        ,exp_state_cpl_nest=EXP_STATE_CPL_NEST             &
                        ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO        &
                        ,mype              =MYPE_LOCAL)
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_RUN succeeded'
      ELSE
        WRITE(0,*)' NMM_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!
      write(0,*) '    END OF NMM_RUN'
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      SUBROUTINE RUN_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE EXECUTES THE DIGITAL FILTERS FOR THE NMM
!***  IF SPECIFIED BY THE USER.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: FILTER_METHOD,HDIFF_ON,MEAN_ON              &
                           ,NDFISTEP,NTIMESTEP
!
      INTEGER(kind=KINT) :: RC
!
      TYPE(ESMF_Clock) :: CLOCK_FILTER
!
      TYPE(ESMF_Time) :: DFITIME,HALFDFITIME,SDFITIME
!
      TYPE(ESMF_TimeInterval) :: HALFDFIINTVAL                          &
                                ,TIMESTEP_FILTER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      method_block: IF(FILTER_METHOD==1)THEN                               !<-- The DFL digital filter.
!
!-----------------------------------------------------------------------
!***  Create a Clock to control the filter's timestepping and then
!***  execute this filter's forward integration.
!-----------------------------------------------------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- This domain's config object
                                    ,value =DFIHR                       &  !<-- The digital filter flag
                                    ,label ='nsecs_dfl:'                &  !<-- Give this label's value to preceding variable
                                    ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,s           =DFIHR                    &
                                 ,rc          =RC)
!
        HALFDFITIME=CURRTIME+HALFDFIINTVAL
        SDFITIME=CURRTIME
        DFITIME=HALFDFITIME+HALFDFIINTVAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clock for the DFL Digital Filter."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        TIMESTEP_FILTER=TIMESTEP(MY_DOMAIN_ID)
!
        CLOCK_FILTER=ESMF_ClockCreate(name     ='CLOCK_DFL'             &  !<-- The Clock for the DFI filter
                                     ,timeStep =TIMESTEP_FILTER         &  !<-- The fundamental timestep in this component
                                     ,startTime=STARTTIME               &  !<-- Start time of filter
                                     ,stopTime =DFITIME                 &  !<-- Stop time of the filter
                                     ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1                                                         !<-- Forward integration so we want horiz diffusion.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- This filter only integrates forward
                          ,domain_grid_comp  =DOMAIN_GRID_COMP          &
                          ,imp_state_domain  =IMP_STATE_DOMAIN          &
                          ,exp_state_domain  =EXP_STATE_DOMAIN          &
                          ,clock_integrate   =CLOCK_FILTER              &
                          ,currtime          =CURRTIME                  &
                          ,starttime         =STARTTIME                 &
                          ,timestep          =TIMESTEP(MY_DOMAIN_ID)    &
                          ,ntimestep         =NTIMESTEP                 &
                          ,dt                =DT(MY_DOMAIN_ID)          &
                          ,filter_method     =FILTER_METHOD             &
                          ,halfdfiintval     =HALFDFIINTVAL             &
                          ,halfdfitime       =HALFDFITIME               &
                          ,npe_print         =NPE_PRINT                 &
                          ,restarted_run     =RESTARTED_RUN             &
                          ,rst_out_00        =RST_OUT_00                &
                          ,i_am_a_fcst_task  =I_AM_A_FCST_TASK          &
                          ,nesting           =NESTING_NMM               &
                          ,i_am_a_nest       =I_AM_A_NEST               &
                          ,my_domain_id      =MY_DOMAIN_ID              &
                          ,comm_to_my_parent =COMM_TO_MY_PARENT         &
                          ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)&
                          ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP &
                          ,imp_state_cpl_nest=IMP_STATE_CPL_NEST        &
                          ,exp_state_cpl_nest=EXP_STATE_CPL_NEST        &
                          ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO   &
                          ,mype              =MYPE)
!
        STARTTIME=CURRTIME                                                 !<-- Start time set to halfway point of filter period
!
        CALL ESMF_ClockSet(clock    =CLOCK_NMM(MY_DOMAIN_ID)            &  !<-- For DFL filter, the starttime of the free forecast
                          ,starttime=STARTTIME                          &  !    moves ahead to the halfway point of the filter
                          ,currtime =CURRTIME                           &  !    interval.
                          ,rc       =RC)
!
        NTIMESTEP_ESMF=NTIMESTEP
!
!-----------------------------------------------------------------------
!
      ELSEIF(FILTER_METHOD==2)THEN  method_block                           !<-- The DDFI digital filter.
!
!-----------------------------------------------------------------------
!
!--------------------------------
!***  The initial backward step.
!--------------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                    ,value = DFIHR                      &  !<-- The digital filter flag
                                    ,label ='nsecs_bckddfi:'            &  !<-- Time duration of this backward part of filter
                                    ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,s           =DFIHR                    &
                                 ,rc          =RC)
!
        HALFDFITIME=STARTTIME-HALFDFIINTVAL
        DFITIME=HALFDFITIME-HALFDFIINTVAL
!
        TIMESTEP_FILTER=-TIMESTEP(MY_DOMAIN_ID)                            !<-- Prepare for backward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clock for the DDFI Digital Filter."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CLOCK_FILTER=ESMF_ClockCreate(name     ='CLOCK_DDFI'            &  !<-- The Clock for the DFI filter
                                     ,timeStep =TIMESTEP_FILTER         &  !<-- The fundamental timestep in this component
                                     ,startTime=STARTTIME               &  !<-- Start time of filter
!!!!!!!!                             ,direction=ESMF_MODE_REVERSE       &  !<-- Reverse the Clock for backward integration
                                     ,stopTime =DFITIME                 &  !<-- Stop time of the filter
                                     ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=0                                                         !<-- Turn off horiz diffusion for backward integration.
        MEAN_ON =1                                                         !<-- Turn off horiz diffusion for backward integration.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Bckward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Bckward'                 &  !<-- The initial backward piece of the filter
                          ,domain_grid_comp  =DOMAIN_GRID_COMP          &
                          ,imp_state_domain  =IMP_STATE_DOMAIN          &
                          ,exp_state_domain  =EXP_STATE_DOMAIN          &
                          ,clock_integrate   =CLOCK_FILTER              &
                          ,currtime          =CURRTIME                  &
                          ,starttime         =STARTTIME                 &
                          ,timestep          =TIMESTEP(MY_DOMAIN_ID)    &
                          ,ntimestep         =NTIMESTEP                 &
                          ,dt                =DT(MY_DOMAIN_ID)          &
                          ,filter_method     =FILTER_METHOD             &
                          ,ndfistep          =NDFISTEP                  &
                          ,npe_print         =NPE_PRINT                 &
                          ,restarted_run     =RESTARTED_RUN             &
                          ,rst_out_00        =RST_OUT_00                &
                          ,i_am_a_fcst_task  =I_AM_A_FCST_TASK          &
                          ,nesting           =NESTING_NMM               &
                          ,i_am_a_nest       =I_AM_A_NEST               &
                          ,my_domain_id      =MY_DOMAIN_ID              &
                          ,comm_to_my_parent =COMM_TO_MY_PARENT         &
                          ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)&
                          ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP &
                          ,imp_state_cpl_nest=IMP_STATE_CPL_NEST        &
                          ,exp_state_cpl_nest=EXP_STATE_CPL_NEST        &
                          ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO   &
                          ,mype              =MYPE)
!
!-----------------------------
!***  The final forward step.
!-----------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object for this domain
                                    ,value = DFIHR                      &  !<-- The digital filter flag
                                    ,label ='nsecs_fwdddfi:'            &  !<-- Time duration of this forward integration
                                    ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,s           =DFIHR                    &
                                 ,rc          =RC)
!
        HALFDFITIME=CURRTIME+HALFDFIINTVAL
        SDFITIME=CURRTIME
        DFITIME=HALFDFITIME+HALFDFIINTVAL
!
        TIMESTEP_FILTER=TIMESTEP(MY_DOMAIN_ID)                            !<-- Prepare for forward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset the Clock for Forward DDFI Digital Filter."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockSet(clock   =CLOCK_FILTER                        &  !<-- Reset the stoptime for the forward part of the filter
                          ,timeStep=TIMESTEP_FILTER                     &  !<-- The fundamental timestep in this component
                          ,stoptime=DFITIME                             &
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1
!
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- The final forward piece of the filter
                          ,domain_grid_comp  =DOMAIN_GRID_COMP          &
                          ,imp_state_domain  =IMP_STATE_DOMAIN          &
                          ,exp_state_domain  =EXP_STATE_DOMAIN          &
                          ,clock_integrate   =CLOCK_FILTER              &
                          ,currtime          =CURRTIME                  &
!!!                       ,starttime         =STARTTIME                 &
                          ,starttime         =CURRTIME                  &  !<-- CURRTIME was set or reset at end of backwward piece
                          ,timestep          =TIMESTEP(MY_DOMAIN_ID)    &
                          ,ntimestep         =NTIMESTEP                 &
                          ,dt                =DT(MY_DOMAIN_ID)          &
                          ,filter_method     =FILTER_METHOD             &
                          ,halfdfiintval     =HALFDFIINTVAL             &
                          ,halfdfitime       =HALFDFITIME               &
                          ,npe_print         =NPE_PRINT                 &
                          ,restarted_run     =RESTARTED_RUN             &
                          ,rst_out_00        =RST_OUT_00                &
                          ,i_am_a_fcst_task  =I_AM_A_FCST_TASK          &
                          ,nesting           =NESTING_NMM               &
                          ,i_am_a_nest       =I_AM_A_NEST               &
                          ,my_domain_id      =MY_DOMAIN_ID              &
                          ,comm_to_my_parent =COMM_TO_MY_PARENT         &
                          ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)&
                          ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP &
                          ,imp_state_cpl_nest=IMP_STATE_CPL_NEST        &
                          ,exp_state_cpl_nest=EXP_STATE_CPL_NEST        &
                          ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO   &
                          ,mype              =MYPE)
!
!-----------------------------------------------------------------------
!
      ELSEIF(FILTER_METHOD==3)THEN  method_block                           !<-- The TDFI digital filter.
!
!-----------------------------------------------------------------------
!
!--------------------------------
!***  The initial backward step.
!--------------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                    ,value = DFIHR                      &  !<-- The digital filter flag
                                    ,label ='nsecs_bcktdfi:'            &  !<-- Give this label's value to preceding variable
                                    ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,s           =DFIHR                    &
                                 ,rc          =RC)
!
        HALFDFITIME=STARTTIME-HALFDFIINTVAL
        DFITIME=HALFDFITIME-HALFDFIINTVAL
!
        TIMESTEP_FILTER=-TIMESTEP(MY_DOMAIN_ID)                            !<-- Prepare for inital backward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clock for the TDFI Digital Filter."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CLOCK_FILTER=ESMF_ClockCreate(name     ='CLOCK_TDFI'            &  !<-- The Clock for the DFI filter
                                     ,timeStep =TIMESTEP(MY_DOMAIN_ID)  &  !<-- The fundamental timestep in this component
                                     ,startTime=STARTTIME               &  !<-- Start time of filter
!!!!!!!!                             ,direction=ESMF_MODE_REVERSE       &  !<-- Reverse the Clock for backward integration
                                     ,stopTime =DFITIME                 &  !<-- Stop time of the filter
                                     ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=0                                                         !<-- Turn off horiz diffusion for backward integration.
        MEAN_ON =0                                                         !<-- Turn off horiz diffusion for backward integration.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Bckward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Bckward'                 &  !<-- The initial backward piece of the filter
                          ,domain_grid_comp  =DOMAIN_GRID_COMP          &
                          ,imp_state_domain  =IMP_STATE_DOMAIN          &
                          ,exp_state_domain  =EXP_STATE_DOMAIN          &
                          ,clock_integrate   =CLOCK_FILTER              &
                          ,currtime          =CURRTIME                  &
                          ,starttime         =STARTTIME                 &
                          ,timestep          =TIMESTEP(MY_DOMAIN_ID)    &
                          ,ntimestep         =NTIMESTEP                 &
                          ,dt                =DT(MY_DOMAIN_ID)          &
                          ,filter_method     =FILTER_METHOD             &
                          ,ndfistep          =NDFISTEP                  &
                          ,npe_print         =NPE_PRINT                 &
                          ,restarted_run     =RESTARTED_RUN             &
                          ,rst_out_00        =RST_OUT_00                &
                          ,i_am_a_fcst_task  =I_AM_A_FCST_TASK          &
                          ,nesting           =NESTING_NMM               &
                          ,i_am_a_nest       =I_AM_A_NEST               &
                          ,my_domain_id      =MY_DOMAIN_ID              &
                          ,comm_to_my_parent =COMM_TO_MY_PARENT         &
                          ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)&
                          ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP &
                          ,imp_state_cpl_nest=IMP_STATE_CPL_NEST        &
                          ,exp_state_cpl_nest=EXP_STATE_CPL_NEST        &
                          ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO   &
                          ,mype              =MYPE)
!
!-----------------------------
!***  The final forward step.
!-----------------------------
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                    ,value =DFIHR                       &  !<-- The digital filter flag
                                    ,label ='nsecs_fwdtdfi:'            &  !<-- Give this label's value to preceding variable
                                    ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,s           =DFIHR                    &
                                 ,rc          =RC)
!
        HALFDFITIME=CURRTIME+HALFDFIINTVAL
        SDFITIME=CURRTIME
        DFITIME=HALFDFITIME+HALFDFIINTVAL
!
        TIMESTEP_FILTER=TIMESTEP(MY_DOMAIN_ID)                            !<-- Prepare for forward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset the Clock for Forward TDFI Digital Filter."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockSet(clock   =CLOCK_FILTER                        &  !<-- Reset the stoptime for the forward part of the filter
                          ,timeStep=TIMESTEP_FILTER                     &  !<-- The fundamental timestep in this component
                          ,stoptime=DFITIME                             &
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1                                                         !<-- Forward integration so we want horiz diffusion.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- The final forward piece of the filter
                          ,domain_grid_comp  =DOMAIN_GRID_COMP          &
                          ,imp_state_domain  =IMP_STATE_DOMAIN          &
                          ,exp_state_domain  =EXP_STATE_DOMAIN          &
                          ,clock_integrate   =CLOCK_FILTER              &
                          ,currtime          =CURRTIME                  &
!!!                       ,starttime         =STARTTIME                 &
                          ,starttime         =CURRTIME                  &  !<-- CURRTIME was set or reset at end of backwward piece
                          ,timestep          =TIMESTEP(MY_DOMAIN_ID)    &
                          ,ntimestep         =NTIMESTEP                 &
                          ,dt                =DT(MY_DOMAIN_ID)          &
                          ,filter_method     =FILTER_METHOD             &
                          ,halfdfiintval     =HALFDFIINTVAL             &
                          ,halfdfitime       =HALFDFITIME               &
                          ,npe_print         =NPE_PRINT                 &
                          ,restarted_run     =RESTARTED_RUN             &
                          ,rst_out_00        =RST_OUT_00                &
                          ,i_am_a_fcst_task  =I_AM_A_FCST_TASK          &
                          ,nesting           =NESTING_NMM               &
                          ,i_am_a_nest       =I_AM_A_NEST               &
                          ,my_domain_id      =MY_DOMAIN_ID              &
                          ,comm_to_my_parent =COMM_TO_MY_PARENT         &
                          ,num_children      =NUM_CHILDREN(MY_DOMAIN_ID)&
                          ,parent_child_cpl  =PARENT_CHILD_COUPLER_COMP &
                          ,imp_state_cpl_nest=IMP_STATE_CPL_NEST        &
                          ,exp_state_cpl_nest=EXP_STATE_CPL_NEST        &
                          ,par_chi_time_ratio=PARENT_CHILD_TIME_RATIO   &
                          ,mype              =MYPE)
!
!-----------------------------------------------------------------------
!
      ENDIF  method_block
!
!-----------------------------------------------------------------------
!***  We do not need the local filter clock any more.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Local Filter Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_FILTER                         &
                             ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Filtering is finished.  Reset the filter flag for the normal
!***  integration that follows.
!-----------------------------------------------------------------------
!
      FILTER_METHOD=0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RUN_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_FINALIZE(NMM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_NMM                                 &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP                   !<-- The NMM component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The NMM import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                       !<-- The NMM export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NMM                       !<-- The NMM component's ESMF Clock
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE                     !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,N
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      write(0,*) "        NMM_FINALIZE"
      RC         =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' NMM_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
      write(0,*) '    END OF NMM_FINALIZE'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      RECURSIVE SUBROUTINE CALL_DOMAIN_INITIALIZE(ID_DOMAIN,CLOCK_NMM)
!
!-----------------------------------------------------------------------
!***  This routine calls DOMAIN_INITIALIZE for all DOMAIN components.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: ID_DOMAIN                                      !<-- ID of the DOMAIN Component to initialize
!
      TYPE(ESMF_Clock),DIMENSION(1:NUM_DOMAINS),INTENT(INOUT) :: CLOCK_NMM !<-- The NMM ESMF Clock
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      INTEGER :: ID_CHILD,IRTN,N,N_CHILDREN
      INTEGER :: RC,RC_CALL_INIT
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_CALL_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      FMT='(I2.2)'
      WRITE(INT_TO_CHAR,FMT)ID_DOMAIN
!
!-----------------------------------------------------------------------
!***  Initialize the DOMAIN component with the ID of ID_DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize DOMAIN Component "//INT_TO_CHAR
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =nmm_int_state%DOMAIN_GRID_COMP(ID_DOMAIN)  &  !<-- The DOMAIN component
                                  ,importState=nmm_int_state%IMP_STATE_DOMAIN(ID_DOMAIN)  &  !<-- The DOMAIN import state
                                  ,exportState=nmm_int_state%EXP_STATE_DOMAIN(ID_DOMAIN)  &  !<-- The DOMAIN export state
                                  ,clock      =CLOCK_NMM(ID_DOMAIN)                       &  !<-- The DOMAIN clock
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CALL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If the domain being initialized has children, force those 
!***  children to wait.  The parent might be generating input for
!***  the children therefore children should not be trying to read
!***  the input in their own initialize steps prematurely.
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(COMM_FULL_DOMAIN,IRTN)
!
!-----------------------------------------------------------------------
!***  If there are children, initialize them.
!-----------------------------------------------------------------------
!
      N_CHILDREN=NUM_CHILDREN(ID_DOMAIN)
!
      IF(N_CHILDREN>0)THEN                                                 !<-- Does the current DOMAIN have any children?
        DO N=1,N_CHILDREN                                                  !<-- If so, loop through the children to Initialize them
          ID_CHILD=ID_CHILDREN(N,ID_DOMAIN)
!         write(0,*)' recursive call  n_children=',n_children,' id_atm=',id_atm,' id_child=',id_child
          CALL CALL_DOMAIN_INITIALIZE(ID_CHILD,CLOCK_NMM)
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CALL_DOMAIN_INITIALIZE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE BOUNDARY_DATA_TO_DOMAIN(EXP_STATE_CPL                  &
                                        ,IMP_STATE_DOMAIN )
!
!-----------------------------------------------------------------------
!***  This routine moves new boundary data for nested domains from the
!***  export state of the Parent-Child coupler to the import state of
!***  the NMM nests' DOMAIN components.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State),INTENT(IN)  :: EXP_STATE_CPL                        !<-- The Parent-Child Coupler's export state
!
      TYPE(ESMF_State),INTENT(OUT) :: IMP_STATE_DOMAIN                     !<-- The nests' DOMAIN import state
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      TYPE SIDES_1D_REAL
        REAL,DIMENSION(:),ALLOCATABLE :: SOUTH
        REAL,DIMENSION(:),ALLOCATABLE :: NORTH
        REAL,DIMENSION(:),ALLOCATABLE :: WEST
        REAL,DIMENSION(:),ALLOCATABLE :: EAST
      END TYPE SIDES_1D_REAL
!
      INTEGER :: ISTAT,KOUNT,RC,RC_BND_MV
!
      TYPE(SIDES_1D_REAL),SAVE :: BOUNDARY_H                            &
                                 ,BOUNDARY_V
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Check each side of the child boundary.  If data is present from
!***  that side in the Parent-Child coupler export state then move it
!***  to the DOMAIN component's import state.
!-----------------------------------------------------------------------
!
!-------------
!***  South H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for South H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='SOUTH_H'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%SOUTH))THEN
          ALLOCATE(BOUNDARY_H%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_h
!
!-------------
!***  South V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for South V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='SOUTH_V'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%SOUTH))THEN
          ALLOCATE(BOUNDARY_V%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_v
!
!-------------
!***  North H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for North H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='NORTH_H'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%NORTH))THEN
          ALLOCATE(BOUNDARY_H%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_h
!
!-------------
!***  North V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for North V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='NORTH_V'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%NORTH))THEN
          ALLOCATE(BOUNDARY_V%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_v
!
!------------
!***  West H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for West H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='WEST_H'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%WEST))THEN
          ALLOCATE(BOUNDARY_H%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='WEST_H'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='WEST_H'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_h
!
!------------
!***  West V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for West V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='WEST_V'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%WEST))THEN
          ALLOCATE(BOUNDARY_V%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='WEST_V'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='WEST_V'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_v
!
!------------
!***  East H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for East H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='EAST_H'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%EAST))THEN
          ALLOCATE(BOUNDARY_H%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='EAST_H'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='EAST_H'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_h
!
!------------
!***  East V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for East V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_CPL                        &   !<-- Look at the Parent-Child Coupler's export state
                            ,name ='EAST_V'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%EAST))THEN
          ALLOCATE(BOUNDARY_V%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='EAST_V'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='EAST_V'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_v
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOUNDARY_DATA_TO_DOMAIN
!
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
