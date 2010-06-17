!-----------------------------------------------------------------------
!
      MODULE module_DOMAIN_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  This is the DOMAIN gridded component module.
!***  It will set up Dynamics, Physics, and coupler subcomponents
!***  and run their Initialize, Run, and Finalize routines.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black - Modified from Wei-yu's version
!   2007-11-20  Black/Henderson - Created an ESMF Clock for the
!                                 ATM Component independent of
!                                 the Main Clock.
!   2007-12-11  Black - Generalized for easier use by any dynamics core.
!   2008-08     Colon - Added conditional checks multiple dynamics cores.
!   2008-10-14  Vasic - Added restart Alarm.
!   2009-05-29  Wang  - Added GFS write grid component
!   2009-07-23  Lu    - Added GOCART grid component
!   2009-08-03  Black - Merging with nesting.
!   2009-08-12  Black - Fixed logic for Physics export when direction of
!                       integration switches from backward to forward.
!   2009-10-05  Wang  - Added GFS ensemble member name and output data at
!                       every nsout timesteps.
!   2009-10-08  W Yang - Ensemble GEFS.
!   2009-11-03  Lu    - Add GOCART and ChemRegistry modules
!   2009-12-17  Lu    - Modify GFS_ATM_INIT routine to create, register,
!                       and initialize GOCART grid component
!   2009-12-22  Lu    - Modify GFS_ATM_INIT routine to create and register
!                       dyn2chem and phy2chem coupler components
!   2009-12-23  Lu    - Modify GFS_INTEGRATE routine to loop thru dyn, phy,
!   2010-02-01  Lu    - Remove dyn2chem coupler component
!   2010-02-05  Wang  - Added restart file for GFS
!   2010-03-04  Lu    - Modify GFS_ATM_INIT (initialization is changed from
!                       DYN-PHY-CPL to DYN-CPL-PHY)
!   2010-03-05  Lu    - Add GOCART_SETUP (to create and register GOCART) and
!                       GOCART_INIT (to initialize GOCART)
!   2010-03-24  Black - Converted to DOMAIN component for NMM-B only.

!
! USAGE: Domain gridded component parts called from subroutines within
!        module_NMM_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_INCLUDE
!
      USE MODULE_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE
!
      USE MODULE_DYNAMICS_GRID_COMP,ONLY: DYN_REGISTER
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &
                                   ,MPI_COMM_COMP                       &
                                   ,MPI_COMM_INTER_ARRAY
!
      USE MODULE_PHYSICS_GRID_COMP,ONLY: PHY_REGISTER
!
      USE MODULE_DYN_PHY_CPL_COMP,ONLY: DYN_PHY_CPL_REGISTER 
!
      USE MODULE_GET_CONFIG_DYN
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_GET_CONFIG_WRITE
!
      USE MODULE_CONTROL,ONLY: TIMEF
      USE MODULE_DIAGNOSE,ONLY: FIELD_STATS
      USE MODULE_NEMSIO
! 
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
      USE MODULE_CLOCKTIMES,ONLY: total_integ_tim
!
      USE MODULE_NMM_CORE_SETUP,ONLY: NMM_SETUP
!
      USE MODULE_NESTING,ONLY: PARENT_DATA_TO_DOMAIN                    &
                              ,PARENT_TO_CHILD_INIT_NMM
!
!-----------------------------------------------------------------------
!***  List other modules with non-generic routines used by DOMAIN.
!-----------------------------------------------------------------------
!
      USE MODULE_WRITE_ROUTINES ,ONLY: WRITE_INIT,WRITE_ASYNC             !<-- These are routines used only when asynchronous
      USE MODULE_WRITE_GRID_COMP,ONLY: WRITE_SETUP                      & !    quilting is specified by the user in the
                                      ,WRITE_DESTROY                      !    configure file for history output.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: DOMAIN_REGISTER 
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: MYPE                                        &  !<-- Each MPI task ID
                           ,NPE_PRINT                                   &
                           ,NUM_TRACERS_CHEM                            &  !<-- Number of chemistry tracer variables
                           ,NUM_TRACERS_MET                             &  !<-- Number of meteorological tracer variables
                           ,WRITE_GROUP_READY_TO_GO                        !<-- The write group to use
      INTEGER(kind=KINT) :: MYPE_GLOBAL                                    !<-- Each MPI task ID
!
      LOGICAL(kind=KLOG) :: QUILTING                                    &  !<-- Is asynchronous quilting specified?
                           ,RESTARTED_RUN                                  !<-- Restarted run logical flag
!
      TYPE(ESMF_Config),SAVE :: CF_1                                       !<-- The principal config object
!
      TYPE(ESMF_VM),SAVE :: VM,VM_LOCAL                                    !<-- The ESMF virtual machine.
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER,SAVE :: DOMAIN_INT_STATE         !<-- The NMM DOMAIN internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE)   ,SAVE :: WRAP                     !<-- The F90 wrap of the NMM DOMAIN internal state
!
      TYPE(ESMF_Time),SAVE :: DFITIME                                   &
                             ,HALFDFITIME 
!
      TYPE(ESMF_TimeInterval),SAVE :: TIMEINTERVAL_CLOCKTIME               !<-- The ESMF time interval between NMM clocktime output
!
      TYPE(ESMF_Logical),SAVE :: PHYSICS_ON                                !<-- Is physics active?
!
!-----------------------------------------------------------------------
!
!---------------------
!***  For NMM Nesting
!---------------------
!
      INTEGER(kind=KINT),SAVE :: MY_DOMAIN_ID                              !<-- Domain IDs; begin with uppermost parent=1
!
      INTEGER(kind=KINT),SAVE :: COMM_FULL_DOMAIN                       &  !<-- Communicator for ALL tasks on domain to be split
                                ,COMM_MY_DOMAIN                         &  !<-- Each domain's local intracommunicator
                                ,NUM_CHILDREN                           &  !<-- Number of (1st generation) children within a domain
                                ,PARENT_CHILD_TIME_RATIO                   !<-- # of child timesteps per parent timestep
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,SAVE :: MY_CHILDREN_ID       !<-- A parent's children's domain IDs
!
      TYPE(ESMF_Logical),SAVE :: I_AM_A_NEST                            &  !<-- Is the domain a nest?
                                ,INPUT_READY                               !<-- If a nest, does its input file already exist?
!
!---------------------------------
!***  For determining clocktimes.
!---------------------------------
!
      REAL(kind=KDBL) :: domain_tim,btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_REGISTER(DOMAIN_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the DOMAIN component's Initialize, Run, and Finalize
!***  routines.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- DOMAIN gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Return code for register
!     
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS   ! Error signal variable
      RC_REG=ESMF_SUCCESS   ! Error signal variable
!
!-----------------------------------------------------------------------
!***  Load the principal configure file to obtain the core name.
!***  We need that here because different numbers of Phases are
!***  needed for the execution of the Run steps of the different 
!***  cores.
!-----------------------------------------------------------------------
!
      CF_1=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_Register: Load principal configure file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigLoadFile(config  =CF_1                            &
                              ,filename='configure_file'                &
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)  
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the DOMAIN INITIALIZE subroutine.  Since it is just one 
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- DOMAIN gridded component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type (Initialize)
                                     ,DOMAIN_INITIALIZE                 &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Run step of the DOMAIN component.
!***  The NMM needs three phases of Run.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 1st Entry Point for the DOMAIN Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The DOMAIN component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,DOMAIN_RUN                        &  !<-- The user's subroutine name for primary integration
                                     ,1                                 &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 2nd Entry Point for the DOMAIN Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The DOMAIN component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,NMM_FILTERING                     &  !<-- Routine to govern digital filtering each timestep
                                     ,2                                 &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 3rd Entry Point for the DOMAIN Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The DOMAIN component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,CALL_WRITE_ASYNC                  &  !<-- Routine to call asynchronous output
                                     ,3                                 &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the DOMAIN FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for DOMAIN Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The DOMAIN component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type (Finalize)
                                     ,DOMAIN_FINALIZE                   &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' DOMAIN_REGISTER succeeded'
      ELSE
        WRITE(0,*)' DOMAIN_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_INITIALIZE(DOMAIN_GRID_COMP                     &
                                  ,IMP_STATE                            &
                                  ,EXP_STATE                            &
                                  ,CLOCK_DOMAIN                         &
                                  ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  This routine sets up fundamental aspects of the model run.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                       &  !<-- The DOMAIN component's import state
                                       ,EXP_STATE                          !<-- The DOMAIN component's export state
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The ESMF Clock from the NMM component.
!
      INTEGER,INTENT(OUT) :: RC_INIT                                       !<-- Return code for Initialize step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: CONFIG_ID,ISTAT,MAX_DOMAINS,N,NFCST,NTSD
!
      INTEGER(kind=KINT) :: IYEAR_FCST                                  &  !<-- Current year from restart file
                           ,IMONTH_FCST                                 &  !<-- Current month from restart file
                           ,IDAY_FCST                                   &  !<-- Current day from restart file
                           ,IHOUR_FCST                                  &  !<-- Current hour from restart file
                           ,IMINUTE_FCST                                &  !<-- Current minute from restart file
                           ,ISECOND_FCST                                   !<-- Current second from restart file
!
      INTEGER(kind=KINT) :: NHOURS_CLOCKTIME                               !<-- Hours between clocktime prints
!
      INTEGER(kind=KINT) :: IERR,IRTN,RC          
!
      INTEGER(ESMF_KIND_I8) :: NTSD_START                                  !<-- Timestep count (>0 for restarted runs)
!
      INTEGER(kind=KINT),DIMENSION(7) :: FCSTDATE
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: DOMAIN_ID_TO_RANK     !<-- Associate configure file IDs with domains
!
      REAL(kind=KFPT) :: SECOND_FCST                                       !<-- Current second from restart file
!
      LOGICAL(kind=KLOG) :: CFILE_EXIST                                 &
                           ,INPUT_READY_FLAG                            &
                           ,INPUT_READY_MY_CHILD                        &
                           ,NEMSIO_INPUT                                &
                           ,OPENED
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT='(I2.2)'
      CHARACTER(64) :: RESTART_FILENAME
      CHARACTER(99) :: CONFIG_FILE_NAME
!
      TYPE(ESMF_Config),DIMENSION(99) :: CF                                !<-- The configure objects for all NMM domains
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STARTTIME                                         !<-- The ESMF start time.
!
      TYPE(ESMF_Grid) :: GRID_DOMAIN                                       !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM DOMAIN component.
      TYPE(ESMF_Grid) :: GRID_DYN                                          !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM dynamics gridded component.
      TYPE(ESMF_Grid) :: GRID_PHY                                          !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM physics gridded component.
!
      TYPE(ESMF_Logical) :: I_AM_A_FCST_TASK,I_AM_A_PARENT
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
      real(kind=8) :: rtim0,rtim1,rtc
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize timing variables.
!-----------------------------------------------------------------------
!
      rtim0=rtc()
      btim0=timef()
      total_integ_tim=0.
!
!-----------------------------------------------------------------------
!***  Allocate the DOMAIN component's internal state.
!-----------------------------------------------------------------------
!
      ALLOCATE(DOMAIN_INT_STATE,stat=RC)
!
      wrap%DOMAIN_INT_STATE=>DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Attach the DOMAIN internal state to the DOMAIN component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach DOMAIN Internal State to Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(DOMAIN_GRID_COMP               &  !<-- The DOMAIN gridded component
                                        ,WRAP                           &  !<-- Pointer to the DOMAIN internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the DOMAIN component.
!***  Call ESMF_GridCompGet to retrieve the VM anywhere you need it.
!***  We need VM now to obtain the MPI task IDs and the local MPI
!***  communicator.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Retrieve VM from DOMAIN Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the DOMAIN component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Obtain Task IDs and Communicator"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine
                     ,localpet       =MYPE                              &  !<-- Each MPI task ID
                     ,mpiCommunicator=COMM_MY_DOMAIN                    &  !<-- This domain's communicator
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the maximum number of domains from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract MAX_DOMAINS from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The DOMAIN import state
                            ,name ='MAX_DOMAINS'                        &  !<-- Name of the attribute to extract
                            ,value=MAX_DOMAINS                          &  !<-- Maximum # of domains
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract this DOMAIN component's domain ID from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain ID from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The DOMAIN import state
                            ,name ='DOMAIN_ID'                          &  !<-- Name of the attribute to extract
                            ,value=MY_DOMAIN_ID                         &  !<-- The ID of this domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the configure file IDs assocaited with each domain.
!-----------------------------------------------------------------------
!
      ALLOCATE(DOMAIN_ID_TO_RANK(1:MAX_DOMAINS),stat=ISTAT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Association of Configure Files with Domains"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The DOMAIN import state
                            ,name     ='DOMAIN_ID_TO_RANK'              &  !<-- Name of the attribute to extract
                            ,count    =MAX_DOMAINS                      &  !<-- Name of the attribute to extract
                            ,valueList=DOMAIN_ID_TO_RANK                &  !<-- The ID of this domain
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now we can load configure files for all domains into memory.
!***  The file name of the uppermost domain is 'configure_file_01'
!***  and is identical to the primary file called 'configure_file'
!***  which is needed in some early parts of the setup.
!-----------------------------------------------------------------------
!
      DO N=1,MAX_DOMAINS                                                   !<-- The number of config files cannot exceed 99
!
        CONFIG_ID=DOMAIN_ID_TO_RANK(N)
        WRITE(INT_TO_CHAR,FMT)CONFIG_ID 
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file names
!
        CFILE_EXIST=.FALSE.
        INQUIRE(file=CONFIG_FILE_NAME,exist=CFILE_EXIST)
!
        IF(CFILE_EXIST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Create the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF(N)=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(N)                       &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ELSE
!
          EXIT
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Will the Write components with asynchronous quilting be used?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Quilting Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =QUILTING                      &  !<-- The quilting flag
                                  ,label ='quilting:'                   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%QUILTING=QUILTING                                   !<-- Save this for the Run step
!
!-----------------------------------------------------------------------
!***  Extract this DOMAIN component's Nest/Not-a-Nest flag
!***  from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Nest/Not-a-Nest Flag from DOMAIN Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The DOMAIN import state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_NEST                          &  !<-- The flag indicating if this domain is a nest
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the ratio of the parent timestep to the child's if
!***  this domain is a nest.
!***  Also extract the flag indicating whether or not the nest's
!***  input file has already been generated by NPS.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_NEST==ESMF_TRUE)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-Child Time Ratio from DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The DOMAIN import state
                              ,name ='Parent-Child Time Ratio'          &  !<-- Name of Attribute
                              ,value=PARENT_CHILD_TIME_RATIO            &  !<-- Ratio of this domain's parent's timestep to its own
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Input Ready Flag from Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                    ,value =INPUT_READY_FLAG              &  !<-- The variable filled (does nest input file exist?
                                    ,label ='input_ready:'                &  !<-- The input datafile for this domain does or does not exist
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(INPUT_READY_FLAG)THEN
          INPUT_READY=ESMF_TRUE
        ELSE
          INPUT_READY=ESMF_FALSE
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Extract the start time from the clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Start Time from Driver Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_DOMAIN                         &  !<-- The ESMF Clock of this domain
                        ,startTime=STARTTIME                            &  !<-- The simulation start time
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CURRTIME=STARTTIME
      NTSD_START=0
!
!-----------------------------------------------------------------------
!***  Extract the NEMSIO_INPUT flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NEMSIO Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                  ,value =NEMSIO_INPUT                &  !<-- The input datafile does or does not have NEMSIO metadata
                                  ,label ='nemsio_input:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the RESTART flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =RESTARTED_RUN                 &  !<-- True => restart; False => cold start
                                  ,label ='restart:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If this is a restarted run then read:
!***    (1) The forecast time that the file was written.
!***    (2) The forecast timestep at which the file was written.
!-----------------------------------------------------------------------
!
      NTSD_START=0
!
      restart: IF(RESTARTED_RUN)THEN                                       !<-- If this is a restarted run, set the current time
!
        WRITE(INT_TO_CHAR,FMT)MY_DOMAIN_ID
!
!----------------------------------------------------------------------
!***  Read the restart data from either pure binary or NEMSIO file.
!-----------------------------------------------------------------------
!
        input: IF(NEMSIO_INPUT)THEN
!
          CALL NEMSIO_INIT()
!
          RESTART_FILENAME='restart_file_'//INT_TO_CHAR//'_nemsio'
          CALL NEMSIO_OPEN(GFILE,RESTART_FILENAME,'read',iret=IRTN)
!
          CALL NEMSIO_GETHEADVAR(GFILE,'FCSTDATE',FCSTDATE,iret=irtn)
!
          IYEAR_FCST  =FCSTDATE(1)
          IMONTH_FCST =FCSTDATE(2)
          IDAY_FCST   =FCSTDATE(3)
          IHOUR_FCST  =FCSTDATE(4)
          IMINUTE_FCST=FCSTDATE(5)
          SECOND_FCST =0.
!
          IF(FCSTDATE(7)/=0)THEN
            SECOND_FCST=FCSTDATE(6)/(FCSTDATE(7)*1.)
          ENDIF
!
          CALL NEMSIO_GETHEADVAR(gfile,'NTIMESTEP',NTSD,iret=irtn)

          CALL NEMSIO_CLOSE(GFILE,iret=IERR)
!
        ELSE                                                                 !<-- Pure binary input
!
          select_unit: DO N=51,59
            INQUIRE(N,OPENED=OPENED)
            IF(.NOT.OPENED)THEN
              NFCST=N
              EXIT select_unit
            ENDIF
          ENDDO select_unit
!
          RESTART_FILENAME='restart_file_'//INT_TO_CHAR
          OPEN(unit=NFCST,file=RESTART_FILENAME,status='old',form='unformatted')
!
          READ(NFCST) IYEAR_FCST                                             !<-- Read time form restart file
          READ(NFCST) IMONTH_FCST                                            !
          READ(NFCST) IDAY_FCST                                              !
          READ(NFCST) IHOUR_FCST                                             !
          READ(NFCST) IMINUTE_FCST                                           !
          READ(NFCST) SECOND_FCST                                            !<--
!
          READ(NFCST) NTSD                                                   !<-- Read timestep from restart file
!
          CLOSE(NFCST)
!
        ENDIF  input
!
!-----------------------------------------------------------------------
!
        ISECOND_FCST=NINT(SECOND_FCST)                                     !<-- ESMF clock needs integer seconds
        NTSD_START=NTSD
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="RESTART: Set the Current Time of the Forecast"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=CURRTIME                                 &  !<-- Current time of the forecast (ESMF)
                         ,yy  =IYEAR_FCST                               &  !<-- Year from restart file
                         ,mm  =IMONTH_FCST                              &  !<-- Month from restart file
                         ,dd  =IDAY_FCST                                &  !<-- Day from restart file
                         ,h   =IHOUR_FCST                               &  !<-- Hour from restart file
                         ,m   =IMINUTE_FCST                             &  !<-- Minute from restart file
                         ,s   =ISECOND_FCST                             &  !<-- Second from restart file
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
      ENDIF restart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  With data from above set the local ESMF Clock
!***  to its correct time and timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the Current Time on the DOMAIN Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockSet(clock       =CLOCK_DOMAIN                      &  !<-- The DOMAIN Component's Clock
                        ,currtime    =CURRTIME                          &  !<-- Current time of simulation
                        ,advanceCount=NTSD_START                        &  !<-- Timestep at this current time
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!***  Create the time interval for printing clocktimes used by model 
!***  sections.  Read in forecast time interval for clocktime output 
!***  as well as the selected task ID that will provide the clocktimes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NHOURS_CLOCKTIME              &  !<-- Fill this variable
                                  ,label ='nhours_clocktime:'           &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read MPI Task ID That Provides Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NPE_PRINT                     &  !<-- Fill this variable
                                  ,label ='npe_print:'                  &  !<-- Give the variable this label's value from the config file
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
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_CLOCKTIME     &  !<-- Time interval between
                               ,h           =NHOURS_CLOCKTIME           &  !<-- Hours between clocktime writes (REAL)
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  How many tracer species are there?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_Init: Extract # of tracers from Config file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =NUM_TRACERS_MET               &  !<-- The variable filled (number of meteorological tracers)
                                  ,label ='num_tracers_met:'            &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =NUM_TRACERS_CHEM              &  !<-- The variable filled (number of chemical tracers)
                                  ,label ='num_tracers_met:'            &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Model-specific routines must be invoked in order to establish
!***  the ESMF Grid.  The different integration grids necessitate
!***  different ways of setting up both the parallelism for
!***  distributed memory runs and the ESMF Grid itself.
!***  When the parallelism is constructed, the local domain limits
!***  need to be inserted into the DOMAIN component's internal state
!***  if quilting is to be used.  See 'IF(QUILTING)THEN' below.
!-----------------------------------------------------------------------
!
      CALL NMM_SETUP(MYPE                                               &
                    ,COMM_MY_DOMAIN                                     &
                    ,CF(MY_DOMAIN_ID)                                   &
                    ,DOMAIN_GRID_COMP                                   &
                    ,DOMAIN_INT_STATE                                   &
                    ,GRID_DOMAIN )
!
!-----------------------------------------------------------------------
!***  Attach the NMM-specific ESMF Grid to the DOMAIN component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the NMM ESMF Grid to the DOMAIN Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=DOMAIN_GRID_COMP                   & !<-- The DOMAIN component
                           ,grid    =GRID_DOMAIN                        & !<-- Attach the ESMF grid to the DOMAIN component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics gridded subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!***  Since there is only a single integration grid, give the
!***  Dynamics the DOMAIN component's grid.
!***  Note that this subcomponent is part of the DOMAIN component's
!***  internal state.  This will be convenient if we need to reach
!***  the Dynamics component via the DOMAIN component such as happens
!***  when Write components are established.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Dynamics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the NMM Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%DYN_GRID_COMP=ESMF_GridCompCreate(               &
                               name      ="Dynamics component"          &  !<-- Name of the new Dynamics gridded component
                              ,config    =CF(MY_DOMAIN_ID)              &  !<-- Attach this configure file to the component
                              ,petList   =domain_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                              ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(domain_int_state%DYN_GRID_COMP      &  ! <-- The Dynamics gridded component
                                   ,DYN_REGISTER                        &  ! <-- The user's subroutineName for Register
                                   ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------
!***  Attach the ESMF Grid to the Dynamics component
!----------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach ESMF Grid to the NMM Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GRID_DYN=GRID_DOMAIN                                                 !<-- For now the Dyn Grid is the same as the DOMAIN Grid
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=domain_int_state%DYN_GRID_COMP     &  !<-- The Dynamics component
                           ,grid    =GRID_DYN                           &  !<-- The Dynamics ESMF grid
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create empty import and export states for the Dynamics component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%IMP_STATE_DYN=ESMF_StateCreate(                  &
                                            stateName="Dynamics Import" &  !<-- The Dynamics import state name
                                           ,statetype=ESMF_STATE_IMPORT &
                                           ,rc       =RC)
!
      domain_int_state%EXP_STATE_DYN=ESMF_StateCreate(                  &
                                            stateName="Dynamics Export" &  !<-- The Dynamics export state name
                                           ,statetype=ESMF_STATE_EXPORT &
                                           ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Physics gridded subcomponent if the user has not
!***  specified a no-physics run.
!***  Register the Initialize, Run, and Finalize steps for it.
!***  Since there is only a single integration grid, give the
!***  Physics component the DOMAIN component's grid.
!***  Note that this subcomponent is part of the DOMAIN component's
!***  internal state.  This will be convenient if we need to reach
!***  the Physics component via the DOMAIN component such as happens
!***  when Write components are established.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the flag from the DOMAIN import state indicating if the
!***  user wants physics to be active.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Physics Flag from DOMAIN Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                      &
                            ,name='PHYSICS_ON'                    &
                            ,value=PHYSICS_ON                     &
                            ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      physics_0: IF(PHYSICS_ON==ESMF_True)THEN
!
!-------------------------------
!***  Create Physics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the NMM Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        domain_int_state%PHY_GRID_COMP=ESMF_GridCompCreate(               &
                                 name      ="Physics component"           &  !<-- Name of the new Physics gridded component
                                ,config    =CF(MY_DOMAIN_ID)              &  !<-- Attach this configure file to the component
                                ,petList   =domain_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register the NMM Physics Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(domain_int_state%PHY_GRID_COMP  &  ! <-- The Physics gridded component
                                     ,PHY_REGISTER                    &  ! <-- The user's subroutineName
                                     ,RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------
!***  Attach the ESMF Grid to the Physics
!-----------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Attach the ESMF Grid to the Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        GRID_PHY=GRID_DOMAIN                                               !<-- For now the Physics Grid is the same as the DOMAIN Grid
!
        CALL ESMF_GridCompSet(gridcomp=domain_int_state%PHY_GRID_COMP   &  !<-- The NMM Physics component
                             ,grid    =GRID_PHY                         &  !<-- The ESMF grid of the Physics component
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!
      ENDIF physics_0
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Physics subcomponent
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%IMP_STATE_PHY=ESMF_StateCreate(                  &
                                            stateName="Physics Import"  &  !<-- The Physics import state
                                           ,statetype=ESMF_STATE_IMPORT &
                                           ,rc       =RC)
!
      domain_int_state%EXP_STATE_PHY=ESMF_StateCreate(                  &
                                            stateName="Physics Export"  &  !<-- The Physics export state
                                           ,statetype=ESMF_STATE_EXPORT &
                                           ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics-Physics coupler subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------------------
!***  Create Dyn-Phy Coupler
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Dynamics-Physics Coupler Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%COUPLER_DYN_PHY_COMP=ESMF_CplCompCreate          &
                                (name   ="Dyn-Phy coupler component"    &
                                ,petList=domain_int_state%PETLIST_FCST  &  !<-- The forecast task IDs
                                ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the Dyn-Phy Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetServices(domain_int_state%COUPLER_DYN_PHY_COMP &  ! <-- The Dyn-Phy coupler component
                                  ,DYN_PHY_CPL_REGISTER                  &  ! <-- The user's subroutineName
                                  ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The Dynamics and Physics import states have already been created.  
!***  Add some key flags to those import states prior to Initialization.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain ID to the Dyn Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_DYN       &  !<-- The Dynamics component import state
                            ,name ='DOMAIN_ID'                          &  !<-- Use this name inside the state
                            ,value=MY_DOMAIN_ID                         &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain ID to the Phy Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_PHY       &  !<-- The Dynamics component import state
                            ,name ='DOMAIN_ID'                          &  !<-- Use this name inside the state
                            ,value=MY_DOMAIN_ID                         &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Insert the flag indicating if the DOMAIN component is a nest.
!***  The Dynamics component needs to know this regarding BC's.
!***  Both Dynamics and Physics need it to properly compute some
!***  fundamental aspects of the nested grids.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Nest Flag to the Dyn Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_DYN       &  !<-- The Dynamics component import state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- Use this name inside the state
                            ,value=I_AM_A_NEST                          &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Nest Flag to the Phy Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_PHY       &  !<-- The Physics component import state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- Use this name inside the state
                            ,value=I_AM_A_NEST                          &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  If this is a nest domain then insert the Parent-Child timestep
!***  ratio into the Dynamics import state since that will be needed
!***  to tell the Dynamics Run step how often to update the boundary
!***  tendencies.
!***  Also insert the flag indicating whether or not the nest domain
!***  alerady has an input file or if one needs to be generated by
!***  its parent.
!------------------------------------------------------------------------
!
      IF(I_AM_A_NEST==ESMF_TRUE)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Parent-Child Time Ratio to the Dyn Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_DYN     &  !<-- The Dynamics component import state
                              ,name ='Parent-Child Time Ratio'          &  !<-- Use this name inside the state
                              ,value=PARENT_CHILD_TIME_RATIO            &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Input-Ready Flag to the Dyn Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_DYN     &  !<-- The Dynamics component import state
                              ,name ='Input Ready'                      &  !<-- Use this name inside the state
                              ,value=INPUT_READY                        &  !<-- Does this nest's input file already exist?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Input-Ready Flag to the Phy Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_PHY     &  !<-- The Physics component import state
                              ,name ='Input Ready'                      &  !<-- Use this name inside the state
                              ,value=INPUT_READY                        &  !<-- Does this nest's input file already exist?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  set up the Write component(s).
!***  This must be done prior to executing the Initialize steps
!***  of the Dynamics and Physics components because the Write
!***  components' import states are required by those steps
!***  when 'quilting' is active and WRITE_SETUP is the routine
!***  in which the Write components' import states are inserted
!***  into the Dynamics and Physics export states which are
!***  in turn part of the DOMAIN internal state.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(QUILTING)THEN
!
        CALL WRITE_SETUP(DOMAIN_GRID_COMP                               &
                        ,DOMAIN_INT_STATE                               &
                        ,CLOCK_DOMAIN )
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the Initialize steps for the gridded subcomponents.
!***  These are the Initialize subroutines specified in the
!***  Register routines called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!--------------
!***  Dynamics
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize the NMM Dynamics Component phase 1"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%DYN_GRID_COMP  &  !<-- The dynamics gridded component
                                  ,importState=domain_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                                  ,exportState=domain_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                                  ,clock      =CLOCK_DOMAIN                    &  !<-- The DOMAIN clock
                                  ,phase      =1                               &
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      physics_1: IF(PHYSICS_ON==ESMF_True)THEN
!
!-----------------------------------------------------------------------
!
!-------------
!***  Physics
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Initialize Physics Component phase 1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%PHY_GRID_COMP  &  !<-- The physics gridded component
                                    ,importState=domain_int_state%IMP_STATE_PHY  &  !<-- The physics import state
                                    ,exportState=domain_int_state%EXP_STATE_PHY  &  !<-- The physics export state
                                    ,clock      =CLOCK_DOMAIN                    &  !<-- The DOMAIN clock
                                    ,phase      =1                               &
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF physics_1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Initialize the Dyn-Phy coupler subcomponent.
!***  The choice of import and export state does not matter
!***  for the initialize step.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dyn-Phy Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompInitialize(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dyn_phy coupler component
                                 ,importState=domain_int_state%EXP_STATE_DYN        &  !<-- The dyn-phy coupler import state
                                 ,exportState=domain_int_state%IMP_STATE_PHY        &  !<-- The dyn-phy coupler export state
                                 ,clock      =CLOCK_DOMAIN                          &  !<-- The DOMAIN Clock
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompInitialize(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dyn_phy coupler component
                                 ,importState=domain_int_state%EXP_STATE_PHY        &  !<-- The dyn-phy coupler import state
                                 ,exportState=domain_int_state%IMP_STATE_DYN        &  !<-- The dyn-phy coupler export state
                                 ,clock      =CLOCK_DOMAIN                          &  !<-- The DOMAIN Clock
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
!--------------
!***  Dynamics
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize the NMM Dynamics Component phase 2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%DYN_GRID_COMP  &  !<-- The dynamics gridded component
                                  ,importState=domain_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                                  ,exportState=domain_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                                  ,clock      =CLOCK_DOMAIN                    &  !<-- The DOMAIN clock
                                  ,phase      =2                               &
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      physics_2: IF(PHYSICS_ON==ESMF_True)THEN
!
!-----------------------------------------------------------------------
!
!-------------
!***  Physics
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Initialize Physics Component phase 2"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%PHY_GRID_COMP  &  !<-- The physics gridded component
                                    ,importState=domain_int_state%IMP_STATE_PHY  &  !<-- The physics import state
                                    ,exportState=domain_int_state%EXP_STATE_PHY  &  !<-- The physics export state
                                    ,clock      =CLOCK_DOMAIN                    &  !<-- The DOMAIN clock
                                    ,phase      =2                               &
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF physics_2
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  execute the Initialize step of the Write component(s).
!***  This must be done after the initialization of the
!***  Dynamics and Physics components because those components'
!***  internal states contain the history output variables.
!***  Pointers to those variables are set during their INITIALIZE
!***  steps and are then loaded into the Write components'
!***  import states which themselves reside in the Dynamics/Physics
!***  export states.
!-----------------------------------------------------------------------
!
      I_AM_A_FCST_TASK=ESMF_TRUE
!
      IF(QUILTING)THEN
!
        CALL WRITE_INIT(DOMAIN_GRID_COMP                                &
                       ,DOMAIN_INT_STATE                                &
                       ,CLOCK_DOMAIN)
!
          IF(MYPE>=domain_int_state%NUM_PES_FCST)THEN
            I_AM_A_FCST_TASK=ESMF_FALSE
          ENDIF
!
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Fcst-or-Write Task Flag to the DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The DOMAIN component export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_FCST_TASK                     &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now extract number of children on this DOMAIN component's domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Extract Number of Children from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The DOMAIN import state
                            ,name ='NUM_CHILDREN'                       &  !<-- Name of the attribute to extract
                            ,value=NUM_CHILDREN                         &  !<-- Put the Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  If the current component/domain is the parent of nests then:
!
!***  (1) Extract the arrays from the Dyn and/or Phy export states
!***      that are required for the children's boundaries and
!***      insert them into the DOMAIN export state since ultimately 
!***      they must be available to the parent in the
!***      Parent-Child Coupler.
!
!***  (2) Check to see if the children have input data ready for them.
!***      If not, do simple nearest neighbor and bilinear  interpolation
!***      from the parent's grid to the children's.  Write out that
!***      interpolated data into files that are waiting for the children 
!***      when they recursively execute DOMAIN_INITIALIZE themselves.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      I_AM_A_PARENT=ESMF_FALSE
!
!-----------------------------------------------------------------------
!***  Extract from the Dynamics and Physics export state the quantities
!***  relevant to the children's boundaries.  Insert them into the 
!***  DOMAIN export state so that DOMAIN_DRIVER can take them and send them 
!***  to the Parent-Child coupler.  Only the forecast tasks participate
!***  in doing this since the Quilt/Write tasks never loaded data into
!***  the Dynamics or Physics export states.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_FCST_TASK==ESMF_TRUE)THEN       
!
        CALL PARENT_DATA_TO_DOMAIN(domain_int_state%EXP_STATE_DYN       &  !<-- The Dynamics export state
                                  ,domain_int_state%EXP_STATE_PHY       &  !<-- The Physics export state
                                  ,EXP_STATE)                              !<-- The DOMAIN export state
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      child_init_block: IF(NUM_CHILDREN>0)THEN                             !<-- Only parents participate                              
!
        I_AM_A_PARENT=ESMF_TRUE
!
!-----------------------------------------------------------------------
!***  Initialize the children's data directly from the parent if
!***  there are no pre-processed input files ready for them.
!***  Files will be written for the children to read in as usual.
!***  Only parent tasks participate.
!-----------------------------------------------------------------------
!
        ALLOCATE(MY_CHILDREN_ID(1:NUM_CHILDREN))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Extract Children's IDs from Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The DOMAIN import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=MY_CHILDREN_ID                 &  !<-- Put the Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        child_init_loop: DO N=1,NUM_CHILDREN                               !<-- Loop through the children
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          MESSAGE_CHECK="Extract Children's Input Flag from Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_CHILDREN_ID(N))     &  !<-- The config object
                                      ,value =INPUT_READY_MY_CHILD      &  !<-- Child's flag for existence of its input file
                                      ,label ='input_ready:'            &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          IF(.NOT.INPUT_READY_MY_CHILD)THEN                                  !<-- INPUT_READY=false -> This child has no input file 
                                                                             !      so parent will generate input.
            CALL PARENT_TO_CHILD_INIT_NMM(MYPE                            &  !<-- This task's rank (in)
                                         ,CF                              &  !<-- Array of configure files (in)
                                         ,MY_DOMAIN_ID                    &  !<-- Each domain's ID (in)
                                         ,MY_CHILDREN_ID(N)               &  !<-- The child's domain ID
                                         ,domain_int_state%DYN_GRID_COMP  &  !<-- The parent's Dynamics Component (inout)
                                         ,domain_int_state%PHY_GRID_COMP  &  !<-- The parent's Physics Component (inout)
                                         ,COMM_MY_DOMAIN )                   !<-- Each domain's intracommunicator
!
          ENDIF
!            
!-----------------------------------------------------------------------
!
        ENDDO child_init_loop
!
!-----------------------------------------------------------------------
!
        DEALLOCATE(MY_CHILDREN_ID)
!
      ENDIF child_init_block
!
!-----------------------------------------------------------------------
!***  Insert into the DOMAIN export state the flag indicating if the
!***  current domain is a parent.  The DOMAIN Driver wants to know this
!***  since most Parent-Child work can be ignored by domains with
!***  no children.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Parent/Not-a-Parent Flag to the DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The DOMAIN component export state
                            ,name ='I-Am-A-Parent Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_PARENT                        &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      rtim1=rtc()
      total_integ_tim=(timef()-btim0)
!
      if(mype==0)write(0,*)' domain_init_tim=',total_integ_tim,(rtim1-rtim0)
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DOMAIN INITIALIZE step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN INITIALIZE step failed RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_RUN(DOMAIN_GRID_COMP                            &
                           ,IMP_STATE                                   &
                           ,EXP_STATE                                   &
                           ,CLOCK_DOMAIN                                &
                           ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The Run step of the DOMAIN component for the NMM.
!***  The forecast tasks execute the Run step of the NMM-B Dynamics.
!***  This is the Run subroutine specified in the Dynamics Register
!***  routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!
      USE MODULE_NESTING,ONLY: BOUNDARY_DATA_STATE_TO_STATE
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN gridded component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                       &  !<-- The DOMAIN Run step's import state
                                       ,EXP_STATE                          !<-- The DOMAIN Run step's export state
!
      TYPE(ESMF_Clock), INTENT(INOUT) :: CLOCK_DOMAIN                      !<-- The DOMAIN ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN                                        !<-- Return code for the Run step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(ESMF_KIND_I4) :: INTEGER_DT,NUMERATOR_DT,IDENOMINATOR_DT     
!
      INTEGER(kind=KINT) :: HDIFF_ON,RC
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF
!
      TYPE(ESMF_Config) :: CF  
!
      real(kind=8) :: rtim0,rtim1,rtc
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rtim0=rtc()
      btim0=timef()
!
      RC    =ESMF_SUCCESS
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!
!     ALLOCATE(DOMAIN_INT_STATE,stat=RC)
!     wrap%DOMAIN_INT_STATE=>DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!
      fcst_pes: IF(MYPE<domain_int_state%NUM_PES_FCST)THEN                 !<-- Only the forecast tasks integrate
!
!-----------------------------------------------------------------------
!***  Extract the timestep from the Clock so that we know the direction
!***  of the integration.  We skip all aspects of Physics if the time
!***  step is negative.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="DOMAIN_Run: Extract the ESMF Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_ClockGet(clock   =CLOCK_DOMAIN                        &
                          ,timeStep=DT_ESMF                             &
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="DOMAIN_Run: Extract Components of the Timestep" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                  &  !<-- the ESMF timestep
                                 ,s           =INTEGER_DT               &  !<-- the integer part of the timestep in seconds
                                 ,sN          =NUMERATOR_DT             &  !<-- the numerator of the fractional second
                                 ,sD          =IDENOMINATOR_DT          &  !<-- the denominator of the fractional second
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  We must transfer the horizontal diffusion flag from the
!***  DOMAIN import state to the Dynamics import state.  The Dynamics
!***  import state is not available outside of the integration time
!***  loop in NMM_INTEGRATE which lies in DOMAIN_DRIVER therefore we
!***  need to perform this transfer each timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Horizontal Diffusion Flag from DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The DOMAIN import state
                              ,name ='HDIFF'                            &  !<-- Name of the attribute to extract
                              ,value=HDIFF_ON                           &  !<-- The ID of this domain
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Horizontal Diffusion Flag to the Dyn Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_DYN     &  !<-- The Dynamics component import state
                              ,name ='HDIFF'                            &  !<-- Use this name inside the state
                              ,value=HDIFF_ON                           &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)  
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If this is a nested domain, new boundary data must first
!***  be moved from the DOMAIN import state to that of the Dynamics
!***  every N timesteps where N is the number of the nest's timesteps 
!***  within one timestep of its parent.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_NEST==ESMF_TRUE)THEN
!         write(6,*)' DOMAIN_RUN calling BOUNDARY_DATA_STATE_TO_STATE'
!         call print_memory()
          CALL BOUNDARY_DATA_STATE_TO_STATE(clock    =CLOCK_DOMAIN                  &  !<-- The DOMAIN Clock
                                           ,ratio    =PARENT_CHILD_TIME_RATIO       &  !<-- # of child timesteps per parent timestep
                                           ,state_in =IMP_STATE                     &  !<-- DOMAIN component's import state
                                           ,state_out=domain_int_state%IMP_STATE_DYN)  !<-- The Dynamics import state
!         write(6,*)' DOMAIN_RUN called BOUNDARY_DATA_STATE_TO_STATE'
!         call print_memory()
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Execute the Run Step for Dynamics"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!       write(6,*)' DOMAIN_RUN calling DYN_RUN'
!       call print_memory()
        CALL ESMF_GridCompRun(gridcomp   =domain_int_state%DYN_GRID_COMP  &  !<-- The dynamics component
                             ,importState=domain_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                             ,exportState=domain_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                             ,clock      =CLOCK_DOMAIN                    &  !<-- The DOMAIN Clock
                             ,rc         =RC)        
!       write(6,*)' DOMAIN_RUN called DYN_RUN'
!       call print_memory()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  If integration is forward and Physics is turned on then proceed
!***  with Physics and the associated coupling to Dynamics.
!-----------------------------------------------------------------------
!
        physics: IF(INTEGER_DT>0.AND.PHYSICS_ON==ESMF_True)THEN            !<-- Physics is active
!
!-----------------------------------------------------------------------
!***  Bring export data from the Dynamics into the Coupler and
!***  point the Physics import state at it.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          MESSAGE_CHECK="Coupler Moves Data from Dynamics to Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          CALL ESMF_CplCompRun(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=domain_int_state%EXP_STATE_DYN        &  !<-- The Coupler import state = Dynamics export state
                              ,exportState=domain_int_state%IMP_STATE_PHY        &  !<-- The Coupler export state = Physics import state
                              ,clock      =CLOCK_DOMAIN                          &  !<-- The DOMAIN clock
                              ,rc         =RC)        
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! - FIX later
 RC=ESMF_SUCCESS
 RC_RUN=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Physics.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          MESSAGE_CHECK="Execute the Run Step for Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          CALL ESMF_GridCompRun(gridcomp   =domain_int_state%PHY_GRID_COMP &  !<-- The physics component
                               ,importState=domain_int_state%IMP_STATE_PHY &  !<-- The physics import state 
                               ,exportState=domain_int_state%EXP_STATE_PHY &  !<-- The physics export state
                               ,clock      =CLOCK_DOMAIN                   &  !<-- The DOMAIN Clock
                               ,rc         =RC)        
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Bring export data from the Physics into the Coupler and
!***  point the Dynamics import state at it.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          MESSAGE_CHECK="DOMAIN_RUN: Coupler Moves Data from Physics to Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!     write(0,*)' DOMAIN_RUN before Run CPL Comp mype=',mype,' integer_dt=',integer_dt
          CALL ESMF_CplCompRun(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=domain_int_state%EXP_STATE_PHY        &  !<-- The Coupler import state = Physics export state
                              ,exportState=domain_int_state%IMP_STATE_DYN        &  !<-- The Coupler export state = Dynamics import state
                              ,clock      =CLOCK_DOMAIN                          &  !<-- The DOMAIN Clock
                              ,rc         =RC)        
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! - FIX later
 RC=ESMF_SUCCESS
 RC_RUN=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!
        ELSE  physics                                                      !<-- Physics is not active
!
!-----------------------------------------------------------------------
!***  If integration is backward or the user has requested that
!***  the forecast be run without Physics then simply redirect
!***  the export data of the Dynamics back to the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          MESSAGE_CHECK="Execute Dyn-Phy Coupler w/o Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          CALL ESMF_CplCompRun(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=domain_int_state%EXP_STATE_DYN        &  !<-- The Coupler import state = Dynamics export state
                              ,exportState=domain_int_state%IMP_STATE_DYN        &  !<-- The Coupler export state = Dynamics import state
                              ,clock      =CLOCK_DOMAIN                          &  !<-- The DOMAIN Clock
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!
        ENDIF  physics
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_pes
!
!-----------------------------------------------------------------------
!
      domain_tim=domain_tim+(timef()-btim0)
!
      rtim1=rtc()
      total_integ_tim=total_integ_tim+(timef()-btim0)
!     write(0,*)'exit DOMAIN_RUN integration time ',total_integ_tim,(rtim1-rtim0)
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DOMAIN RUN step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN RUN step failed RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_FINALIZE(DOMAIN_GRID_COMP                       &
                                ,IMP_STATE                              &
                                ,EXP_STATE                              &
                                ,CLOCK_DOMAIN                           &
                                ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  This routine Finalizes the DOMAIN gridded component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN gridded component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                       &  !<-- The DOMAIN finalize step's import state
                                       ,EXP_STATE                          !<-- The DOMAIN finalize step's export state
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The DOMAIN ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FINALIZE                                   !<-- Return code for the Finalize step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J
      INTEGER :: RC                                                        ! The final error signal variables.
!
      CHARACTER(50):: MODE
!
      TYPE(ESMF_Logical) :: PHYSICS_ON     
!
      TYPE(ESMF_Config) :: CF                                              !<-- The config object
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP                             !<-- The F90 wrap of the DOMAIN internal state
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The DOMAIN internal state pointer
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Retrieve the config object CF from the DOMAIN component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from DOMAIN Component"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the diabatic/adiabatic flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Adiabatic Flag from Config Object"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(TRIM(MODE)=='TRUE')THEN
        PHYSICS_ON=ESMF_False
        write(0,*)' Finalize without physics coupling. '
      ELSE
        PHYSICS_ON=ESMF_True 
        write(0,*)' Finalize with physics coupling. '
      ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve the DOMAIN component's internal state.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The DOMAIN component
                                        ,WRAP                           &  !<-- The F90 wrap of the DOMAIN internal state
                                        ,RC)
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Finalize each of the subcomponents.
!-----------------------------------------------------------------------
!
!--------------
!***  Dynamics
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =domain_int_state%DYN_GRID_COMP &
                                ,importState=domain_int_state%IMP_STATE_DYN &
                                ,exportState=domain_int_state%EXP_STATE_DYN &
                                ,clock      =CLOCK_DOMAIN                &
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINALIZE=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(PHYSICS_ON==ESMF_True)THEN
!
!-------------
!***  Physics
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompFinalize(gridcomp   =domain_int_state%PHY_GRID_COMP &
                                  ,importState=domain_int_state%IMP_STATE_PHY &
                                  ,exportState=domain_int_state%EXP_STATE_PHY &
                                  ,clock      =CLOCK_DOMAIN                   &
                                  ,rc         =RC)
      ENDIF 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINALIZE=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------
!***  Dynamics-Physics coupler
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompFinalize(cplcomp    =domain_int_state%COUPLER_DYN_PHY_COMP &
                               ,importState=domain_int_state%EXP_STATE_DYN        &
                               ,exportState=domain_int_state%IMP_STATE_PHY        &
                               ,clock      =CLOCK_DOMAIN                          &
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy all States.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(state=domain_int_state%IMP_STATE_DYN       &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=domain_int_state%EXP_STATE_DYN       &
                            ,rc   =RC)
!
      IF(PHYSICS_ON==ESMF_True)THEN
        CALL ESMF_StateDestroy(state=domain_int_state%IMP_STATE_PHY     &
                              ,rc   =RC)
!
        CALL ESMF_StateDestroy(state=domain_int_state%EXP_STATE_PHY     &
                              ,rc   =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  finalize and destroy objects related to the Write components.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(domain_int_state%QUILTING)THEN
        CALL WRITE_DESTROY(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Destroy the DOMAIN Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy DOMAIN Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_DOMAIN                         &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy all subcomponents.
!-----------------------------------------------------------------------
!
!--------------
!***  Dynamics
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompDestroy(gridcomp=domain_int_state%DYN_GRID_COMP & 
                               ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  Physics
!-------------
!
      IF(PHYSICS_ON==ESMF_TRUE)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Destroy Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompDestroy(gridcomp=domain_int_state%PHY_GRID_COMP  &
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!------------------------------
!***  Dynamics-Physics coupler
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompDestroy(cplcomp=domain_int_state%COUPLER_DYN_PHY_COMP &
                              ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
        WRITE(0,*)'DOMAIN FINALIZE step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN FINALIZE step failed'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_FILTERING(DOMAIN_GRID_COMP                         &
                              ,IMP_STATE                                &
                              ,EXP_STATE                                &
                              ,CLOCK_DOMAIN                             &
                              ,RC_FILT)
!
!-----------------------------------------------------------------------
!***  Phase 2 of the Run step of the DOMAIN component.
!***  This phase is only relevant when digital filtering is
!***  in effect and executes at the end of each timestep
!***  after the Dynamics and Physics (in phase 1). 
!
!***  Called from subroutine NMM_INTEGRATE.
!-----------------------------------------------------------------------
!
      USE module_DIGITAL_FILTER_NMM
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN gridded component
!
      TYPE(ESMF_State),INTENT(IN) :: IMP_STATE                          &  !<-- The DOMAIN import state
                                    ,EXP_STATE                             !<-- The DOMAIN export state
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The DOMAIN ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FILT                                       !<-- Return code for this step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: MEAN_ON,NDFISTEP                            &
                           ,NUM_TRACERS_CHEM,NUM_TRACERS_MET
!
      LOGICAL(kind=KLOG),SAVE :: FIRST_PASS=.TRUE.
!
      CHARACTER(7) :: CLOCK_DIRECTION
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The current time of Clock_DOMAIN
                        ,STARTTIME                                         !<-- The start time of Clock_DOMAIN
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER,SAVE :: DOMAIN_INT_STATE         !<-- The DOMAIN internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE),SAVE :: WRAP                        !<-- The F90 wrap of the DOMAIN internal state
!
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_FILT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  We need to extract the DOMAIN internal state only once.
!-----------------------------------------------------------------------
!
      IF(FIRST_PASS)THEN
!
        FIRST_PASS=.FALSE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="NMM_Filtering: Extract the DOMAIN Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP             &  !<-- The DOMAIN component
                                          ,WRAP                         &  !<-- The F90 wrap of the DOMAIN internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  What are the start time and the current time?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract StartTime,CurrentTime"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_DOMAIN                        &
                        ,startTime=STARTTIME                           &
                        ,currTime =CURRTIME                            &
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is the Clock direction?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract Clock Direction."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- Extract the direction of the Clock from the import state
                            ,name ='Clock_Direction'                    &
                            ,value=CLOCK_DIRECTION                      &
                            ,rc   =RC )
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- Extract MEAN_ON flag from import state
                            ,name ='MEAN_ON'                            &
                            ,value=MEAN_ON                              &
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<domain_int_state%NUM_PES_FCST)THEN               !<-- Only forecast tasks deal with the integration
!
!-----------------------------------------------------------------------
!
!-----------------------
!***  The initial stage
!-----------------------
!
        IF(CURRTIME==STARTTIME)THEN
!
          CALL ESMF_AttributeGet(state=IMP_STATE                        &  !<-- Extract the filter value NDFISTEP from the import state
                                ,name ='NDFISTEP'                       &
                                ,value=NDFISTEP                         &
                                ,rc   =RC )
!
          CALL DIGITAL_FILTER_DYN_INIT_NMM(domain_int_state%IMP_STATE_DYN &
                                          ,NDFISTEP                       &
                                          ,NUM_TRACERS_MET                &
                                          ,NUM_TRACERS_CHEM)
!
          CALL DIGITAL_FILTER_PHY_INIT_NMM(domain_int_state%IMP_STATE_PHY)
!
        ENDIF
!
!-----------------------------------------------------------------------
        direction: IF(CLOCK_DIRECTION=='Forward')THEN
!-----------------------------------------------------------------------
!
!-------------------------
!***  The summation stage
!-------------------------
!
          IF(CURRTIME>=STARTTIME)THEN
!
            CALL DIGITAL_FILTER_DYN_SUM_NMM(domain_int_state%IMP_STATE_DYN &
                                           ,MEAN_ON                        &
                                           ,NUM_TRACERS_MET                &
                                           ,NUM_TRACERS_CHEM)
          ENDIF
!
!---------------------
!
          IF(CURRTIME==HALFDFITIME)THEN
!
            CALL DIGITAL_FILTER_PHY_SAVE_NMM(domain_int_state%IMP_STATE_PHY)
!
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          IF(CURRTIME==DFITIME)THEN
            write(0,*)' DFI at final dfitime'
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(domain_int_state%IMP_STATE_DYN &
                                               ,NUM_TRACERS_MET                &
                                               ,NUM_TRACERS_CHEM)
!
            CALL DIGITAL_FILTER_PHY_RESTORE_NMM(domain_int_state%IMP_STATE_PHY)
!
            CALL ESMF_ClockPrint(clock  =CLOCK_DOMAIN                   &
                                ,options="currtime string"              &
                                ,rc     =RC)
          ENDIF
!
!-----------------------------------------------------------------------
        ELSEIF(CLOCK_DIRECTION=='Bckward')THEN
!-----------------------------------------------------------------------
!
!-------------------------
!***  The summation stage
!-------------------------
!
          IF(CURRTIME<=STARTTIME)THEN
            CALL DIGITAL_FILTER_DYN_SUM_NMM(domain_int_state%IMP_STATE_DYN &
                                           ,MEAN_ON                        &
                                           ,NUM_TRACERS_MET                &
                                           ,NUM_TRACERS_CHEM)
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          IF(CURRTIME==DFITIME)THEN
            write(0,*)' DFI at final dfitime '
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(domain_int_state%IMP_STATE_DYN &
                                               ,NUM_TRACERS_MET                &
                                               ,NUM_TRACERS_CHEM)
!
! ----------------------------------------------------------------------
!
            CALL ESMF_ClockPrint(clock  =CLOCK_DOMAIN                   &
                                ,options="currtime string"              &
                                ,rc     =RC)
          ENDIF
!-----------------------------------------------------------------------
!
        ENDIF direction
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_FILTERING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CALL_WRITE_ASYNC(DOMAIN_GRID_COMP                      &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_DOMAIN                          &
                                 ,RC_RUN2)
!
!-----------------------------------------------------------------------
!***  Phase 3 of the Run step of the NMM DOMAIN component.
!***  It initiates the writing of history/restart files
!***  from each DOMAIN component.
!-----------------------------------------------------------------------
!
      USE module_NMM_INTEGRATE,ONLY : ALARM_HISTORY                     &
                                     ,ALARM_RESTART
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN gridded component
!
      TYPE(ESMF_State),INTENT(IN)    :: IMP_STATE                          !<-- The DOMAIN Run step's import state
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE                          !<-- The DOMAIN Run step's export state
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The DOMAIN ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN2                                       !<-- Return code for the Run step 
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The DOMAIN internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP                             !<-- The F90 wrap of the DOMAIN internal state
!
      CHARACTER(ESMF_MAXSTR) :: CWRT
!
      INTEGER(kind=KINT) :: RC                                             !<-- Error signal variable.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_RUN2=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Write a history file at the end of the appropriate timesteps.
!-----------------------------------------------------------------------
!
      IF(ESMF_AlarmIsRinging(alarm=ALARM_HISTORY                        &  !<-- The history output alarm
                            ,rc   =RC))THEN
!
!-----------------------------------------------------------------------
!***  Retrieve the DOMAIN component's internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Run2: Retrieve DOMAIN Component's Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the DOMAIN internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Execute the writing of a history file.
!-----------------------------------------------------------------------
!
!     write(6,*)' calling WRITE_ASYNC'
!     call print_memory()
        IF(domain_int_state%QUILTING)THEN
          CWRT='History'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN                                 &
                          ,MYPE                                         &
                          ,CWRT)
        ENDIF
!     write(6,*)' called WRITE_ASYNC'
!     call print_memory()
!
      ENDIF 
!
!-----------------------------------------------------------------------
!***  Write a restart file at the end of the appropriate timesteps.
!-----------------------------------------------------------------------
!
      IF(ESMF_AlarmIsRinging(alarm=ALARM_RESTART                        &  !<-- The restart output alarm
                            ,rc   =RC))THEN
!
!-----------------------------------------------------------------------
!***  Retrieve the DOMAIN component's internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Run2: Retrieve DOMAIN Component's Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the DOMAIN internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Execute the writing of a restart file.
!-----------------------------------------------------------------------
!
!     write(6,*)' calling WRITE_ASYNC for restart'
!     call print_memory()
        IF(domain_int_state%QUILTING)THEN
          CWRT='Restart'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN                                 &
                          ,MYPE                                         &
                          ,CWRT)
        ENDIF
!     write(6,*)' called WRITE_ASYNC for restart'
!     call print_memory()
!
      ENDIF 
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CALL_WRITE_ASYNC
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE module_DOMAIN_GRID_COMP
!
!-----------------------------------------------------------------------
