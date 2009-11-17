!-----------------------------------------------------------------------
!
      MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE ATM (Atmosphere) GRIDDED COMPONENT MODULE.
!***  IT WILL SET UP DYNAMICS, PHYSICS, AND COUPLER SUBCOMPONENTS
!***  AND RUN THEIR INITIALIZE, RUN, AND FINALIZE ROUTINES.
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
!   2009-11-03  Lu    - Add GOCART and ChemRegistry modules
!
! USAGE: ATM Gridded component parts called from subroutines within
!        module_ATM_DRIVER_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_INCLUDE
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
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
      USE GFS_DYNAMICS_GRID_COMP_MOD,ONLY: GFS_DYN_SETSERVICES
      USE GFS_PHYSICS_GRID_COMP_MOD ,ONLY: GFS_PHY_SETSERVICES
      USE MODULE_GFS_INTEGRATE      ,ONLY: GFS_INTEGRATE
!
      USE ATMOS_DYN_PHY_CPL_COMP_MOD,ONLY: ATM_CPL_SETSERVICES
!
      USE MODULE_GFS_MPI_DEF        ,ONLY: ENSMEM_NAME                  &
                                          ,PETLIST_FCST                 &
                                          ,WRITE_GROUPS
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
      USE MODULE_CLOCKTIMES,ONLY: total_integ_tim
!
!-----------------------------------------------------------------------
!***  LIST MODULES FOR GSFC CHEMISTRY PACKAGE
!-----------------------------------------------------------------------
!
      USE GOCART_GridCompMod    , ONLY: GOCART_SETSERVICES => SETSERVICES 
!
      USE Chem_RegistryMod
!
!-----------------------------------------------------------------------
!***  LIST OTHER MODULES WITH NON-GENERIC ROUTINES USED BY ATM.
!-----------------------------------------------------------------------
!
      USE MODULE_WRITE_ROUTINES ,ONLY: WRITE_INIT,WRITE_ASYNC             !<-- These are routines used only when asynchronous
      USE MODULE_WRITE_GRID_COMP,ONLY: WRITE_SETUP                      & !    quilting is specified by the user in the
                                      ,WRITE_DESTROY                      !    configure file for NMM output.
!
      USE MODULE_WRITE_GRID_COMP_GFS,ONLY: WRITE_SETUP_GFS              & !<-- These are routines used only when asynchronous
                                          ,WRITE_DESTROY_GFS              !    quilting is specified by the user in the
      USE MODULE_WRITE_ROUTINES_GFS ,ONLY: WRITE_INIT_GFS                 !    configure file for GFS output.

!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
#include "../../../inc/ESMF_LogMacros.inc"
!-----------------------------------------------------------------------
! 
      PRIVATE
!
      PUBLIC :: ATM_REGISTER 
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: INPES,JNPES                                 &  !<-- MPI tasks in I and J directions
                           ,MYPE                                        &  !<-- Each MPI task ID
                           ,NPE_PRINT                                   &
                           ,NUM_TRACERS_CHEM                            &  !<-- Number of chemistry tracer variables
                           ,NUM_TRACERS_MET                             &  !<-- Number of meteorological tracer variables
                           ,WRITE_GROUP_READY_TO_GO                        !<-- The write group to use
!
      CHARACTER(3),SAVE :: CORE                                            !<-- The name of the selected dynamic core
!
      CHARACTER(ESMF_MAXSTR) :: INFILE="restart_file"
!
      LOGICAL(kind=KLOG) :: QUILTING                                    &  !<-- Is asynchronous quilting specified?
                           ,RESTARTED_RUN                               &  !<-- Original/restarted run logical flag
                           ,STANDALONE_POST                                !<-- Logical flag for running standalone post
!
      TYPE(ESMF_Config),SAVE :: CF_1                                       !<-- The principal config object
!
      TYPE(ESMF_VM),SAVE :: VM                                             !<-- The ESMF virtual machine.
!
      TYPE(ATM_INTERNAL_STATE),POINTER,SAVE :: ATM_INT_STATE               !<-- The NMM ATM internal state pointer
      TYPE(WRAP_ATM_INTERNAL_STATE)   ,SAVE :: WRAP                        !<-- The F90 wrap of the NMM ATM internal state
!
      TYPE(ESMF_Time),SAVE :: DFITIME                                   &
                             ,HALFDFITIME 
!
      TYPE(ESMF_TimeInterval),SAVE :: HALFDFIINTVAL_BCK                 &
                                     ,HALFDFIINTVAL_FWD
!
      TYPE(ESMF_TimeInterval),SAVE :: TIMEINTERVAL_CLOCKTIME               !<-- The ESMF time interval between NMM clocktime output
!
      TYPE(ESMF_Logical),SAVE :: PHYSICS_ON                                !<-- Is physics active?
!
!-----------------------------------------------------------------------
!***  FOR GSFC CHEMISTRY PACKAGE
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),SAVE :: GC_GFS_CHEM                              !<-- The GFS chemistry component 
      TYPE(ESMF_State),   SAVE :: IMP_GFS_CHEM,EXP_GFS_CHEM                !<-- Import/export states for GFS Chemistry
!
      TYPE(Chem_Registry),SAVE :: REG                                      !<-- The GOCART Chem_Registry

      TYPE(ESMF_Logical), SAVE :: CHEMISTRY_ON                             !<-- Is chemistry active?

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
!-----------------------
!***  For GFS Ensembles
!-----------------------
!
      TYPE(ESMF_GridComp),SAVE :: GC_GFS_DYN                            &
                                 ,GC_GFS_PHY
!
      TYPE(ESMF_GridComp),DIMENSION(:),POINTER,SAVE :: WRT_COMPS
!
      TYPE(ESMF_CplComp), SAVE :: GC_ATM_CPL
!
      TYPE(ESMF_State),SAVE :: IMP_GFS_DYN,EXP_GFS_DYN                  &  !<-- Import/export states for GFS Dynamics
                              ,IMP_GFS_PHY,EXP_GFS_PHY                  &  !<-- Import/export states for GFS Physics
                              ,IMP_GFS_WRT,EXP_GFS_WRT                     !<-- Import/export states for GFS Write
!
      TYPE(ESMF_TimeInterval),SAVE :: TIMEINTERVAL_GFS_OUTPUT              !<-- The ESMF time interval between GFS history output
!
!---------------------------------
!***  For determining clocktimes.
!---------------------------------
!
      REAL(kind=KDBL) :: atm_tim,btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_REGISTER(ATM_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the ATM gridded component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- ATM gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Return code for register
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
      RC=ESMF_SUCCESS       ! Error signal variable
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
      MESSAGE_CHECK="ATM_Register: Load principal configure file"
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
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Dynamic Core Name"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_1                          &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the ATM INITIALIZE subroutine.  Since it is just one 
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
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,ATM_INITIALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Run step of the ATM component.
!***  The GFS needs one phase while the NMM needs three.
!-----------------------------------------------------------------------
!
!-------------------------
      IF(CORE=='gfs')THEN         
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Set 1st Entry Point for ATM Run"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                       ,ESMF_SETRUN                       &  !<-- Subroutine type
                                       ,ATM_RUN                           &  !<-- The primary Dynamics / Physics /Coupler sequence
                                       ,ESMF_SINGLEPHASE                  &
                                       ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!
!-------------------------
      ELSEIF(CORE=='nmm')THEN    
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Set 1st Entry Point for the ATM Run Step"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                       ,ESMF_SETRUN                       &  !<-- Subroutine type
                                       ,ATM_RUN                           &  !<-- The primary Dynamics /Physics /Coupler sequence
                                       ,1                                 &
                                       ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Set 2nd Entry Point for the ATM Run Step"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                   &  !<-- ATM gridded component
                                       ,ESMF_SETRUN                     &  !<-- Subroutine type
                                       ,NMM_ATM_FILTERING               &  !<-- Routine to govern digital filtering each timestep
                                       ,2                               &
                                       ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Set 3rd Entry Point for the ATM Run Step"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                   &  !<-- ATM gridded component
                                       ,ESMF_SETRUN                     &  !<-- Subroutine type
                                       ,CALL_WRITE_ASYNC                &  !<-- Routine to call asynchronous output
                                       ,3                               &
                                       ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------
      ENDIF
!-----------
!
!-----------------------------------------------------------------------
!***  Register the ATM FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type
                                     ,ATM_FINALIZE                      &  !<-- User's subroutine name
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
!       WRITE(0,*)' ATM_SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)' ATM_SET_SERVICES FAILED  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_INITIALIZE(ATM_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_ATM                               &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE SETS UP FUNDAMENTAL ASPECTS OF THE MODEL RUN.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                    &  !<-- The ATM component's import state
                                          ,EXP_STATE                       !<-- The ATM component's export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component.
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT                         !<-- Return code for Initialize step
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: RC,RC_FINAL
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
!***  The following subroutines used by the different cores to 
!***  perform the Initialize step of the ATM component are internal
!***  to ATM_INITIALIZE.
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
!
        CALL NMM_ATM_INIT(ATM_GRID_COMP                                 &
                         ,IMP_STATE                                     &
                         ,EXP_STATE                                     &
                         ,CLOCK_ATM )
!
      ELSEIF(CORE=='gfs')THEN
!
        CALL GFS_ATM_INIT(ATM_GRID_COMP                                 &
                         ,CLOCK_ATM )
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'ATM INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      rtim1=rtc()
      total_integ_tim=(timef()-btim0)
!
      if(mype==0)write(0,*)' atm_init_tim=',total_integ_tim,(rtim1-rtim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_INITIALIZE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_ATM_INIT(ATM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_ATM )
!
!-----------------------------------------------------------------------
!***  THIS SUBROUTINE PERFORMS THE INITIALIZE STEP OF THE ATM COMPONENT
!***  FOR THE NMM.
!
!***  CALLED FROM SUBROUTINE ATM_INITIALIZE.
!-----------------------------------------------------------------------
!
      USE MODULE_NMM_CORE_SETUP,ONLY: NMM_SETUP 
!
      USE MODULE_NESTING,ONLY: PARENT_DATA_TO_ATM                       &
                              ,PARENT_TO_CHILD_INIT_NMM
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
!
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                    &  !<-- The ATM component's import state
                                          ,EXP_STATE                       !<-- The ATM component's export state
!
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component.
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,IRTN,N,NFCST,NTSD
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
      INTEGER(kind=KINT) :: RC,RC_INIT
!
      INTEGER(kind=KINT),DIMENSION(7) :: FCSTDATE
!
      INTEGER(ESMF_KIND_I8) :: NTSD_START                                  !<-- Timestep count (>0 for restarted runs)
!
      REAL(kind=KFPT) :: SECOND_FCST                                       !<-- Current second from restart file
!
      LOGICAL(kind=KLOG) :: INPUT_READY_FLAG                            &
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
      TYPE(ESMF_Grid) :: GRID_NMM_ATM                                      !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM ATM gridded component.
      TYPE(ESMF_Grid) :: GRID_NMM_DYN                                      !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM dynamics gridded component.
      TYPE(ESMF_Grid) :: GRID_NMM_PHY                                      !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM physics gridded component.
!
      TYPE(ESMF_Logical) :: I_AM_A_FCST_TASK,I_AM_A_PARENT
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Load configure files for all domains into memory.
!***  The file name of the uppermost domain is 'configure_file_01'
!***  and is identical to the primary file called 'configure_file'
!***  which is needed in some early parts of the setup.
!-----------------------------------------------------------------------
!
      DO N=1,99                                                            !<-- The number of config files cannot exceed 99
        CF(N)=ESMF_ConfigCreate(rc=RC)
!
        WRITE(INT_TO_CHAR,FMT)N
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file names
!
        CALL ESMF_ConfigLoadFile(config  =CF(N)                         &
                                ,filename=CONFIG_FILE_NAME              &
                                ,rc      =RC)
        IF(RC/=0)EXIT                                                      !<-- Exit loop after running out of config files
      ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate the NMM ATM component's internal state.
!-----------------------------------------------------------------------
!
      ALLOCATE(ATM_INT_STATE,stat=RC)
!
      wrap%ATM_INT_STATE=>ATM_INT_STATE
!
!-----------------------------------------------------------------------
!***  Attach the NMM ATM internal state to the ATM gridded component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach ATM Internal State to Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- Pointer to the ATM internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the ATM gridded component.
!***  Call ESMF_GridCompGet to retrieve the VM anywhere you need it.
!***  We need VM now to obtain the MPI task IDs and the local MPI
!***  communicator.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_INIT: Retrieve VM from ATM Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_INIT: Obtain Task IDs and Communicator"
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
!***  Extract this ATM component's domain ID from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain ID from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The ATM import state
                            ,name ='DOMAIN_ID'                          &  !<-- Name of the attribute to extract
                            ,value=MY_DOMAIN_ID                         &  !<-- The ID of this domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
      atm_int_state%QUILTING=QUILTING                                      !<-- Save this for the Run step
!
!-----------------------------------------------------------------------
!***  Extract this ATM component's Nest/Not-a-Nest flag
!***  from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Nest/Not-a-Nest Flag from ATM Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The ATM import state
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
        MESSAGE_CHECK="Extract Parent-Child Time Ratio from ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The ATM import state
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
      MESSAGE_CHECK="NMM_ATM_INIT: Start Time from Driver Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_ATM                            &  !<-- The ESMF Clock of this domain
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
      MESSAGE_CHECK="Set the Current Time on the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockSet(clock       =CLOCK_ATM                         &  !<-- The ATM Component's Clock
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
      MESSAGE_CHECK="ATM_Init: Extract # of tracers from Config file"
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
!***  need to be inserted into the ATM component's internal state
!***  if quilting is to be used.  See 'IF(QUILTING)THEN' below.
!-----------------------------------------------------------------------
!
      CALL NMM_SETUP(MYPE                                               &
                    ,COMM_MY_DOMAIN                                     &
                    ,CF(MY_DOMAIN_ID)                                   &
                    ,ATM_GRID_COMP                                      &
                    ,ATM_INT_STATE                                      &
                    ,GRID_NMM_ATM)
!
!-----------------------------------------------------------------------
!***  Attach the NMM-specific ESMF Grid to the ATM gridded component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the NMM ESMF Grid to the ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=ATM_GRID_COMP                      & !<-- The ATM gridded component
                           ,grid    =GRID_NMM_ATM                       & !<-- Attach the ESMF grid to the ATM component
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
!***  Dynamics the ATM component's grid.
!***  Note that this subcomponent is part of the ATM component's
!***  internal state.  This will be convenient if we need to reach
!***  the Dynamics component via the ATM component such as happens
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
      atm_int_state%DYN_GRID_COMP=ESMF_GridCompCreate(                  &
                                  name      ="Dynamics component"       &  !<-- Name of the new Dynamics gridded component
                                 ,config    =CF(MY_DOMAIN_ID)           &  !<-- Attach this configure file to the component
                                 ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
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
      CALL ESMF_GridCompSetServices(atm_int_state%DYN_GRID_COMP         &  ! <-- The Dynamics gridded component
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
      GRID_NMM_DYN=GRID_NMM_ATM                                            !<-- For now the Dyn Grid is the same as the ATM Grid
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=atm_int_state%DYN_GRID_COMP        &  !<-- The Dynamics component
                           ,grid    =GRID_NMM_DYN                       &  !<-- The Dynamics ESMF grid
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
      atm_int_state%IMP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Import" &  !<-- The Dynamics import state name
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
      atm_int_state%EXP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Export" &  !<-- The Dynamics export state name
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
!***  Physics component the ATM component's grid.
!***  Note that this subcomponent is part of the ATM component's
!***  internal state.  This will be convenient if we need to reach
!***  the Physics component via the ATM component such as happens
!***  when Write components are established.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the flag from the ATM import state indicating if the
!***  user wants physics to be active.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Physics Flag from ATM Import State"
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
      physics_1: IF(PHYSICS_ON==ESMF_True)THEN
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
        atm_int_state%PHY_GRID_COMP=ESMF_GridCompCreate(                  &
                                    name      ="Physics component"        &  !<-- Name of the new Physics gridded component
                                   ,config    =CF(MY_DOMAIN_ID)           &  !<-- Attach this configure file to the component
                                   ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
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
        CALL ESMF_GridCompSetServices(atm_int_state%PHY_GRID_COMP     &  ! <-- The Physics gridded component
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
        GRID_NMM_PHY=GRID_NMM_ATM                                          !<-- For now the Physics Grid is the same as the ATM Grid
!
        CALL ESMF_GridCompSet(gridcomp=atm_int_state%PHY_GRID_COMP      &  !<-- The NMM Physics component
                             ,grid    =GRID_NMM_PHY                     &  !<-- The ESMF grid of the Physics component
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!
      ENDIF physics_1
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
      atm_int_state%IMP_STATE_PHY=ESMF_StateCreate(stateName="Physics Import"  &  !<-- The Physics import state
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
      atm_int_state%EXP_STATE_PHY=ESMF_StateCreate(stateName="Physics Export"  &  !<-- The Physics export state
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
      atm_int_state%COUPLER_DYN_PHY_COMP=ESMF_CplCompCreate                   &
                                         (name   ="Dyn-Phy coupler component" &
                                         ,petList=atm_int_state%PETLIST_FCST  &  !<-- The forecast task IDs
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
      CALL ESMF_CplCompSetServices(atm_int_state%COUPLER_DYN_PHY_COMP &  ! <-- The Dyn-Phy coupler component
                                  ,DYN_PHY_CPL_REGISTER               &  ! <-- The user's subroutineName
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
      CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_DYN          &  !<-- The Dynamics component import state
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
      CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_PHY          &  !<-- The Dynamics component import state
                            ,name ='DOMAIN_ID'                          &  !<-- Use this name inside the state
                            ,value=MY_DOMAIN_ID                         &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Insert the flag indicating if the ATM component is a nest.
!***  The Dynamics component needs to know this regarding BC's.
!***  Both Dynamics and Physics need it to properly compute some
!***  fundamental aspects of the nested grids.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="ATM_INIT: Add Nest Flag to the Dyn Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_DYN          &  !<-- The Dynamics component import state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- Use this name inside the state
                            ,value=I_AM_A_NEST                          &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="ATM_INIT: Add Nest Flag to the Phy Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_PHY          &  !<-- The Physics component import state
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
        MESSAGE_CHECK="ATM_INIT: Add Parent-Child Time Ratio to the Dyn Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_DYN        &  !<-- The Dynamics component import state
                              ,name ='Parent-Child Time Ratio'          &  !<-- Use this name inside the state
                              ,value=PARENT_CHILD_TIME_RATIO            &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="ATM_INIT: Add Input-Ready Flag to the Dyn Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_DYN        &  !<-- The Dynamics component import state
                              ,name ='Input Ready'                      &  !<-- Use this name inside the state
                              ,value=INPUT_READY                        &  !<-- Does this nest's input file already exist?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="ATM_INIT: Add Input-Ready Flag to the Phy Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_PHY        &  !<-- The Physics component import state
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
!***  in turn part of the ATM internal state.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(QUILTING)THEN
!
        CALL WRITE_SETUP(ATM_GRID_COMP                                  &
                        ,ATM_INT_STATE                                  &
                        ,CLOCK_ATM)
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
      MESSAGE_CHECK="Initialize the NMM Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%DYN_GRID_COMP  &  !<-- The dynamics gridded component
                                  ,importState=atm_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                                  ,exportState=atm_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                                  ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                  ,phase      =ESMF_SINGLEPHASE             &
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
        MESSAGE_CHECK="Initialize Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%PHY_GRID_COMP  &  !<-- The physics gridded component
                                    ,importState=atm_int_state%IMP_STATE_PHY  &  !<-- The physics import state
                                    ,exportState=atm_int_state%EXP_STATE_PHY  &  !<-- The physics export state
                                    ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                    ,phase      =ESMF_SINGLEPHASE             &
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
      CALL ESMF_CplCompInitialize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dyn_phy coupler component
                                 ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The dyn-phy coupler import state
                                 ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The dyn-phy coupler export state
                                 ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
        CALL WRITE_INIT(ATM_GRID_COMP                                   &
                       ,ATM_INT_STATE                                   &
                       ,CLOCK_ATM)
!
          IF(MYPE>=atm_int_state%NUM_PES_FCST)THEN
            I_AM_A_FCST_TASK=ESMF_FALSE
          ENDIF
!
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Fcst-or-Write Task Flag to the ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The ATM component export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_FCST_TASK                     &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now extract number of children on this ATM component's domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Extract Number of Children from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The ATM import state
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
!***      insert them into the ATM export state since ultimately 
!***      they must be available to the parent in the
!***      Parent-Child Coupler.
!
!***  (2) Check to see if the children have input data ready for them.
!***      If not, do simple nearest neighbor and bilinear  interpolation
!***      from the parent's grid to the children's.  Write out that
!***      interpolated data into files that are waiting for the children 
!***      when they recursively execute ATM_INITIALIZE themselves.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      I_AM_A_PARENT=ESMF_FALSE
!
!-----------------------------------------------------------------------
!***  Extract from the Dynamics and Physics export state the quantities
!***  relevant to the children's boundaries.  Insert them into the 
!***  ATM export state so that ATM_DRIVER can take them and send them 
!***  to the Parent-Child coupler.  Only the forecast tasks participate
!***  in doing this since the Quilt/Write tasks never loaded data into
!***  the Dynamics or Physics export states.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_FCST_TASK==ESMF_TRUE)THEN       
!
        CALL PARENT_DATA_TO_ATM(atm_int_state%EXP_STATE_DYN             &  !<-- The Dynamics export state
                               ,atm_int_state%EXP_STATE_PHY             &  !<-- The Physics export state
                               ,EXP_STATE)                                 !<-- The ATM export state
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      child_init_block: IF(NUM_CHILDREN>0.AND..NOT.RESTARTED_RUN)THEN      !<-- Only parents participate                              
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
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The ATM import state
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
          IF(.NOT.INPUT_READY_MY_CHILD)THEN                                !<-- INPUT_READY=false -> This child has no input file 
                                                                           !      so parent will generate input.
            CALL PARENT_TO_CHILD_INIT_NMM(MYPE                          &  !<-- This task's rank (in)
                                         ,CF                            &  !<-- Array of configure files (in)
                                         ,MY_DOMAIN_ID                  &  !<-- Each domain's ID (in)
                                         ,MY_CHILDREN_ID(N)             &  !<-- The child's domain ID
                                         ,atm_int_state%DYN_GRID_COMP   &  !<-- The parent's Dynamics Component (inout)
                                         ,atm_int_state%PHY_GRID_COMP   &  !<-- The parent's Physics Component (inout)
                                         ,COMM_MY_DOMAIN )                 !<-- Each domain's intracommunicator
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
!***  Insert into the ATM export state the flag indicating if the
!***  current domain is a parent.  The ATM Driver wants to know this
!***  since most Parent-Child work can be ignored by domains with
!***  no children.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Parent/Not-a-Parent Flag to the ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The ATM component export state
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
      END SUBROUTINE NMM_ATM_INIT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_ATM_INIT(ATM_GRID_COMP                             &
                             ,CLOCK_ATM )
!
!-----------------------------------------------------------------------
!***  THIS SUBROUTINE PERFORMS THE INITIALIZE STEP OF THE ATM COMPONENT
!***  FOR THE GFS.
!
!***  CALLED FROM SUBROUTINE ATM_INITIALIZE.
!-----------------------------------------------------------------------
!
      USE MODULE_GFS_CORE_SETUP,ONLY: GFS_SETUP
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
!
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component.
!
!---------------------
!***  Local Variables
!---------------------
!
!-----------------------------------------------------------------------
      INTEGER(kind=KINT) :: MEMBER_ID,MODENS,NENS                       &
                           ,NFHOUT,NFMOUT,NFSOUT,NSOUT                  &
                           ,TOTAL_MEMBER,TOTAL_TASKS
!
      INTEGER(kind=KINT) :: RC,RC_INIT,IERR
!
      REAL(kind=KFPT) :: DELTIM
!
      CHARACTER(50) :: MODE
!
      TYPE(ESMF_Config) :: CF                                              !<-- The configure object for all GFS members
!
      TYPE(ESMF_Grid) :: GRID_GFS_DYN                                      !<-- The ESMF grid for the integration attached to
                                                                           !     the GFS dynamics gridded component.
      TYPE(ESMF_Grid) :: GRID_GFS_PHY                                      !<-- The ESMF grid for the integration attached to
                                                                           !     the GFS physics gridded component.
      TYPE(ESMF_Grid) :: GRID_GFS_ATM                                      !<-- The ESMF grid for the integration attached to
                                                                           !     the GFS ATM gridded component.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Retrieve the ESMF configure object.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve GFS Config Object from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve global VM then the total number of tasks for
!***  then entire system.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve global VM for GFS"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_VMGetGlobal(vm=VM                                       &  !<-- The virtual machine
                           ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GFS_ATM_INIT: Obtain MPI Task IDs from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm     =VM                                        &  !<-- The virtual machine
                     ,pecount=TOTAL_TASKS                               &  !<-- # of MPI tasks for entire GFS system
                     ,localpet=MYPE                                     &  !<-- Each MPI task ID
                     ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Obtain the total number of GFS ensemble members from the 
!***  configure file then create the members' IDs and names.
!-----------------------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =TOTAL_MEMBER                  &  !<-- Fill this variable 
                                  ,label ='total_member:'               &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
      NENS=TOTAL_TASKS/TOTAL_MEMBER
      MODENS=TOTAL_TASKS-NENS*TOTAL_MEMBER
!
      IF(MYPE<MODENS*(NENS+1)) THEN
        MEMBER_ID=MYPE/(NENS+1)+1
      ELSE
        MEMBER_ID=MODENS+(MYPE-MODENS*NENS)/NENS+1
      ENDIF
!
      IF(TOTAL_MEMBER==1) THEN
        ENSMEM_NAME=' '
      ELSE
        WRITE(ENSMEM_NAME,'("_",i2.2)') MEMBER_ID
      ENDIF
!
!-----------------------------------------------------------------------
!***  Model-specific routines must be invoked in order to establish
!***  the ESMF Grid.  The different integration grids necessitate
!***  different ways of setting up both the parallelism for
!***  distributed memory runs and the ESMF Grid itself.
!***  When the parallelism is constructed, the local domain limits
!***  need to be inserted into the ATM component's internal state
!***  if quilting is to be used.
!-----------------------------------------------------------------------
!
      CALL GFS_SETUP(ATM_GRID_COMP                                      &
                    ,GRID_GFS_ATM)
!
!-----------------------------------------------------------------------
!***  Now establish the frequency of forecast output.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract History Output Interval from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NFHOUT                        &  !<-- Hours between GFS history output
                                  ,label ='nfhout:'                     &  !<-- Give the variable this label's value
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NSOUT                         &  !<-- Fill this variable bel's value from the config file
                                  ,label ='nsout:'                      &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =DELTIM                        &  !<-- Fill this variable
                                  ,label ='deltim:'                     &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(NSOUT>0) THEN
        NFHOUT=INT(NSOUT*DELTIM/3600.)
        NFMOUT=INT((NSOUT*DELTIM-NFHOUT*3600.)/60.)
        NFSOUT=INT(NSOUT*DELTIM-NFHOUT*3600.-NFMOUT*60)
      ELSE
        NFMOUT=0
        NFSOUT=0
      ENDIF
      write(0,*)'nfhout=',nfhout,'nfmout=',nfmout,'nfsout=',nfsout
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set GFS History Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_GFS_OUTPUT    &  !<-- ESMF time interval between GFS history output
                               ,h           =NFHOUT                     &  !<-- Hours between GFS history output
                               ,m           =NFMOUT                     &  !<-- Minutes between GFS history output
                               ,s           =NFSOUT                     &  !<-- Seconds between GFS history output
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Read and print Chem_Registry
!-----------------------------------------------------------------------
!
      REG = Chem_RegistryCreate ( IERR )

      IF(REG%doing_gocart)THEN                                             !<-- GOCART => Chemistry on
        CHEMISTRY_ON=ESMF_True
        write(0,*)' Initialize with gocart coupling '
      ELSE                                                                 !<-- no GOCART => Chemistry off
        CHEMISTRY_ON=ESMF_False
        write(0,*)' Initialize without gocart coupling '
      ENDIF

      CALL Chem_RegistryPrint ( REG )

      CALL Chem_RegistryDestroy ( REG, IERR )
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics gridded subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Dynamics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the GFS Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GC_GFS_DYN=ESMF_GridCompCreate(name      ="dynamics component"    &
                                    ,configFile='dyn_namelist.rc'       &
                                    ,petList   =PETLIST_FCST            &
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
      MESSAGE_CHECK="Register GFS Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(GC_GFS_DYN                          &  !<-- The GFS Dynamics gridded component
                                   ,GFS_DYN_SETSERVICES                 &  !<-- The user's subroutineName for Register
                                   ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create empty Import and Export states for the Dynamics component
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_GFS_DYN=ESMF_StateCreate(statename="dynamics import"          &
                                  ,statetype=esmf_state_import          &
                                  ,rc       =RC)
!
      EXP_GFS_DYN=ESMF_StateCreate(statename="dynamics export"          &
                                  ,statetype=esmf_state_export          &
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Is this an adiabatic (no physics) run?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Physics On/Off Switch from GFS Config Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(MODE)=='.true.')THEN                                         !<-- Adiabatic => Physics off
        PHYSICS_ON=ESMF_False
        write(0,*)' Initialize without physics coupling '
      ELSE                                                                 !<-- Not adiabatic => Physics on
        PHYSICS_ON=ESMF_True
        write(0,*)' Initialize with physics coupling '
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Physics gridded subcomponent if physics is turned on.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(PHYSICS_ON==ESMF_True)THEN
!
!-------------------------------
!***  Create Physics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the GFS Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        GC_GFS_PHY=ESMF_GridCompCreate(name      ="physics component"   &
                                      ,configfile='phy_namelist.rc'     &
                                      ,petList   =PETLIST_FCST          &
                                      ,rc        =RC)
        write(0,*)'in GFS_ATM_INIT after phys comp created, petlist_fcst=',petlist_fcst
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------
!***  Register the Init, Run, and Finalize steps.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Physics Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(GC_GFS_PHY                        &
                                     ,GFS_PHY_SETSERVICES               &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Physics subcomponent.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for GFS Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_GFS_PHY=ESMF_StateCreate(statename="physics import"           &
                                  ,statetype=ESMF_STATE_IMPORT          &
                                  ,rc       =RC)
!
      EXP_GFS_PHY=ESMF_StateCreate(statename="physics export"           &
                                  ,statetype=ESMF_STATE_EXPORT          &
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
      MESSAGE_CHECK="Create the GFS Dynamics-Physics Coupler Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GC_ATM_CPL=ESMF_CplCompCreate(name   ="coupler component"         &
                                   ,petList=PETLIST_FCST                &
                                   ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------
!***  Register the Init, Run, and Finalize steps.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the Dyn-Phy Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetServices(GC_ATM_CPL                           &  !<-- The GFS Dynamics/Physics coupler component
                                  ,ATM_CPL_SETSERVICES                  &  !<-- The user's subroutine name for Register
                                  ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Will the Write components with asynchronous quilting be used?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Quilting Flag from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The GFS config object
                                  ,value =QUILTING                      &  !<-- The quilting flag
                                  ,label ='quilting:'                   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Setup the Write component(s) (which may run without quilting).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      write(0,*)'before write_setup_gfs, allocate,write_groups=',write_groups
      ALLOCATE(WRT_COMPS(WRITE_GROUPS))
      CALL WRITE_SETUP_GFS(ATM_GRID_COMP                                &
                          ,WRT_COMPS                                    &
                          ,EXP_GFS_DYN                                  &
                          ,EXP_GFS_PHY                                  &
                          ,IMP_GFS_WRT                                  &
                          ,EXP_GFS_WRT)
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
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize GFS Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =GC_GFS_DYN               &
                                  ,importstate=IMP_GFS_DYN              &
                                  ,exportstate=EXP_GFS_DYN              &
                                  ,clock      =CLOCK_ATM                &  
                                  ,phase      =ESMF_SINGLEPHASE         &
                                  ,rc         =RC)
!
      GRID_GFS_DYN=GRID_GFS_ATM                                            !<-- Use the ATM Grid for the Dynamics
!
      CALL ESMF_GridCompSet(gridcomp=GC_GFS_DYN                         &
                           ,grid    =GRID_GFS_DYN                       &
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  PHYSICS
!-------------
!
      IF(PHYSICS_ON==ESMF_True)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Initialize GFS Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompInitialize(gridcomp   =GC_GFS_PHY             &
                                    ,importstate=IMP_GFS_PHY            &
                                    ,exportstate=EXP_GFS_PHY            &
                                    ,clock      =CLOCK_ATM              &
                                    ,phase      =ESMF_SINGLEPHASE       &
                                    ,rc         =RC)
!
        GRID_GFS_PHY=GRID_GFS_ATM                                          !<-- Use the ATM Grid for the Physics
!
        CALL ESMF_GridCompSet(gridcomp=GC_GFS_PHY                       &
                             ,grid    =GRID_GFS_PHY                     &
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Initialize the Dyn-Phy Coupler subcomponent.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dyn-Phy Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompInitialize(cplcomp    =GC_ATM_CPL                &
                                 ,importstate=EXP_GFS_DYN               &
                                 ,exportstate=IMP_GFS_PHY               &
                                 ,clock      =CLOCK_ATM                 &
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the Write component(s).
!-----------------------------------------------------------------------
!
      CALL WRITE_INIT_GFS(ATM_GRID_COMP                                 &
                         ,WRT_COMPS                                     &
                         ,IMP_GFS_WRT                                   &
                         ,EXP_GFS_WRT                                   &
                         ,CLOCK_ATM                                     &
                         ,WRITE_GROUP_READY_TO_GO)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_ATM_INIT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_RUN(ATM_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  RUN THE ATM (Atmosphere) GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM Run step's import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM Run step's export
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_RUN                          !<-- Return code for the Run step
!
!---------------------
!***  Local variables
!---------------------
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
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
!
        CALL NMM_ATM_RUN(ATM_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
      ELSEIF(CORE=='gfs')THEN
!
        CALL GFS_ATM_RUN(ATM_GRID_COMP                                  &
                        ,CLOCK_ATM)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      rtim1=rtc()
      total_integ_tim=total_integ_tim+(timef()-btim0)
!     write(0,*)'exit ATM_RUN integration time ',total_integ_tim,(rtim1-rtim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_ATM_RUN(ATM_GRID_COMP                              &
                            ,IMP_STATE                                  &
                            ,CLOCK_ATM                                  &
                            ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  THE RUN STEP OF THE ATM COMPONENT FOR THE NMM.
!***  THE FORECAST TASKS EXECUTE THE RUN STEP OF THE NMM-B DYNAMICS.
!***  THIS IS THE RUN SUBROUTINE SPECIFIED IN THE DYNAMICS REGISTER
!***  ROUTINE CALLED IN ESMF_GridCompSetServices ABOVE.
!-----------------------------------------------------------------------
!
      USE MODULE_NESTING,ONLY: BOUNDARY_DATA_STATE_TO_STATE
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- ATM gridded component
!
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The ATM import state
!
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM ESMF Clock
!
      INTEGER(kind=KINT) ,INTENT(INOUT) :: RC_RUN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(ESMF_KIND_I4) :: INTEGER_DT,NUMERATOR_DT,IDENOMINATOR_DT     
!
      INTEGER(kind=KINT) :: HDIFF_ON,RC
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
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
!     ALLOCATE(ATM_INT_STATE,stat=RC)
!     wrap%ATM_INT_STATE=>ATM_INT_STATE
!
!-----------------------------------------------------------------------
!
      fcst_pes: IF(MYPE<atm_int_state%NUM_PES_FCST)THEN                    !<-- Only the forecast tasks integrate
!
!-----------------------------------------------------------------------
!***  Extract the timestep from the Clock so that we know the direction
!***  of the integration.  We skip all aspects of Physics if the time
!***  step is negative.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="ATM_Run: Extract the ESMF Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &
                          ,timeStep=DT_ESMF                             &
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="ATM_Run: Extract Components of the Timestep" 
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
!***  ATM import state to the Dynamics import state.  The Dynamics
!***  import state is not available outside of the integration time
!***  loop in NMM_INTEGRATE which lies in ATM_DRIVER therefore we
!***  need to perform this transfer each timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Horizontal Diffusion Flag from ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The ATM import state
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
        CALL ESMF_AttributeSet(state=atm_int_state%IMP_STATE_DYN        &  !<-- The Dynamics component import state
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
!***  be moved from the ATM import state to that of the Dynamics
!***  every N timesteps where N is the number of the nest's timesteps 
!***  within one timestep of its parent.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_NEST==ESMF_TRUE)THEN
!         write(6,*)' ATM_RUN calling BOUNDARY_DATA_STATE_TO_STATE'
!         call print_memory()
          CALL BOUNDARY_DATA_STATE_TO_STATE(clock    =CLOCK_ATM                    &  !<-- The ATM Clock
                                           ,ratio    =PARENT_CHILD_TIME_RATIO      &  !<-- # of child timesteps per parent timestep
                                           ,state_in =IMP_STATE                    &  !<-- ATM component's import state
                                           ,state_out=atm_int_state%IMP_STATE_DYN)    !<-- The Dynamics import state
!         write(6,*)' ATM_RUN called BOUNDARY_DATA_STATE_TO_STATE'
!         call print_memory()
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Execute the Run Step for Dynamics"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!       write(6,*)' ATM_RUN calling DYN_RUN'
!       call print_memory()
        CALL ESMF_GridCompRun(gridcomp   =atm_int_state%DYN_GRID_COMP   &  !<-- The dynamics component
                             ,importState=atm_int_state%IMP_STATE_DYN   &  !<-- The dynamics import state
                             ,exportState=atm_int_state%EXP_STATE_DYN   &  !<-- The dynamics export state
                             ,clock      =CLOCK_ATM                     &  !<-- The ATM clock
                             ,rc         =RC)        
!       write(6,*)' ATM_RUN called DYN_RUN'
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
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The Coupler import state = Dynamics export state
                              ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The Coupler export state = Physics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM clock
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
          CALL ESMF_GridCompRun(gridcomp   =atm_int_state%PHY_GRID_COMP &  !<-- The physics component
                               ,importState=atm_int_state%IMP_STATE_PHY &  !<-- The physics import state 
                               ,exportState=atm_int_state%EXP_STATE_PHY &  !<-- The physics export state
                               ,clock      =CLOCK_ATM                   &  !<-- The ATM Clock
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
          MESSAGE_CHECK="NMM_ATM_RUN: Coupler Moves Data from Physics to Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=atm_int_state%EXP_STATE_PHY        &  !<-- The Coupler import state = Physics export state
                              ,exportState=atm_int_state%IMP_STATE_DYN        &  !<-- The Coupler export state = Dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
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
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The Dynamics-Physics coupler component
                              ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The Coupler import state = Dynamics export state
                              ,exportState=atm_int_state%IMP_STATE_DYN        &  !<-- The Coupler export state = Dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
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
      atm_tim=atm_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'ATM RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_ATM_RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_ATM_RUN(ATM_GRID_COMP                              &
                            ,CLOCK_ATM)
!
!-----------------------------------------------------------------------
!***  THE RUN STEP OF THE ATM COMPONENT FOR THE GFS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- ATM gridded component
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_ATM                          !<-- The ATM ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: DFIHR                                       &
                           ,NTIMESTEP                                      !<-- The current forecast timestep (integer)
!
      INTEGER(kind=KINT) :: RC,RC_RUN                                      !<-- Error signal variables.
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF) (integer)
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION                            &  !<-- The forecast length (ESMF)
                                ,TIMESTEP                                  !<-- The fundamental timestep, seconds (ESMF)
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STARTTIME                                         !<-- The ESMF start time.
!
      TYPE(ESMF_Config) :: CF                                              !<-- The configure object for all GFS members
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the fundamental time information from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve GFS Timestep from the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &  !<-- The ESMF Clock
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,timestep    =TIMESTEP                          &  !<-- The model's timestep length
                        ,starttime   =STARTTIME                         &  !<-- The forecast start time
                        ,currtime    =CURRTIME                          &  !<-- The Clock's current time
                        ,runduration =RUNDURATION                       &  !<-- The length of the forecast
                        ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF                                             !<-- Convert timestep from ESMF to integer
!
!-----------------------------------------------------------------------
!***  We need the DFI filter duration.  Extract the configure file 
!***  from the ATM component and then obtain the filter duration.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- Tha ATM component
                           ,config  =CF                                 &  !<-- The configure object
                           ,rc      =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =DFIHR                         &  !<-- The DFI filter duration
                                  ,label ='nhours_dfini:'               &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
!-----------------------------------------------------------------------
!***  Execute the GFS forecast runstream.
!-----------------------------------------------------------------------
!
      CALL GFS_INTEGRATE(GC_GFS_DYN                                     &
                        ,GC_GFS_PHY                                     &
                        ,GC_ATM_CPL                                     &
                        ,WRT_COMPS                                      &
                        ,IMP_GFS_DYN                                    &
                        ,EXP_GFS_DYN                                    &
                        ,IMP_GFS_PHY                                    &
                        ,EXP_GFS_PHY                                    &
                        ,IMP_GFS_WRT                                    &
                        ,EXP_GFS_WRT                                    &
                        ,CLOCK_ATM                                      &
                        ,TIMEINTERVAL_GFS_OUTPUT                        &
                        ,QUILTING                                       &
                        ,WRITE_GROUP_READY_TO_GO                        &
                        ,CURRTIME                                       &
                        ,STARTTIME                                      &
                        ,NTIMESTEP                                      &
                        ,TIMESTEP                                       &
                        ,DFIHR                                          &
                        ,MYPE                                           &
                        ,PHYSICS_ON)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_ATM_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_FINALIZE(ATM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_ATM                                 &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The main ESMF Clock
      TYPE(ESMF_Config)                 :: CF                               !<-- The config object
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_FINALIZE                     !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                             !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
!
      INTEGER :: I,J
      INTEGER :: RC,RC_FINAL                                                ! The final error signal variables.
!
      CHARACTER(50):: MODE
!
      TYPE(ESMF_Logical) :: PHYSICS_ON     
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINAL   =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from ATM Component"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Dynamic Core Name"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
      IF(CORE=='nmm')THEN
        IF(TRIM(MODE)=='TRUE')THEN
          PHYSICS_ON=ESMF_False
          write(0,*)' Finalize without physics coupling. '
        ELSE
          PHYSICS_ON=ESMF_True 
          write(0,*)' Finalize with physics coupling. '
        ENDIF
      ELSEIF(CORE=='gfs')THEN
        IF(TRIM(MODE)=='.TRUE.')THEN
          PHYSICS_ON=ESMF_False
          write(0,*)' Finalize without physics coupling. '
        ELSE
          PHYSICS_ON=ESMF_True
          write(0,*)' Finalize with physics coupling. '
        ENDIF
      ENDIF


!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM GRIDDED COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
        CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                &  !<-- The ATM gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the ATM internal state
                                          ,RC)
!
        ATM_INT_STATE=>wrap%ATM_INT_STATE
      
      ENDIF 
!
!-----------------------------------------------------------------------
!***  FINALIZE EACH OF THE SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm')THEN
        CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%DYN_GRID_COMP &
                                  ,importState=atm_int_state%IMP_STATE_DYN &
                                  ,exportState=atm_int_state%EXP_STATE_DYN &
                                  ,clock      =CLOCK_ATM                   &
                                  ,rc         =RC)
!
      ELSEIF(CORE=='gfs')THEN
        CALL ESMF_GridCompFinalize(gridcomp   =gc_gfs_dyn                  &
                                 ,importstate=imp_gfs_dyn                  &
                                 ,exportstate=exp_gfs_dyn                  &
                                 ,clock      =CLOCK_ATM                    &
                                 ,rc         =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINAL=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(PHYSICS_ON==ESMF_True)THEN
!
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(CORE=='nmm')THEN
          CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%PHY_GRID_COMP &
                                    ,importState=atm_int_state%IMP_STATE_PHY &
                                    ,exportState=atm_int_state%EXP_STATE_PHY &
                                    ,clock      =CLOCK_ATM                   &
                                    ,rc         =RC)
        ELSEIF(CORE=='gfs')THEN
          CALL ESMF_GridCompFinalize(gridcomp   =GC_GFS_PHY                  &
                                    ,importstate=IMP_GFS_PHY                 &
                                    ,exportstate=EXP_GFS_PHY                 &
                                    ,clock      =CLOCK_ATM                   &
                                    ,rc         =RC)
        ENDIF
!
      ENDIF 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINAL=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  DYNAMICS-PHYSICS COUPLER
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm')THEN
        CALL ESMF_CplCompFinalize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &
                                 ,importState=atm_int_state%EXP_STATE_DYN        &
                                 ,exportState=atm_int_state%IMP_STATE_PHY        &
                                 ,clock      =CLOCK_ATM                          &
                                 ,rc         =RC)
      ELSEIF(CORE=='gfs')THEN
        CALL ESMF_CplCompFinalize(cplcomp    =GC_ATM_CPL                &
                                 ,importstate=EXP_GFS_DYN               &
                                 ,exportstate=IMP_GFS_PHY               &
                                 ,clock      =CLOCK_ATM                 &
                                 ,rc         =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY ALL STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm')THEN
        CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_DYN        &
                              ,rc   =RC)
!
        CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_DYN        &
                              ,rc   =RC)
!
        IF(PHYSICS_ON==ESMF_True)THEN
          CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_PHY      &
                                ,rc   =RC)
!
          CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_PHY      &
                                ,rc   =RC)
        ENDIF
!
      ELSEIF(CORE=='gfs')THEN
        CALL ESMF_StateDestroy(state=IMP_GFS_DYN, rc=RC)
        CALL ESMF_StateDestroy(state=EXP_GFS_DYN, rc=RC)
        IF(PHYSICS_ON==ESMF_True)THEN
          CALL ESMF_StateDestroy(state=IMP_GFS_PHY, rc=RC)
          CALL ESMF_StateDestroy(state=EXP_GFS_PHY, rc=RC)
        ENDIF
!
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  FINALIZE AND DESTROY OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
!
        IF(atm_int_state%QUILTING)THEN
          CALL WRITE_DESTROY(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
        ENDIF
!
!-----------------------------------------------------------------------
!***  DESTROY THE ATM CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Destroy ATM Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockDestroy(clock=CLOCK_ATM                          &
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!jw for GFS destroy WRT grid component
      IF (CORE=='gfs') THEN
        IF(QUILTING) THEN
          CALL WRITE_DESTROY_GFS(ATM_GRID_COMP,WRT_COMPS,               &
            IMP_GFS_WRT,EXP_GFS_WRT,CLOCK_ATM) 
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  DESTROY ALL SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm')THEN
        CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%DYN_GRID_COMP  &
                                 ,rc      =RC)
      ELSEIF(CORE=='gfs')THEN
        CALL ESMF_GridCompDestroy(gridcomp=GC_GFS_DYN                   &
                                 ,rc      =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON==ESMF_True) THEN
!
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%PHY_GRID_COMP    &
                               ,rc      =RC)
      ELSE
        CALL ESMF_GridCompDestroy(gridcomp=GC_GFS_PHY                   &
                                 ,rc      =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!------------------------------
!***  DYNAMICS-PHYSICS COUPLER
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm')THEN
        CALL ESMF_CplCompDestroy(cplcomp=atm_int_state%COUPLER_DYN_PHY_COMP &
                                ,rc     =RC)
      ELSE
        CALL ESMF_CplCompDestroy(cplcomp=GC_ATM_CPL                      &
                                ,rc     =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM FINALIZE STEP FAILED'
      ENDIF
!
      IF(PRESENT(RC_FINALIZE))THEN
        RC_FINALIZE=RC_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_ATM_FILTERING(ATM_GRID_COMP                        &
                                  ,IMP_STATE                            &
                                  ,EXP_STATE                            &
                                  ,CLOCK_ATM                            &
                                  ,RC_FILT)
!
!-----------------------------------------------------------------------
!***  PHASE 2 OF THE RUN STEP OF THE ATM (Atmosphere) GRIDDED COMPONENT.
!***  THIS PHASE IS ONLY RELEVANT WHEN DIGITAL FILTERING IS
!***  IN EFFECT AND EXECUTES AT THE END OF EACH TIMESTEP
!***  AFTER THE DYNAMICS AND PHYSICS (IN PHASE 1). 
!
!***  CALLED FROM SUBROUTINE NMM_INTEGRATE.
!-----------------------------------------------------------------------
!
      USE MODULE_DIGITAL_FILTER_NMM
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(IN)    :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_FILT                         !<-- Return code for this step
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
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The current time of Clock_ATM
                        ,STARTTIME                                         !<-- The start time of CLock_ATM
!
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
!
      TYPE(WRAP_ATM_INTERNAL_STATE) :: WRAP                                !<-- The F90 wrap of the ATM internal state
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
!***  We need to extract the ATM internal state only once.
!-----------------------------------------------------------------------
!
      IF(FIRST_PASS)THEN
!
        FIRST_PASS=.FALSE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="NMM_ATM_Filtering: Extract the ATM Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                &  !<-- The ATM component
                                          ,WRAP                         &  !<-- The F90 wrap of the ATM internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        ATM_INT_STATE=>wrap%ATM_INT_STATE
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  What are the start time and the current time?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_FILTERING: Extract StartTime,CurrentTime"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_ATM                           &
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
      MESSAGE_CHECK="NMM_ATM_FILTERING: Extract Clock Direction."
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
      fcst_tasks: IF(MYPE<atm_int_state%NUM_PES_FCST)THEN                  !<-- Only forecast tasks deal with the integration
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
          CALL DIGITAL_FILTER_DYN_INIT_NMM(atm_int_state%IMP_STATE_DYN  &
                                          ,NDFISTEP                     &
                                          ,NUM_TRACERS_MET              &
                                          ,NUM_TRACERS_CHEM)
!
          CALL DIGITAL_FILTER_PHY_INIT_NMM(atm_int_state%IMP_STATE_PHY)
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
            CALL DIGITAL_FILTER_DYN_SUM_NMM(atm_int_state%IMP_STATE_DYN &
                                           ,MEAN_ON                     &
                                           ,NUM_TRACERS_MET             &
                                           ,NUM_TRACERS_CHEM)
          ENDIF
!
!---------------------
!
          IF(CURRTIME==HALFDFITIME)THEN
!
            CALL DIGITAL_FILTER_PHY_SAVE_NMM(atm_int_state%IMP_STATE_PHY)
!
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          IF(CURRTIME==DFITIME)THEN
            write(0,*)' DFI at final dfitime'
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(atm_int_state%IMP_STATE_DYN &
                                               ,NUM_TRACERS_MET             &
                                               ,NUM_TRACERS_CHEM)
!
            CALL DIGITAL_FILTER_PHY_RESTORE_NMM(atm_int_state%IMP_STATE_PHY)
!
            CALL ESMF_ClockPrint(clock  =CLOCK_ATM                          &
                                ,options="currtime string"                  &
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
            CALL DIGITAL_FILTER_DYN_SUM_NMM(atm_int_state%IMP_STATE_DYN &
                                           ,MEAN_ON                     &
                                           ,NUM_TRACERS_MET             &
                                           ,NUM_TRACERS_CHEM)
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          IF(CURRTIME==DFITIME)THEN
            write(0,*)' DFI at final dfitime '
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(atm_int_state%IMP_STATE_DYN &
                                               ,NUM_TRACERS_MET             &
                                               ,NUM_TRACERS_CHEM)
!
! ----------------------------------------------------------------------
!
            CALL ESMF_ClockPrint(clock  =CLOCK_ATM                      &
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
      END SUBROUTINE NMM_ATM_FILTERING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CALL_WRITE_ASYNC(ATM_GRID_COMP                         &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_ATM                             &
                                 ,RC_RUN2)
!
!-----------------------------------------------------------------------
!***  PHASE 3 OF THE RUN STEP OF THE NMM ATM COMPONENT.
!***  IT INITIATES THE WRITING OF HISTORY/RESTART FILES
!***  FROM EACH ATM COMPONENT.
!-----------------------------------------------------------------------
!
      USE module_NMM_INTEGRATE,ONLY : ALARM_HISTORY                     &
                                     ,ALARM_RESTART
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(IN)    :: IMP_STATE                       !<-- The ATM Run step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM Run step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_RUN2                         !<-- Return code for the Run step 
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
!
      TYPE(WRAP_ATM_INTERNAL_STATE) :: WRAP                                !<-- The F90 wrap of the ATM internal state
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
!***  WRITE A HISTORY FILE AT THE END OF THE APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!
      IF(ESMF_AlarmIsRinging(alarm=ALARM_HISTORY                        &  !<-- The history output alarm
                            ,rc   =RC))THEN
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Run2: Retrieve ATM Component's Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                &  !<-- The ATM gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the ATM internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        ATM_INT_STATE=>wrap%ATM_INT_STATE
!
!-----------------------------------------------------------------------
!***  EXECUTE THE WRITING OF A HISTORY FILE.
!-----------------------------------------------------------------------
!
!     write(6,*)' calling WRITE_ASYNC'
!     call print_memory()
        IF(atm_int_state%QUILTING)THEN
          CWRT='History'
          CALL WRITE_ASYNC(ATM_GRID_COMP                                &
                          ,ATM_INT_STATE                                &
                          ,CLOCK_ATM                                    &
                          ,MYPE                                         &
                          ,CWRT)
        ENDIF
!     write(6,*)' called WRITE_ASYNC'
!     call print_memory()
!
      ENDIF 
!
!-----------------------------------------------------------------------
!***  WRITE A RESTART FILE AT THE END OF THE APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!
      IF(ESMF_AlarmIsRinging(alarm=ALARM_RESTART                        &  !<-- The restart output alarm
                            ,rc   =RC))THEN
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Run2: Retrieve ATM Component's Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                &  !<-- The ATM gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the ATM internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        ATM_INT_STATE=>wrap%ATM_INT_STATE
!
!-----------------------------------------------------------------------
!***  EXECUTE THE WRITING OF A RESTART FILE.
!-----------------------------------------------------------------------
!
!     write(6,*)' calling WRITE_ASYNC for restart'
!     call print_memory()
        IF(atm_int_state%QUILTING)THEN
          CWRT='Restart'
          CALL WRITE_ASYNC(ATM_GRID_COMP                                &
                          ,ATM_INT_STATE                                &
                          ,CLOCK_ATM                                    &
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
      END MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
