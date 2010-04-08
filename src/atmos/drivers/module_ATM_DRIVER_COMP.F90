!-----------------------------------------------------------------------
!
      MODULE MODULE_ATM_DRIVER_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE ATM (Atmosphere) DRIVER MODULE.
!***  IT WILL SET UP THE ATM GRIDDED COMPONENTS THEN
!***  RUN THEIR INITIALIZE, RUN, AND FINALIZE ROUTINES.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!               Black - Original setup code for the NMM-B
!               Juang - Original setup code for the GFS
!   2008-04-22  Black/Henderson - Created module to more easily handle
!                                 nested domains.
!   2008-08-27  Black - Incorporated nested domain boundary data.
!   2009-07-30  Black - Merged into NEMS trunk code.
!   2009-11-03  Weiyu - Modified for the ensemble GEFS.
!   2009-11-30  Wang  - removed writing "core" into atm_namelist.rc
!
! USAGE: ATM_DRIVER parts called from MAIN_ESMF.F90
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_INCLUDE
!
      USE MODULE_ATM_DRIVER_INTERNAL_STATE,ONLY:                        &
                                     ATM_DRIVER_INTERNAL_STATE          &
                                    ,WRAP_ATM_DRIVER_INTERNAL_STATE
!
      USE MODULE_NMM_INTEGRATE,ONLY: NMM_INTEGRATE
!
      USE MODULE_ATM_GRID_COMP,ONLY: ATM_REGISTER                          !<-- The Register (or SetServices) routine for ATM_GRID_COMP
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!***  MODULES USED BY ATM_DRIVER ONLY FOR THE NMM.
!-----------------------------------------------------------------------
!
      USE MODULE_NESTING,ONLY: PARENT_CHILD_COMMS
!
      USE MODULE_PARENT_CHILD_CPL_COMP,ONLY: PARENT_CHILD_CPL_REGISTER  &  !<-- The Register routine for PARENT_CHILD Coupler
                                            ,PARENT_CHILD_COUPLER_SETUP
!
      USE MODULE_CONTROL,ONLY: TIMEF
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
      PUBLIC :: ATM_DRIVER_REGISTER
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
      CHARACTER(3) :: CORE                                                 !<-- The name of the selected dynamic core
!
      CHARACTER(ESMF_MAXSTR) :: CLOCK_ATM_DRV_NAME                         !<-- Name of the ESMF ATM Driver Clock
!
!
      LOGICAL(kind=KLOG) :: RESTARTED_RUN                                  !<-- Flag indicating if this is a restarted run
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
      TYPE(ESMF_Clock),DIMENSION(:),ALLOCATABLE :: CLOCK_ATM_DRV           !<-- The ATM_DRIVER ESMF Clocks
!
      TYPE(ESMF_Config),SAVE :: CF_DRIVER                                  !<-- The lead configure object for the ATM Driver
!
      TYPE(WRAP_ATM_DRIVER_INTERNAL_STATE) :: WRAP                         !<-- The F90 wrap of the ATM internal state
!
      TYPE(ATM_DRIVER_INTERNAL_STATE),POINTER :: ATM_DRV_INT_STATE         !<-- The ATM_DRIVER internal state pointer
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
                                                  ,PETLIST_ATM             !<-- List of task IDs for each domain (ATM Component)
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
!------------------------------
!***  For GFS Ensemble Members
!------------------------------
!
      INTEGER(kind=KINT)  :: MEMBER_ID, TOTAL_MEMBER, MYPE_GLOBAL
      TYPE(ESMF_VM), SAVE :: VM_GLOBAL
!
      INTEGER(kind=KINT),DIMENSION(:),  ALLOCATABLE :: PE_MEMBER           !<-- Tasks for each member
      INTEGER(kind=KINT),DIMENSION(:,:),ALLOCATABLE :: PETLIST             !<-- Task list for each member
!
      CHARACTER(ESMF_MAXSTR) :: IMPSTATENAME                            &  !<-- Import state name of the ATM components
                               ,EXPSTATENAME                               !<-- Export state name of the ATM components
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: GRIDCOMPNAME      !<-- Names of each member's ATM component
!
      TYPE(ESMF_GridComp),DIMENSION(:),ALLOCATABLE :: GC_ATM               !<-- ATM components for each member
!
      TYPE(ESMF_State),SAVE :: IMP_ATM                                  &  ! Import state for GFS ATM components
                              ,EXP_ATM                                     ! Export state for GFS ATM components
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
      SUBROUTINE ATM_DRIVER_REGISTER(ATM_DRIVER_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the ATM Driver component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_DRIVER_COMP              !<-- ATM Driver component
!
      INTEGER,INTENT(OUT)               :: RC_REG                       !<-- Return code for register
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
      RC    =ESMF_SUCCESS   ! Error signal variable
      RC_REG=ESMF_SUCCESS   ! Error signal variable

!-----------------------------------------------------------------------
!***  Register the ATM Driver Initialize subroutine.  Since it is just
!***  one subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM_Driver Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_DRIVER_COMP                   &  !<-- ATM Driver component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,ATM_DRIVER_INITIALIZE             &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the ATM_DRIVER RUN subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM_Driver Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_DRIVER_COMP                   &  !<-- ATM Driver component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type
                                     ,ATM_DRIVER_RUN                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the ATM FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM_DRIVER Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_DRIVER_COMP                   &  !<-- ATM Driver component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type
                                     ,ATM_DRIVER_FINALIZE               &  !<-- User's subroutine name
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
!       WRITE(0,*)' ATM_DRIVER_REGISTER SUCCEEDED'
      ELSE
        WRITE(0,*)' ATM_DRIVER_REGISTER FAILED  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_DRIVER_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_DRIVER_INITIALIZE(ATM_DRIVER_COMP                  &
                                      ,IMP_STATE                        &
                                      ,EXP_STATE                        &
                                      ,CLOCK_MAIN                       &
                                      ,RC_INIT_DRV)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CREATES THE ATM GRIDDED COMPONENTS AND EXECUTES
!***  THEIR INITIALIZE STEP.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_DRIVER_COMP                !<-- The ATM_DRIVER component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE,EXP_STATE            !<-- The ATM Driver import/export states
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_MAIN                     !<-- The main program's ESMF Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT_DRV                    !<-- Return code for ATM_DRIVER Initialize step
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(ESMF_LOGICAL)      :: Cpl_flag  
      INTEGER                 :: RC,RC_FINAL
!
      LOGICAL :: ATM_INIT_SUCCESS=.FALSE.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
      RC         =ESMF_SUCCESS
      RC_INIT_DRV=ESMF_SUCCESS
!
!     ntimestep=0
!
!-----------------------------------------------------------------------
!***  Allocate the ATM_DRIVER component's internal state.
!-----------------------------------------------------------------------
!
      ALLOCATE(ATM_DRV_INT_STATE,stat=RC)
!
      wrap%ATM_DRV_INT_STATE=>ATM_DRV_INT_STATE
!
!-----------------------------------------------------------------------
!***  Attach the ATM_DRIVER internal state to the ATM_DRIVER component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Attach ATM Driver Internal State to the Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetInternalState(ATM_DRIVER_COMP                &  !<-- The ATM_DRIVER component
                                        ,WRAP                           &  !<-- Pointer to the ATM_DRIVER internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the ATM_DRIVER component.
!***  We need VM now to obtain the MPI task IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Retrieve VM from ATM_DRIVER Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompGet(gridcomp=ATM_DRIVER_COMP                    &  !<-- The ATM_DRIVER component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  The different cores handle configure files differently
!***  but to begin we will read the core name from the leading
!***  configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Lead Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_DRIVER=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config  =CF_DRIVER                       & !<-- The lead configure object for ATM Driver
                              ,filename='configure_file'                & !<-- The name of the configure file
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Core Name from the Lead Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_DRIVER                     &  !<-- The lead config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now each dynamic core performs its general preliminary work
!***  including the creation of the ATM gridded components.
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Obtain MPI Task IDs from VM"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                       ,localpet=MYPE                                     &  !<-- Each MPI task ID
                       ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL NMM_ATM_DRIVER_INIT(CLOCK_MAIN)
!
      ELSEIF(CORE=='gfs')THEN
!
        CALL GFS_ATM_DRIVER_INIT(CLOCK_MAIN,     &
                                 IMP_STATE,      &
                                 EXP_STATE)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Write the final error signal.
!-----------------------------------------------------------------------
!
      IF(RC_INIT_DRV==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM_DRIVER INITIALIZE STEP SUCCEEDED'
        ATM_INIT_SUCCESS=.TRUE.
      ELSE
        WRITE(0,*)'ATM_DRIVER INITIALIZE STEP FAILED RC_INIT_DRV=',RC_INIT_DRV
      ENDIF
!
!-----------------------------------------------------------------------
!
      atm_drv_init=(timef()-btim0)
      write(0,*)' Clocktime: ATM_DRIVER_INITIALIZE=',atm_drv_init*1.e-3
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_DRIVER_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_DRIVER_RUN(ATM_DRIVER_COMP                         &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_MAIN                              &
                               ,RC_RUN_DRV)
!
!-----------------------------------------------------------------------
!***  RUN THE ATM_DRIVER COMPONENT.  THIS ROUTINE EXECUTES THE
!***  RUN STEPS OF ALL ATM COMPONENTS.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_DRIVER_COMP                 !<-- The ATM_DRIVER component
      TYPE(ESMF_State),   INTENT(IN)    :: IMP_STATE                       !<-- The ATM_DRIVER Run step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM_DRIVER Run step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_RUN_DRV                      !<-- Return code for the ATM_DRIVER Run step 
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: RC                                             !<-- Error signal variable.
      INTEGER(kind=KINT) :: ID_ATM,N
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF_INT)
!
      CHARACTER(2) :: INT_TO_CHAR
      CHARACTER(6) :: FMT='(I2.2)'
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_RUN_DRV=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Call the appropriate subroutine to execute the given core's
!***  Run step of the ATM_DRIVER component.
!
!***  The main integration runstream of the NMM appears in subroutine
!***  NMM_INTEGRATE which is called by NMM_ATM_DRIVER_RUN in the
!***  ATM_DRIVER component.
!
!***  The main integration runstream of the GFS appears in subroutine
!***  GFS_INTEGRATE which is called by GFS_ATM_RUN in the ATM component.
!-----------------------------------------------------------------------
!
      IF(CORE=='nmm')THEN
!
        CALL NMM_ATM_DRIVER_RUN
!
      ELSEIF(CORE=='gfs')THEN
!
        CALL GFS_ATM_DRIVER_RUN(CLOCK_MAIN)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_DRIVER_RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_DRIVER_FINALIZE(ATM_DRIVER_COMP                    &
                                    ,IMP_STATE                          &
                                    ,EXP_STATE                          &
                                    ,CLOCK_MAIN                         &
                                    ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_DRIVER_COMP                 !<-- The ATM_DRIVER component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_MAIN                      !<-- The main ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_FINALIZE                     !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,N
      INTEGER :: RC,RC_FINAL_DRV                                           ! The final error signal variables.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_FINAL_DRV=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  We want everyone present before finalizing and destroying objects.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   MESSAGE_CHECK=" Initial Barrier in ATM_DRIVER_FINALIZE"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!  
!!!   CALL ESMF_VMBarrier(vm=VM                                         &
!!!                      ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Retrieve the ATM_DRIVER component's internal state.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(ATM_DRIVER_COMP                &  !<-- The ATM_DRIVER component
                                        ,WRAP                           &  !<-- The F90 wrap of the ATM internal state
                                        ,RC)
!
      ATM_DRV_INT_STATE=>wrap%ATM_DRV_INT_STATE
!
!-----------------------------------------------------------------------
!***  Finalize the ATM subcomponents.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!     MESSAGE_CHECK="Finalize ATM Components"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!  
!!!     CALL ESMF_GridCompFinalize(gridcomp   =atm_drv_int_state%ATM_GRID_COMP(N) &
!!!                               ,importState=atm_drv_int_state%IMP_STATE_ATM(N) &
!!!                               ,exportState=atm_drv_int_state%EXP_STATE_ATM(N) & 
!!!                               ,clock      =CLOCK_ATM_DRV(N)                   &
!!!                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Destroy ATM subcomponents.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!     MESSAGE_CHECK="Destroy ATM Components"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!!!     CALL ESMF_GridCompDestroy(gridcomp=atm_drv_int_state%ATM_GRID_COMP(N) &
!!!                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Destroy the ATM_DRIVER clocks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!     MESSAGE_CHECK="Destroy ATM_DRIVER Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!     CALL ESMF_ClockDestroy(clock=CLOCK_ATM_DRV(N)                   &
!!!                           ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Finalize the Parent-Child coupler if it exists.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   MESSAGE_CHECK="Finalize Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!  
!!!   CALL ESMF_CplCompFinalize(cplcomp     =PARENT_CHILD_COUPLER_COMP    &  !<-- The Parent-Child Coupler Component
!!!                             ,importState=IMP_STATE_CPL_NEST           &  !<-- The Nesting coupler import state
!!!                             ,exportState=EXP_STATE_CPL_NEST           &  !<-- The Nesting coupler export state
!!!                             ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!
!-----------------------------------------------------------------------
!***  Destroy all states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   MESSAGE_CHECK="Destroy ATM_DRIVER States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!!!   CALL ESMF_StateDestroy(state=atm_drv_int_state%IMP_STATE_ATM(N)   &
!!!                         ,rc   =RC)
!
!!!   CALL ESMF_StateDestroy(state=atm_drv_int_state%EXP_STATE_ATM(N)   &
!!!                         ,rc   =RC)
!
!!!   CALL ESMF_StateDestroy(state=IMP_STATE_CPL_NEST                   &
!!!                         ,rc   =RC)
!
!!!   CALL ESMF_StateDestroy(state=EXP_STATE_CPL_NEST                   &
!!!                         ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_FINAL_DRV==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM_DRIVER FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM_DRIVER FINALIZE STEP FAILED RC_FINAL_DRC=',RC_FINAL_DRV
      ENDIF
!
!     IF(PRESENT(RC_FINALIZE))THEN
        RC_FINALIZE=RC_FINAL_DRV
!     ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_DRIVER_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_ATM_DRIVER_INIT(CLOCK_MAIN)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CREATES CONFIGURE OBJECTS, CLOCKS, AND THE
!***  INDIVIDUAL ATM COMPONENTS FOR ALL NMM DOMAINS.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_MAIN                         !<-- The main ESMF Clock
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_DOM,ID_X,ISTAT,N
!
      INTEGER(kind=KINT) :: NHOURS_FCST                                 &  !<-- Length of forecast in hours
                           ,NHOURS_HISTORY                              &  !<-- Hours between history output
                           ,NHOURS_RESTART                              &  !<-- Hours between restart output
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
      CHARACTER(ESMF_MAXSTR) :: ATM_GRID_COMP_BASE='ATM Gridded Component ' &
                               ,ATM_GRID_COMP_NAME 
!
      TYPE(ESMF_TimeInterval) :: TIMEINTERVAL_RECV_FROM_PARENT             !<-- ESMF time interval between Recv times from parent
!
      TYPE(ESMF_Logical) :: PHYSICS_ON                                     !<-- Does the integration include physics?
!
      TYPE(ESMF_Config) :: CF_X                                            !<-- The config objects (one per domain)
!
      INTEGER(kind=KINT) :: RC,RC_NMM_DRV_INIT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC             =ESMF_SUCCESS
      RC_NMM_DRV_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Create and load all of the configure objects.  All domains
!***  are functionally equivalent thus each has its own configure
!***  file.  We are counting the configure files as we create the
!***  ESMF configure objects thus we will know how many different
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
          MESSAGE_CHECK="Create Temporary Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF_X=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Load the Temp Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF_X                        &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Domain ID From Temp Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_X                      &  !<-- The config object
                                      ,value =ID_X                      &  !<-- The domain's ID
                                      ,label ='my_domain_id:'           &  !<-- Take value from this config labelious variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF(ID_X)=ESMF_ConfigCreate(rc=RC)                                !<-- Domain's ID is its element in the CF array
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Destroy Temporary Config Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigDestroy(config=CF_X                           &  !<-- The temporary config object
                                 ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Load the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(ID_X)                    &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DOMAIN_ID_TO_RANK(ID_X)=N                                         !<-- The configure file rank for a given domain ID 
          RANK_TO_DOMAIN_ID(N)=ID_X                                         !<-- The domain ID for a given configure file rank
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
!***    (2) Create an ATM subcomponent for all domains;
!***    (3) Call ATM_INIT recursively for all domains.
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
                               ,PETLIST_ATM )                              !<-- List of task IDs for each domain (ATM Component) (out)
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
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract INPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =INPES                       &  !<-- The domain's fcst tasks in I
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract JNPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =JNPES                       &  !<-- The domain's fcst tasks in J
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Write_Groups From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_GROUPS                &  !<-- The number of Write groups on this domain
                                    ,label ='write_groups:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Write Tasks Per Group From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_TASKS_PER_GROUP       &  !<-- The number of tasks per Write group
                                    ,label ='write_tasks_per_group:'    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
        ALLOCATE(PETLIST_ATM(1:N_TASKS,1))
!
        DO N=1,N_TASKS
          PETLIST_ATM(N,1)=N-1                                             !<-- The list of task IDs for the ATM Component
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
!***  Allocate the ATM import/export states.
!-----------------------------------------------------------------------
!
      ALLOCATE(atm_drv_int_state%IMP_STATE_ATM(1:NUM_DOMAINS),stat=ISTAT)
      ALLOCATE(atm_drv_int_state%EXP_STATE_ATM(1:NUM_DOMAINS),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  Create the ATM import/export states.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        atm_drv_int_state%IMP_STATE_ATM(ID_DOM)=ESMF_StateCreate(        &  !<-- ATM import state
                                           statename='ATM Import State'  &  !<-- ATM import state name
                                          ,statetype= ESMF_STATE_IMPORT  &
                                          ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        atm_drv_int_state%EXP_STATE_ATM(N)=ESMF_StateCreate(             &  !<-- ATM export state
                                           statename='ATM Export State'  &  !<-- ATM export state name
                                          ,statetype= ESMF_STATE_EXPORT  &
                                          ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Restart Flag from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                  ,value =RESTARTED_RUN               &  !<-- Logical flag indicating if this is a restarted run
                                  ,label ='restart:'                  &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Each task will create Clocks for all domains for simplicity
!***  in executing the major DO loops over the ATM components.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create the domains' clocks with their timesteps, start times,
!***  and run durations.  
!-----------------------------------------------------------------------
!
      ALLOCATE(CLOCK_ATM_DRV(1:NUM_DOMAINS))
      ALLOCATE(TIMESTEP     (1:NUM_DOMAINS))
      ALLOCATE(DT           (1:NUM_DOMAINS))
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
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Timestep from Config File"
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Establish the timesteps for all of the domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Set Timestep Interval"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP(ID_DOM)         &  !<-- The fundamental timestep on domain N (sec) (ESMF)
                                 ,s           =TIMESTEP_SEC_WHOLE       &
                                 ,sn          =TIMESTEP_SEC_NUMERATOR   &
                                 ,sd          =TIMESTEP_SEC_DENOMINATOR &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
                                    ,value =NHOURS_HISTORY              &  !<-- Fill this variable
                                    ,label ='nhours_history:'           &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
                                 ,h           =NHOURS_HISTORY           &  !<-- Hours between history
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
                                    ,value =NHOURS_RESTART              &  !<-- Fill this variable
                                    ,label ='nhours_restart:'           &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
                                 ,h           =NHOURS_RESTART           &  !<-- Hours between restart output (integer)
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO timeinfo_loop
!
!-----------------------------------------------------------------------
!***  Obtain the forecast start time from the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Start Time from Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_MAIN                         &  !<-- The Main ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time (ESMF)
!!!                     ,runDuration=RUNDURATION                        &  !<-- The simulation run duration (ESMF)
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Forecast Length from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &
                                    ,value =NHOURS_FCST                 &
                                    ,label ='nhours_fcst:'              &
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NSECONDS_FCST=NHOURS_FCST*3600                                     !<-- The forecast length (sec) (REAL)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Set the Forecast Length"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION              &  !<-- The forecast length (sec) (ESMF)
                                 ,s           =NSECONDS_FCST            &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  With data from above, create the ESMF Clocks to control
!***  the timestepping within the ATM Driver component(s).
!***  Each domain will set its own clock in the initialize 
!***  step of ATM_GRID_COMP.
!-----------------------------------------------------------------------
!
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        CLOCK_ATM_DRV_NAME='CLOCK_ATM_DRV_'//INT_TO_CHAR 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clock for the ATM DRIVER Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CLOCK_ATM_DRV(ID_DOM)=ESMF_ClockCreate(name       =CLOCK_ATM_DRV_NAME  &  !<-- The ATM_DRIVER Clock's name
                                              ,timeStep   =TIMESTEP(ID_DOM)    &  !<-- The fundamental timestep in this component 
                                              ,startTime  =STARTTIME           &  !<-- Start time of simulation
                                              ,runDuration=RUNDURATION         &  !<-- Duration of simulation
                                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDDO clock_loop
!
!-----------------------------------------------------------------------
!***  Allocate the ATM gridded component(s).
!-----------------------------------------------------------------------
!
      ALLOCATE(atm_drv_int_state%ATM_GRID_COMP(1:NUM_DOMAINS),stat=ISTAT)
!
      IF(ISTAT/=0)THEN
        WRITE(0,*)' ERROR: Failed to allocate ATM_GRID_COMP'
        WRITE(6,*)' ERROR: Failed to allocate ATM_GRID_COMP'
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the ATM gridded component(s).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      atm_comp_create: DO N=1,NUM_DOMAINS
!
!-----------------------------------------------------------------------
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
        WRITE(INT_TO_CHAR,FMT)ID_DOM 
        ATM_GRID_COMP_NAME=ATM_GRID_COMP_BASE//INT_TO_CHAR                 !<-- Append domain ID to ATM Grid Comp name
        N_TASKS=NTASKS_DOMAIN(ID_DOM)                                      !<-- # of tasks on this domain
        PETLIST=>PETLIST_ATM(1:N_TASKS,N)                                  !<-- The PETlist for this domain
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="ATM_DRIVER: Create ATM_GRID_COMP"//INT_TO_CHAR
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        atm_drv_int_state%ATM_GRID_COMP(ID_DOM)=ESMF_GridCompCreate(    &  !<-- The ATM Component for this domain
                                            name   =ATM_GRID_COMP_NAME  &  !<-- Name of the new ATM gridded component
                                           ,config =CF(N)               &  !<-- This domain's configure file
                                           ,petList=PETLIST             &  !<-- The IDs of tasks that will run on this domain
                                           ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Register the ATM components' Init, Run, Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register ATM Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(atm_drv_int_state%ATM_GRID_COMP(ID_DOM)  &  !<-- The ATM gridded component
                                     ,ATM_REGISTER                             &  !<-- User's subroutineName
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(ID_DOM/=MY_DOMAIN_ID)CYCLE                                      !<-- Only need to load Import State properly for my domain
!
!-----------------------------------------------------------------------
!***  Check the configure flag indicating whether or not to run
!***  adiabatically (i.e., with no physics).  Insert the flag
!***  into the ATM import state.
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
        MESSAGE_CHECK="Add Physics flag to the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component's import state
                              ,name ='PHYSICS_ON'                                  &  !<-- The flag indicating if physics is active
                              ,value=PHYSICS_ON                                    &  !<-- The value being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the maximum number of domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add MAX_DOMAINS to the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component's import state
                              ,name ='MAX_DOMAINS'                                 &  !<-- Maximum # of domains
                              ,value=MAX_DOMAINS                                   &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the domain IDs into the ATM import state(s) along with
!***  the number of children and the children's domain IDs.
!***  Also insert a flag as to whether the ATM component is a nest.
!
!***  Note that all tasks are aware of all domains' IDs, 
!***  number of children, and those children's domain IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Domain IDs to the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component's import state
                              ,name ='DOMAIN_ID'                                   &  !<-- This ATM Component's domain ID
                              ,value=N                                             &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
        CALL ESMF_AttributeSet(state    =atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component's import state
                              ,name     ='DOMAIN_ID_TO_RANK'                           &  !<-- Adding Attribute with this name
                              ,count    =MAX_DOMAINS                                   &  !<-- Total # of domains
                              ,valueList=DOMAIN_ID_TO_RANK                             &  !<-- Configure file IDs linked to each domain
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Number of Children to the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component's import state
                              ,name ='NUM_CHILDREN'                                &  !<-- This ATM Component's # of children
                              ,value=NUM_CHILDREN(MY_DOMAIN_ID)                    &  !<-- Insert this into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(COMM_TO_MY_PARENT<0)THEN
          I_AM_A_NEST=ESMF_FALSE
        ELSE
          I_AM_A_NEST=ESMF_TRUE
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Nest/Not-a-Nest Flag to the ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component import state
                              ,name ='I-Am-A-Nest Flag'                            &  !<-- Name of Attribute
                              ,value=I_AM_A_NEST                                   &  !<-- Logical nest flag
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nesting_block_2: IF(NESTING_NMM)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Domain IDs of Children to the ATM Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          LENGTH=MAX(1,NUM_CHILDREN(N))
          CHILD_ID=>ID_CHILDREN(1:LENGTH,N)                                !<-- Select only the IDs of this Component's children
!
          CALL ESMF_AttributeSet(state    =atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component import state
                                ,name     ='CHILD_IDs'                                   &  !<-- The children's IDs of this ATM Component
                                ,count    =LENGTH                                        &  !<-- Length of inserted array
                                ,valueList=CHILD_ID                                      &  !<-- Insert this into the import state
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I_AM_A_NEST==ESMF_TRUE)THEN
!
            PARENT_CHILD_TIME_RATIO=NINT(DT(ID_PARENTS(ID_DOM))/DT(ID_DOM))      !<-- Ratio of parent's timestep to this nest's
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Add Parent-Child Time Ratio to ATM Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID) &  !<-- This ATM component import state
                                  ,name ='Parent-Child Time Ratio'                     &  !<-- Name of Attribute
                                  ,value=PARENT_CHILD_TIME_RATIO                       &  !<-- # of child timesteps per parent timestep
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
          ENDIF
!
        ENDIF nesting_block_2
!
!-----------------------------------------------------------------------
!
      ENDDO atm_comp_create
!
!-----------------------------------------------------------------------
!***  At this point, ATM components for each domain have been created 
!***  and registered.  Now they need to be initialized. 
!
!***  The following call will initialize ATM_GRID_COMP for domain #1.
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
      CALL CALL_ATM_INITIALIZE(1,CLOCK_ATM_DRV)                            !<-- Initiate cascade of ATM Initialize calls for all domains
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
      MESSAGE_CHECK="ATM_DRIVER_INIT: Extract Fcst-or-Write Flag from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=atm_drv_int_state%EXP_STATE_ATM(MY_DOMAIN_ID)   &  !<-- The ATM component export state
                            ,name ='Fcst-or-Write Flag'                            &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_FCST_TASK                                &  !<-- Am I a forecast task?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
        CALL PARENT_CHILD_COUPLER_SETUP(NUM_DOMAINS                                   &  !      
                                       ,MY_DOMAIN_ID                                  &  !      
                                       ,NUM_CHILDREN(MY_DOMAIN_ID)                    &  !      
                                       ,COMM_TO_MY_CHILDREN                           &  !      
                                       ,COMM_TO_MY_PARENT                             &  !
                                       ,COMM_MY_DOMAIN                                &  !
                                       ,DT                                            &  !  
                                       ,CHILD_ID                                      &  !     ^
                                       ,atm_drv_int_state%EXP_STATE_ATM(MY_DOMAIN_ID) &  !     |  
                                       ,FTASKS_DOMAIN                                 &  !     |  
                                       ,ID_PARENTS                                    &  !     |
                                       ,DOMAIN_ID_TO_RANK                             &  !     |
                                       ,MAX_DOMAINS                                   &  !   Input
!                                                                                          ----------
                                       ,IMP_STATE_CPL_NEST                            &  !   Output
                                       ,EXP_STATE_CPL_NEST                            &  !     |
                                       ,PARENT_CHILD_COUPLER_COMP )                      !     v
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
                                     ,clock      =CLOCK_ATM_DRV(MY_DOMAIN_ID) &  !<-- The ATM Clock
                                     ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
!!!     IF(MY_DOMAIN_ID==1)THEN
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
                       ,clock            =CLOCK_ATM_DRV(MY_DOMAIN_ID)    &  !<-- Each domain's ATM Driver Clock
                       ,ringTime         =STARTTIME                      &  !<-- First time the Alarm rings (ESMF)
                       ,ringInterval     =TIMEINTERVAL_RECV_FROM_PARENT  &  !<-- Recv from my parent at this frequency (ESMF)
                       ,ringTimeStepCount=1                              &  !<-- The Alarm rings for this many timesteps
                       ,sticky           =.false.                        &  !<-- Alarm does not ring until turned off
                       ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Clocktime Output Alarm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALARM_CLOCKTIME=ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'           &
                                      ,clock            =CLOCK_ATM_DRV(MY_DOMAIN_ID) &  !<-- ATM Clock
                                      ,ringInterval     =INTERVAL_CLOCKTIME          &  !<-- Time interval between clocktime prints (ESMF)
                                      ,ringTimeStepCount=1                           &  !<-- The Alarm rings for this many timesteps
                                      ,sticky           =.false.                     &  !<-- Alarm does not ring until turned off
                                      ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_DRV_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_ATM_DRIVER_INIT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_ATM_DRIVER_RUN
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE EXECUTES THE INTEGRATION TIMELOOP FOR THE NMM
!***  THROUGH A CALL TO SUBROUTINE NMM_INTEGRATE.
!***  THAT IS PRECEDED BY DIGITAL FILTERING IF IT IS REQUESTED.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: FILTER_METHOD,HDIFF_ON,MYPE_LOCAL,NTIMESTEP
!
      INTEGER(kind=KINT) :: RC,RC_NMM_ATM_DRIVER_RUN
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      TYPE(ESMF_Time) :: CURRTIME
!
      TYPE(ESMF_State) :: IMP_STATE_ATM                                 &
                         ,EXP_STATE_ATM 
!
      TYPE(ESMF_GridComp) :: ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
      RC                   =ESMF_SUCCESS
      RC_NMM_ATM_DRIVER_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      ATM_GRID_COMP=atm_drv_int_state%ATM_GRID_COMP(MY_DOMAIN_ID)
      IMP_STATE_ATM=atm_drv_int_state%IMP_STATE_ATM(MY_DOMAIN_ID)
      EXP_STATE_ATM=atm_drv_int_state%EXP_STATE_ATM(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!***  Obtain current information from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Driver_Run: Get current time info from the Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM_DRV(MY_DOMAIN_ID)       &
                        ,starttime   =STARTTIME                         &
                        ,currtime    =CURRTIME                          &
                        ,advanceCount=NTIMESTEP_ESMF                    &
                        ,runduration =RUNDURATION                       &
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CURRTIME=STARTTIME
      NTIMESTEP=NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***  We need the local MPI task ID on the given NMM domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_DRIVER_RUN: Retrieve VM from ATM Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_ATM_DRIVER_RUN: Obtain the Local Task ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine for current ATM component
                     ,localpet       =MYPE_LOCAL                        &  !<-- Each MPI task ID
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Driver_Run: Put Filter Method into ATM import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_ATM                        &  !<-- This ATM component's import state for filter
                            ,name ='Filter_Method'                      &  !<-- Flag for type of digital filter
                            ,value=FILTER_METHOD                        &  !<-- Value of digital filter flag
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set default value for horizontal diffusion flag (1-->ON).
!-----------------------------------------------------------------------
!
      HDIFF_ON=1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Driver_Run: Put Horizontal Diffusion Flag intp ATM import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_ATM                        &  !<-- This ATM component's import state for horiz diff
                            ,name ='HDIFF'                              &  !<-- Flag for diffusion on/off
                            ,value=HDIFF_ON                             &  !<-- Value of horizontal diffusion flag
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
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
      CALL NMM_INTEGRATE(atm_grid_comp     =ATM_GRID_COMP                  &
                        ,imp_state_atm     =IMP_STATE_ATM                  &
                        ,exp_state_atm     =EXP_STATE_ATM                  &
                        ,clock_integrate   =CLOCK_ATM_DRV(MY_DOMAIN_ID)    &
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
      INTEGER(kind=KINT) :: RC,RC_NMM_ATM_DRIVER_RUN
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1                                                         !<-- Forward integration so we want horiz diffusion.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- This filter only integrates forward
                          ,atm_grid_comp     =ATM_GRID_COMP             &
                          ,imp_state_atm     =IMP_STATE_ATM             &
                          ,exp_state_atm     =EXP_STATE_ATM             &
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
        CALL ESMF_ClockSet(clock    =CLOCK_ATM_DRV(MY_DOMAIN_ID)        &  !<-- For DFL filter, the starttime of the free forecast
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=0                                                         !<-- Turn off horiz diffusion for backward integration.
        MEAN_ON =1                                                         !<-- Turn off horiz diffusion for backward integration.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Bckward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Bckward'                 &  !<-- The initial backward piece of the filter
                          ,atm_grid_comp     =ATM_GRID_COMP             &
                          ,imp_state_atm     =IMP_STATE_ATM             &
                          ,exp_state_atm     =EXP_STATE_ATM             &
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1                                                         
!
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- The final forward piece of the filter
                          ,atm_grid_comp     =ATM_GRID_COMP             &
                          ,imp_state_atm     =IMP_STATE_ATM             &
                          ,exp_state_atm     =EXP_STATE_ATM             &
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=0                                                         !<-- Turn off horiz diffusion for backward integration.
        MEAN_ON =0                                                         !<-- Turn off horiz diffusion for backward integration.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Bckward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Bckward'                 &  !<-- The initial backward piece of the filter
                          ,atm_grid_comp     =ATM_GRID_COMP             &
                          ,imp_state_atm     =IMP_STATE_ATM             &
                          ,exp_state_atm     =EXP_STATE_ATM             &
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
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        HDIFF_ON=1                                                         !<-- Forward integration so we want horiz diffusion.
        MEAN_ON =1                                                         !<-- Forward integration so we want horiz diffusion.
        NDFISTEP=HALFDFIINTVAL/TIMESTEP(MY_DOMAIN_ID)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='Clock_Direction'                  &
                              ,value='Forward'                          &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                              ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='MEAN_ON'                          &
                              ,value=MEAN_ON                            &
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=IMP_STATE_ATM                      &  !<-- This ATM component's import state for filter
                              ,name ='NDFISTEP'                         &
                              ,value=NDFISTEP                           &
                              ,rc   =RC)
!
        CALL NMM_INTEGRATE(clock_direction   ='Forward'                 &  !<-- The final forward piece of the filter
                          ,atm_grid_comp     =ATM_GRID_COMP             &
                          ,imp_state_atm     =IMP_STATE_ATM             &
                          ,exp_state_atm     =EXP_STATE_ATM             &
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NMM_ATM_DRIVER_RUN)
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
      END SUBROUTINE NMM_ATM_DRIVER_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_ATM_DRIVER_INIT(CLOCK_MAIN,     &
                                     IMP_STATE,      &
                                     EXP_STATE)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CREATES THE CONFIGURE OBJECT AND THE
!***  INDIVIDUAL ATM COMPONENTS FOR ALL GFS MEMBERS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Clock), INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main program's ESMF Clock
      TYPE(ESMF_State), INTENT(INOUT) :: IMP_STATE, EXP_STATE            !<-- The ATM Driver import/export sta
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      INTEGER            :: I,IJ,J,ME,PE_MAX,TASKS
!
      CHARACTER(20)      :: PELAB
      TYPE(ESMF_LOGICAL) :: Cpl_flag
!
      INTEGER            :: RC,RC_INIT_DRV
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_INIT_DRV=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      CALL ESMF_ConfigDestroy(config=CF_DRIVER                          &  !<-- We no longer need the original CF from ATM_DRIVER_INIT
                             ,rc    =RC)
!
      CF_DRIVER=ESMF_ConfigCreate(rc=RC)
!
!-----------------------------------------------------------------------
!***  Load the GFS configure file into the ATM Driver configure object.
!-----------------------------------------------------------------------
!
      CALL ESMF_ConfigLoadFile(config  =CF_DRIVER                       & !<-- The ATM Driver configure object
                              ,filename='atm_namelist.rc'               & !<-- The GFS configure file
                              ,rc      =RC)
!
!-----------------------------------------------------------------------
!***  Extract the total number of GFS ensemble members.
!-----------------------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(config=CF_DRIVER                     &
                                  ,value =TOTAL_MEMBER                  &
                                  ,label ='total_member:'               &
                                  ,rc    =RC)
!
!-----------------------------------------------------------------------
!***  Allocate a standard set of arrays for each GFS member.
!-----------------------------------------------------------------------
!
      ALLOCATE(GC_ATM      (TOTAL_MEMBER))
      ALLOCATE(GRIDCOMPNAME(TOTAL_MEMBER))
      ALLOCATE(PE_MEMBER   (TOTAL_MEMBER))
!
!-----------------------------------------------------------------------
!***  Obtain the total task count and local ID from the VM.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGetGlobal(VM_GLOBAL, rc = RC)
      CALL ESMF_VmGet(VM_GLOBAL,              &
                      pecount  = TASKS,       &
                      localpet = MYPE_GLOBAL, &
                      rc       = RC)
!
!-----------------------------------------------------------------------
!***  For each member create the ATM component and state names 
!***  and fill in the task information.
!-----------------------------------------------------------------------
!
      WRITE(IMPSTATENAME   , '("atm import state")')
      WRITE(EXPSTATENAME   , '("atm export state")')

      PE_MEMBER = 0
      DO I=1,TOTAL_MEMBER
          WRITE(GRIDCOMPNAME(I), '("atm main grid component",I2.2)') I
!
          WRITE(PELAB, '("PE_MEMBER", I2.2, ":")') I
          CALL ESMF_ConfigGetAttribute(CF_DRIVER, PE_MEMBER(I), label = PELAB, rc = RC)
          IF(PE_MEMBER(I) == 0) PE_MEMBER(I) = TASKS / TOTAL_MEMBER
      ENDDO
!
      PE_MAX = 1
      DO I=1,TOTAL_MEMBER
        PE_MAX = MAX(PE_MAX,PE_MEMBER(I))
      ENDDO
!
!-----------------------------------------------------------------------
!***  Set up the PE list.
!-----------------------------------------------------------------------
!
      ALLOCATE(PETLIST(PE_MAX, TOTAL_MEMBER))
!
      IJ = 0
      DO J = 1, TOTAL_MEMBER
        DO I = 1, PE_MEMBER(J)
          PETLIST(I, J) = IJ
          IF(MYPE_GLOBAL == IJ) THEN
            MEMBER_ID = J
          END IF
          IJ = IJ+1
        END DO
      END DO
!
!-----------------------------------------------------------------------
!***  Now we can create the ATM components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create GFS ATM Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1,TOTAL_MEMBER
        GC_ATM(I) = ESMF_GridCompCreate (                               &
                              name         = GRIDCOMPNAME(I)            &
                             ,gridcomptype = ESMF_ATM                   &
                             ,petlist      = PETLIST(1:PE_MEMBER(I),I)  &
                             ,config       = CF_DRIVER                  &
                             ,rc           = RC)
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Init, Run, and Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register ATM Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
        CALL ESMF_GridCompSetServices(GC_ATM(I)                         &  !<-- The GFS ATM gridded components
                                     ,ATM_REGISTER                      &  !<-- User's subroutineName
                                     ,RC)
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the import and export states for the ATM components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_ATM = ESMF_StateCreate(statename = IMPSTATENAME               &
                                ,statetype = ESMF_STATE_IMPORT          &
                                ,rc        = RC)
!
      EXP_ATM = ESMF_StateCreate(statename = EXPSTATENAME               &
                                ,statetype = ESMF_STATE_EXPORT          &
                                ,rc        = RC)
      CALL ESMF_StateAdd(IMP_STATE, IMP_ATM, rc = RC)
      CALL ESMF_StateAdd(EXP_STATE, EXP_ATM, rc = RC)

      CALL ESMF_AttributeGet(IMP_STATE, 'Cpl_flag', Cpl_flag, rc = rc)
      CALL ESMF_AttributeSet(IMP_ATM,   'Cpl_flag', Cpl_flag, rc = rc)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Obtain the forecast times from the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GFS_ATM_DRIVER_INIT: Start Time from Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_MAIN                         &  !<-- The Main ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time (ESMF)
                        ,runDuration=RUNDURATION                        &  !<-- The simulation run duration (ESMF)
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract fundamental timestep information from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Timestep Information from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_DRIVER                     &  !<-- The config object
                                  ,value =TIMESTEP_SEC_WHOLE            &  !<-- The variable filled (integer part of timestep (sec))
                                  ,label ='dt_int:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF_DRIVER                     &  !<-- The config object
                                  ,value =TIMESTEP_SEC_NUMERATOR        &  !<-- The variable filled (numerator of timestep fraction)
                                  ,label ='dt_num:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF_DRIVER                     &  !<-- The config object
                                  ,value =TIMESTEP_SEC_DENOMINATOR      &  !<-- The variable filled (denominator of timestep fraction)
                                  ,label ='dt_den:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Establish the timestep for the Driver Clock.
!-----------------------------------------------------------------------
!
      ALLOCATE(TIMESTEP(1))                                                !    and timestep
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Time Step Interval in Driver Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP(1)                &  !<-- Driver clock's fundamental timestep (sec) (ESMF)
                               ,s           =TIMESTEP_SEC_WHOLE         &  !<-- Whole part of timestep
                               ,sn          =TIMESTEP_SEC_NUMERATOR     &  !<-- Numerator of fractional part
                               ,sd          =TIMESTEP_SEC_DENOMINATOR   &  !<-- Denominator of fractional part
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Driver Clock.
!-----------------------------------------------------------------------
!
!      ALLOCATE(CLOCK_ATM_DRV(1))                                           !<-- GFS ensemble members need only one Clock
!
!      CLOCK_ATM_DRV_NAME='CLOCK_ATM_DRV'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      MESSAGE_CHECK="Create the Clock for the ATM DRIVER Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!      CLOCK_ATM_DRV(1)=ESMF_ClockCreate(name       =CLOCK_ATM_DRV_NAME  &  !<-- The ATM_DRIVER Clock's name
!                                       ,timeStep   =TIMESTEP(1)         &  !<-- The fundamental timestep in this component
!                                       ,startTime  =STARTTIME           &  !<-- Start time of simulation
!                                       ,runDuration=RUNDURATION         &  !<-- Duration of simulation
!                                       ,rc         =RC)
       CALL ESMF_ClockSet(CLOCK_MAIN, timeStep = TIMESTEP(1), rc = rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT_DRV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
!-----------------------------------------------------------------------
!***  Initialize the ATM components of the GFS members.
!-----------------------------------------------------------------------
!
      DO I=1,TOTAL_MEMBER
        IF(MEMBER_ID==I)THEN
          CALL ESMF_GridCompInitialize(gridcomp   =GC_ATM(I)            &  !<-- The ATM gridded component
                                      ,importstate=IMP_ATM              &  !<-- The ATM component's import state
                                      ,exportstate=EXP_ATM              &  !<-- The ATM component's export state
                                      ,clock      =CLOCK_MAIN           &  !<-- The ESMF clock
                                      ,phase      =ESMF_SINGLEPHASE     &
                                      ,rc         =RC)
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_ATM_DRIVER_INIT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_ATM_DRIVER_RUN(CLOCK_MAIN)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE EXECUTES THE INTEGRATION OF THE GFS MEMBERS
!***  THROUGH AN ESMF CALL TO SUBROUTINE ATM_RUN.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_MAIN                         !<-- The main program's ESMF Clock
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: I,RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Call the Run step of the ATM components
!***  for all GFS ensemble members.
!-----------------------------------------------------------------------
!
      DO I = 1, TOTAL_MEMBER
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompRun(gridcomp   =GC_ATM(I)                   &  !<-- The ATM gridded component of this GFS member
                               ,importstate=IMP_ATM                     &  !<-- The ATM component's import state
                               ,exportstate=EXP_ATM                     &  !<-- The ATM component's export state
                               ,clock      =CLOCK_MAIN                  &  !<-- The ESMF clock
                               ,phase      =ESMF_SINGLEPHASE            &
                               ,rc         =RC)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_ATM_DRIVER_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      RECURSIVE SUBROUTINE CALL_ATM_INITIALIZE(ID_ATM,CLOCK_ATM_DRV)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CALLS ATM_INITIALIZE FOR ALL NMM ATM COMPONENTS.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: ID_ATM                                        !<-- ID of the ATM Component (domain) to initialize
!
      TYPE(ESMF_Clock),DIMENSION(1:NUM_DOMAINS),INTENT(INOUT) :: CLOCK_ATM_DRV  !<-- The ATM_DRIVER's ESMF Clock
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
      WRITE(INT_TO_CHAR,FMT)ID_ATM
!
!-----------------------------------------------------------------------
!***  Initialize the ATM component with the ID of ID_ATM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize ATM Component "//INT_TO_CHAR
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =atm_drv_int_state%ATM_GRID_COMP(ID_ATM)  &  !<-- The ATM gridded component
                                  ,importState=atm_drv_int_state%IMP_STATE_ATM(ID_ATM)  &  !<-- The ATM import state
                                  ,exportState=atm_drv_int_state%EXP_STATE_ATM(ID_ATM)  &  !<-- The ATM export state
                                  ,clock      =CLOCK_ATM_DRV(ID_ATM)                    &  !<-- The ATM clock
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
      N_CHILDREN=NUM_CHILDREN(ID_ATM)
!
      IF(N_CHILDREN>0)THEN                                                 !<-- Does the current ATM domain have any children?
        DO N=1,N_CHILDREN                                                  !<-- If so, loop through the children to Initialize them
          ID_CHILD=ID_CHILDREN(N,ID_ATM)
!         write(0,*)' recursive call  n_children=',n_children,' id_atm=',id_atm,' id_child=',id_child
          CALL CALL_ATM_INITIALIZE(ID_CHILD,CLOCK_ATM_DRV)
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CALL_ATM_INITIALIZE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE BOUNDARY_DATA_TO_ATM(EXP_STATE_CPL                     &
                                     ,IMP_STATE_ATM )
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE MOVES NEW BOUNDARY DATA FOR NESTED DOMAINS FROM THE
!***  EXPORT STATE OF THE PARENT-CHILD COUPLER TO THE IMPORT STATE OF
!***  THE NMM NESTS' ATM COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State),INTENT(IN)  :: EXP_STATE_CPL                        !<-- The Parent-Child Coupler's export state
!
      TYPE(ESMF_State),INTENT(OUT) :: IMP_STATE_ATM                        !<-- The nests' ATM component import state
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
!***  to the ATM component's import state.
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
        MESSAGE_CHECK="Insert South H Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert South V Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert North H Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert North V Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert West H Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert West V Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert East H Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
        MESSAGE_CHECK="Insert East V Data into Nest ATM Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_ATM                  &   !<-- Insert data into nest's ATM import state
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
      END SUBROUTINE BOUNDARY_DATA_TO_ATM
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_ATM_DRIVER_COMP
!
!-----------------------------------------------------------------------
