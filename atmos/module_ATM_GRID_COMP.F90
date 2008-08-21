!-----------------------------------------------------------------------
!
      MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE ATM (Atmosphere) GRIDDED COMPONENT MODULE.
!***  IT WILL SET UP GRIDDED AND COUPLER SUBCOMPONENTS
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
!
! USAGE: ATM Gridded component parts CALLed from MAIN_ESMF.F90
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_INCLUDE
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
      
!
      USE MODULE_DYNAMICS_GRID_COMP  
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &
                                   ,MPI_COMM_COMP                       &
                                   ,MPI_COMM_INTER_ARRAY                &
                                   ,MYPE_SHARE

!
      USE MODULE_PHYSICS_GRID_COMP,ONLY: PHY_REGISTER
!
      USE MODULE_DYN_PHY_CPL_COMP,ONLY: DYN_PHY_CPL_REGISTER 
!
      USE MODULE_GET_CONFIG_DYN
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_GET_CONFIG_WRITE
      USE MODULE_CONTROL,ONLY: TIMEF
      USE MODULE_DIAGNOSE,ONLY : FIELD_STATS
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
      
      
      USE gfs_dynamics_grid_comp_mod  ,only: gfs_dyn_setservices

      USE gfs_physics_grid_comp_mod   ,only: gfs_phy_setservices

      USE atmos_dyn_phy_cpl_comp_mod  ,only: atm_cpl_setservices
 
      USE atmos_grid_comp_mod, only: atmos_setservices, atmos_initialize, atmos_run, atmos_finalize        
                                  
                                  

!
!-----------------------------------------------------------------------
!***  LIST OTHER MODULES WITH NON_GENERIC ROUTINES USED BY ATM.
!-----------------------------------------------------------------------
!
      USE MODULE_WRITE_ROUTINES ,ONLY: WRITE_INIT,WRITE_ASYNC             !<-- These are routines used only when asynchronous
      USE MODULE_WRITE_GRID_COMP,ONLY: WRITE_SETUP                      & !    quilting is specified by the user in the
                                      ,WRITE_DESTROY                      !    configure file for history output.
!
      USE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
#include "ESMF_LogMacros.inc"
!-----------------------------------------------------------------------
! 
      PRIVATE
!
      PUBLIC :: ATM_REGISTER

!
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT)       :: DT                                         !<-- The fundamental timestep (sec)
!
      INTEGER(KIND=KINT)    :: MYPE                                       !<-- Each MPI task ID

      INTEGER(KIND=KINT)    :: NPE_PRINT  
!
      CHARACTER(3)          :: CORE                                       !<-- The name of the selected dynamic core
!
      TYPE(ESMF_Clock),SAVE :: CLOCK_ATM                                  !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_VM),SAVE    :: VM                                         !<-- The ESMF virtual machine.
!
      TYPE(ESMF_Alarm),SAVE :: ALARM_CLOCKTIME                            !<-- The ESMF Alarm for clocktime prints
      TYPE(ESMF_Alarm),SAVE :: ALARM_HISTORY                              !<-- The ESMF Alarm for history output
!
!
      TYPE(ESMF_GridComp),save :: gc_gfs_dyn
      TYPE(ESMF_GridComp),save :: gc_gfs_phy
      TYPE(ESMF_CplComp), save :: gc_atm_cpl

!
      TYPE(ESMF_State),save :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),save :: imp_gfs_phy,exp_gfs_phy
      INTEGER :: inpes,jnpes  ! mpi tasks in i and j

!      INTEGER(KIND=KINT),PUBLIC :: IM,JM,LM
!      INTEGER(KIND=KINT)        :: NUM_PES
!
!
!      LOGICAL :: ADVECT_TRACERS                                         & !<-- Flag for advecting tracers
!                ,OLD_PASSIVE                                            & !<-- Flag for old passive advection
!                ,OPERATIONAL_PHYSICS                                      !<-- Flag to designate use of operational physics suite
!
!-----------------------------------------------------------------------
!***  FOR DETERMINING CLOCKTIMES OF VARIOUS PIECES OF THE DYNAMICS.
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT) :: btim,btim0
!
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
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                !<-- ATM gridded component
!
      INTEGER,INTENT(OUT)               :: RC_REG                       !<-- Return code for register
!     
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)                :: CF                            !<-- The config object
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC=ESMF_SUCCESS       ! Error signal variable
      RC_REG=ESMF_SUCCESS   ! Error signal variable

!-----------------------------------------------------------------------
!***  Register the ATM INITIALIZE subroutine.  Since it is just one 
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
!***  Register the ATM RUN subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM Run"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type
                                     ,ATM_RUN                           &  !<-- User's subroutine name
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
      MESSAGE_CHECK="Set Entry Point for ATM Finalize"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_INITIALIZE(ATM_GRID_COMP                           &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_MAIN                              &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE SETS UP FUNDAMENTAL ASPECTS OF THE MODEL RUN.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                !<-- The ATM gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE,EXP_STATE          !<-- The ATM component's import/export states
      TYPE(ESMF_Clock)   ,INTENT(INOUT)    :: CLOCK_MAIN                   !<-- The main program's ESMF Clock!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT                      !<-- Return code for Initialize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                          !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                 !<-- The ATM internal state pointer
!
      TYPE(ESMF_Config)                :: CF                            !<-- The config object
!
      TYPE(ESMF_Grid)                  :: GRID_ATM                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the ATM gridded component.
      TYPE(ESMF_Grid)                  :: GRID_DYN                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the dynamics gridded component.
      TYPE(ESMF_Grid)                  :: GRID_PHY                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the physics gridded component.
!
      INTEGER                          :: TIMESTEP_SEC_WHOLE            !<-- Integer part of timestep in seconds
      INTEGER                          :: TIMESTEP_SEC_NUMERATOR        !<-- Numerator of timestep fraction
      INTEGER                          :: TIMESTEP_SEC_DENOMINATOR      !<-- Denominator of timestep fraction
      TYPE(ESMF_TimeInterval)          :: TIMESTEP                      !<-- The ESMF timestep (s)
      TYPE(ESMF_TimeInterval)          :: RUNDURATION                   !<-- The ESMF simulation length (h)
      TYPE(ESMF_Time)                  :: CURRTIME                      !<-- The ESMF current time.
      TYPE(ESMF_Time)                  :: STARTTIME                     !<-- The ESMF start time.
      TYPE(ESMF_Grid)  :: grid_gfs_dyn     ! the ESMF grid for the integration attached to
                                                   ! the dynamics gridded component.
      TYPE(ESMF_Grid)  :: grid_gfs_phy     ! the ESMF grid for the integration attached to
!
      TYPE(ESMF_Grid)                  :: grid_atmos    ! the esmf grid for the integration attached to

      INTEGER                          :: NHOURS_CLOCKTIME              !<-- Hours between clocktime prints
      INTEGER                          :: NHOURS_HISTORY                !<-- Hours between history output
      TYPE(ESMF_TimeInterval)          :: TIMEINTERVAL_CLOCKTIME        !<-- ESMF time interval between clocktime prints (h)
      TYPE(ESMF_TimeInterval)          :: TIMEINTERVAL_HISTORY          !<-- ESMF time interval between history output (h)
!
      LOGICAL                          :: DIST_MEM,PHYSICS_ON           !<-- Logical flag for distributed
      LOGICAL                          :: QUILTING                      !<-- Logical flag for quilting in
!
      INTEGER                          :: I,IERR,RC,IRTN,J,K,N,NN
      INTEGER                          :: RC_FINAL
      CHARACTER(50)                    :: MODE
!
      type(ESMF_DistGrid)          :: DistGrid_atmos

!
      integer, dimension(2)         :: ncounts     ! parameter array to set up the
                                                   ! size of the 2-d ESMF grid.
      integer, dimension(2)         :: min,max     ! parameter arrays to set up the
                                                   ! start number and the end number of
                                                   ! the ESMF grid in each dimension.
      integer,dimension(ESMF_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      integer                      :: im,jm,lm                  ! full grid dimensions
      integer                      :: num_pes,num_pes_fcst,num_pes_tot
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator!
      logical                      :: global
      RC     =ESMF_success


!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  INITIALIZE TIMING VARIABLES.
!-----------------------------------------------------------------------
!
      btim0=timef()
      total_tim=0.
      totalsum_tim=0.
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from ATM Component"
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
!***  EXTRACT THE DYNAMIC CORE NAME.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Dynamic Core Name"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!***  ALLOCATE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
      ALLOCATE(ATM_INT_STATE,stat=RC)
!
      wrap%ATM_INT_STATE=>ATM_INT_STATE
!
!-----------------------------------------------------------------------
!***  ATTACH THE ATM INTERNAL STATE TO THE ATM GRIDDED COMPONENT.
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
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO OBTAIN THE MPI TASK IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from ATM Gridded Component"
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
      ENDIF
!-----------------------------------------------------------------------
      IF (CORE=='nmm') THEN
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
!-----------------------------------------------------------------------
!***  WILL THE WRITE COMPONENTS WITH ASYNCHRONOUS QUILTING BE USED?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Quilting Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =QUILTING                      &  !<-- The variable filled (flag
                                  ,label ='quilting:'                   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      atm_int_state%QUILTING=QUILTING                                      !<-- Save this for the ATM's Run step
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  EXTRACT FUNDAMENTAL TIMESTEP INFORMATION FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Timestep Information from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_WHOLE            &  !<-- The variable filled (integer part of timestep (sec))
                                  ,label ='dt_int:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_NUMERATOR        &  !<-- The variable filled (numerator of timestep fraction)
                                  ,label ='dt_num:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_DENOMINATOR      &  !<-- The variable filled (denominator of timestep fraction)
                                  ,label ='dt_den:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ESTABLISH THE TIMESTEP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set up Time Step Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                   &  !<-- The model's fundamental timestep (sec) (ESMF)
                               ,s           =TIMESTEP_SEC_WHOLE         &
                               ,sn          =TIMESTEP_SEC_NUMERATOR     &
                               ,sd          =TIMESTEP_SEC_DENOMINATOR   &
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DT=TIMESTEP_SEC_WHOLE+                                            &  !<-- The model's fundamental timestep (sec) (REAL)
         REAL(TIMESTEP_SEC_NUMERATOR)/REAL(TIMESTEP_SEC_DENOMINATOR)
!
!-----------------------------------------------------------------------
!***  OBTAIN THE FORECAST START/END TIMES FROM THE MAIN CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Start Time and Fcst Duration from Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_MAIN                         &  !<-- The Main ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time
                        ,runDuration=RUNDURATION                        &  !<-- The simulation run duration
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CURRTIME=STARTTIME                                                   !<-- Set the current forecast time
!
!-----------------------------------------------------------------------
!***  WITH DATA FROM ABOVE, CREATE THE LOCAL ESMF CLOCK
!***  TO CONTROL THE TIMESTEPPING IN THE ATM COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Local Clock for the ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CLOCK_ATM=ESMF_ClockCreate(name       ='CLOCK_ATM'                &  !<-- The ATM Clock's name
                                ,timestep   =TIMESTEP                   &  !<-- The fundamental timestep in this component
                                ,starttime  =STARTTIME                  &  !<-- Start time of simulation
                                ,runduration=RUNDURATION                &  !<-- Duration of simulation
                                ,rc         =RC)
!
      CALL ESMF_ClockSet(clock   =CLOCK_ATM                             &  !<-- The ATM Component's Clock
                        ,currtime=CURRTIME                              &  !<-- Current time of simulation
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  GET THE HISTORY OUTPUT INTERVAL (HOURS) FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain History Interval from the Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NHOURS_HISTORY                &  !<-- Fill this variable
                                  ,label ='nhours_history:'             &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE HISTORY FILE OUTPUT ALARM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the History Output Alarm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_HISTORY           &  !<-- Time interval between
                               ,h           =NHOURS_HISTORY                 &  !<-- Hours between history
                               ,rc          =RC)
!
      ALARM_HISTORY=ESMF_AlarmCreate(name             ='ALARM_HISTORY'      &
                                    ,clock            =CLOCK_ATM            &  !<-- ATM Clock
                                    ,ringTime         =STARTTIME            &  !<-- Forecast start time (ESMF)
                                    ,ringInterval     =TIMEINTERVAL_HISTORY &  !<-- Time interval between
                                    ,ringTimeStepCount=1                    &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.              &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE ALARM FOR PRINTING CLOCKTIMES USED BY MODEL PIECES.
!***  READ IN FORECAST TIME INTERVAL FOR CLOCKTIME OUTPUT AS WELL AS
!***  THE SELECTED TASK ID THAT WILL PROVIDE THE CLOCKTIMES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
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
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NPE_PRINT                     &  !<-- Fill this variable
                                  ,label ='npe_print:'                  &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create ESMF Clocktime Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_CLOCKTIME         &  !<-- Time interval between
                               ,h           =NHOURS_CLOCKTIME               &  !<-- Hours between clocktime writes (REAL)
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
      ALARM_CLOCKTIME=ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'      &
                                    ,clock              =CLOCK_ATM              &  !<-- ATM Clock
                                    ,ringTime           =STARTTIME              &  !<-- Forecast start time (ESMF)
                                    ,ringInterval       =TIMEINTERVAL_CLOCKTIME &  !<-- Time interval between clocktime prints (ESMF)
                                    ,ringTimeStepCount  =1                      &  !<-- The Alarm rings for this many timesteps
                                    ,sticky             =.false.                &  !<-- Alarm does not ring until turned off
                                    ,rc                 =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF

      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  MODEL-SPECIFIC ROUTINES MUST BE INVOKED IN ORDER TO ESTABLISH
!***  THE ESMF Grid.  THE DIFFERENT INTEGRATION GRIDS NECESSITATE
!***  DIFFERENT WAYS OF SETTING UP BOTH THE PARALLELISM FOR
!***  DISTRIBUTED MEMORY RUNS AND THE ESMF Grid ITSELF.
!***  WHEN THE PARALLELISM IS CONSTRUCTED, THE LOCAL DOMAIN LIMITS
!***  NEED TO BE INSERTED INTO THE ATM COMPONENT'S INTERNAL STATE
!***  IF QUILTING IS TO BE USED.  SEE 'IF(QUILTING)THEN' BELOW.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(CORE=='nmm')THEN
        CALL NMM_SETUP(ATM_GRID_COMP,ATM_INT_STATE,GRID_ATM)
!-----------------------------------------------------------------------
!***  ATTACH THE CORE-SPECIFIC ESMF GRID TO THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the ESMF Grid to the ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=ATM_GRID_COMP                      & !<-- The ATM gridded component
                           ,grid    =GRID_ATM                           & !<-- Attach the ESMF grid to the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!-----------------------------------------------------------------------

      ELSE IF (CORE=='gfs') THEN
        CALL GFS_SETUP(ATM_GRID_COMP,grid_atmos)
      ENDIF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE DYNAMICS GRIDDED SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
!***  SINCE THERE IS ONLY A SINGLE INTEGRATION GRID, GIVE THE
!***  DYNAMICS THE ATM COMPONENT'S GRID.
!***  NOTE THAT THIS SUBCOMPONENT IS PART OF THE ATM COMPONENT'S
!***  INTERNAL STATE.  THIS WILL BE CONVENIENT IF WE NEED TO REACH
!***  THE Dynamics COMPONENT VIA THE ATM COMPONENT SUCH AS HAPPENS
!***  WHEN WRITE COMPONENTS ARE ESTABLISHED.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Dynamics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%DYN_GRID_COMP=ESMF_GridCompCreate(                  &
                                  name      ="Dynamics component"       &  !<-- Name of the new Dynamics gridded component
                                 ,configFile='configure_file'           &  !<-- Attach this configure file to the component
                                 ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                                 ,rc        =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_gfs_dyn=esmf_gridcompcreate( name      ="dynamics component"   &
                                     ,configFile='dyn_namelist.rc'      &
                                     ,rc        =RC)
      ENDIF

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
      MESSAGE_CHECK="Register Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompSetServices(atm_int_state%DYN_GRID_COMP         &  ! <-- The Dynamics gridded component
                                   ,DYN_REGISTER                        &  ! <-- The user's subroutineName for Register
                                   ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompsetservices(gc_gfs_dyn                   &
                                   ,gfs_dyn_setservices             &
                                   ,RC)
      ENDIF

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
      MESSAGE_CHECK="Attach ESMF Grid to the Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      
      IF (CORE=='nmm') THEN
      GRID_DYN=GRID_ATM
      CALL ESMF_GridCompSet(gridcomp=atm_int_state%DYN_GRID_COMP        &  !<-- The Dynamics component
                           ,grid    =GRID_DYN                           &  !<-- The Dynamics ESMF grid
                           ,rc      =RC)
      ENDIF
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
      IF (CORE=='nmm') THEN
      atm_int_state%IMP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Import" &  !<-- The Dynamics import state name
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
      atm_int_state%EXP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Export" &  !<-- The Dynamics export state name
                                                  ,statetype=ESMF_STATE_EXPORT &
                                                  ,rc       =RC)
      ELSE IF (CORE=='gfs') THEN
      imp_gfs_dyn=esmf_statecreate(statename="dynamics import"        &
                                    ,statetype=esmf_state_import      &
                                    ,rc       =RC)
!
      exp_gfs_dyn=esmf_statecreate(statename="dynamics export"        &
                                    ,statetype=esmf_state_export      &
                                    ,rc       =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------------------
!***  CHECK CONDITION TO RUN ADIABATIC CASE (NO PHYSICS) 
!---------------------------------------------------------
!
      call ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' initialize without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' initialize with physics coupling '
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' initialize without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' initialize with physics coupling '
      endif
      ENDIF
!
      IF (PHYSICS_ON) THEN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE PHYSICS GRIDDED SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
!***  SINCE THERE IS ONLY A SINGLE INTEGRATION GRID, GIVE THE
!***  PHYSICS THE ATM COMPONENT'S GRID.
!***  NOTE THAT THIS SUBCOMPONENT IS PART OF THE ATM COMPONENT'S
!***  INTERNAL STATE.  THIS WILL BE CONVENIENT IF WE NEED TO REACH
!***  THE Physics COMPONENT VIA THE ATM COMPONENT SUCH AS HAPPENS
!***  WHEN WRITE COMPONENTS ARE ESTABLISHED.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Physics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%PHY_GRID_COMP=ESMF_GridCompCreate(                  &
                                  name      ="Physics component"        &  !<-- Name of the new Physics gridded component
                                 ,configFile='configure_file'           &  !<-- Attach this configure file to the component
                                 ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                                 ,rc        =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_gfs_phy=esmf_gridcompcreate( name      ="physics component"    &
                                     ,configfile='phy_namelist.rc'      &
                                     ,rc        =rc)
      ENDIF
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
      MESSAGE_CHECK="Register Physics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompSetServices(atm_int_state%PHY_GRID_COMP         &  ! <-- The Physics gridded component
                                   ,PHY_REGISTER                        &  ! <-- The user's subroutineName
                                   ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompsetservices(gc_gfs_phy                   &
                                   ,gfs_phy_setservices             &
                                   ,rc)
      ENDIF
!
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      GRID_PHY=GRID_ATM
      CALL ESMF_GridCompSet(gridcomp=atm_int_state%PHY_GRID_COMP        &  !<-- The Physics component
                           ,grid    =GRID_PHY                           &  !<-- The Physics ESMF grid
                           ,rc      =RC)
!      ELSE IF (CORE=='gfs') THEN
!      grid_gfs_phy=grid_atmos
!      call esmf_gridcompset(gridcomp= gc_gfs_phy                         &
!                            ,grid    =grid_gfs_phy                       &
!                            ,rc	     =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF !PHYSICS_ON
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
      IF (CORE=='nmm') THEN
      atm_int_state%IMP_STATE_PHY=ESMF_StateCreate(stateName="Physics Import"  &  !<-- The Physics import
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
 
      atm_int_state%EXP_STATE_PHY=ESMF_StateCreate(stateName="Physics Export"  &  !<-- The Physics export
                                                  ,statetype=ESMF_STATE_EXPORT &
                                                  ,rc       =RC)
      
      ELSE IF (CORE=='gfs') THEN
      imp_gfs_phy=esmf_statecreate(statename="physics import"         &
                                    ,statetype=esmf_state_import      &
                                    ,rc       =rc)
!
      exp_gfs_phy=esmf_statecreate(statename="physics export"         &
                                    ,statetype=esmf_state_export      &
                                    ,rc       =rc)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE DYNAMICS-PHYSICS COUPLER SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
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
      IF (CORE=='nmm') THEN
      atm_int_state%COUPLER_DYN_PHY_COMP=ESMF_CplCompCreate                   &
                                         (name   ="Dyn-Phy coupler component" &
                                         ,petList=atm_int_state%PETLIST_FCST  &  !<-- The forecast task IDs
                                         ,rc     =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_atm_cpl=esmf_cplcompcreate(name="coupler component"            &
                                     ,rc  =RC)
      ENDIF

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
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompSetServices(atm_int_state%COUPLER_DYN_PHY_COMP   &  ! <-- The Dyn-Phy coupler component
                                  ,DYN_PHY_CPL_REGISTER                 &  ! <-- The user's subroutineName
                                  ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_cplcompsetservices(gc_atm_cpl                &
                                  ,atm_cpl_setservices              &
                                  ,RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  SETUP THE WRITE COMPONENT(S).
!***  THIS MUST BE DONE PRIOR TO EXECUTING THE INITIALIZE STEPS
!***  OF THE DYNAMICS AND PHYSICS COMPONENTS BECAUSE THE WRITE
!***  COMPONENTS' IMPORT STATES ARE REQUIRED BY THOSE STEPS
!***  WHEN 'quilting' IS ACTIVE AND WRITE_SETUP IS THE ROUTINE
!***  IN WHICH THE WRITE COMPONENTS' IMPORT STATES ARE INSERTED
!***  INTO THE DYNAMICS AND PHYSICS EXPORT STATES WHICH ARE
!***  IN TURN PART OF THE ATM INTERNAL STATE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF (CORE=='nmm') THEN
      IF(QUILTING)THEN
        CALL WRITE_SETUP(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
      ENDIF
      ENDIF
      
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEPS FOR THE GRIDDED SUBCOMPONENTS.
!***  THESE ARE THE INITIALIZE SUBROUTINES SPECIFIED IN THE
!***  REGISTER ROUTINES CALLED IN ESMF_GridCompSetServices ABOVE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%DYN_GRID_COMP  &  !<-- The dynamics gridded component
                                  ,importState=atm_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                                  ,exportState=atm_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                                  ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                  ,phase      =ESMF_SINGLEPHASE             &
                                  ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call ESMF_GridCompInitialize(gridcomp   =gc_gfs_dyn            &
                                  ,importstate=imp_gfs_dyn           &
                                  ,exportstate=exp_gfs_dyn           &
                                  ,clock      =CLOCK_MAIN            &
                                  ,phase      =esmf_singlephase      &
                                  ,rc         =RC)
      grid_gfs_dyn=grid_atmos
      call esmf_gridcompset(         gc_gfs_dyn      &
                           ,grid    =grid_gfs_dyn    &
                           ,rc      =rc)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON) THEN
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%PHY_GRID_COMP  &  !<-- The physics gridded component
                                  ,importState=atm_int_state%IMP_STATE_PHY  &  !<-- The physics import state
                                  ,exportState=atm_int_state%EXP_STATE_PHY  &  !<-- The physics export state
                                  ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                  ,phase      =ESMF_SINGLEPHASE             &
                                  ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompinitialize(gridcomp   =gc_gfs_phy            &
                                  ,importstate=imp_gfs_phy           &
                                  ,exportstate=exp_gfs_phy           &
                                  ,clock      =CLOCK_MAIN            &
                                  ,phase      =esmf_singlephase      &
                                  ,rc         =RC)
      grid_gfs_phy=grid_atmos
      call esmf_gridcompset(         gc_gfs_phy      &
                           ,grid    =grid_gfs_phy    &
                           ,rc      =RC)

      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF !PHYSICS_ON
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  INITIALIZE THE DYN-PHY COUPLER SUBCOMPONENT.
!***  THE CHOICE OF IMPORT AND EXPORT STATE DOES NOT MATTER
!***  FOR THE INITIALIZE STEP.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dyn-Phy Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompInitialize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dyn_phy coupler component
                                 ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The dyn-phy coupler import state
                                 ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The dyn-phy coupler export state
                                 ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                                 ,rc         =RC)
      ELSE If (CORE=='gfs') THEN
      call esmf_cplcompinitialize(cplcomp    =gc_atm_cpl              &
                                 ,importstate=exp_gfs_dyn             &
                                 ,exportstate=imp_gfs_phy             &
                                 ,clock      =CLOCK_MAIN              &
                                 ,rc         =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENT(S).
!***  THIS MUST BE DONE AFTER THE INITIALIZATION OF THE
!***  DYNAMICS AND PHYSICS COMPONENTS BECAUSE THOSE COMPONENTS'
!***  INTERNAL STATES CONTAIN THE HISTORY OUTPUT VARIABLES.
!***  POINTERS TO THOSE VARIABLES ARE SET DURING THEIR INITIALIZE
!***  STEPS AND ARE THEN LOADED INTO THE WRITE COMPONENTS'
!***  IMPORT STATES WHICH THEMSELVES RESIDE IN THE DYNAMICS/PHYSICS
!***  EXPORT STATES.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
      IF(QUILTING)THEN
        CALL WRITE_INIT(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
      ENDIF
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
      total_tim=timef()-btim0
!
      if(mype==0)write(0,*)' atm_init_tim=',total_tim*1.e-3
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_INITIALIZE
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_RUN(ATM_GRID_COMP                                  &
                        ,IMP_STATE,EXP_STATE                            &
                        ,CLOCK_MAIN                                     &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  RUN THE ATM (Atmosphere) GRIDDED COMPONENT.
!-----------------------------------------------------------------------
      USE MODULE_DIGITAL_FILTER_NMM
      USE MODULE_DIGITAL_FILTER_GFS 
!      USE MODULE_ALARMS

!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(IN)    :: IMP_STATE                       !<-- The ATM Run step's import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM Run step's export
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_RUN                          !<-- Return code for the Run step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)         :: CF  
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                             !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
      TYPE(ESMF_Time)                  :: CURRTIME                         !<-- The current forecast time
      TYPE(ESMF_Time)                  :: STARTTIME                        !<-- The current forecast time
      TYPE(ESMF_Time)                  :: SDFITIME
      TYPE(ESMF_Time)                  :: DFITIME                          !<-- Digital filter time interval
      TYPE(ESMF_Time)                  :: HALFDFITIME                      !<-- Digital filter time interval
      TYPE(ESMF_TimeInterval)          :: HALFDFIINTVAL                    !<-- Digital filter time interval
      TYPE(ESMF_TimeInterval)          :: TIMESTEP                         !<-- Digital filter time interval
      TYPE(ESMF_Alarm), save           :: alarm(2)
!
      INTEGER(KIND=KINT)         :: RC,NUM_PES_FCST                        !<-- Error signal variable.
      INTEGER(KIND=KINT)         :: I,IER,J
      INTEGER(KIND=KINT)         :: YY,MM,DD,H,M,S,NDFISTEP
      INTEGER(KIND=KINT),SAVE    :: NTIMESTEP                              !<-- The current forecast timestep (INT)
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF_INT)
      INTEGER(KIND=KINT),SAVE    :: DFIHR                                  !<-- Digital filter time interval
      CHARACTER(50)              :: MODE
      LOGICAL                    :: PHYSICS_ON
      
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
      DFIHR=0
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component                     
                            ,config  =CF                                 &  !<-- The config object (~namelist)
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve ATM Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- The F90 wrap of the ATM internal state
                                        ,RC)
!
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ATM_INT_STATE=>wrap%ATM_INT_STATE
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,rc          =RC)
!
      NTIMESTEP=NTIMESTEP_ESMF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF 
      
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- the ATM gridded component
                           ,vm      =VM                                 &
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF (CORE=='gfs') THEN
      MESSAGE_CHECK="Retrieve VM from ATM Component"
      call esmf_gridcompget(gridcomp=ATM_GRID_COMP                      &
                           ,vm      =VM                                 &
                           ,config  =CF                                 &
                           ,rc      =RC)
!
      call esmf_vmget(     vm       =VM                                 &
                           ,localpet=MYPE                               &
                     	   ,rc      =RC)
      ENDIF

!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' run without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' run with physics coupling.'
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' run without physics coupling.'
      else
        PHYSICS_ON=.true.
        print *,' run with physics coupling.'
      endif
      ENDIF

!
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =DFIHR                         &
                                  ,label ='nhours_dfini:'               &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ------------------------------------------------------------------
!
           IF (DFIHR .GT. 0) THEN
           

! -------------------- initial stage -------------------------------
!
           CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,h=DFIHR,rc=RC)
           IF (CORE=='nmm') THEN
           CALL ESMF_ClockGet(clock =CLOCK_ATM                          &
                          ,starttime =STARTTIME                         &
                          ,timestep =TIMESTEP                           &
                          ,rc =RC)
           ELSE IF (CORE=='gfs') THEN
           CALL ESMF_ClockGet(clock =CLOCK_MAIN                          &
                          ,starttime =STARTTIME                         &
                          ,timestep =TIMESTEP                           &
                          ,rc =RC)
           ENDIF
           NDFISTEP = HALFDFIINTVAL / TIMESTEP
           print *,'NDFISTEP=',NDFISTEP
           HALFDFITIME = STARTTIME + HALFDFIINTVAL
           SDFITIME = STARTTIME + HALFDFIINTVAL/NDFISTEP
           DFITIME = HALFDFITIME + HALFDFIINTVAL
           IF (CORE=='nmm') THEN
           alarm(1) = ESMF_AlarmCreate("Half of digital filter", CLOCK_ATM, &
                                  ringTime=HALFDFITIME,                 &
                                  ringTimeStepCount=1, sticky=.false.,  rc=RC)
           alarm(2) = ESMF_AlarmCreate("Fulldigital filter", CLOCK_ATM,     &
                                  ringTime=DFITIME,                     &
                                  ringTimeStepCount=1, sticky=.false.,  rc=RC)
           ELSE IF (CORE=='gfs') THEN
           alarm(1) = ESMF_AlarmCreate("Half of digital filter", CLOCK_MAIN, &
                                  ringTime=HALFDFITIME,                 &
                                  ringTimeStepCount=1, sticky=.false.,  rc=RC)
           alarm(2) = ESMF_AlarmCreate("Fulldigital filter", CLOCK_MAIN,     &
                                  ringTime=DFITIME,                     &
                                  ringTimeStepCount=1, sticky=.false.,  rc=RC)
           ENDIF
           ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE INTEGRATION TIME LOOP OF THE ATMOSPHERE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      total_tim=total_tim+timef()-btim0
      IF (CORE=='nmm') THEN
!
      timeloop: DO WHILE (.NOT.ESMF_ClockIsStopTime(CLOCK_ATM,RC))
!
!-----------------------------------------------------------------------
!
        btim0=timef()
!-----------------------------------------------------------------------
!***  WRITE A HISTORY FILE AT THE START OF THE FIRST TIMESTEP
!***  OTHERWISE WRITE IT AT THE END OF THE APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!  
       output_0: IF(NTIMESTEP==0)THEN
!
          IF(atm_int_state%QUILTING)THEN
            CALL WRITE_ASYNC(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM,MYPE)
          ENDIF
!
        ENDIF output_0
!
!-----------------------------------------------------------------------
!***  THE FORECAST TASKS EXECUTE THE RUN STEP OF THE DYNAMICS.
!***  THIS IS THE RUN SUBROUTINE SPECIFIED IN
!***  THE DYNAMICS REGISTER ROUTINE CALLED IN
!***  ESMF_GridCompSetServices ABOVE.
!-----------------------------------------------------------------------
!
        fcst_pes: IF(MYPE<atm_int_state%NUM_PES_FCST)THEN                  !<-- Only the forecast tasks integrate
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run Step for Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompRun(gridcomp   =atm_int_state%DYN_GRID_COMP &  !<-- The dynamics component
                               ,importState=atm_int_state%IMP_STATE_DYN &  !<-- The dynamics import state
                               ,exportState=atm_int_state%EXP_STATE_DYN &  !<-- The dynamics export state
                               ,clock      =CLOCK_ATM                   &  !<-- The ATM clock
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  BRING EXPORT DATA FROM THE DYNAMICS INTO THE COUPLER
!***  AND EXPORT IT TO THE PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Coupler moves Information from Dynamics to Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The coupler import state = dynamics export state
                              ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The coupler export state = physics import state
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
!***  EXECUTE THE RUN STEP OF THE PHYSICS.
!-----------------------------------------------------------------------
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run Step for Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
         IF (PHYSICS_ON) THEN
!
          CALL ESMF_GridCompRun(gridcomp   =atm_int_state%PHY_GRID_COMP &  !<-- The physics component
                               ,importState=atm_int_state%IMP_STATE_PHY &  !<-- The physics import state
                               ,exportState=atm_int_state%EXP_STATE_PHY &  !<-- The physics export state
                               ,clock      =CLOCK_ATM                   &  !<-- The ATM Clock
                               ,rc         =RC)

         ELSE


          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%IMP_STATE_PHY        &  !<-- The coupler import state = physics export state
                              ,exportState=atm_int_state%EXP_STATE_PHY        &  !<-- The coupler export state = dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                              ,rc         =RC)

          ENDIF

          
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  BRING EXPORT DATA FROM THE PHYSICS INTO THE COUPLER
!***  AND EXPORT IT TO THE DYNAMICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Coupler moves Information from Physics to Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%EXP_STATE_PHY        &  !<-- The coupler import state = physics export state
                              ,exportState=atm_int_state%IMP_STATE_DYN        &  !<-- The coupler export state = dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                              ,rc         =RC)
!
!-----------------------------------------------------------------------
!
        ENDIF fcst_pes

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! - FIX later
 RC=ESMF_SUCCESS
 RC_RUN=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_ATM                          &
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK AND PRINT FORECAST TIME.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        IF(MYPE==0)THEN
          WRITE(0,25)NTIMESTEP-1,NTIMESTEP*DT/3600.
   25     FORMAT(' Finished Timestep ',i5,' ending at ',f7.3,' hours')
        ENDIF
        call ESMF_ClockGet( clock       =CLOCK_ATM                              &
                            ,currtime=CURRTIME                  &
                            ,rc      =RC)

!
        IF (DFIHR .GT. 0) THEN
 
        IF(MYPE<atm_int_state%NUM_PES_FCST)THEN
        print *,'NTIMESTEP=',NTIMESTEP

         IF ( CURRTIME .eq. SDFITIME ) THEN
           CALL DIGITAL_FILTER_DYN_INIT_NMM(atm_int_state%IMP_STATE_DYN,NDFISTEP)

! -------------------- inital summation  ------------------------------
!
                IF ( PHYSICS_ON ) THEN
                        CALL DIGITAL_FILTER_PHY_INIT_NMM(atm_int_state%IMP_STATE_PHY)
                ENDIF
          ENDIF


!
! -------------------- summation stage ---------------------------------
!
           
      	      CALL DIGITAL_FILTER_DYN_SUM_NMM(atm_int_state%IMP_STATE_DYN)

!
! ----------------------------------------------------------------------
!
              IF( PHYSICS_ON ) then
              IF( CURRTIME .eq. HALFDFITIME ) THEN 
              
               IF (ESMF_AlarmIsRinging(alarm(1), RC)) then
                CALL DIGITAL_FILTER_PHY_SAVE_NMM(atm_int_state%IMP_STATE_PHY)
               ENDIF
              ENDIF
              ENDIF
!
! ----------------------------------------------------------------------
!
! --------------------- final stage ------------------------------------
!
           IF ( CURRTIME .EQ. DFITIME ) THEN
             IF (ESMF_AlarmIsRinging(alarm(2), RC)) then
              print *,' dfi at finaldfitime '
                CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(atm_int_state%IMP_STATE_DYN)
              IF( PHYSICS_ON ) then
                CALL DIGITAL_FILTER_PHY_RESTORE_NMM(atm_int_state%IMP_STATE_PHY)
              ENDIF
!
! ----------------------------------------------------------------------
!
!              CALL ESMF_ClockSet(clock   =CLOCK_ATM                     &
!                                ,currtime=HALFDFITIME                   &
!                                ,rc      =RC)
              DFITIME = STARTTIME
              DFIHR = 0
              CALL ESMF_ClockPrint(clock=CLOCK_ATM                      &
                              ,options="currtime string"                &
                              ,rc=RC)
             ENDIF
            ENDIF
!
! -----------------------------------------------------------------------
!
        ENDIF
        ENDIF !< -- DFIHR



!-----------------------------------------------------------------------
!***  WRITE A HISTORY FILE IF THE ALARM INDICATES IT IS TIME TO DO SO
!***  AFTER TIMESTEP 0.
!-----------------------------------------------------------------------
!
        output: IF(ESMF_AlarmIsRinging(alarm=ALARM_HISTORY              &  !<-- The history output alarm
                                      ,rc   =RC))THEN
!
          IF(atm_int_state%QUILTING)THEN
            CALL WRITE_ASYNC(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM,MYPE)
          ENDIF
!
        ENDIF output
!
!
!-----------------------------------------------------------------------
!
        total_tim=total_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!***  PRINT CLOCKTIMES OF INTEGRATION SECTIONS
!***  ON MPI TASK OF CHOICE.
!-----------------------------------------------------------------------
!
        clocktimes: IF(ESMF_AlarmIsRinging(alarm=ALARM_CLOCKTIME        &  !<-- The alarm to print clocktimes used by model parts
                                          ,rc   =RC))THEN

             CALL PRINT_CLOCKTIMES(NTIMESTEP,MYPE,NPE_PRINT)

        ENDIF clocktimes
!-----------------------------------------------------------------------
!
      ENDDO timeloop

      ELSE IF (CORE=='gfs') THEN

          do while (.not.esmf_clockisstoptime(CLOCK_MAIN,rc))
          call esmf_logwrite("execute dynamics",esmf_log_info,rc=rc)
          call esmf_gridcomprun(gridcomp   =gc_gfs_dyn          &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_MAIN           &
                               ,rc         =RC)
!
          call err_msg(RC,'execute dynamics',RC_RUN)
          call esmf_logwrite("couple dyn_exp-to-phy_imp",      &
                             esmf_log_info,rc=rc)
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl          &
                              ,importstate=exp_gfs_dyn         &
                              ,exportstate=imp_gfs_phy         &
                              ,clock      =CLOCK_MAIN           &
                              ,rc         =RC)
!
          call err_msg(RC,'couple dyn-to-phy',RC_RUN)
!
          IF (PHYSICS_ON) THEN
            call esmf_logwrite("execute physics",esmf_log_info,rc=rc)
            call esmf_gridcomprun(gridcomp   =gc_gfs_phy            &
                                 ,importstate=imp_gfs_phy           &
                                 ,exportstate=exp_gfs_phy           &
                                 ,clock      =CLOCK_MAIN             &
                                 ,rc         =RC)
           call err_msg(RC,'execute physics',RC_RUN)
          ELSE
           call esmf_logwrite("pass phy_imp to phy_exp ",       &
                               esmf_log_info,rc=rc)
!
            call esmf_cplcomprun(            gc_atm_cpl          &
                                ,importstate=imp_gfs_phy         &
                                ,exportstate=exp_gfs_phy         &
                                ,clock      =CLOCK_MAIN          &
                                ,rc         =RC)
!
            call err_msg(RC,'pass phy_imp-to-phy_exp',RC_RUN)
          ENDIF

 !
!
!----------------
!-----------------------------------------------------------------------
!***  bring export data from the physics into the coupler
!***  and export it to the dynamics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple phy_exp-to-dyn_imp",         &
                             esmf_log_info,rc=RC)
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl             &
                              ,importstate=exp_gfs_phy            &
                              ,exportstate=imp_gfs_dyn            &
                              ,clock      =CLOCK_MAIN             &
                              ,rc         =RC)
!
          call err_msg(RC,'couple phy_exp-to-dyn_imp',RC_RUN)
 
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_MAIN                          &
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK AND PRINT FORECAST TIME.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK_MAIN                       &
                          ,currtime=CURRTIME                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
!
        IF (DFIHR .GT. 0) THEN
! --------------------- first stage -----------------------------------
           IF ( CURRTIME .eq. SDFITIME ) THEN 

           CALL DIGITAL_FILTER_DYN_INIT_GFS(imp_gfs_dyn,NDFISTEP)

! -------------------- inital summation  ------------------------------
!
                IF ( PHYSICS_ON ) THEN
                        CALL DIGITAL_FILTER_PHY_INIT_GFS(imp_gfs_phy)
                ENDIF
             ENDIF

!
! -------------------- summation stage ---------------------------------
!
              CALL DIGITAL_FILTER_DYN_SUM_GFS(imp_gfs_dyn)
!
! ----------------------------------------------------------------------
!
              IF( PHYSICS_ON ) then
              IF( CURRTIME .eq. HALFDFITIME ) THEN
               IF (ESMF_AlarmIsRinging(alarm(1), RC)) then
                CALL DIGITAL_FILTER_PHY_SAVE_GFS(imp_gfs_phy)
               ENDIF
              ENDIF
              ENDIF
!
! ----------------------------------------------------------------------
!
! --------------------- final stage ------------------------------------
!
             IF (ESMF_AlarmIsRinging(alarm(2), RC)) then
              print *,' dfi at finaldfitime '
                CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(imp_gfs_dyn)
              IF( PHYSICS_ON ) then
                CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
              ENDIF
!
! ----------------------------------------------------------------------
!
              CALL ESMF_ClockSet(clock   =CLOCK_MAIN                     &
                                ,currtime=HALFDFITIME                   &
                                ,rc      =RC)
              DFITIME = STARTTIME
              DFIHR = 0
              CALL ESMF_ClockPrint(clock=CLOCK_MAIN                      &
                              ,options="currtime string"                &
                              ,rc=RC)
             ENDIF
        ENDIF !< -- DFIHR

        total_tim=total_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!
       ENDDO
        MESSAGE_CHECK="last step dynamics"
                call esmf_gridcomprun(gridcomp=gc_gfs_dyn       &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_MAIN          &
                               ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN) 
      ENDIF
      CALL ESMF_ClockPrint(   clock=CLOCK_MAIN                      &
                              ,options="stoptime string"            &
                              ,rc=rc)


!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_FINALIZE(ATM_GRID_COMP                             &
                             ,IMP_STATE,EXP_STATE                       &
                             ,CLOCK_MAIN                                &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main ESMF Clock
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
      CHARACTER(50):: MODE
      LOGICAL    :: PHYSICS_ON     

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINAL   =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
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
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' finalize without physics coupling. '
      else
        PHYSICS_ON=.true.
        print *,' finalize with physics coupling. '
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' finalize without physics coupling. '
      else
        PHYSICS_ON=.true.
        print *,' finalize with physics coupling. '
      endif
      ENDIF


!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM GRIDDED COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- The F90 wrap of the ATM internal state
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
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%DYN_GRID_COMP &
                                ,importState=atm_int_state%IMP_STATE_DYN &
                                ,exportState=atm_int_state%EXP_STATE_DYN &
                                ,clock      =CLOCK_ATM                   &
                                ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompfinalize(gridcomp   =gc_gfs_dyn                  &
                               ,importstate=imp_gfs_dyn                 &
                               ,exportstate=exp_gfs_dyn                 &
                               ,clock      =CLOCK_MAIN                  &
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
      IF (PHYSICS_ON) THEN
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%PHY_GRID_COMP &
                                ,importState=atm_int_state%IMP_STATE_PHY &
                                ,exportState=atm_int_state%EXP_STATE_PHY &
                                ,clock      =CLOCK_ATM                   &
                                ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompfinalize(gridcomp   =gc_gfs_phy                 &
                                ,importstate=imp_gfs_phy                &
                                ,exportstate=exp_gfs_phy                &
                                ,clock      =CLOCK_MAIN                 &
                                ,rc         =RC)
      ENDIF
      ENDIF !PHYSICS_ON
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
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompFinalize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &
                               ,importState=atm_int_state%EXP_STATE_DYN        &
                               ,exportState=atm_int_state%IMP_STATE_PHY        &
                               ,clock      =CLOCK_ATM                          &
                               ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_cplcompfinalize(cplcomp    = gc_atm_cpl                  &
                               ,importstate=exp_gfs_dyn                 &
                               ,exportstate=imp_gfs_phy                 &
                               ,clock      =CLOCK_MAIN                  &
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
      IF (CORE=='nmm') THEN
      CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_DYN          &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_DYN          &
                            ,rc   =RC)
!
      IF (PHYSICS_ON) THEN
      CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_PHY          &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_PHY          &
                            ,rc   =RC)
      ENDIF
      ELSE IF (CORE=='gfs') THEN
      call esmf_statedestroy(state=imp_gfs_dyn, rc=RC)
      call esmf_statedestroy(state=exp_gfs_dyn, rc=RC)
      IF (PHYSICS_ON) THEN
      call esmf_statedestroy(state=imp_gfs_phy, rc=RC)
      call esmf_statedestroy(state=exp_gfs_phy, rc=RC)
      ENDIF
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
      IF (CORE=='nmm') THEN
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_ATM                            &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       ENDIF
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
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%DYN_GRID_COMP      &
                               ,rc      =RC)
      ELSE 
      call esmf_gridcompdestroy(gridcomp=gc_gfs_dyn, rc=RC)
      ENDIF
    
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON) THEN
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
      CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%PHY_GRID_COMP      &
                               ,rc      =RC)
      ELSE
      call esmf_gridcompdestroy(gridcomp=gc_gfs_phy, rc=RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!------------------------------
!***  DYNAMICS-PHYSICS COUPLER
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompDestroy(cplcomp=atm_int_state%COUPLER_DYN_PHY_COMP &
                              ,rc     =RC)
      ELSE
      call esmf_cplcompdestroy(cplcomp=gc_atm_cpl, rc=RC)
      ENDIF
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
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_SETUP(ATM_GRID_COMP,ATM_INT_STATE,GRID_ATM)
! 
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
      USE MODULE_ATM_INTERNAL_STATE
!
      USE MODULE_DM_PARALLEL,ONLY : DECOMP                              &
                                   ,LOCAL_ISTART,LOCAL_IEND             &
                                   ,LOCAL_JSTART,LOCAL_JEND             &
                                   ,SETUP_SERVERS
!
      USE MODULE_INCLUDE
!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE 
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Grid)         ,INTENT(OUT)   :: GRID_ATM                  !<-- The ESMF GRID for the NMM integration grid
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)   :: CF                                           !<-- The config object
      TYPE(ESMF_VM)       :: VM                                           !<-- The ESMF virtual machine.
      TYPE(ESMF_DELayout) :: MY_DE_LAYOUT                                 !<-- The ESMF layout type array (for tasks).
!
      INTEGER             :: IM,JM                                        !<-- Horizontal dimensions of the full integration grid
      INTEGER(KIND=KINT)  :: INPES,JNPES                                  !<-- MPI tasks in I and J directions
      INTEGER             :: LM                                           !<-- Number of atmospheric model layers
      INTEGER             :: MPI_INTRA,MPI_INTRA_B                        !<-- The MPI intra-communicator
      INTEGER(KIND=KINT)  :: MYPE                                         !<-- My MPI task ID
      INTEGER(KIND=KINT)  :: NUM_PES_FCST                                 !<-- Number of MPI tasks applied to the forecast
      INTEGER             :: NUM_PES_TOT                                  !<-- Total # of MPI tasks in the job
      INTEGER             :: WRITE_GROUPS                               & !<-- Number of groups of write tasks
                            ,WRITE_TASKS_PER_GROUP                        !<-- #of tasks in each write group
!
      INTEGER,DIMENSION(2)             :: NCOUNTS                         !<-- Parameter array to set up the
                                                                          !    size of the 2-D ESMF grid.
      INTEGER,DIMENSION(2)             :: I1                              !<-- # of I and J points in each fcst task's subdomain
!
      INTEGER,DIMENSION(2)  :: MIN,MAX                                    !<-- Parameter arrays to set up the
                                                                          !    start number and the end number of
                                                                          !    the ESMF grid in each dimension.
!
      CHARACTER(50)      :: MODE                                          !<-- Flag for global or regional run
!
      LOGICAL            :: GLOBAL                                        !<-- .TRUE. => global ; .FALSE. => regional
!
      INTEGER :: I,J,K,N,NUM_PES,RC,RC_CORE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO SET UP THE DE LAYOUT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve VM from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SET UP PARAMETERS FOR MPI COMMUNICATIONS.
!***  USE ESMF UTILITY TO GET PE IDENTIFICATION AND TOTAL NUMBER OF PEs
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get Task IDs and Count"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,petcount=NUM_PES_TOT                              &  !<-- Total # of tasks (fcst + quilt)
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES=NUM_PES_TOT
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIGURE OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve Config Object from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ESTABLISH THE TASK LAYOUT INCLUDING THE WRITE TASKS.
!***  GET THE GLOBAL MPI COMMUNICATOR AND THE NUMBER OF
!***  OF FORECAST TASKS IN THE I AND J DIRECTIONS AND
!***  GIVE THOSE TO SETUP_SERVERS WHICH WILL SPLIT THE
!***  COMMUNICATOR BETWEEN FORECAST AND QUILT/WRITE TASKS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve Global Communicator"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &
                     ,mpiCommunicator=MPI_INTRA                         &  !<-- The global communicator
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get INPES/JNPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL MPI_COMM_DUP(MPI_INTRA,MPI_INTRA_B,RC)                          !<-- Use a duplicate of the communicator for safety
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SET UP QUILT/WRITE TASK SPECIFICATIONS.
!***  FIRST RETRIEVE THE TASK AND GROUP COUNTS FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get Write Task/Group Info from Config File"
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SEGREGATE THE FORECAST TASKS FROM THE QUILT/WRITE TASKS.
!-----------------------------------------------------------------------
!
      CALL SETUP_SERVERS(MYPE,INPES,JNPES,NUM_PES                       &
                        ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP             &
                        ,MPI_INTRA_B)
!
!***
!***  NOTE: At this point, NUM_PES is the number of Forecast tasks only.
!***
!-----------------------------------------------------------------------
!
      NUM_PES_FCST=INPES*JNPES                                           !<-- Number of forecast tasks
      atm_int_state%NUM_PES_FCST=NUM_PES_FCST                            !<-- Save this for the ATM's Run step
!
!-----------------------------------------------------------------------
!***  CREATE DE LAYOUT BASED ON THE I TASKS BY J TASKS SPECIFIED IN
!***  THE CONFIG FILE.
!***  THIS REFERS ONLY TO FORECAST TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        MY_DE_LAYOUT=ESMF_DELayoutCreate(            VM                 &  !<-- The ESMF virtual machine
                                        ,deCountList=(/INPES,JNPES/)    &  !<-- User-specified I-task by J-task layout
                                        ,rc         =RC)
      ENDIF
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF GRID.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DIMENSIONS OF THE DOMAIN FROM THE CONFIGURE FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get IM,JM,LM from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =IM                            &  !<-- I dimension of full domain
                                  ,label ='im:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =JM                            &  !<-- J dimension of full domain
                                  ,label ='jm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =LM                            &  !<-- Vertical dimension of full domain
                                  ,label ='lm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------
!***  RETRIEVE THE FORECAST DOMAIN MODE FROM THE CONFIG FILE.
!------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get GLOBAL/REGIONAL Mode from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =MODE                          &  !<-- Flag for global (true) or regional (false) run
                                  ,label ='global:'                     &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(MODE)=='true')THEN
        GLOBAL=.true.
        print *,'Initialized as a global run.'
      ELSE
        GLOBAL=.false.
        print *,'Initialized as a regional run.'
      ENDIF
!
!-----------------------------------------------------------------------
!***  IF THIS IS A GLOBAL MODE FORECAST, EXTEND IM AND JM.
!***  THE FIRST DIMENSION OF NCOUNTS IS THE I DIMENSION FOR PARALLELIZATION.
!***  THE SECOND DIMENSION OF NCOUNTS IS THE J DIMENSION.
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN      !<-- Global mode horizontal dimensions.
        NCOUNTS(1)=IM+2
        NCOUNTS(2)=JM+2
      ELSE                !<-- Regional mode horizontal dimensions.
        NCOUNTS(1)=IM
        NCOUNTS(2)=JM
      ENDIF
!
      MAX(1)=NCOUNTS(1)
      MAX(2)=NCOUNTS(2)
!
      MIN(1)=1
      MIN(2)=1
!
!-----------------------------------------------------------------------
!***  NOW CREATE THE ATM GRIDDED COMPONENT's ESMF GRID
!***  FOR THE NMM's INTEGRATION GRID.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Create the ESMF Grid"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GRID_ATM=ESMF_GridCreateShapeTile(regDecomp     =(/INPES,JNPES/)    &  !<-- I x J task layout
                                       ,minIndex      =(/MIN(1),MIN(2)/)  &  !<-- Min indices in I and J
                                       ,maxIndex      =(/MAX(1),MAX(2)/)  &  !<-- Max indices in I and J
                                       ,gridEdgeLWidth=(/0,0/)            &  !<-- Padding, lower edges for noncentered stagger
                                       ,gridEdgeUWidth=(/0,0/)            &  !<-- Padding, upper edges for noncentered stagger
                                       ,name          ="GRID"             &  !<-- Name of the Grid
                                       ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  GET THE LOCAL ARRAY SIZES FOR THE ATM GRID.
!***  ONLY FORECAST TASKS ARE RELEVANT HERE.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                              !<-- Select only fcst tasks
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="CORE_SETUP: Get EMSF Sizes of Local Subdomains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridGet(grid              =GRID_ATM                   &
                         ,localDe           =0                          &
                         ,staggerloc        =ESMF_STAGGERLOC_CENTER     &
                         ,computationalCount=I1                         & !<-- # of local points in I and J on each task
                         ,rc                =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  USING 'computationalCount' FROM ARRAY I1 OBTAINED IN THE
!***  PREVIOUS CALL, GENERATE ALL OF THE LOCAL TASK INDEX LIMITS
!***  FOR ALL FORECAST TASKS.  
!***  THE USER, NOT ESMF, DOES THIS WORK.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        CALL DECOMP(MYPE,INPES,JNPES,NUM_PES_FCST,IM,JM,LM,GLOBAL,I1)
!
        ALLOCATE(atm_int_state%LOCAL_ISTART(0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_IEND  (0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_JSTART(0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_JEND  (0:NUM_PES_FCST-1))
!
        DO N=0,NUM_PES_FCST-1
          atm_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)                    !<-- Starting I for all forecasts' subdomains
          atm_int_state%LOCAL_IEND  (N)=LOCAL_IEND  (N)                    !<-- Ending I for all forecasts' subdomains
          atm_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)                    !<-- Starting J for all forecasts' subdomains
          atm_int_state%LOCAL_JEND  (N)=LOCAL_JEND  (N)                    !<-- Ending J for all forecasts' subdomains
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_SETUP

!-----------------------------------------------------------------------
      SUBROUTINE GFS_SETUP(gc_atm,grid_atmos)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
!
      USE MODULE_INCLUDE
!


!
      type(ESMF_gridcomp),intent(inout) :: gc_atm
      type(ESMF_grid),intent(out)  :: grid_atmos    ! the ESMF grid for the integration attached to

            

!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      type(ESMF_config)            :: cf           ! the config object
!
      type(ESMF_DistGrid)          :: DistGrid_atmos

!
      integer, dimension(2)         :: ncounts     ! parameter array to set up the
                                                   ! size of the 2-d ESMF grid.
      integer, dimension(2)         :: min,max     ! parameter arrays to set up the
                                                   ! start number and the end number of
                                                   ! the ESMF grid in each dimension.
      real(ESMF_kind_r8),dimension(ESMF_maxgriddim) :: mincoords,maxcoords
      integer,dimension(ESMF_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      integer                      :: im,jm,lm                  ! full grid dimensions
      integer                      :: mype,num_pes,num_pes_fcst,num_pes_tot
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator
!
      integer                      :: rc,irtn
      integer                      :: RC_RUN
      logical                      :: global
      character(50)                :: mode 
      rc     =ESMF_success
      RC_RUN=ESMF_success


!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet(gridcomp=gc_atm                      &  !<-- The ATM gridded component
                           ,config  =cf                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)

!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO SET UP THE DE LAYOUT.
!-----------------------------------------------------------------------
!
      CALL ESMF_logwrite("retrieve the config object and vm ",          &
                         ESMF_log_info,rc=RC)
      
!
      CALL ESMF_gridcompget(         gc_atm                             &
                           ,config  =cf                                 &
                           ,vm      =vm                                 &
                           ,rc      =RC)
!
!
!-----------------------------------------------------------------------
!***  set up parameters for mpi communications.
!***  use ESMF utility to get pe identification and total number of pes
!***  (referred to here as nodes.)
!-----------------------------------------------------------------------
!
      CALL ESMF_logwrite("get mype and nodes from vm",ESMF_log_info,rc=RC)
!
      CALL ESMF_vmget(vm                        &  !<-- the virtual machine
                     ,localpet=mype             &  !<-- local pe rank
                     ,petcount=num_pes          &  !<-- total # of tasks
                     ,rc      =RC)
!
      num_pes_tot=num_pes

! Allocate the local index array i2 to store the local size information of the
! ditributed grid.  Information is based per dimension and per De.
!-----------------------------------------------------------------------------
      ALLOCATE(i2(2, num_pes))

!
!***  note: at this point, num_pes is the total number of mpi tasks,
!***        i.e., forecast tasks + quilt tasks.
!
!
!-----------------------------------------------------------------------
!***  establish the task layout including the quilt servers
!***  here in the main gridded component.  get the global
!***  mpi communicator and give it to setup_servers who will
!***  split it between forecast and quilt tasks.
!-----------------------------------------------------------------------
!
      CALL ESMF_vmget(vm                                  &
                     ,mpicommunicator=mpi_intra           &  !<-- the global communicator
                     ,rc             =RC)
!
      CALL mpi_comm_dup(mpi_intra,mpi_intra_b,rc)
!
      CALL ESMF_configgetattribute(cf                      &
                                  ,value =inpes            &  !<-- # of fcst tasks in i direction
                                  ,label ='inpes:'         &
                                  ,rc    =RC)
!
      CALL ESMF_configgetattribute(cf                      &
                                  ,value =jnpes            &  !<-- # of fcst tasks in j direction
                                  ,label ='jnpes:'         &
                                  ,rc    =RC)

      num_pes_fcst=num_pes
!-----------------------------------------------------------------------
!***  create the ESMF grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  the first dimension of ncounts is the i dimension for parallelization.
!***  the second dimension of ncounts is the j dimension.
!-----------------------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =im                            &
                                  ,label ='im:'                         &
                                  ,rc    =RC)
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =jm                            &
                                  ,label ='jm:'                         &
                                  ,rc    =RC)
!
!---------------------------------------------------------
!***  if this is a global mode forecast, extend im and jm.
!***  retrieve the mode from the config file.
!---------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='global:'                     &
                                  ,rc    =RC)
!
      if(trim(mode)=='.true.')then
        GLOBAL=.true.
        print *,'Run as a global model.'
      else
        GLOBAL=.false.
        print *,'Run as a regional model.'
      endif
!
      if(GLOBAL)then      !<-- global mode horizontal dimensions.
        ncounts(1)=im+2
        ncounts(2)=jm+2
      else                !<-- regional mode horizontal dimensions.
        ncounts(1)=im
        ncounts(2)=jm
      endif
!
!
      max(1)=ncounts(1)
      max(2)=ncounts(2)
!
      min(1)=1
      min(2)=1
!
!
!
!-----------------------------------------------------------------------
!***  now create the main gridded component's ESMF grid.
!-----------------------------------------------------------------------
!
! Create the ESMF DistGrid_atmos.
!--------------------------------
!      CALL ESMF_LogWrite("Create DistGrid_atmos", ESMF_LOG_INFO, rc = rc)

      DistGrid_atmos = ESMF_DistGridCreate(minIndex  = min,              &
                                           maxIndex  = max,              &
                                           regDecomp = (/inpes, jnpes/), &
                                           rc        = rc)


! Create the ESMF grid_atmos based on the created ESMF DistGrid_atmos information.
!---------------------------------------------------------------------------------
      CALL ESMF_logwrite("create grid_atmos",ESMF_log_info,rc=RC)
!
      grid_atmos = ESMF_GridCreate(name     = "grid_atmos",   &
                                   distgrid = DistGrid_atmos, &
                                   rc       = rc)

!
!-----------------------------------------------------------------------
!***  attach the ESMF grid to the main gridded component.
!-----------------------------------------------------------------------
!
!      CALL ESMF_gridcompset(         gc_atm                          &
!                           ,grid    =grid_atmos                      &
!                           ,rc      =RC)
!
!-----------------------------------------------------------------------
!***  get the local array sizes for the main grid.
!***  again, only forecast tasks are relevant here.
!-----------------------------------------------------------------------
!
      if(mype<num_pes_fcst)then
          i2 = 0
          CALL ESMF_DistGridGet(DistGrid_atmos, indexCountPDimPDe = i2, rc = rc)
      endif
!
!-----------------------------------------------------------------------
!***  USING 'computationalCount' FROM ARRAY I1 OBTAINED IN THE
!***  PREVIOUS CALL, GENERATE ALL OF THE LOCAL TASK INDEX LIMITS
!***  FOR ALL FORECAST TASKS.
!***  THE USER, NOT ESMF, DOES THIS WORK.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_SETUP

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
