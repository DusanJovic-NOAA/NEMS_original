!-----------------------------------------------------------------------
!
      PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------
!***  Main Program for NEMS.
!***  Define ESMF data types and procedures.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black   - Modified from Wei-yu's version
!   2007-09     Black   - Create the Clock here.
!   2009-08     Colon   - Unified NEM-NMM & NEMS-GFS
!   2009-06-29  Black   - Modified for addition of NMM nesting;
!                         added new ATM Driver Component.
!   2009-09     Lu      - Add w3tage calls for resource statistics
!   2009-08     W. Yang - Ensemble GEFS Concurrency Code.
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
!
!-----------------------------------------------------------------------
!***  USE the model's ATM gridded component module.  Although it
!***  contains the calls to Register and the top level Initialize,
!***  Run, and Finalize, only the Register routine is public.
!-----------------------------------------------------------------------
!
       USE module_ATM_DRIVER_COMP, ONLY: ATM_DRIVER_REGISTER
       USE GEFS_CplComp_ESMFMod,   ONLY: GEFS_CplCompSetServices
!
!-----------------------------------------------------------------------
!***  The following module contains error-checking.
!-----------------------------------------------------------------------
!
       USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!***  This is some information required by the ESMF Log error utility.
!-----------------------------------------------------------------------
!
#include "../../inc/ESMF_LogMacros.inc"
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  Local variables.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      TYPE(ESMF_VM)       :: VM                                          !<-- The ESMF virtual machine,
                                                                         !    which contains and manages
                                                                         !    the computer CPU resource
                                                                         !    for the ESMF grid components.
      TYPE(ESMF_GridComp) :: ATM_DRIVER_COMP                             !<-- The ATM Driver gridded component.
!
      TYPE(ESMF_State) :: ATM_DRIVER_EXP_STATE                        &  !<-- The ATM Driver's export state
                         ,ATM_DRIVER_IMP_STATE                           !<-- The ATM Driver's import state
                                                                         !    as both import and export
                                                                         !    for the ATM component.
      
      TYPE(ESMF_CplComp)               :: CplGEFS       ! the ESMF GFS ensemble run coupler gridded
                                                        ! component.
      TYPE(ESMF_State)                 :: impGEFS       ! the ESMF Coupler import state. ! Right now not used.
      TYPE(ESMF_State)                 :: expGEFS       ! the ESMF Coupler export state. ! Right now not used.
!
      TYPE(ESMF_Clock)    :: CLOCK_MAIN                                  !<-- The ESMF time management clock
!
      TYPE(ESMF_Config)   :: CF_MAIN                                     !<-- The Configure object
!
!-----------------------------------------------------------------------
!
      INTEGER              :: MYPE                                      &  !<-- The local MPI task ID
                           ,YY,MM,DD                                    &  !<-- Time variables for date
                           ,HH,MNS,SEC                                  &  !<-- Time variables for time of day
                           ,NHOURS_FCST                                 &  !<-- Length of forecast in hours
                           ,NSECONDS_FCST                               &  !<-- Length of forecast in seconds
                           ,TIMESTEP_SEC_WHOLE                          &  !<-- Integer part of timestep
                           ,TIMESTEP_SEC_NUMERATOR                      &  !<-- Numerator of fractional part
                           ,TIMESTEP_SEC_DENOMINATOR                       !<-- Denominator of fractional part

      TYPE(ESMF_TimeInterval) :: RUNDURATION                             !<-- The ESMF time. The total forecast hours.
!
      TYPE(ESMF_TimeInterval) :: TIMESTEP                                !<-- The ESMF timestep length (we only need a dummy here)
!
      TYPE(ESMF_Time)         :: STARTTIME                               !<-- The ESMF start time.
!
      INTEGER                 :: RC=ESMF_SUCCESS                         !<-- The running error signal
!
      INTEGER                 :: RC_MAIN                                 !<-- The final value of the
      INTEGER                 :: i

      INTEGER                 :: hh_increase
      INTEGER                 :: hh_start
      INTEGER                 :: hh_final
      INTEGER                 :: Number_start
      INTEGER                 :: Number_final
      LOGICAL                 :: Ens_sps       ! control of stochastic perturbation scheme (sps)
      TYPE(ESMF_LOGICAL)      :: Cpl_flag      ! control of atm running, start from ensemble coupler or not.

!
!-----------------------------------------------------------------------
!***  Declare the ATM gridded component and state names.
!-----------------------------------------------------------------------
!
      CHARACTER(3)                       :: CORE

      CHARACTER(ESMF_MAXSTR)             :: CplCompName
      CHARACTER(ESMF_MAXSTR)             :: impGEFSName
      CHARACTER(ESMF_MAXSTR)             :: expGEFSName
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      include 'fexcp.h'
      call signal(11, xl__trce)
!
!-----------------------------------------------------------------------
!
      CALL w3tagb('nems     ',0000,0000,0000,'np23   ')

!-----------------------------------------------------------------------
!***  Initialize the final error signal.
!-----------------------------------------------------------------------
!
      RC_MAIN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize the ESMF framework. 
!-----------------------------------------------------------------------
! 
      CALL ESMF_Initialize(VM             =VM                           & !<-- The ESMF Virtual Machine
                          ,defaultCalendar=ESMF_CAL_GREGORIAN           & !<-- Set up the default calendar.
                          ,defaultlogtype =ESMF_LOG_MULTI               & !<-- Define multiple log error output file;
                                                                          !    each task has its own log error output file.
                          ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Obtain the local task ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VmGet(vm      =VM                                       &  !<-- The ESMF Virtual Machine
                     ,localpet=MYPE                                     &  !<-- The local MPI task ID
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set up the default log.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Up ESMF Log"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_LogSet(verbose    =ESMF_TRUE                            &
                      ,flush      =ESMF_TRUE                            &
                      ,rootOnly   =ESMF_FALSE                           &
                      ,halt       =ESMF_LOG_HALTNEVER                   & !<-- The job will not stop automatically
                                                                          !    when an ESMF error occurs.
                      ,maxElements=1                                    & ! Maximum number of elements in the log
                                                                          ! before printing them to the log file.
                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_MAIN=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config  =CF_MAIN            & !<-- The Configure object
                              ,filename='configure_file'   & !<-- The name of the configure file 
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the Dynamic core from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Dynamic core from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &  !<-- The config object for the MAIN program
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Label in config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      Ens_sps = .false.

      IF(CORE=='gfs') THEN
!---------------------------------------------------------------------
! Get some parameters of VM and time from config file.
!---------------------------------------------------------------------
          CALL GetTotalMember_EnsembleRunTime(Ens_sps, hh_increase, &
              hh_start, hh_final, rc)
          IF(Ens_sps) THEN
              WRITE(CplCompName, '("GEFS Coupler Grid Component Name")')
              WRITE(impGEFSName, '("GEFS Import State")')
              WRITE(expGEFSName, '("GEFS Export State")')
          END IF
      END IF
!
!-----------------------------------------------------------------------
!***  Create the ATM Driver gridded component which will create and
!***  control the individual ATM forecast components.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM Driver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ATM_DRIVER_COMP=ESMF_GridCompCreate(name        ='ATM DRIVER COMP'  & !<-- ATM Driver component name
                                         ,gridcomptype=ESMF_ATM           & !<-- A Type that the user can query.
                                         ,configFile  ='configure_file'   & !<-- Link the user-created configure file.
                                         ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!CREATE THE ESMF GFS COUPLER GRID COMPONNT,
!-------------------------------------------
      IF(Ens_sps .AND. CORE == 'gfs') THEN
          MESSAGE_CHECK= "Create the GEFS Coupler Grid Component"
!          CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
    
          CplGEFS = ESMF_CplCompCreate(name = CplCompName, rc = rc)

          CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)
      END IF 
!
!-----------------------------------------------------------------------
!***  Register the ATM Driver component's initialize, run and finalize
!***  routines.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register ATM_DRIVER Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(ATM_DRIVER_COMP                       &  !<-- The ATM Driver component
                                   ,ATM_DRIVER_REGISTER                   &  !<-- User's subroutineName
                                   ,RC)

      IF(Ens_sps .AND. CORE == 'gfs') THEN
          CALL ESMF_CplCompSetServices(CplGEFS, GEFS_CplCompSetServices, rc)
      END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the main ESMF Clock.
!***  The Clock is needed for all calls to Init, Run, and Finalize.
!
!***  A timestep is needed to create the Clock but actual timesteps
!***  are handled by ATM_DRIVER because multiple Clocks may be needed
!***  therefore a dummy value is used here.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      TIMESTEP_SEC_WHOLE      =1                                           !<-- Dummy timestep values
      TIMESTEP_SEC_NUMERATOR  =0                                           !
      TIMESTEP_SEC_DENOMINATOR=1                                           !<--
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set up Time Step Interval in Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                   &  !<-- Main Clock's timestep
                               ,s           =TIMESTEP_SEC_WHOLE         &  !<-- Whole part of timestep
                               ,sn          =TIMESTEP_SEC_NUMERATOR     &  !<-- Numerator of fractional part
                               ,sd          =TIMESTEP_SEC_DENOMINATOR   &  !<-- Denominator of fractional part
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the start time in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Year from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =YY                            &
                                  ,label ='start_year:'                 &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Month from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MM                            &
                                  ,LABEL ='start_month:'                &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Day from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =DD                            &
                                  ,label ='start_day:'                  &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Hour from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =HH                            &
                                  ,label ='start_hour:'                 &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Minute from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MNS                           &
                                  ,label ='start_minute:'               &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Second from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =SEC                           &
                                  ,label ='start_second:'               &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Start Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeSet(time=STARTTIME                                  &  !<-- The start time of the forecast (ESMF)
                       ,yy  =YY                                         &  !<-- Year from config file
                       ,mm  =MM                                         &  !<-- Month from config file
                       ,dd  =DD                                         &  !<-- Day from config file
                       ,h   =HH                                         &  !<-- Hour from config file
                       ,m   =MNS                                        &  !<-- Minute from config file
                       ,s   =SEC                                        &  !<-- Second from config file
                       ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the run duration in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Forecast Length from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(Ens_sps) THEN
          CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                      ,value =NHOURS_FCST                   &
                                      ,label ='nhours_fcst1:'               &
                                      ,rc    =RC)
      ELSE
          CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                      ,value =NHOURS_FCST                   &
                                      ,label ='nhours_fcst:'                &
                                      ,rc    =RC)
      END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NSECONDS_FCST=NHOURS_FCST*3600                                       !<-- The forecast length (sec) (REAL)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Length"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION                &  !<-- The forecast length (s) (ESMF)
                               ,s           =NSECONDS_FCST              &
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now the Main Clock can be created.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CLOCK_MAIN=ESMF_ClockCreate(name       ='CLOCK_MAIN'              &  !<-- The top-level ESMF Clock
                                 ,timestep   =TIMESTEP                  &  !<-- Timestep needed by the Clock (ESMF)
                                 ,starttime  =STARTTIME                 &  !<-- The integration start time (ESMF)
                                 ,runduration=RUNDURATION               &  !<-- The integration duration (ESMF)
                                 ,rc         =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the ATM Driver component's import/export states.
!***  Currently they are not required to perform an actual function.
!***  When ATM Driver is coupled to another component then the
!***  import/export states will actually be utilized.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM Driver Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ATM_DRIVER_IMP_STATE=ESMF_StateCreate(statename='ATM Driver Import State' &
                                           ,rc       =RC)
!
      ATM_DRIVER_EXP_STATE=ESMF_StateCreate(statename='ATM Driver Export State' &
                                           ,rc       =RC)
      
      IF(CORE == 'gfs') THEN
          Cpl_flag = ESMF_FALSE
          CALL ESMF_AttributeSet(ATM_DRIVER_IMP_STATE,                     &
                                 'Cpl_flag', Cpl_flag, rc = rc)

          IF(Ens_sps) THEN
              impGEFS = ESMF_StateCreate(statename = impGEFSName,              &
                                    statetype      = ESMF_STATE_IMPORT,        &
                                    rc             = rc)

              expGEFS = ESMF_StateCreate(statename = expGEFSName,              &
                                    statetype      = ESMF_STATE_EXPORT,        &
                                    rc             = rc)
! Add the GFS ESMF states as the nested states into the GEFS states.
!-------------------------------------------------------------------
              MESSAGE_CHECK= "Add the GFS states into the GEFS states"
!              CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

              CALL ESMF_StateAdd(impGEFS, ATM_DRIVER_EXP_STATE, rc = rc)
              CALL ESMF_StateAdd(expGEFS, ATM_DRIVER_IMP_STATE, rc = rc)

              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
          END IF
      END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the INITIALIZE step for the ATM Driver component.
!***  The initialize routine that is called here as well as the
!***  run and finalize routines invoked below are those specified
!***  in the Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Driver Component Initialize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =ATM_DRIVER_COMP          &  !<-- The ATM Driver component
                                  ,importstate=ATM_DRIVER_IMP_STATE     &  !<-- The ATM Driver's import state
                                  ,exportstate=ATM_DRIVER_EXP_STATE     &  !<-- The ATM Driver's export state
                                  ,clock      =CLOCK_MAIN               &  !<-- The ESMF clock
                                  ,phase      =ESMF_SINGLEPHASE         &
                                  ,rc         =RC)
!
!!INITIALIZE THE GFS COUPLER GRID COMPONENT:
!
      IF(Ens_sps .AND. CORE == 'gfs') THEN
          MESSAGE_CHECK = "Calling the GEFS COUPLER Initialize"
!          CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

          CALL ESMF_CplCompInitialize(CplGEFS,                        &
                                      importstate = impGEFS,          &
                                      exportstate = expGEFS,          &
                                      clock       = CLOCK_MAIN,       &
                                      phase       = ESMF_SINGLEPHASE, &
                                      rc          = rc)
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
      END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the RUN step for the ATM Driver component.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Driver Component Run Step"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompRun(gridcomp   =ATM_DRIVER_COMP                 &  !<-- The ATM Driver component
                           ,importstate=ATM_DRIVER_IMP_STATE            &  !<-- The ATM Driver's import state
                           ,exportstate=ATM_DRIVER_EXP_STATE            &  !<-- The ATM Driver's export state
                           ,clock      =CLOCK_MAIN                      &  !<-- The ESMF clock
                           ,phase      =ESMF_SINGLEPHASE                &
                           ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      IF(Ens_sps .AND. CORE == 'gfs') THEN
         Number_start = hh_start / hh_increase + 1
         Number_final = hh_final / hh_increase - 1
         PRINT *, 'DHOUCoup ', Number_start, Number_final, hh_start, hh_final, hh_increase

         DO i = Number_start, Number_final
!
!!RUNNING THE COUPLER GRID COMPONENT:
!
             PRINT *, 'DHOU CPL NO i=',i
             MESSAGE_CHECK = "Calling the GFS Coupler Run"
             CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

             CALL ESMF_CplCompRun(CplGEFS,                        &
                                  importstate = impGEFS,          &
                                  exportstate = expGEFS,          &
                                  clock       = CLOCK_MAIN,       &
                                  phase       = ESMF_SINGLEPHASE, &
                                  rc          = rc)
             CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)
!
!!RUNNING THE GFS GRID COMPONENT AGAIN:
!
             MESSAGE_CHECK = "Get runDuration from clock"
!             CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

             CALL ESMF_VMBarrier(vm, rc = rc)

             CALL ESMF_ClockGet(CLOCK_MAIN, runDuration = runDuration, rc = rc)

             CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)

             MESSAGE_CHECK = "Adjust clock"
!             CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

             CALL ESMF_TimeIntervalGet(runDuration, h = hh,       rc = rc)
             hh = hh + hh_increase

             CALL ESMF_TimeIntervalSet(runDuration, h = hh,       rc = rc)
             CALL ESMF_ClockSet(CLOCK_MAIN, runDuration = runDuration, rc = rc)

             CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)

             MESSAGE_CHECK = "Run the GFS RUN"
!             CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

             CALL ESMF_GridCompRun(gridcomp   =ATM_DRIVER_COMP                 &  !<-- The ATM Driver component
                                  ,importstate=ATM_DRIVER_IMP_STATE            &  !<-- The ATM Driver's import state
                                  ,exportstate=ATM_DRIVER_EXP_STATE            &  !<-- The ATM Driver's export state
                                  ,clock      =CLOCK_MAIN                      &  !<-- The ESMF clock
                                  ,phase      =ESMF_SINGLEPHASE                &
                                  ,rc         =RC)
             CALL ERR_MSG(RC,MESSAGE_CHECK, RC_MAIN)
             IF(i < Number_final) THEN
                 PRINT*, 'Complete GFS Run Cycle = ', i + 1
             END IF
         END DO
     END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the FINALIZE step for the ATM Driver component.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Driver Component Finalize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =ATM_DRIVER_COMP            &  !<-- The ATM Driver component
                                ,importstate=ATM_DRIVER_IMP_STATE       &  !<-- The ATM Driver's import state
                                ,exportstate=ATM_DRIVER_EXP_STATE       &  !<-- The ATM Driver's export state
                                ,clock      =CLOCK_MAIN                 &  !<-- The ESMF clock
                                ,phase      =ESMF_SINGLEPHASE           &
                                ,rc         =RC)
!
!!AFTER RUNNING, FINALIZE THE COUPLER COMPONENT:
!
      IF(Ens_sps .AND. CORE == 'gfs') THEN
          MESSAGE_CHECK = "Calling the GEFS Coupler Finalize"
!          CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

          CALL ESMF_CplCompFinalize  (CplGEFS,                        &
                                      importstate = impGEFS,          &
                                      exportstate = expGEFS,          &
                                      clock       = CLOCK_MAIN,       &
                                      phase       = ESMF_SINGLEPHASE, &
                                      rc          = rc)

          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
      END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_MAIN                           &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Configure object.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigDestroy(config=CF_MAIN                            &
                             ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      MESSAGE_CHECK = "Destroy the ESMF states"
!      CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

      CALL ESMF_StateDestroy(ATM_DRIVER_IMP_STATE, rc = rc)
      CALL ESMF_StateDestroy(ATM_DRIVER_EXP_STATE, rc = rc)

      IF(Ens_sps .AND. CORE == 'gfs') THEN
          CALL ESMF_StateDestroy(impGEFS, rc = rc)
          CALL ESMF_StateDestroy(expGEFS, rc = rc)
      END IF

      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)

      MESSAGE_CHECK = "Destroy ESMF Grid Comp and Cpl Comp"
!      CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)

      CALL ESMF_GridCompDestroy(ATM_DRIVER_COMP, rc = rc)

      IF(Ens_sps .AND. CORE == 'gfs') THEN
          CALL ESMF_CplCompDestroy(CplGEFS, rc= rc)
      END IF

      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_MAIN)

!-----------------------------------------------------------------------
!***  Shut down the ESMF system.
!-----------------------------------------------------------------------
!
      CALL ESMF_Finalize()
!
      IF(RC_MAIN==ESMF_SUCCESS)THEN
        WRITE(0,*)'MODEL FORECAST RUN SUCCEEDED'
      ELSE
        WRITE(0,*)'MODEL FORECAST RUN FAILED  RC_MAIN=',RC_MAIN
      ENDIF
!
!-----------------------------------------------------------------------
!
      if(mype==0)call w3tage('nems     ')
!
!-----------------------------------------------------------------------
!
      END PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------





 SUBROUTINE GetTotalMember_EnsembleRunTime(Ens_sps,      &
                                           hh_increase,  &
                                           hh_start,     &
                                           hh_final,     &
                                           rc)

 USE ESMF_Mod

 IMPLICIT none

 LOGICAL,                        INTENT(out) :: Ens_sps
 INTEGER,                        INTENT(out) :: hh_increase
 INTEGER,                        INTENT(out) :: hh_start
 INTEGER,                        INTENT(out) :: hh_final
 INTEGER,                        INTENT(out) :: rc
 TYPE(ESMF_Config)                           :: Cf
 CHARACTER(ESMF_MAXSTR)                      :: Cf_fname

 rc       = ESMF_SUCCESS
 Cf       = ESMF_ConfigCreate(rc = rc)
 Cf_fname = 'configure_file'

 CALL ESMF_ConfigLoadFile(Cf, Cf_fname, rc = rc)

 CALL ESMF_ConfigGetAttribute(Cf,                      &
                              Ens_sps,                 &
                              label = 'ENS_SPS:',      &
                              rc    = rc)

 IF(Ens_sps) THEN
     CALL ESMF_ConfigGetAttribute(Cf,                      &
                                  hh_increase,             &
                                  label = 'HH_INCREASE:',  &
                                  rc    = rc)

     CALL ESMF_ConfigGetAttribute(Cf,                      &
                                  hh_start,                &
                                  label = 'HH_START:',     &
                                  rc    = rc)

     CALL ESMF_ConfigGetAttribute(Cf,                      &
                                  hh_final,                &
                                  label = 'HH_FINAL:',     &
                                  rc    = rc)
 END IF

 END SUBROUTINE GetTotalMember_EnsembleRunTime
