#include "./ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE module_NEMS_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the NEMS component.
!-----------------------------------------------------------------------
!
!***  The NEMS component lies in the hierarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!              /|\
!             / | \
!          ATM/OCN/ICE/WAV/LND/IPM/HYD .. components
!          |    |
!          |    |
!          |    |
!          |    (MOM5, HYCOM, etc.)
!          |
!          CORE component (GSM, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!  2011-05-11  Theurich & Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10-04  Yang  - Modified for using the ESMF 5.2.0r library.
!  2013-07     Theurich - Macro based ESMF error handling
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE module_NEMS_INTERNAL_STATE,ONLY: NEMS_INTERNAL_STATE          &
                                          ,WRAP_NEMS_INTERNAL_STATE
!
      USE ENS_CplComp_ESMFMod,ONLY: ENS_CplCompSetServices
!
      USE module_EARTH_GRID_COMP
!
      USE module_NEMS_UTILS,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NEMS_REGISTER
!
!-----------------------------------------------------------------------
!
      INTEGER :: MEMBER_ID                                              &
                ,TOTAL_MEMBER
!
      INTEGER :: HH_INCREASE                                            &
                ,HH_START                                               &
                ,HH_FINAL
!
      INTEGER :: NUMBER_START                                           &
                ,NUMBER_FINAL
!
      INTEGER,DIMENSION(:),   ALLOCATABLE :: PE_MEMBER                     !<-- Tasks for each member
      INTEGER,DIMENSION(:, :),ALLOCATABLE :: PETLIST                       !<-- Task list for each member
!
      LOGICAL :: ENS_SPS                                                   !<-- Control of Stochastic Perturbation Scheme (SPS)
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: IMP_EARTH_NAME    !<-- Import state name of the EARTH components
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: EXP_EARTH_NAME    !<-- Export state name of the EARTH components
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: GC_EARTH_NAME     !<-- Name of the EARTH component
!
      TYPE(NEMS_INTERNAL_STATE),POINTER,SAVE :: NEMS_INT_STATE
      TYPE(WRAP_NEMS_INTERNAL_STATE)   ,SAVE :: WRAP
!
      TYPE(ESMF_Clock), SAVE :: CLOCK_NEMS                                 !<-- The ESMF Clock of the NEMS component
      TYPE(ESMF_Config),SAVE :: CF_NEMS                                    !<-- The configure object of the NEMS component
      TYPE(ESMF_VM),    SAVE :: VM_GLOBAL
      TYPE(ESMF_TIME),  SAVE :: STARTTIME, CURRTIME
!
      TYPE(ESMF_CplComp),SAVE :: ENS_CPL_COMP                              !<-- Ensemble Coupler components
      TYPE(ESMF_State),  SAVE :: ENS_CPL_IMP_STATE                         !<-- Import state of the Ensemble Coupler component
      TYPE(ESMF_State),  SAVE :: ENS_CPL_EXP_STATE                         !<-- Export state of the Ensemble Coupler component
!
      TYPE(ESMF_GridComp),DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_GRID_COMP !<-- EARTH components for each member
      TYPE(ESMF_State),   DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_IMP_STATE !<-- Import state of the EARTH component
      TYPE(ESMF_State),   DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_EXP_STATE !<-- Export state of the EARTH component
!
!-----------------------------------------------------------------------
!
      CONTAINS

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_REGISTER(NEMS_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NEMS_GRID_COMP                                !<-- The NEMS gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Error return code
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
      RC_REG = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Register the NEMS Initialize, Run, and Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Subroutine type (Initialize)
                                     ,NEMS_INITIALIZE                   &  !<-- User's subroutine name
                                     ,phase=1                           &
                                     ,rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_METHOD_RUN                   &  !<-- Subroutine type (Run)
                                     ,NEMS_RUN                          &  !<-- User's subroutine name
                                     ,phase=1                           &
                                     ,rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_METHOD_FINALIZE              &  !<-- Subroutine type (Finalize)
                                     ,NEMS_FINALIZE                     &  !<-- User's subroutine name
                                     ,phase=1                           &
                                     ,rc=RC)
      ESMF_ERR_RETURN(RC,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END SUBROUTINE NEMS_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_INITIALIZE(NEMS_GRID_COMP                         &
                                ,IMP_STATE                              &
                                ,EXP_STATE                              &
                                ,CLOCK_MAIN                             &
                                ,RC_INIT)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NEMS_GRID_COMP                                !<-- The NEMS component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NEMS import state
                         ,EXP_STATE                                        !<-- The NEMS export state
!
      TYPE(ESMF_Clock) :: CLOCK_MAIN                                       !<-- The main Clock
!
      INTEGER,INTENT(OUT) :: RC_INIT                                       !<-- Error return code
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION, restartOffset
!
      CHARACTER(20) :: PELAB
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: EARTH_COMP_NAME &  !<-- Names of each member's EARTH component
                                                        ,IMP_EARTH_NAME  &  !<-- Import state name of the EARTH components
                                                        ,EXP_EARTH_NAME     !<-- Export state name of the EARTH components
!
      INTEGER :: I,IJ,J,RC,RC_USER
!
      INTEGER :: MYPE_GLOBAL                                            &
                ,NHOURS_FCST                                            &
                ,NSECONDS_FCST                                          &
                ,PE_MAX                                                 &
                ,TASKS                                                  &
                ,fhrot
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: PETLIST                        !<-- Task list for each ensemble member
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_INIT = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the main Clock within
!***  the NEMS component.
!-----------------------------------------------------------------------
!
      CLOCK_NEMS=CLOCK_MAIN
!
!-----------------------------------------------------------------------
!***  What is the start time on the NEMS clock?
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract the start time of the NEMS clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_ClockGet(clock     = CLOCK_NEMS                         &
                        ,startTime = STARTTIME                          &
                        ,rc = RC)

      ESMF_ERR_RETURN(RC,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Allocate the NEMS component's internal state, point at it,
!***  and attach it to the NEMS component.
!-----------------------------------------------------------------------
!
      ALLOCATE(NEMS_INT_STATE,stat=RC)
      wrap%NEMS_INT_STATE=>NEMS_INT_STATE
!
      CALL ESMF_GridCompSetInternalState(NEMS_GRID_COMP                 &  !<--The NEMS component
                                        ,WRAP                           &  !<-- Pointer to the NEMS internal state
                                        ,RC)     
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Get the global VM (Virtual Machine).
!***  Obtain the total task count and the local task ID.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGetGlobal(vm = VM_GLOBAL                              &  !<-- The ESMF global Virtual Machine
                           ,rc = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)
!
      CALL ESMF_VMGet(vm       = VM_GLOBAL                              &  !<-- The ESMF global Virtual Machine
                     ,pecount  = TASKS                                  &  !<-- Total # of MPI tasks
                     ,localpet = MYPE_GLOBAL                            &  !<-- This task's global rank
                     ,rc       = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)

!jm 20140818
!jm   call rsl_lite_error_dup1( mype_global )                              !<-- Enable this for per-MPI task output
!jm                                                                        !<-- obtain other bits from JM

!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the NEMS configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the NEMS Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_NEMS=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config   = CF_NEMS                       &  !<-- The configure object
                              ,filename = 'model_configure'             &  !<-- The name of the configure file
                              ,rc       = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)
!-----------------------------------------------------------------------
!***  Get the ensemble stochastic coupling flag from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the Ensemble Stochastic Coupling Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config = CF_NEMS                     &  !<-- The NEMS configure object
                                  ,value  = ENS_SPS                     &  !<-- Value of control flag for 
                                                                           !    stochastic perturbation scheme
                                  ,label  = 'ENS_SPS:'                  &  !<-- Flag's label in configure file
                                  ,rc     = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Reset the run duration in the NEMS Clock for the ensemble 
!***  stochastic perturbation case.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NEMS: Extract Forecast Length from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_NEMS                     &
                                    ,value =NHOURS_FCST                 &
                                    ,label ='nhours_fcst1:'             &
                                    ,rc    =RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NSECONDS_FCST=NHOURS_FCST*3600                                     !<-- The forecast length (sec) (Integer)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NEMS: Set the Forecast Length"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION              &  !<-- The forecast length (s) (ESMF)
                                 ,s           =NSECONDS_FCST            &  !<-- The forecast length (s) (Integer)
                                 ,rc          =RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
        CALL ESMF_ClockSet(clock       = CLOCK_NEMS                     &  !<-- The NEMS Clock
                          ,runDuration = RUNDURATION                    &  !<-- The forecast length (s) (ESMF)
                          ,rc          = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
      END IF
!
!-----------------------------------------------------------------------
!***  Extract the total number of EARTH ensemble members.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the Total Number of the EARTH Members from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config = CF_NEMS                     &  !<-- The NEMS configure object
                                  ,value  = TOTAL_MEMBER                &  !<-- Total # of ensemble members
                                  ,label  = 'total_member:'             &  !<-- Flag in configure file 
                                  ,rc     = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Allocate a standard set of arrays for each ensemble member.
!-----------------------------------------------------------------------
!
      ALLOCATE(EARTH_GRID_COMP (TOTAL_MEMBER))
      ALLOCATE(EARTH_IMP_STATE (TOTAL_MEMBER))
      ALLOCATE(EARTH_EXP_STATE (TOTAL_MEMBER))
      ALLOCATE(EARTH_COMP_NAME (TOTAL_MEMBER))
      ALLOCATE(IMP_EARTH_NAME  (TOTAL_MEMBER))
      ALLOCATE(EXP_EARTH_NAME  (TOTAL_MEMBER))
      ALLOCATE(PE_MEMBER       (TOTAL_MEMBER))
!
!-----------------------------------------------------------------------
!***  For each member create the names of the EARTH components and 
!***  ESMF states then fill in the task information.
!-----------------------------------------------------------------------
!
      PE_MEMBER = 0
!
      DO I = 1, TOTAL_MEMBER
!
        WRITE(EARTH_COMP_NAME(I), '("EARTH grid component", I2.2)') I
        WRITE(IMP_EARTH_NAME (I), '("EARTH import state",   I2.2)') I
        WRITE(EXP_EARTH_NAME (I), '("EARTH export state",   I2.2)') I
!
        WRITE(PELAB, '("PE_MEMBER", I2.2, ":")') I
! 
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = PE_MEMBER(I)              &
                                    ,label  = TRIM(PELAB)               &
                                    ,rc     = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
        IF(PE_MEMBER(I) == 0) PE_MEMBER(i) = TASKS / TOTAL_MEMBER
!
      END DO
!
      PE_MAX = 1
!
      DO I = 1, TOTAL_MEMBER
        PE_MAX = MAX(PE_MAX, PE_MEMBER(I))
      END DO
!
!-----------------------------------------------------------------------
!***  Set up the PE list.
!-----------------------------------------------------------------------
!
      ALLOCATE(PETLIST(1:PE_MAX, 1:TOTAL_MEMBER))
!
      IJ = 0
!
      DO J = 1, TOTAL_MEMBER
!
        DO I = 1, PE_MEMBER(J)
          PETLIST(I, J) = IJ
!
          IF(MYPE_GLOBAL == IJ) THEN
            MEMBER_ID = J
          END IF
!
          IJ = IJ + 1
        END DO
!
      END DO
!
!-----------------------------------------------------------------------
!***  Get the clock parameters for the ensemble stochastic perturbation 
!***  cycles from the config file.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the Ensemble Clock Parameters from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_increase               &
                                    ,label  = 'HH_INCREASE:'            &
                                    ,rc     = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_start                  &
                                    ,label  = 'HH_START:'               &
                                    ,rc     = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_final                  &
                                    ,label  = 'HH_FINAL:'               &
                                    ,rc     = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
        NUMBER_START = HH_START / HH_INCREASE + 1
        NUMBER_FINAL = HH_FINAL / HH_INCREASE - 1
!
        PRINT *, 'ENS Coup ', NUMBER_START, NUMBER_FINAL                &
               , HH_START, HH_FINAL, HH_INCREASE
!
      END IF
!
!-----------------------------------------------------------------------
!***  Create the EARTH grid components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create EARTH grid Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
        EARTH_GRID_COMP(i) = ESMF_GridCompCreate (                      &
                             name         = EARTH_COMP_NAME(I)          &  !<-- Name of element I of the EARTH component array
                            ,petlist      = PETLIST(1:PE_MEMBER(I), I)  &  !<-- Element I's PE list
                            ,config       = CF_NEMS                     &  !<-- Associate the NEMS config object with this element
                            ,rc           = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      END DO
!
!-----------------------------------------------------------------------
!***  Create the Ensemble Coupler component inside of
!***  the NEMS internal state.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
        ENS_CPL_COMP=ESMF_CplCompCreate(name = "ENS Cpl component"      &  
                                       ,rc   = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      END IF
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines of
!***  each element in the EARTH component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register EARTH Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
        CALL ESMF_GridCompSetServices(EARTH_GRID_COMP(I)                &  !<-- The EARTH gridded components
                                     ,EARTH_REGISTER                    &  !<-- User's name for the Register routine
                                     ,rc=RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      END DO
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines of
!***  the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Ensemble Coupler Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompSetServices(ENS_CPL_COMP                       &
                                    ,ENS_CplCompSetServices             &  !<-- The user's name for the Register routine
                                    ,rc=RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
      END IF
!
!-----------------------------------------------------------------------
!***  Create the EARTH import and export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the EARTH import states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1,TOTAL_MEMBER
        EARTH_IMP_STATE(I) = ESMF_StateCreate(                          &
                                         name = IMP_EARTH_NAME(I)  &
                                        ,stateintent = ESMF_STATEINTENT_IMPORT  &
                                        ,rc        = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      END DO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the EARTH export states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1,TOTAL_MEMBER
        EARTH_EXP_STATE(I) = ESMF_StateCreate(                                 &
                                         name   = EXP_EARTH_NAME(I)       &
                                        ,stateintent = ESMF_STATEINTENT_EXPORT &
                                        ,rc          = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      END DO
!
!-----------------------------------------------------------------------
!***  Create the Ensemble Coupler import and export states.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Ensemble Coupler import state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENS_CPL_IMP_STATE=ESMF_StateCreate(name   = "ENS_CPL_Import"         &
                                          ,stateintent = ESMF_STATEINTENT_IMPORT  &
                                          ,rc          = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Ensemble Coupler export state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENS_CPL_EXP_STATE=ESMF_StateCreate(name = "ENS_CPL_Export"           &
                                          ,stateintent = ESMF_STATEINTENT_EXPORT  &
                                          ,rc        = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Nest the EARTH export/import states into the import/export states
!***  of the ensemble coupler.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK= "Add the EARTH states into the ENS_CPL states"
!       CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
        DO I = 1, TOTAL_MEMBER
          IF(MEMBER_ID == I) THEN
            CALL ESMF_StateAddReplace(ENS_CPL_IMP_STATE,(/EARTH_EXP_STATE(I)/), rc = RC)
            ESMF_ERR_RETURN(RC,RC_INIT)
            CALL ESMF_StateAddReplace(ENS_CPL_EXP_STATE,(/EARTH_IMP_STATE(I)/), rc = RC)
            ESMF_ERR_RETURN(RC,RC_INIT)
          END IF
        END DO
!
      END IF

!
!-----------------------------------------------------------------------
!***  Adjust the currTime of the NEMS clock: CLOCK_NEMS
!***  if the fhrot is > 0
!***  This will correctly set the NEMS and EARTH clocks in case of
!***  Restart-From-History.
!-----------------------------------------------------------------------

      CALL ESMF_ConfigGetAttribute(config = CF_NEMS &
                                   ,value  = fhrot &
                                   ,label  = 'fhrot:' &
                                   ,default = 0 &
                                   ,rc     = RC)
      ESMF_ERR_RETURN(RC,RC_INIT)
      if (fhrot > 0) then
        CALL ESMF_TimeIntervalSet(restartOffset, h=fhrot,rc=RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
        CURRTIME = STARTTIME + restartOffset
        call ESMF_ClockSet(CLOCK_NEMS, currTime=CURRTIME, rc=RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
      endif
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of each element of the EARTH
!***  component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Initialize step of the EARTH component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompInitialize(gridcomp    = EARTH_GRID_COMP(I)  &
                                      ,importState = EARTH_IMP_STATE(I)  &
                                      ,exportState = EARTH_EXP_STATE(I)  &
                                      ,clock       = CLOCK_NEMS          &
                                      ,phase       = 1                   &
                                      ,userRc      = RC_USER             &
                                      ,rc          = RC)
          ESMF_ERR_RETURN(RC,RC_INIT)
          ESMF_ERR_RETURN(RC_USER,RC_INIT)
        END IF
!
      END DO
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Execute the Initialize step of the Ensemble Coupler component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompInitialize(cplcomp    =ENS_CPL_COMP            &
                                   ,importState=ENS_CPL_IMP_STATE       &
                                   ,exportState=ENS_CPL_EXP_STATE       &
                                   ,clock      =CLOCK_NEMS              &
                                   ,phase      =1                       &
                                   ,userRc     =RC_USER                 &
                                   ,rc         =RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
        ESMF_ERR_RETURN(RC_USER,RC_INIT)
!
      END IF
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!=-=-=-=-=-=-=-=-=-=-=-= ADDED BY SWD 27Nov2017  =-=-=-=-=-=-=-=-=-=-=-=
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!-----------------------------------------------------------------------
!***    Reconcile 'EARTH_exp_state' across all pets. [It was created only on src_petlist.]
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        !MESSAGE_CHECK = "  DEMO_init:  Reconciling 'src_exp_state' across global VM"
        MESSAGE_CHECK = "  NEMS_init:  Reconciling 'earth_exp_state' across global VM"
        CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO i = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == i) THEN
        CALL ESMF_StateReconcile(state=earth_exp_state(i)     &
                                ,rc=rc)
!
          ESMF_ERR_RETURN(RC,RC_INIT)
!
          ESMF_ERR_RETURN(RC_USER,RC_INIT)
!
        END IF
!
      END DO
!
        CALL ESMF_VMBarrier(vm=VM_GLOBAL)
!
!-----------------------------------------------------------------------
!***  Get the global VM (Virtual Machine).
!***  Obtain the total task count and the local task ID.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !MESSAGE_CHECK="  DEMO_init:  Get the VM of the demo Grid Component"
      MESSAGE_CHECK="  NEMS_init:  Get the VM of the NEMS Grid Component"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp = nems_grid_comp          &
                           ,vm       = vm_global               &
                           ,rc       = rc)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !MESSAGE_CHECK="  DEMO_init:  Get the #_global_tasks and mype_global from the global VM"
      MESSAGE_CHECK="  NEMS_init:  Get the #_global_tasks and mype_global from the global VM"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_VMGet(vm       = vm_global       &  !<-- The ESMF global Virtual Machine
                     ,pecount  = global_tasks    &  !<-- Total # of MPI tasks
                     ,localpet = mype_global     &  !<-- This task's global rank
                     ,rc       = RC)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***    Get Attributes from the configure file to create the PETlist for the dst_gridcomp
!-----------------------------------------------------------------------
!
      CALL ESMF_VMBarrier(vm=vm_global)
!
!-----------------------------------------------------------------------
!***  Extract the total number of (independent) nested domains
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      !MESSAGE_CHECK="  DEMO_init:  Extract the Total Number of Storms (nested domains) from Config File"
      MESSAGE_CHECK="  NEMS_init:  Extract the Total Number of Storms (nested domains) from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      cf_main=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config   = cf_main         &  !<-- The configure object
                              ,filename = 'demo.conf'     &  !<-- The name of the configure file
                              ,rc       = RC)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
!
      CALL ESMF_ConfigGetAttribute(config = cf_main                     &  !<-- The DEMO configure object
                                  ,value  = nstorms                     &  !<-- Total # of ensemble members
                                  ,label  = 'total_storms:'             &  !<-- Flag in configure file
                                  ,rc     = RC)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Allocate a set of arrays for each storm.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        !MESSAGE_CHECK="  DEMO_init:  Allocate array to create unique 'dst gridded component' for each storm"
        MESSAGE_CHECK="  NEMS_init:  Allocate array to create unique 'dst gridded component' for each storm"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(dst_gridcomp  (nstorms))
      ALLOCATE(dst_imp_state (nstorms))
      ALLOCATE(dst_exp_state (nstorms))
      ALLOCATE(storm_name    (nstorms))
      ALLOCATE(configfile    (nstorms))
      ALLOCATE(cf_dst        (nstorms))
!
      ALLOCATE(src2dst_cplcomp   (nstorms))
      ALLOCATE(src2dst_imp_state (nstorms))
      ALLOCATE(src2dst_exp_state (nstorms))
!
      do i=1,nstorms
        write(storm_name(i), '("storm0", I2.2)') i
        write(configfile(i), '("dst.conf_0", I2.2)') i
      end do
!
        CALL ESMF_VMBarrier(vm=vm_global)
!
!-----------------------------------------------------------------------
!***  Create the configure object for (each of?) the 'dst gridded components'.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="  NEMS_init:  Create the dst.conf object"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        cf_dst(i)=ESMF_ConfigCreate(rc=RC)
      end do
!
        CALL ESMF_VMBarrier(vm=vm_global)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="  NEMS_init:  Load the dst.conf object"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        CALL ESMF_ConfigLoadFile(config   = cf_dst(i)           &  !<-- The configure object
                                ,filename = 'dst.conf_001'      &  !<-- The name of the configure file
                                ,rc       = RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="  NEMS_init:  Extract Number of compute tasks in the I direction for each dst gridcomp"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=cf_dst(i)             &
                                    ,value = dst_inpes            &
                                    ,label ='dst_inpes:'          &
                                    ,rc    =RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="  DEMO_init:  Extract Number of compute tasks in the J direction for each dst gridcomp"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=cf_dst(i)             &
                                    ,value = dst_jnpes            &
                                    ,label ='dst_jnpes:'          &
                                    ,rc    =RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Create petlist for the DST Gridded Component
!-----------------------------------------------------------------------
!
      dst_pets=dst_inpes*dst_jnpes
!
      allocate(dst_petlist(nstorms,dst_pets))
!
      !*define petlist for each storm
      do i=1,nstorms
        do j=1,dst_pets
          dst_petlist(i,j)=(src_petlist_count-1)+((i-1)*dst_pets)+j
        end do
      end do
!
!-----------------------------------------------------------------------
!***  Create the 'dst' Gridded Components
!***
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Create dst Gridded Components"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        dst_gridcomp(i) = ESMF_GridCompCreate(name=storm_name(i)         &  !<-- Name of the DST component
                                             ,petlist=dst_petlist(i,:)   &  !<-- Element I's PE list
                                             ,configFile=configfile(i)   &  !<-- Associate the name of the config $
                                             ,rc=RC)
!
        if(rc.ne.0)write(0,*)"Check total processors in run script vs. # domains."
        ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines
!***  of the 'dst' gridded component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Register dst Init, Run, Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        CALL ESMF_GridCompSetServices(gridcomp=dst_gridcomp(i)         &  !<-- The ith dst gridded component
                                     ,userRoutine=dst_register        &  !<-- User's name for the Register routine
                                     ,rc=RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Create the DST component's import/export states.
!***  Currently they are not required to perform an actual function.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="  NEMS_init:  Create the (array of) DST import & export states"
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        write(imp_state_name, '("storm0",I2.2,"_imp_state")') i
!
        dst_imp_state(i) = ESMF_StateCreate(name=imp_state_name     &
                                           ,stateintent=ESMF_STATEINTENT_IMPORT  &
                                           ,rc=RC)
!
        write(exp_state_name, '("storm0",I2.2,"_exp_state")') i
!
        dst_exp_state(i) = ESMF_StateCreate(name=exp_state_name     &
                                           ,stateintent=ESMF_STATEINTENT_EXPORT  &
                                           ,rc=RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***    Associate a config file with each 'dst' gridded component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        MESSAGE_CHECK = "  DEMO: Associate a config file with each 'dst' gridded component"
!        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompGet(dst_gridcomp(i)                &
                               ,name=myname                    &
                               ,config=cf_dst(i)               &
                               ,configFileIsPresent=cfpresent  &
                               ,configFile=cfname              &     !<== import the config file with the name
                               ,rc=rc)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
!
!-----------------------------------------------------------------------
!***    Get Attributes from the configure file to create the clock for the dst_gridcomp
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        MESSAGE_CHECK="  DEMO_init:  Extract timestep from configure file for i-th dst gridded component"
!        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=cf_dst(i)          &
                                    ,value = dt_int            &
                                    ,label ='dt_int:'          &
                                    ,rc    =rc)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        MESSAGE_CHECK="  DEMO_init:  Extract Number of compute tasks in the I direction for the dst gridded compo$
!        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=cf_dst(i)          &
                                    ,value = dt_num            &
                                    ,label ='dt_num:'          &
                                    ,rc    =rc)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        MESSAGE_CHECK="  DEMO_init:  Extract Number of compute tasks in the I direction for the dst gridded compo$
!        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=cf_dst(i)          &
                                    ,value = dt_den            &
                                    ,label ='dt_den:'          &
                                    ,rc    =rc)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Set time interval in Destination Clock"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=timestep_dst   &  !<-- Destination Clock's timestep
                               ,s           =dt_int     &  !<-- Whole part of timestep
                               ,sn          =dt_num     &  !<-- Numerator of fractional part
                               ,sd          =dt_den     &  !<-- Denominator of fractional part
                               ,rc          =RC)
!
      ESMF_ERR_ABORT(RC)
!
!-----------------------------------------------------------------------
!***  Create the DST component's clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Creating clock_dst"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      clock_dst=ESMF_ClockCreate(name       ='CLOCK_DST'           &  !<-- The top-level ESMF Clock
                                ,startTime  =starttime_demo        &  !<-- use same starttime as demo!!
                                ,runDuration=runduration_demo      &  !<-- use same runduration as demo!!
                                ,timeStep   =timestep_dst          &  !<--  use newly defined timestep!!
                                ,rc         =RC)
!
      ESMF_ERR_ABORT(RC)
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the DST component.
!***  The Initialize routine that is called here as well as the
!***  Run and Finalize routines (invoked below) are those which are specified
!***  in the Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Execute (array of) 'dst_initialize'."
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
       CALL ESMF_GridCompInitialize(gridcomp    = dst_gridcomp(i)    &
                                   ,importState = dst_imp_state(i)   &
                                   ,exportState = dst_exp_state(i)   &
                                   ,clock       = clock_dst          &
                                   ,phase       = 1                  &
                                   ,userRc      = RC_USER            &
                                   ,rc          = RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
        ESMF_ERR_RETURN(RC_USER,RC_INIT)
      end do
!
      CALL ESMF_VMBarrier(vm=vm_global)
!
!-----------------------------------------------------------------------
!***    Reconcile 'dst_exp_state(i)' across all pets. [It was created only on dst_petlist(i).]
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "  DEMO_init:  Reconciling 'dst_exp_state(i)' across global VM"
        CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        do i=1,nstorms
          CALL ESMF_StateReconcile(state=dst_exp_state(i)     &
                                  ,rc=rc)
!
          ESMF_ERR_RETURN(rc,rc_init)
        end do
!
        CALL ESMF_VMBarrier(vm=VM_global)
!
!-----------------------------------------------------------------------
!***  Create the SRC2DST Coupler Component
!***
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Create SRC2DST Coupler Component for each storm"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      !Use all PETs to create the SRC2dst component.
      do i=1,nstorms
      src2dst_cplcomp(i) = ESMF_CplCompCreate(name=storm_name(i)      &  !<-- Name of the SRC component
                                             !,petlist=src2dst_petlist(i,:)  &  !<-- i-th element's PETlist
                                             ,configFile=configfile(i)   &  !<-- Associate the name of the config $
                                             ,rc=RC)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines
!***  of the SRC2DST coupler component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Register SRC2DST Init, Run, Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        call ESMF_CPLCompSetServices(cplcomp    =src2dst_cplcomp(i)       &  !<-- The SRC gridded component
                                    ,userRoutine=src2dst_register      &  !<-- User's name for the Register routine
                                    ,rc=RC)
!
        ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Create the SRC2DST component's import/export states.
!***  Currently they are not required to perform an actual function.  <-- They're not?? 18Jul2017
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Create the SRC2DST import & export states"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        write(imp_state_name, '("storm0",I2.2,"_cpl_imp_state")') i
        src2dst_imp_state(i) = ESMF_StateCreate(name=imp_state_name                  &
                                               ,stateintent=ESMF_STATEINTENT_IMPORT  &
                                               ,rc=RC)
!
        write(exp_state_name, '("storm0",I2.2,"_cpl_exp_state")') i
        src2dst_exp_state(i) = ESMF_StateCreate(name=exp_state_name                  &
                                               ,stateintent=ESMF_STATEINTENT_EXPORT  &
                                               ,rc=RC)
!
      ESMF_ERR_RETURN(RC,RC_INIT)
      end do
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the SRC2DST component.
!***  The Initialize routine that is called here as well as the
!***  Run and Finalize routines invoked below are those specified
!***  in the Register routine called in ESMF_CplCompSetServices above.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create the SRC2DST component's clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Creating clock_src2dst"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      clock_src2dst=ESMF_ClockCreate(name    = 'CLOCK_SRC2DST'       &  !<-- The coupler Clock
                                ,startTime   = starttime_demo        &  !<-- use same starttime as dst!!
                                ,runDuration = runduration_demo      &  !<-- use same runduration as dst!!
                                ,timeStep    = timestep_dst          &  !<-- use timestep as dst!!
                                ,rc          = RC)
!
      ESMF_ERR_ABORT(RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="  NEMS_init:  Execute (array of) 'src2dst_initialize'"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      do i=1,nstorms
        call ESMF_CplCompInitialize(cplcomp     = src2dst_cplcomp(i)    &
                                   ,importState = src2dst_imp_state(i)  &
                                   ,exportState = src2dst_exp_state(i)  &
                                   ,clock       = clock_src2dst         &
                                   ,phase       = 1                     &
                                   ,userRc      = RC_USER               &
                                   ,rc          = RC)
        ESMF_ERR_RETURN(RC,RC_INIT)
        ESMF_ERR_RETURN(RC_USER,RC_INIT)
      end do
!
      CALL ESMF_VMBarrier(VM_global)
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!=-=-=-=-=-=-=-=-=-=-=-= END OF STEVE'S DOING  =-=-=-=-=-=-=-=-=-=-=-=-=
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_RUN(NEMS_GRID_COMP                                &
                         ,IMP_STATE                                     &
                         ,EXP_STATE                                     &
                         ,CLOCK_MAIN                                    &
                         ,RC_RUN)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NEMS_GRID_COMP                                !<-- The NEMS component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NEMS import state
                         ,EXP_STATE                                        !<-- The NEMS export state
!
      TYPE(ESMF_Clock) :: CLOCK_MAIN                                       !<-- The main Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN                                        !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: HH,I,J,RC,RC_USER
!
      TYPE(ESMF_Time) :: CURRTIME
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
      type(ESMF_TimeInterval) :: timeStep     !swd
!
!=-=-=-=-= Steve's doing -=-=-=-=-=
      !* local integers for printing clock information
      integer(ESMF_KIND_I4) :: yy, h, m, s         !<-- why these are ESMF_KIND...
      integer               :: mm, dd              !<-- and these are not, is beyond me....
!=-=-= END Steve's doing -=-=-=-=-=
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_RUN = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the main Clock within
!***  the NEMS component.
!-----------------------------------------------------------------------
!
      CLOCK_NEMS=CLOCK_MAIN
!
!-----------------------------------------------------------------------
!***  Execute the Run step of each element in the EARTH component
!***  array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Run step of the EARTH components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompRun(gridcomp    = EARTH_GRID_COMP(I)        &
                               ,importState = EARTH_IMP_STATE(I)        &
                               ,exportState = EARTH_EXP_STATE(I)        &
                               ,clock       = CLOCK_NEMS                &
                               ,phase       = 1                         &
                               ,userRc      = RC_USER                   &
                               ,rc          = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
          ESMF_ERR_RETURN(RC_USER,RC_RUN)
!
!=-=-=-=  Steve's doing  =-=-=-=-=-=-=-=-=-=--=
          CALL ESMF_GridCompRun(gridcomp    = EARTH_GRID_COMP(I)        &
                               ,importState = EARTH_IMP_STATE(I)        &
                               ,exportState = EARTH_EXP_STATE(I)        &
                               ,clock       = CLOCK_NEMS                &
                               ,phase       = 2                         &
                               ,userRc      = RC_USER                   &
                               ,rc          = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
          ESMF_ERR_RETURN(RC_USER,RC_RUN)
!=-=-= END Steve's doing  =-=-=-=-=-=-=-=-=-=--=
!
        END IF
!
      END DO
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      sps: IF(ENS_SPS) THEN
!
!-----------------------------------------------------------------------
!
        DO I = NUMBER_START, NUMBER_FINAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run step of the Ensemble Coupler component"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =ENS_CPL_COMP                 &
                              ,importState=ENS_CPL_IMP_STATE            &
                              ,exportState=ENS_CPL_EXP_STATE            &
                              ,clock      =CLOCK_NEMS                   &
                              ,phase      =1                            &
                              ,userRc     =RC_USER                      &
                              ,rc         =RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
          ESMF_ERR_RETURN(RC_USER,RC_RUN)
!
          CALL ESMF_VMBarrier(vm = VM_GLOBAL, rc = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
!-----------------------------------------------------------------------
!***  Adjust the ESMF clock for the next run cycle.
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Update the current time of the NEMS clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock       = CLOCK_NEMS                   &
                            ,runDuration = RUNDURATION                  &
                            ,rc          = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
          CURRTIME = STARTTIME + RUNDURATION
!
          CALL ESMF_ClockSet(clock    = CLOCK_NEMS                      &
                            ,currTime = CURRTIME                        &
                            ,rc       = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Adjust clock - add one more cycle run duration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeInterval = RUNDURATION          &
                                   ,h            = HH                   &
                                   ,rc           = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
          HH = HH + HH_INCREASE
!
          CALL ESMF_TimeIntervalSet(timeInterval = RUNDURATION          &
                                   ,h            = hh                   &
                                   ,rc           = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
          CALL ESMF_ClockSet(clock       = CLOCK_NEMS                   &
                            ,runDuration = RUNDURATION                  &
                            ,rc = RC)
          ESMF_ERR_RETURN(RC,RC_RUN)
!
!-----------------------------------------------------------------------
!***  Execute the Run step of each element in the EARTH component
!***  array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run step of the EARTH component"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO J = 1, TOTAL_MEMBER
!
            IF(MEMBER_ID == J) THEN
                CALL ESMF_GridCompRun(gridcomp    = EARTH_GRID_COMP(J)  &
                                     ,importState = EARTH_IMP_STATE(J)  &
                                     ,exportState = EARTH_EXP_STATE(J)  &
                                     ,clock       = CLOCK_NEMS          &
                                     ,phase       = 1                   &
                                     ,userRc      = RC_USER             &
                                     ,rc          = RC)
                ESMF_ERR_RETURN(RC,RC_RUN)
                ESMF_ERR_RETURN(RC_USER,RC_RUN)
            END IF
!
          END DO
!
          PRINT*, 'Complete EARTH Run Cycle ', I + 1
        END DO
!
!-----------------------------------------------------------------------
!
      ELSE sps
!
         CALL ESMF_ClockGet(clock       = CLOCK_NEMS                    &
                           ,runDuration = RUNDURATION                   &
                           ,rc          = RC)
         ESMF_ERR_RETURN(RC,RC_RUN)
!
!-----------------------------------------------------------------------
!
      END IF sps
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_FINALIZE(NEMS_GRID_COMP                           &
                              ,IMP_STATE                                &
                              ,EXP_STATE                                &
                              ,CLOCK_MAIN                               &
                              ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NEMS_GRID_COMP                                !<-- The NEMS component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NEMS import state
                         ,EXP_STATE                                        !<-- The NEMS export state
!
      TYPE(ESMF_Clock) :: CLOCK_MAIN                                       !<-- The main Clock
!
      INTEGER,INTENT(OUT) :: RC_FINALIZE                                   !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,RC,RC_USER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_FINALIZE = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Execute the Finalize step of each element of the 
!***  EARTH component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Finalize step of the EARTH component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompFinalize(gridcomp    = EARTH_GRID_COMP(I)   &
                                    ,importState = EARTH_IMP_STATE(I)   &
                                    ,exportState = EARTH_EXP_STATE(I)   &
                                    ,clock       = CLOCK_NEMS           &
                                    ,phase       = 1                    &
                                    ,userRc      = RC_USER              &
                                    ,rc          = RC)
          ESMF_ERR_RETURN(RC,RC_FINALIZE)
          ESMF_ERR_RETURN(RC_USER,RC_FINALIZE)
        END IF
!
      END DO
!
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Execute the Finalize step of the Ensemble Coupler component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompFinalize(cplcomp    =ENS_CPL_COMP              &
                                 ,importState=ENS_CPL_IMP_STATE         &
                                 ,exportState=ENS_CPL_EXP_STATE         &
                                 ,clock      =CLOCK_NEMS                &
                                 ,phase      =1                         &
                                 ,userRc     =RC_USER                   &
                                 ,rc         =RC)
        ESMF_ERR_RETURN(RC,RC_FINALIZE)
        ESMF_ERR_RETURN(RC_USER,RC_FINALIZE)
!
      END IF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NEMS_GRID_COMP
!
!-----------------------------------------------------------------------
