#include "../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#else
#define ESMF_520r
#endif

!-----------------------------------------------------------------------
!
      MODULE module_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the ATM component.
!-----------------------------------------------------------------------
!
!***  The ATM component lies in the heirarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!               |
!               |
!          ATM/OCEAN/ICE components
!               |
!               |
!          CORE component (GFS, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!  2011-05-11  Theurich & Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10/04  Yang  - Modified for using the ESMF 5.2.0r library.
!-----------------------------------------------------------------------
!
      USE esmf_mod

#ifdef WITH_NUOPC
      use NUOPC
      use NUOPC_Model, only: &
        model_routine_SS    => routine_SetServices, &
        model_label_Advance => label_Advance
#endif

!
      USE module_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
!
      USE module_NMM_GRID_COMP,ONLY: NMM_REGISTER
      USE module_GFS_GRID_COMP,ONLY: GFS_REGISTER
      USE module_FIM_GRID_COMP,ONLY: FIM_REGISTER
      USE module_GEN_GRID_COMP,ONLY: GEN_REGISTER   ! For the "Generic Core" gridded component.
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
      PUBLIC :: ATM_REGISTER
!
!-----------------------------------------------------------------------
!
      TYPE(ATM_INTERNAL_STATE),POINTER,SAVE :: ATM_INT_STATE
      TYPE(WRAP_ATM_INTERNAL_STATE)   ,SAVE :: WRAP
!
      TYPE(ESMF_Clock),SAVE :: CLOCK_ATM                                   !<-- The Clock of the ATM component
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
!***  Register the Init, Run, and Finalize routines of 
!***  the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      INTEGER            ,INTENT(OUT)   :: RC_REG                          !<-- Error return code
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

#ifdef WITH_NUOPC
      ! the NUOPC model component will register the generic methods
      call model_routine_SS(ATM_GRID_COMP, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! NUOPC_Model requires InitP2, even if it is NOOP
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        NOOP, phase=2, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! NUOPC_Model requires InitP4, even if it is NOOP
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        NOOP, phase=4, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! attach specializing method(s)
      call ESMF_MethodAdd(ATM_GRID_COMP, label=model_label_Advance, &
        userRoutine=ATM_ADVANCE, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for ATM Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- The ATM component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Subroutine type (Initialize)
                                     ,ATM_INITIALIZE                    &  !<-- User's subroutine name
#ifdef ESMF_3
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
#ifdef ESMF_520r
                                     ,phase=1                           &
                                     ,rc=RC)
#else
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for ATM Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifndef WITH_NUOPC

      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- The ATM component
                                     ,ESMF_METHOD_RUN                   &  !<-- Subroutine type (Initialize)
                                     ,ATM_RUN                           &  !<-- User's subroutine name
#ifdef ESMF_3
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
#ifdef ESMF_520r
                                     ,phase=1                           &
                                     ,rc=RC)
#else
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
#endif
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for ATM Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- The ATM component
                                     ,ESMF_METHOD_FINALIZE              &  !<-- Subroutine type (Initialize)
                                     ,ATM_FINALIZE                      &  !<-- User's subroutine name
#ifdef ESMF_3
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
#ifdef ESMF_520r
                                     ,phase=1                           &
                                     ,rc=RC)
#else
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_REGISTER succeeded'
      ELSE
        WRITE(0,*)' ATM_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_REGISTER
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
#ifdef WITH_NUOPC
  subroutine Noop(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
#endif

      SUBROUTINE ATM_INITIALIZE(ATM_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_EARTH                             &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  The Initialize step of the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM export state
      TYPE(ESMF_Clock)                  :: CLOCK_EARTH                     !<-- The Clock of the EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
      TYPE(ESMF_Config) :: CF
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_INIT = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
#ifdef WITH_NUOPC
      call NUOPC_ClockPrintCurrTime(CLOCK_EARTH, &
        string="entering ATM_INITIALIZE with CLOCK_EARTH current: ")
      call NUOPC_ClockPrintStartTime(CLOCK_EARTH, &
        string="entering ATM_INITIALIZE with CLOCK_EARTH start:   ")
      call NUOPC_ClockPrintStopTime(CLOCK_EARTH, &
        string="entering ATM_INITIALIZE with CLOCK_EARTH stop:    ")
#endif

!
!-----------------------------------------------------------------------
!***  Allocate the ATM component's internal state, point at it,
!***  and attach it to the ATM component.
!-----------------------------------------------------------------------
!
      ALLOCATE(ATM_INT_STATE,stat=RC)
      wrap%ATM_INT_STATE=>ATM_INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the ATM Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(ATM_GRID_COMP                  &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the EARTH Clock within
!***  the ATM component.
!-----------------------------------------------------------------------
!
#ifdef WITH_NUOPC
      call NUOPC_GridCompSetClock(ATM_GRID_COMP, CLOCK_EARTH, rc=RC_INIT)
      if (ESMF_LogFoundError(rcToCheck=RC_INIT, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      atm_int_state%CLOCK_ATM = ESMF_ClockCreate(CLOCK_EARTH, rc=RC_INIT)
      if (ESMF_LogFoundError(rcToCheck=RC_INIT, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

#else

      atm_int_state%CLOCK_ATM=CLOCK_EARTH

#endif

!
!-----------------------------------------------------------------------
!***  Create the configure object for the ATM configure file which
!***  specifies the dynamic core.
!-----------------------------------------------------------------------
!
      CF=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Load the ATM configure file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigLoadFile(config=CF ,filename='atmos.configure' ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the configure object to the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the configure file to the ATM component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM component
                           ,config  =CF                                 &  !<-- The associated configure object
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the dynamic core name from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract dynamic core from the ATM configure file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ATM configure object
                                  ,value =atm_int_state%CORE            &  !<-- The dynamic core name
                                  ,label ='core:'                       &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ATM subcomponent and its associated import/export
!***  states for the core name that was extracted.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_GRID_COMP=ESMF_GridCompCreate(name=TRIM(atm_int_state%CORE)//' component' &
                                                      ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the configure object to the CORE component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     MESSAGE_CHECK="Attach the configure file to the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     CALL ESMF_GridCompSet(gridcomp=atm_int_state%CORE_GRID_COMP       &  !<-- The ATM component
!                          ,config  =CF                                 &  !<-- The associated configure object
!                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the subcomponent's Init, Run, and Finalize subroutines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the CORE component's Init, Run, and Finalize steps"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      SELECT CASE(atm_int_state%CORE)
!
#ifdef ESMF_3
        CASE('nmm')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,NMM_REGISTER                   &
                                        ,RC)
!
        CASE('gfs')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GFS_REGISTER                   &
                                        ,RC)
!
        CASE('fim')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,FIM_REGISTER                   &
                                        ,RC)

        CASE('gen')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GEN_REGISTER                   &
                                        ,RC)

#else
        CASE('nmm')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,NMM_REGISTER                   &
                                        ,rc=RC)
!
        CASE('gfs')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GFS_REGISTER                   &
                                        ,rc=RC)
!
        CASE('fim')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,FIM_REGISTER                   &
                                        ,rc=RC)

        CASE('gen')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GEN_REGISTER                   &
                                        ,rc=RC)
#endif
        CASE DEFAULT
          write(0,*)' ATM_INITIALIZE requires unknown core: ',TRIM(atm_int_state%CORE)                      
!
      END SELECT
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Core component's import/export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_IMP_STATE=ESMF_StateCreate(STATENAME   = "CORE Import"            &
                                                   ,stateintent = ESMF_STATEINTENT_IMPORT  &
                                                   ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE export state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_EXP_STATE=ESMF_StateCreate(STATENAME   = "CORE Export"            &
                                                   ,stateintent = ESMF_STATEINTENT_EXPORT  &
                                                   ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Nest the import/export states of the CORE component into the
!***  analgous states of the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK= "Add the CORE states into the ATMOS states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifndef WITH_NUOPC
      CALL ESMF_StateAdd(IMP_STATE,LISTWRAPPER(atm_int_state%CORE_IMP_STATE),rc = RC)
      CALL ESMF_StateAdd(EXP_STATE,LISTWRAPPER(atm_int_state%CORE_EXP_STATE),rc = RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Initialize the CORE component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                  ,importState=atm_int_state%CORE_IMP_STATE &
                                  ,exportState=atm_int_state%CORE_EXP_STATE &
                                  ,clock      =atm_int_state%CLOCK_ATM      &
                                  ,rc         =RC)
#else
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                  ,importState=atm_int_state%CORE_IMP_STATE &
                                  ,exportState=atm_int_state%CORE_EXP_STATE &
                                  ,clock      =atm_int_state%CLOCK_ATM      &
                                  ,phase      =ESMF_SINGLEPHASE             &
                                  ,rc         =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' ATM_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
#ifdef WITH_NUOPC
      call NUOPC_ClockPrintCurrTime(CLOCK_EARTH, &
        string="leaving  ATM_INITIALIZE with CLOCK_EARTH current: ")
      call NUOPC_ClockPrintStartTime(CLOCK_EARTH, &
        string="leaving  ATM_INITIALIZE with CLOCK_EARTH start:   ")
      call NUOPC_ClockPrintStopTime(CLOCK_EARTH, &
        string="leaving  ATM_INITIALIZE with CLOCK_EARTH stop:    ")

      call NUOPC_ClockPrintCurrTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_INITIALIZE with CLOCK_ATM current: ")
      call NUOPC_ClockPrintStartTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_INITIALIZE with CLOCK_ATM start:   ")
      call NUOPC_ClockPrintStopTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_INITIALIZE with CLOCK_ATM stop:    ")
#endif
!
!-----------------------------------------------------------------------
!

      END SUBROUTINE ATM_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
#ifndef WITH_NUOPC

      SUBROUTINE ATM_RUN(ATM_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_EARTH                                    &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The Run step of the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
! 
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM export state
      TYPE(ESMF_Clock)                  :: CLOCK_EARTH                     !<-- The Clock of the EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
      TYPE(ESMF_Time) :: CURRTIME                                       &
                        ,STARTTIME
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the EARTH Clock within
!***  the ATM component.
!-----------------------------------------------------------------------
!
      atm_int_state%CLOCK_ATM=CLOCK_EARTH
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the selected dynamic core.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Run step of the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
      CALL ESMF_GridCompRun(gridcomp   =atm_int_state%CORE_GRID_COMP    &
                           ,importState=atm_int_state%CORE_IMP_STATE    &
                           ,exportState=atm_int_state%CORE_EXP_STATE    &
                           ,clock      =atm_int_state%CLOCK_ATM         &
                           ,rc         =RC)
#else
      CALL ESMF_GridCompRun(gridcomp   =atm_int_state%CORE_GRID_COMP    &
                           ,importState=atm_int_state%CORE_IMP_STATE    &
                           ,exportState=atm_int_state%CORE_EXP_STATE    &
                           ,clock      =atm_int_state%CLOCK_ATM         &
                           ,phase      =ESMF_SINGLEPHASE                & 
                           ,rc         =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Update the ATMOS clock.
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Update the current time of the ATMOS clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =atm_int_state%CLOCK_ATM            &
                        ,startTime  =STARTTIME                          &
                        ,runDuration=RUNDURATION                        &
                        ,rc         =RC)
!
      CURRTIME=STARTTIME+RUNDURATION
!
      CALL ESMF_ClockSet(clock   =atm_int_state%CLOCK_ATM               &
                        ,currTime=CURRTIME                              &
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_RUN succeeded'
      ELSE
        WRITE(0,*)' ATM_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

#else

      SUBROUTINE ATM_ADVANCE(ATM_GRID_COMP                                  &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The Run step of the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
      TYPE(ESMF_Time) :: CURRTIME                                       &
                        ,STARTTIME
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
      type(ESMF_Pointer)      :: this
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Use the internal Clock set by NUOPC layer for ATM
!-----------------------------------------------------------------------
!
      call NUOPC_ClockPrintCurrTime(atm_int_state%CLOCK_ATM, &
        string="entering ATM_ADVANCE with CLOCK_ATM current: ")
      call NUOPC_ClockPrintStartTime(atm_int_state%CLOCK_ATM, &
        string="entering ATM_ADVANCE with CLOCK_ATM start:   ")
      call NUOPC_ClockPrintStopTime(atm_int_state%CLOCK_ATM, &
        string="entering ATM_ADVANCE with CLOCK_ATM stop:    ")

!-----------------------------------------------------------------------
!***  Execute the Run step of the selected dynamic core.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Run step of the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompRun(gridcomp   =atm_int_state%CORE_GRID_COMP    &
                           ,importState=atm_int_state%CORE_IMP_STATE    &
                           ,exportState=atm_int_state%CORE_EXP_STATE    &
                           ,clock      =atm_int_state%CLOCK_ATM         &
                           ,phase      =1                               &
                           ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_RUN succeeded'
      ELSE
        WRITE(0,*)' ATM_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!-----------------------------------------------------------------------
!
      call NUOPC_ClockPrintCurrTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_ADVANCE with CLOCK_ATM current: ")
      call NUOPC_ClockPrintStartTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_ADVANCE with CLOCK_ATM start:   ")
      call NUOPC_ClockPrintStopTime(atm_int_state%CLOCK_ATM, &
        string="leaving  ATM_ADVANCE with CLOCK_ATM stop:    ")

!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_ADVANCE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
#endif

      SUBROUTINE ATM_FINALIZE(ATM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_EARTH                               &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  Finalize the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM import state
      TYPE(ESMF_Clock)                  :: CLOCK_EARTH                     !<-- The Clock of the EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE                     !<-- Error return code
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
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Finalize step of the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                ,importState=atm_int_state%CORE_IMP_STATE &
                                ,exportState=atm_int_state%CORE_EXP_STATE &
                                ,clock      =atm_int_state%CLOCK_ATM      &
                                ,rc         =RC)
#else
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                ,importState=atm_int_state%CORE_IMP_STATE &
                                ,exportState=atm_int_state%CORE_EXP_STATE &
                                ,clock      =atm_int_state%CLOCK_ATM      &
                                ,phase      =ESMF_SINGLEPHASE             &
                                ,rc         =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
#ifdef WITH_NUOPC
!-----------------------------------------------------------------------
!
      call ESMF_ClockDestroy(atm_int_state%CLOCK_ATM, rc=RC_FINALIZE)
      if (ESMF_LogFoundError(rcToCheck=RC_FINALIZE, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' ATM_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
