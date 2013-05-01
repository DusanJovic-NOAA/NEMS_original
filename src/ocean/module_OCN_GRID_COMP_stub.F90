#include "../ESMFVersionDefine.h"
#include "../NUOPC_Switch.h"

!-----------------------------------------------------------------------
!
      MODULE module_OCN_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the OCN component.
!-----------------------------------------------------------------------
!
!***  The OCN component lies in the heirarchy seen here:
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
!          OCN/OCN/ICE components
!               |
!               |
!          CORE component (GFS, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!  2011-05-11  Theurich & Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!
      USE ESMF_MOD
!
#ifdef WITH_NUOPC
      use NUOPC
      use NUOPC_Model, only: &
        model_routine_SS    => routine_SetServices, &
        model_label_Advance => label_Advance
#endif
!
      USE module_OCN_INTERNAL_STATE,ONLY: OCN_INTERNAL_STATE            &
                                         ,WRAP_OCN_INTERNAL_STATE

!-----------------------------------------------------------------------
! The different ocean module core registers go here
!-----------------------------------------------------------------------
      USE module_HYCOM_GRID_COMP,ONLY: HYCOM_REGISTER

      USE module_GENOCN_GRID_COMP,ONLY: GENOCN_REGISTER   ! For the "Generic Ocen Core" gridded component.

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
      PUBLIC :: OCN_REGISTER
!
!-----------------------------------------------------------------------
!
      TYPE(OCN_INTERNAL_STATE),POINTER,SAVE :: OCN_INT_STATE
      TYPE(WRAP_OCN_INTERNAL_STATE)   ,SAVE :: WRAP
!
      TYPE(ESMF_Clock),SAVE :: CLOCK_OCN                                   !<-- The Clock of the OCN component
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE OCN_REGISTER(OCN_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  Register the Init, Run, and Finalize routines of 
!***  the OCN component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: OCN_GRID_COMP                   !<-- The OCN component
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
      call model_routine_SS(OCN_GRID_COMP, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! NUOPC_Model requires InitP0, even if it is NOOP
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        NOOP, phase=0, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! attach specializing method(s)
      call ESMF_MethodAdd(OCN_GRID_COMP, label=model_label_Advance, &
        userRoutine=OCN_ADVANCE, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif


! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for OCN Initialize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_INITIALIZE                      &  !<-- Subroutine type (Initialize)
                                     ,OCN_INITIALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_INITIALIZE                      &  !<-- Subroutine type (Initialize)
                                     ,OCN_INITIALIZE                    &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

#ifndef WITH_NUOPC
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for OCN Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_RUN                       &  !<-- Subroutine type (Run)
                                     ,OCN_RUN                           &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_RUN                       &  !<-- Subroutine type (Run)
                                     ,OCN_RUN                           &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
#endif   /* WITH_NUOPC  */

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for OCN Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_FINALIZE                     &  !<-- Subroutine type (Finalize)
                                     ,OCN_FINALIZE                      &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(OCN_GRID_COMP                     &  !<-- The OCN component
                                     ,ESMF_METHOD_FINALIZE                     &  !<-- Subroutine type (Finalize)
                                     ,OCN_FINALIZE                      &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
        WRITE(0,*)' OCN_REGISTER succeeded'
      ELSE
        WRITE(0,*)' OCN_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE OCN_REGISTER
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
#endif
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE OCN_INITIALIZE(OCN_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_EARTH                             &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  The Initialize step of the OCN component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: OCN_GRID_COMP                   !<-- The OCN component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The OCN import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The OCN export state
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
      RC_INIT = ESMF_SUCCESS

!-----------------------------------------------------------------------

      END SUBROUTINE OCN_INITIALIZE
!-----------------------------------------------------------------------


#ifndef WITH_NUOPC
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE OCN_RUN(OCN_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_EARTH                                    &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The Run step of the OCN component.
!-----------------------------------------------------------------------

!------------------------
!***  Argument Variables
!------------------------
! 
      TYPE(ESMF_GridComp)               :: OCN_GRID_COMP                   !<-- The OCN component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The OCN import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The OCN export state
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
      RC_RUN = 0


      END SUBROUTINE OCN_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

#else

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

      SUBROUTINE OCN_ADVANCE(OCN_GRID_COMP,RC_RUN)

!-----------------------------------------------------------------------
!***  The NUOPC Run step of the OCN component.
!-----------------------------------------------------------------------

!------------------------
!***  Argument Variables
!------------------------

      TYPE(ESMF_GridComp)               :: OCN_GRID_COMP                   !<-- The OCN component
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code

!---------------------
!***  Local Variables
!---------------------

      INTEGER :: RC

      TYPE(ESMF_Time) :: CURRTIME                                       &
                        ,STARTTIME

      TYPE(ESMF_TimeInterval) :: RUNDURATION
      type(ESMF_Pointer)      :: this
!
!-----------------------------------------------------------------------
      RC_RUN = 0

!-----------------------------------------------------------------------
      END SUBROUTINE OCN_ADVANCE
!-----------------------------------------------------------------------
#endif



!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

!
      SUBROUTINE OCN_FINALIZE(OCN_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_EARTH                               &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  Finalize the OCN component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: OCN_GRID_COMP                   !<-- The OCN component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The OCN import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The OCN import state
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
      RC_FINALIZE = 0

!-----------------------------------------------------------------------
      END SUBROUTINE OCN_FINALIZE
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      END MODULE module_OCN_GRID_COMP
!-----------------------------------------------------------------------

