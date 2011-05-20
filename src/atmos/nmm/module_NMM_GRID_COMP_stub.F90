#include "../../ESMFVersionDefine.h"

!  2011-05-11  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!--------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
      MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
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
      INTEGER :: DUMMY
!
!-----------------------------------------------------------------------

      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_REGISTER(NMM_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP
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
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETINIT                      &
                                     ,NMM_INITIALIZE                    &
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETRUN                       &
                                     ,NMM_RUN                           &
                                     ,1                                 &
                                     ,RC)
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETFINAL                     &
                                     ,NMM_FINALIZE                      &
                                     ,1                                 &
                                     ,RC)
!
#else
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETINIT                      &
                                     ,NMM_INITIALIZE                    &
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETRUN                       &
                                     ,NMM_RUN                           &
                                     ,phase=1                           &
                                     ,rc=RC)
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_SETFINAL                     &
                                     ,NMM_FINALIZE                      &
                                     ,phase=1                           &
                                     ,rc=RC)
!
#endif

!-----------------------------------------------------------------------
!
      RC_REG = ESMF_SUCCESS
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
                               ,CLOCK_NMM                               &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NMM
      INTEGER            ,INTENT(OUT)   :: RC_INIT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      write(0,*) "        NMM_INITIALIZE stub"
!
      RC_INIT = ESMF_SUCCESS
!
      write(0,*) "        END OF NMM_INITIALIZE stub"
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
                        ,CLOCK_NMM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NMM
      INTEGER            ,INTENT(OUT)   :: RC_RUN
!
!-----------------------------------------------------------------------
!
      write(0,*) "        NMM_RUN stub"
!
      RC_RUN=ESMF_SUCCESS
!
      write(0,*) "        END OF NMM_RUN stub"
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
      TYPE(ESMF_GridComp),INTENT(INOUT) :: NMM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_NMM
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE
!
!-----------------------------------------------------------------------
!
      write(0,*) "        NMM_FINALIZE stub"
!
      RC_FINALIZE=ESMF_SUCCESS
!
      write(0,*) "        END OF NMM_FINALIZE stub"
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
