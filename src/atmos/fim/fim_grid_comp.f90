      MODULE module_FIM_GRID_COMP

      USE ESMF_MOD
      USE FIM_INTERNAL_STATE_MOD ,ONLY: FIM_INTERNAL_STATE            &
                                       ,WRAP_FIM_INTERNAL_STATE

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: FIM_REGISTER

      TYPE(FIM_INTERNAL_STATE),POINTER,SAVE :: FIM_INT_STATE
      TYPE(WRAP_FIM_INTERNAL_STATE)   ,SAVE :: WRAP

      CONTAINS

!#######################################################################

      SUBROUTINE FIM_REGISTER(FIM_GRID_COMP,RC_REG)
      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      INTEGER            ,INTENT(OUT)   :: RC_REG

      INTEGER :: RC

      write(0,*) "    FIM_REGISTER"
      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETINIT ,FIM_INITIALIZE ,ESMF_SINGLEPHASE ,RC)
      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETRUN  ,FIM_RUN        ,1                ,RC)
      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETFINAL,FIM_FINALIZE   ,ESMF_SINGLEPHASE ,RC)

      RC_REG = ESMF_SUCCESS
      write(0,*) "    END OF FIM_REGISTER"

      END SUBROUTINE FIM_REGISTER

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_INITIALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_INIT)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_INIT

      INTEGER :: RC
      TYPE(ESMF_Config) :: CF

      write(0,*) "        FIM_INITIALIZE"
      RC_INIT = ESMF_SUCCESS

      ALLOCATE(FIM_INT_STATE,stat=RC)
      WRAP%FIM_INT_STATE=>FIM_INT_STATE

      CALL ESMF_GridCompSetInternalState(FIM_GRID_COMP ,WRAP ,RC)

      CF=ESMF_ConfigCreate(rc=RC)
      CALL ESMF_ConfigLoadFile(config=CF ,filename='fim.configure' ,rc=RC)
      CALL ESMF_GridCompSet(FIM_GRID_COMP, config=CF, rc=RC)

      write(0,*) "        END OF FIM_INITIALIZE"
      END SUBROUTINE FIM_INITIALIZE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_RUN(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_RUN)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_RUN

      write(0,*) "        FIM_RUN"
      RC_RUN=ESMF_SUCCESS

      write(0,*) "        END OF FIM_RUN"
      END SUBROUTINE FIM_RUN

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_FINALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_FINALIZE)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE

      write(0,*) "        FIM_FINALIZE"
      RC_FINALIZE=ESMF_SUCCESS

      write(0,*) "        END OF FIM_FINALIZE"
      END SUBROUTINE FIM_FINALIZE

!#######################################################################

      END MODULE module_FIM_GRID_COMP
