      MODULE module_GFS_GRID_COMP

      USE ESMF_MOD

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: GFS_REGISTER

      INTEGER :: DUMMY

      CONTAINS

!#######################################################################

      SUBROUTINE GFS_REGISTER(GFS_GRID_COMP,RC_REG)
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GFS_GRID_COMP
      INTEGER            ,INTENT(OUT)   :: RC_REG

      INTEGER :: RC

      write(0,*) "    GFS_REGISTER stub"
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP ,ESMF_SETINIT ,GFS_INITIALIZE ,ESMF_SINGLEPHASE ,RC)
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP ,ESMF_SETRUN  ,GFS_RUN        ,1                ,RC)
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP ,ESMF_SETFINAL,GFS_FINALIZE   ,ESMF_SINGLEPHASE ,RC)

      RC_REG = ESMF_SUCCESS
      write(0,*) "    END OF GFS_REGISTER stub"

      END SUBROUTINE GFS_REGISTER

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GFS_INITIALIZE(GFS_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_GFS ,RC_INIT)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: GFS_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_GFS
      INTEGER            ,INTENT(OUT)   :: RC_INIT

      write(0,*) "        GFS_INITIALIZE stub"
      RC_INIT = ESMF_SUCCESS
      write(0,*) "        END OF GFS_INITIALIZE stub"

      END SUBROUTINE GFS_INITIALIZE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GFS_RUN(GFS_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_GFS ,RC_RUN)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: GFS_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_GFS
      INTEGER            ,INTENT(OUT)   :: RC_RUN

      write(0,*) "        GFS_RUN stub"
      RC_RUN=ESMF_SUCCESS
      write(0,*) "        END OF GFS_RUN stub"

      END SUBROUTINE GFS_RUN

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GFS_FINALIZE(GFS_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_GFS ,RC_FINALIZE)

      TYPE(ESMF_GridComp),INTENT(INOUT) :: GFS_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_GFS
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE

      write(0,*) "        GFS_FINALIZE stub"
      RC_FINALIZE=ESMF_SUCCESS
      write(0,*) "        END OF GFS_FINALIZE stub"

      END SUBROUTINE GFS_FINALIZE

!#######################################################################

      END MODULE module_GFS_GRID_COMP
