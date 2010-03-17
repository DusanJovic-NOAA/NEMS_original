!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PHYSICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE MAIN GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS PHYSICS INITIALIZE, ETC.) 
!***  IN MODULE_ATM_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-28  Vasic - Removed counters (computed in SET_INTERNAL_STATE_PHY)
!   2008-10-    Vasic - Restart capability
!   2008-08-23  Janjic - Removed uz0h, vz0h
!   2008-08-23  Janjic - General hybrid coordinate
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: PHY_REGISTER
!
!-----------------------------------------------------------------------
!!!!  INCLUDE 'kind.inc'
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Physics Gridded Component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Register
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_INITIALIZE(GRID_COMP,IMP_STATE,EXP_STATE,CLOCK     &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_RUN(GRID_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_RUN)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_FINALIZE(GRID_COMP                                 &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK_ATM                                 &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE PHYSICS COMPONENT.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Physics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE                 !<-- The Physics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE_WRITE                 !<-- The Physics export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM component's ESMF Clock.
!
      INTEGER            ,INTENT(OUT)   :: RCFINAL
!      
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_FINALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
