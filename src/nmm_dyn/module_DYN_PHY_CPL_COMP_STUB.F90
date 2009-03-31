!-----------------------------------------------------------------------
!
      MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE COUPLER'S REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  IN MODULE_ATM_GRID_COMP.F.
!
!***  THE COUPLER PROVIDES 2-WAY COUPLING BETWEEN THE DYNAMICS AND
!***  PHYSICS GRIDDED COMPONENTS BY TRANSFERING THEIR EXPORT AND
!***  IMPORT STATES BETWEEN THE TWO.
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
      PUBLIC :: DYN_PHY_CPL_REGISTER
!
!-----------------------------------------------------------------------
      INCLUDE '../../inc/kind.inc'
!-----------------------------------------------------------------------
!
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_PHY_CPL_REGISTER(CPL_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- Coupler component
!
      INTEGER,INTENT(OUT)              :: RC_REG                          !<-- Return code for register
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_PHY_CPL_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_INITIALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK      &
                               ,RC_CPL)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SET UP THE COUPLER.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,  INTENT(OUT)   :: RC_CPL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_RUN(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_CPL)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  RUN THE COUPLER TO TRANSFER DATA BETWEEN THE GRIDDED COMPONENTS.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
      
!
      INTEGER,OPTIONAL,INTENT(OUT)     :: RC_CPL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_FINALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_CPL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE COUPLER.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                         !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                        !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                        !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                            !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,   INTENT(OUT)  :: RC_CPL
!      
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_FINALIZE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
