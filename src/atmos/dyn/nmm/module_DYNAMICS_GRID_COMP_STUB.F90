!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.) 
!***  IN MODULE_MAIN_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-30  Janjic - Add CONVECTION='none' to OPERATIONAL_PHYSICS.
!               Janjic - Fix lower J limit in FFTFHN(WATER).
!   2008-08-23  Janjic - General pressure-sigma hybrid
!               Janjic - Consistent nonhydrostatic correction in the
!                        first term of the pressure gradient force
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
      PUBLIC :: DYN_REGISTER                                            
!
!
!
!-----------------------------------------------------------------------
      INCLUDE '../../../../inc/kind.inc'
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
      SUBROUTINE DYN_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Dynamics gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Dyn register
!
!
      END SUBROUTINE DYN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_INITIALIZE(GRID_COMP                               &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_ATM                               &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CARRY OUT ALL NECESSARY SETUPS FOR THE MODEL DYNAMICS.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
      INCLUDE 'mpif.h'
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_ATM                       !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL,    INTENT(OUT)  :: RC_INIT
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_RUN(GRID_COMP                                      &
                        ,IMP_STATE,EXP_STATE                            &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                     !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                     !<-- The Dynamics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                     !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM                     !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_FINALIZE(GRID_COMP                                 &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK_ATM                                 &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE DYNAMICS COMPONENT.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE                 !<-- The Dynamics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE_WRITE                 !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM component's ESMF Clock.
!
      INTEGER            ,INTENT(OUT)   :: RCFINAL
!      
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_FINALIZE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
