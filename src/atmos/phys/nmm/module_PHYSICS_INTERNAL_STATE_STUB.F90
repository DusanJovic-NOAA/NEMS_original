!-----------------------------------------------------------------------
      MODULE MODULE_PHYSICS_INTERNAL_STATE
!-----------------------------------------------------------------------
!***  This module declares the derived datatype called INTERNAL_STATE.
!***  For now the components of this datatype will be everything needed
!***  to advance the model integration, i.e. everything that would be
!***  part of a restart file.  Specifically this will include those
!***  quantities that evolve during the integration, the namelist
!***  variables, and the grid decomposition variables.
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!                   
!   2008-07-28  Vasic - Precompute accumulation counters.
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
!      USE MODULE_INCLUDE
!      USE MODULE_ERR_MSG     ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      TYPE INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END TYPE INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  THE INTERNAL_STATE TYPE IS SUPPORTED BY A C POINTER (NOT AN F90
!***  POINTER) AND THEREFORE THE FOLLOWING TYPE IS NEEDED.
!-----------------------------------------------------------------------
!
      TYPE WRAP_INTERNAL_STATE
        TYPE(INTERNAL_STATE),POINTER :: INT_STATE
      END TYPE WRAP_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
