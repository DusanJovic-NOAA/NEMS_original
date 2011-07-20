!-----------------------------------------------------------------------
!
      MODULE module_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  The NMM component's ESMF internal state.
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
      PUBLIC :: NMM_INTERNAL_STATE                                      &
               ,WRAP_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------

      TYPE NMM_INTERNAL_STATE

        TYPE(ESMF_GridComp),DIMENSION(:),POINTER :: DOMAIN_GRID_COMP           !<-- Gridded components of all domains
!
        TYPE(ESMF_State),DIMENSION(:),POINTER :: IMP_STATE_DOMAIN              !<-- The import state of the DOMAIN components
        TYPE(ESMF_State),DIMENSION(:),POINTER :: EXP_STATE_DOMAIN              !<-- The export state of the DOMAIN components
!
      END TYPE NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_NMM_INTERNAL_STATE
!
        TYPE(NMM_INTERNAL_STATE),POINTER :: NMM_INT_STATE
!
      END TYPE WRAP_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
