#include "../../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#define ESMF_LogFoundError ESMF_LogMsgFoundError
#else
#define ESMF_520r
#endif

!-----------------------------------------------------------------------
!
      MODULE module_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  The NMM component's ESMF internal state.
!-----------------------------------------------------------------------
!
#ifdef ESMF_520r
      USE esmf
#else
      USE esmf_mod
#endif
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
