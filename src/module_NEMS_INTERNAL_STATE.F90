#include "./ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#define ESMF_LogFoundError ESMF_LogMsgFoundError
#else
#define ESMF_520r
#endif

!-----------------------------------------------------------------------
!
      MODULE module_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Contents of the ESMF internal state of the NEMS component.
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
      PUBLIC :: NEMS_INTERNAL_STATE                                     &
               ,WRAP_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE NEMS_INTERNAL_STATE
!
        REAL :: DUMMY1
!
      END TYPE NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_NEMS_INTERNAL_STATE
!
        REAL :: DUMMY2
!
        TYPE(NEMS_INTERNAL_STATE),POINTER :: NEMS_INT_STATE
!
      END TYPE WRAP_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
