!-----------------------------------------------------------------------
      MODULE MODULE_WRITE_INTERNAL_STATE_GFS
!-----------------------------------------------------------------------
!***  THE INTERNAL STATE OF THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!***
!***  HISTORY
!***
!       xx Feb 2007:  W. Yang - Originator
!       14 Jun 2007:  T. Black - Name revisions
!       14 Aug 2007:  T. Black - Some pointers changed to arrays
!       11 Sep 2007:  T. Black - Updates for quilting
!       15 Aug 2008:  J. Wang  - Add NEMSIO variables
!       16 Sep 2008:  J. Wang  - 3-D output arrays revert to 2-D
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      TYPE WRITE_INTERNAL_STATE_GFS

      INTEGER :: IDUMMY

      END TYPE WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------------------
!***  THIS STATE IS SUPPORTED BY C POINTERS BUT NOT F90 POINTERS
!***  THEREFORE WE NEED THIS WRAP.
!-----------------------------------------------------------
!
      TYPE WRITE_WRAP_GFS
        TYPE(WRITE_INTERNAL_STATE_GFS),POINTER :: WRITE_INT_STATE
      END TYPE WRITE_WRAP_GFS

!-----------------------------------------------------------
!
      END MODULE MODULE_WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------
