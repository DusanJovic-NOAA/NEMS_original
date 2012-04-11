#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE module_GENOCN_INTERNAL_STATE
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
      PUBLIC :: GENOCN_INTERNAL_STATE                                    &
               ,WRAP_GENOCN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE GENOCN_INTERNAL_STATE
          TYPE(ESMF_GridComp) :: GENOCN_GRID_COMP

          TYPE(ESMF_State   ) :: GENOCN_IMP_STATE                          &  !<-- Import/export states for GENOCN
                                ,GENOCN_EXP_STATE

          TYPE(ESMF_Clock   ) :: CLOCK_GENOCN

#ifdef ESMF_3
          TYPE(ESMF_LOGICAL)  :: Cpl_flag
#else
          LOGICAL             :: Cpl_flag
#endif

          INTEGER             :: MYPE                                  &  !<-- Each MPI task ID
                                ,WRITE_GROUP_READY_TO_GO                  !<-- The write group to use
 
          LOGICAL             :: QUILTING                                 !<-- Is asynchronous quilting specified?
 
          TYPE(ESMF_GridComp), DIMENSION(:), POINTER :: WRT_COMPS
 
          TYPE(ESMF_TimeInterval) :: TIMEINTERVAL_GENOCN_OUTPUT              !<-- The ESMF time interval between GENOCN history output
      END TYPE GENOCN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_GENOCN_INTERNAL_STATE
!
      TYPE(GENOCN_INTERNAL_STATE),POINTER :: GENOCN_INT_STATE
!
      END TYPE WRAP_GENOCN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_GENOCN_INTERNAL_STATE
!
!-----------------------------------------------------------------------

