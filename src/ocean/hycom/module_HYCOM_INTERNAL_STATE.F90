#include "../../ESMFVersionDefine.h"

! February 2011    Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                              ESMF 5 library and the the ESMF 3.1.0rp2 library.
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
      MODULE module_HYCOM_INTERNAL_STATE
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
      PUBLIC :: HYCOM_INTERNAL_STATE                                    &
               ,WRAP_HYCOM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE HYCOM_INTERNAL_STATE
          TYPE(ESMF_GridComp) :: HYCOM_GRID_COMP

          TYPE(ESMF_State   ) :: HYCOM_IMP_STATE                          &  !<-- Import/export states for HYCOM
                                ,HYCOM_EXP_STATE

          TYPE(ESMF_Clock   ) :: CLOCK_HYCOM

#ifdef ESMF_3
          TYPE(ESMF_LOGICAL)  :: Cpl_flag
#else
          LOGICAL             :: Cpl_flag
#endif

          INTEGER             :: MYPE                                  &  !<-- Each MPI task ID
                                ,WRITE_GROUP_READY_TO_GO                  !<-- The write group to use
 
          LOGICAL             :: QUILTING                                 !<-- Is asynchronous quilting specified?
 
          TYPE(ESMF_GridComp), DIMENSION(:), POINTER :: WRT_COMPS
 
          TYPE(ESMF_TimeInterval) :: TIMEINTERVAL_HYCOM_OUTPUT              !<-- The ESMF time interval between HYCOM history output
      END TYPE HYCOM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_HYCOM_INTERNAL_STATE
!
      TYPE(HYCOM_INTERNAL_STATE),POINTER :: HYCOM_INT_STATE
!
      END TYPE WRAP_HYCOM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_HYCOM_INTERNAL_STATE
!
!-----------------------------------------------------------------------

