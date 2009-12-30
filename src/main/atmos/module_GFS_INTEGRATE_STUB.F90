!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_INTEGRATE
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
!
      PUBLIC :: GFS_INTEGRATE
!
!-----------------------------------------------------------------------
      INCLUDE '../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GFS_INTEGRATE(GC_GFS_DYN                               &
                              ,GC_GFS_PHY                               &
                              ,GC_GFS_CHEM                              &
                              ,GC_ATM_CPL                               &
                              ,GC_DYN2CHEM_CPL                          &
                              ,GC_PHY2CHEM_CPL                          &
                              ,WRT_COMPS                                &
                              ,IMP_GFS_DYN                              &
                              ,EXP_GFS_DYN                              &
                              ,IMP_GFS_PHY                              &
                              ,EXP_GFS_PHY                              &
                              ,IMP_GFS_CHEM                             &
                              ,EXP_GFS_CHEM                             &
                              ,IMP_GFS_WRT                              &
                              ,EXP_GFS_WRT                              &
                              ,CLOCK_MAIN                               &
                              ,OUTPUT_INTERVAL                          &
                              ,QUILTING                                 &
                              ,WRITE_GROUP_READY_TO_GO                  &
                              ,CURRTIME                                 &
                              ,STARTTIME                                &
                              ,NTIMESTEP                                &
                              ,TIMESTEP                                 &
                              ,DFIHR                                    &
                              ,MYPE                                     &
                              ,PHYSICS_ON                               &
                              ,CHEMISTRY_ON)

!
!-----------------------------------------------------------------------
!
!

      TYPE(ESMF_GridComp),INTENT(INOUT)      :: GC_GFS_DYN
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: GC_GFS_PHY
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: GC_GFS_CHEM
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: GC_ATM_CPL
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: GC_DYN2CHEM_CPL
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: GC_PHY2CHEM_CPL
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: WRT_COMPS(:)
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_GFS_DYN,EXP_GFS_DYN
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_GFS_PHY,EXP_GFS_PHY
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_GFS_CHEM,EXP_GFS_CHEM
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_GFS_WRT,EXP_GFS_WRT
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_MAIN                  !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)          :: CURRTIME                    !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      INTEGER(kind=KINT),INTENT(INOUT)       :: DFIHR, NTIMESTEP
      INTEGER(kind=KINT),INTENT(IN)          :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(IN)     :: TIMESTEP                    !<-- The ESMF timestep (s)
      TYPE(ESMF_Logical),INTENT(IN)          :: PHYSICS_ON
      TYPE(ESMF_Logical),INTENT(IN)          :: CHEMISTRY_ON
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: OUTPUT_INTERVAL
      LOGICAL(kind=KLOG),INTENT(IN)          :: QUILTING
      INTEGER(kind=KINT),INTENT(INOUT)       :: WRITE_GROUP_READY_TO_GO
!
      END SUBROUTINE GFS_INTEGRATE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_INTEGRATE
