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

      SUBROUTINE GFS_INTEGRATE(gc_gfs_dyn                               &
                                 ,gc_gfs_phy                            &
                                 ,gc_atm_cpl                            &
                                 ,wrt_comps                             &
                                 ,imp_gfs_dyn                           &
                                 ,exp_gfs_dyn                           &
                                 ,imp_gfs_phy                           &
                                 ,exp_gfs_phy                           &
                                 ,imp_gfs_wrt                           &
                                 ,exp_gfs_wrt                           &
                                 ,CLOCK_MAIN                            &
                                 ,OUTPUT_INTERVAL                       &
                                 ,quilting                              &
                                 ,WRITE_GROUP_READY_TO_GO               &
                                 ,CURRTIME                              &
                                 ,STARTTIME                             &
                                 ,NTIMESTEP                             &
                                 ,TIMESTEP                              &
                                 ,DFIHR                                 &
                                 ,MYPE                                  &
                                 ,PHYSICS_ON)

!
!-----------------------------------------------------------------------
!
!

      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: gc_gfs_phy
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: gc_atm_cpl
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: wrt_comps(:)
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_phy,exp_gfs_phy
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_wrt,exp_gfs_wrt
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_MAIN                         !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)             :: CURRTIME                           !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)             :: STARTTIME
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR, NTIMESTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(IN)          :: TIMESTEP                      !<-- The ESMF timestep (s)
      LOGICAL,INTENT(IN)                     :: PHYSICS_ON
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: output_interval
      LOGICAL,INTENT(IN)                     :: QUILTING
      INTEGER(KIND=KINT),INTENT(INOUT)       :: WRITE_GROUP_READY_TO_GO
!
      END SUBROUTINE GFS_INTEGRATE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_INTEGRATE
