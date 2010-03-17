!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_CORE_SETUP
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
      PUBLIC :: GFS_SETUP       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!-----------------------------------------------------------------------
!
      CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE GFS_SETUP(gc_atm,grid_atmos)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
!
      type(ESMF_gridcomp),intent(inout) :: gc_atm
      type(ESMF_grid),intent(out)  :: grid_atmos    ! the ESMF grid for the integration attached to
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_SETUP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_CORE_SETUP
!
!-----------------------------------------------------------------------

