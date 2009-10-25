      module module_digital_filter_nmm
!
! a generic digital filter for any model under ESMF 
!
! March 2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
!
      use esmf_mod
      implicit none

      contains

! ---------------------------------------------------------------
! subroutine for dynamics
! ---------------------------------------------------------------
      subroutine digital_filter_dyn_init_nmm(dyn_state,ndfistep,NUM_WATER,NUM_TRACERS)
      type(esmf_state), intent(in) :: dyn_state   
      INTEGER, intent(in)          :: ndfistep
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS

      end subroutine digital_filter_dyn_init_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_sum_nmm(dyn_state,MEAN_ON,NUM_WATER,NUM_TRACERS)
      USE MODULE_DM_PARALLEL
!
      type(esmf_state), intent(in)  :: dyn_state 
      INTEGER, intent(in)          :: MEAN_ON,NUM_WATER,NUM_TRACERS

      end subroutine digital_filter_dyn_sum_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_average_nmm(dyn_state,NUM_WATER,NUM_TRACERS)
!
      type(esmf_state), intent(inout) :: dyn_state
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS


      end subroutine digital_filter_dyn_average_nmm
!----------------------------------------------------------------------------
      subroutine digital_filter_phy_init_nmm(phy_state)
!
      type(esmf_state), intent(in) :: phy_state
!
      end subroutine digital_filter_phy_init_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_phy_save_nmm(phy_state)
!
      type(esmf_state), intent(in) :: phy_state
      end subroutine digital_filter_phy_save_nmm

! ---------------------------------------------------------------
      subroutine digital_filter_phy_restore_nmm(phy_state)
!
      type(esmf_state), intent(inout) :: phy_state
      
      end subroutine digital_filter_phy_restore_nmm

      end module module_digital_filter_nmm
