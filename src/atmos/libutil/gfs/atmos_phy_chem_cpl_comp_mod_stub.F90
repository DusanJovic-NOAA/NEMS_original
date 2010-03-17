!
      module atmos_phy_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!
!**  this module holds the phy_to_chem coupler's register and run routines
!**  they are called from the main gc in module_main_grid_comp
!**  
!
!! Code Revision:
!! 11Nov 2009     Sarah Lu, First Crack
!! 18Nov 2009     Sarah Lu, Revise coupler run to do data copy
!! 23Dec 2009     Sarah Lu, The stub version
!-----------------------------------------------------------------------

      use ESMF_MOD

!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private

      public :: setservices

      contains

!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine setservices(GC, RC_REG)

!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      implicit none
      type(ESMF_cplcomp),intent(inout) :: gc         ! coupler component
!
      integer,intent(out) :: rc_reg                  ! return code for register
!
      end subroutine setservices


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine run(GC, PHY_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      implicit none
      type(ESMF_cplcomp),intent(inout) :: GC
      type(ESMF_state),  intent(inout) :: PHY_EXP_STATE    ! coupler import state
      type(ESMF_state),  intent(inout) :: CHEM_IMP_STATE   ! coupler export state
      type(ESMF_clock),  intent(in)    :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

      END subroutine run


      END module atmos_phy_chem_cpl_comp_mod



