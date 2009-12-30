!
      module atmos_dyn_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!
!**  this module holds the dyn_to_chem coupler's register and run routines
!**  
!
!! Code Revision:
!! 18Nov 2009     Sarah Lu, First Crack
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
!
      implicit none
!
      type(ESMF_cplcomp),intent(inout) :: gc         ! coupler component
!
      integer,intent(out) :: rc_reg                  ! return code for register
!
      end subroutine setservices


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine run(GC, DYN_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
      implicit none
!
      type(ESMF_cplcomp),intent(inout) :: GC
      type(ESMF_state),  intent(inout) :: DYN_EXP_STATE    ! coupler import state
      type(ESMF_state),  intent(inout) :: CHEM_IMP_STATE   ! coupler export state
      type(ESMF_clock),  intent(in)    :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

      END subroutine run

      END module atmos_dyn_chem_cpl_comp_mod



