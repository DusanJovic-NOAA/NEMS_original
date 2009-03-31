#include "../../inc/options.h"
!-----------------------------------------------------------------------
!
      module atmos_dyn_phy_cpl_comp_mod
!
!-----------------------------------------------------------------------
!
!***  this module holds the coupler's register, init, run, and finalize 
!***  routines.  they are called from the main gridded component
!***  in module_main_grid_comp.f.
!
!***  the coupler provides 2-way coupling between the dynamics and
!***  physics gridded components by transfering their export and
!***  import states between the two.
!
!-----------------------------------------------------------------------
!
      use esmf_mod
!      use module_dm_parallel,only : ids,ide,jds,jde                     &
!                                   ,ims,ime,jms,jme                     &
!                                   ,its,ite,jts,jte                     &
!                                   ,mype_share
!      use module_control,only : lm
!      use module_export_import_data
!      use atmos_err_msg_mod
!
!-----------------------------------------------------------------------
!
      implicit none
!      integer 	lm
!
!-----------------------------------------------------------------------
!
      private
!
      public :: atm_cpl_setservices
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine atm_cpl_setservices(gc_atm_cpl,rc_reg)
!
!-----------------------------------------------------------------------
!***  register the coupler component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp),intent(inout) :: gc_atm_cpl 	! coupler component
!
      integer,intent(out) :: rc_reg               	! return code for register
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      end subroutine atm_cpl_setservices
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine atm_cpl_initialize(gc_atm_cpl,imp_state,exp_state         &
                               ,clock,rc_cpl)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  set up the coupler.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  argument variables.
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp),intent(inout) :: gc_atm_cpl
      type(esmf_state),  intent(inout) :: imp_state
      type(esmf_state),  intent(inout) :: exp_state
      type(esmf_clock),  intent(in)    :: clock
!
      integer,           intent(out)   :: rc_cpl
!
!
      end subroutine atm_cpl_initialize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine atm_cpl_run(gc_atm_cpl,imp_state,exp_state                &
                             ,clock,rc_cpl)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  run the coupler to transfer data between the gridded components.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp),intent(inout) :: gc_atm_cpl
      type(esmf_state),  intent(inout) :: imp_state
      type(esmf_state),  intent(inout) :: exp_state
      type(esmf_clock),  intent(in)    :: clock
!
      integer,           intent(out)   :: rc_cpl
!
!
      end subroutine atm_cpl_run
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine atm_cpl_finalize(gc_atm_cpl,imp_state,exp_state           &
                                  ,clock,rc_cpl)
!
!-----------------------------------------------------------------------
!***  finalize the coupler.
!-----------------------------------------------------------------------
!
!
      type(esmf_cplcomp),intent(inout) :: gc_atm_cpl
      type(esmf_state),  intent(inout) :: imp_state
      type(esmf_state),  intent(inout) :: exp_state
      type(esmf_clock),  intent(in)    :: clock
!
      integer,           intent(out)   :: rc_cpl
!      
!-----------------------------------------------------------------------
!
      end subroutine atm_cpl_finalize
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      end module atmos_dyn_phy_cpl_comp_mod
!
!-----------------------------------------------------------------------
