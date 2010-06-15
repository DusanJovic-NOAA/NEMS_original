! !module: gfs_physics_grid_comp_mod --- 
!                       esmf gridded component of gfs physics
!
! !description: gfs physics gridded component main module.
!
! !revision history:
!
!  july     2007     shrinivas moorthi
!  november 2007     hann-ming henry juang 
!                           
!
! !interface:
!
      module gfs_physics_grid_comp_mod
 
!!uses:
!------
      use esmf_mod

!      use gfs_physics_err_msg_mod,        ONLY: gfs_physics_err_msg,        &
!                                                gfs_physics_err_msg_final
!      use gfs_physics_initialize_mod,     ONLY: gfs_physics_initialize
!      use gfs_physics_run_mod,            ONLY: gfs_physics_run
!      use gfs_physics_finalize_mod,       ONLY: gfs_physics_finalize
!      USE gfs_physics_getcf_mod,          ONLY: gfs_physics_getcf
!      USE gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state, &
!                                                gfs_phy_wrap
!      USE mpi_def,                        ONLY: mpi_comm_all
!      USE layout1,                        ONLY: me
!      USE date_def,                       ONLY: idate, fhour
!      USE namelist_physics_def,           ONLY: fhini, fhmax

!      implicit none

!      private   ! by default, data is private to this module

      public gfs_phy_setservices	! only set service is public

!eop
!-------------------------------------------------------------------------


      contains


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_setservices --- 
!           set services for gfs physics gridded component.
! 
! !interface:
!
      subroutine gfs_phy_setservices (gc_gfs_phy, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp), intent(in)  :: gc_gfs_phy 	! gridded component
      integer,             intent(out) :: rc    	! return code
     

      end subroutine gfs_phy_setservices





!----------------------------------------------------------------------
!bop
! !routine:  gfs_phy_initialize --- initialize routine to initialize 
!                                   and set up the gfs running job.
!
! !description: this subroutine initializes the gfs running before
!               the main running loop.
!
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     h.-m. h. juang
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_phy_initialize(gc_gfs_phy, 			&
                                   imp_gfs_phy, exp_gfs_phy, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
!      use gfs_physics_states_mod
!      use gfs_physics_grid_create_mod
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp), intent(inout) :: gc_gfs_phy 
      type(esmf_state),    intent(inout) :: imp_gfs_phy
      type(esmf_state),    intent(inout) :: exp_gfs_phy
      type(esmf_clock),    intent(inout) :: clock

!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc  

      end subroutine gfs_phy_initialize





!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_run --- 
!           main grid component routine to run the gfs physics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs physics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  december 2007     juang
!
! !interface:
!

      subroutine gfs_phy_run(gc_gfs_phy, 				&
                            imp_gfs_phy, exp_gfs_phy, clock, rc)

!      use gfs_physics_states_mod
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout) :: gc_gfs_phy   
      type(esmf_state),    intent(in)    :: imp_gfs_phy 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock),    intent(inout) :: clock
      type(esmf_state),    intent(inout) :: exp_gfs_phy
      integer,             intent(out)   :: rc   
      
      end subroutine gfs_phy_run


!----------------------------------------------------------------------
!bop
!
! !routine: finalize --- finalizing routine to finish the 
!                        gfs running job.
!
! !description: this subroutine will finish the gfs computations,
! !             and will release the memory space.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     juang for dynamics only
!  july     2007     juang for physics only
!
! !interface:

      subroutine gfs_phy_finalize(gc_gfs_phy, 				&
                                 imp_gfs_phy, exp_gfs_phy, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout)  :: gc_gfs_phy
      type(esmf_state),    intent(inout)  :: imp_gfs_phy
      type(esmf_state),    intent(inout)  :: exp_gfs_phy
      type(esmf_clock),    intent(inout)  :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

      end subroutine gfs_phy_finalize
      end module gfs_physics_grid_comp_mod
