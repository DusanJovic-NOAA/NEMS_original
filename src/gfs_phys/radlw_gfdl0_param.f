!!!!!  ==========================================================  !!!!!
!!!!!            gfdl0 radiation package description               !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the gfdl0 package includes these parts:                           !
!                                                                      !
!       'radlw_gfdl0_param.f'                                          !
!       'radlw_gfdl0_datatb.f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    the 'radlw_gfdl0_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_gfdl0_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_avplank'     -- plank flux data                  !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!                                                                      !
!    the 'radlw_gfdl0_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main rrtm lw radiation routine                  !
!       'rlwinit'   -- initialization routine                          !
!                                                                      !
!    all the lw radiation subprograms become contained subprograms     !
!    in module 'module_radlw_main' and many of them are not directly   !
!    accessable from places outside the module.                        !
!                                                                      !
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_gfdl0_param.f'                                          !
!       'radlw_gfdl0_datatb.f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radlw_cntr_para      !
!........................................!
!
        implicit   none
!
        integer :: ilwrate, iaerlw, iflagcld, ico2tran

!
!  ---  set up control parameters for lw radiation
!
        parameter ( ilwrate=2 )     !===> ... lw heating rate unit selection
                                    ! =1: output in k/day
                        !(default)  ! =2: output in k/second

        parameter ( iaerlw=0 )      !===> ... control flag for aerosols ** not yet **
                        !(default)  ! =0: do not include aerosol effect
!not yet                !           ! >0: include aerosol effect calc for one broad band

        parameter ( iflagcld=1 )    !===> ... control flag for cloud optical proerties
                                    ! =0: input cloud opt depth, ignor iflagice setting
                        !(default)  ! =1: input cwp, cip, use ncar ccm3 method

!not    parameter ( ico2tran=0 )    !===> ... control flag for co2 transm data
                        !(default)  ! =0: read in pre-calculated co2 transmission data
                                    ! =1: compute co2 transm data during initialization

!
!........................................!
      end module module_radlw_cntr_para  !
!========================================!




!========================================!
      module module_radlw_parameters     !
!........................................!

      use machine,                 only : kind_phys

      implicit none
!
      public
!
!  ---  define type construct for radiation fluxes at toa
!
      type :: topflw_type
        real (kind=kind_phys) :: upfxc         ! total sky upward flux at toa
        real (kind=kind_phys) :: upfx0         ! clear sky upward flux at toa
      end type
!
!  ---  define type construct for radiation fluxes at surface
!
      type :: sfcflw_type
        real (kind=kind_phys) :: upfxc         ! total sky upward flux at sfc
        real (kind=kind_phys) :: dnfxc         ! total sky downward flux at sfc
        real (kind=kind_phys) :: dnfx0         ! clear sky downward flux at sfc
      end type
!
!  ---  define type construct for optional radiation flux profiles
!
      type :: proflw_type
        real (kind=kind_phys) :: netfc         ! level net flux for total sky
        real (kind=kind_phys) :: netf0         ! level net flux for clear sky
      end type
!
!  ---  parameter constants for lw band structures
!
      integer, parameter :: NBLX    = 47
      integer, parameter :: NBLY    = 15
      integer, parameter :: NBLW    = 163
      integer, parameter :: N5040   = 5040

      integer, parameter :: NBDLW   = 1        ! set as 1 to simplify lw
                                               ! aerosol calculations

!  ---  band spectrum structures (boundaries of wavenumber in cm**-1)
!  ***  to avoid complication in computation, treat the quantities
!       such as aerosol optical properties as in one single broad band
      real (kind=kind_phys) :: wvnlw1(NBDLW), wvnlw2(NBDLW)
!     data  wvnlw1 / 250.0 /        ! corresponding to 40 mu
      data  wvnlw1 / 400.0 /        ! corresponding to 25 mu
      data  wvnlw2 / 2500. /        ! corresponding to 4  mu

!........................................!
      end module module_radlw_parameters !
!========================================!

