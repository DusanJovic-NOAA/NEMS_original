!!!!!  ==========================================================  !!!!!
!!!!!             gfdl1 radiation package description              !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the gfdl1 package includes these parts:                           !
!                                                                      !
!       'radlw_gfdl1_param.f'                                          !
!       'radlw_gfdl1_datatb.f'                                         !
!       'radlw_gfdl1_main.f'                                           !
!                                                                      !
!    the 'radlw_gfdl1_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_gfdl1_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_avplank'     -- plank flux data                  !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!                                                                      !
!    the 'radlw_gfdl88_main.f' contains the main module:               !
!                                                                      !
!       'module_radlw_main'                                            !
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
!       'radlw_gfdl1_param,f'                                          !
!       'radlw_gfdl1_datatb,f'                                         !
!       'radlw_gfdl1_main.f'                                           !
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
        integer :: ilwrate, iaerlw, ithkcld, iflgco2, ico2tfs,          &
     &             iflgcld, ifch4n2o, ich4n2otf, iflgcfc

!
!  ---  set up control parameters for lw radiation
!
        parameter ( ilwrate=2 )     !===> ... lw heating rate unit selection
                                    ! =1: output in k/day
                        !(default)  ! =2: output in k/second

        parameter ( iaerlw=1 )      !===> ... control flag for aerosols
                                    ! =0: do not include aerosol effect
                        !(default)  ! =1: aeros opt prop are calc for each spectral band
                                    ! =2: broad band aeros opt prop are used for all bands

        parameter ( iflgcld=1 )     !===> ... control flag for cloud optical proerties
                                    ! =0: input cloud optical depth
                        !(default)  ! =1: input cloud condensates: cwp, cip, etc...

        parameter ( ithkcld=0 )     !===> ... control flag for thick cloud adjustment
                        !(default)  ! =0: no peseudo-conv adj for maxi ovlp clouds
                                    ! =1: do peseudo-conv adj for maxi ovlp clouds

        parameter ( iflgco2=1 )     !===> ... control flag for co2 gas effect
                        !(default)  ! =1: logarithmic press interp co2 coeff
                                    ! =2: linear press interp co2 coeff

        parameter ( ico2tfs=2 )     !===> ... control flag for co2 transmission func
                                    ! =1: calc and write out co2 transmission functions
                        !(default)  ! =2: read in co2 transm funcs from pre calc table

        parameter ( ifch4n2o=1 )    !===> ... control flag for ch4 and n2o gases
                                    ! =0: do not include ch4 and n2o gases
                        !(default)  ! =1: ch4,n2o gases without lbl temp interpolation
                                    ! =2: ch4,n2o gases with lbl temp interpolation

        parameter ( ich4n2otf=2 )   !===> ... control flag for ch4 transmission func
                                    ! =1: calc and write out ch4 and n2o transm funcs
                        !(default)  ! =2: read in ch4 and n2o transm funcs from pre calc tables

        parameter ( iflgcfc=1 )     !===> ... control flag for cfc gases effect
                                    ! =0: do not include cfc gases
                        !(default)  ! =1: include 4 types of cfc effect

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
!  ---  define type constructs for table-lookup
!
      type :: tab1_type
        real (kind=kind_phys), dimension(:,:), pointer :: vae,td,md,cd
      end type tab1_type

      type :: tab3_type
        real (kind=kind_phys), dimension(:,:), pointer :: vae, td
      end type tab3_type

      type :: axis_type
        integer :: first_col
        real (kind=kind_phys) :: min_val, max_val, tab_inc
      end type axis_type
!
!  ---  parameter constants for lw band structures
!
      integer, parameter :: NBCO215     = 3
      integer, parameter :: NSTDCO2LVLS = 496  ! num of levs at which lbl co2
                                               ! transmission functions computed
      integer, parameter :: IOFFSET     = 32   ! for ckd h2o continuum
      integer, parameter :: NBLY        = 48   ! for ckd h2o continuum
      integer, parameter :: NBLX        = 48
      integer, parameter :: NBLW        = 300
      integer, parameter :: NBLWCFC     = 8    ! num of bands with cfc included
      integer, parameter :: IOFFH2O     = 16   ! offset from the absorption tables for
                                               ! freq calc
      integer, parameter :: NBANDS      = 8
      integer, parameter :: NBCNTN      = 40   ! bands for ckd h2o continuum

      integer, parameter :: NTTABH2O    = 28
      integer, parameter :: NUTABH2O    = 180

      integer, parameter :: NBFL        = 12   ! num of bands in snow water paramterization
      integer, parameter :: NLWCLDB     = 7    ! num of freq bands for cloud prop

!! ---  due to complicity of the R-T scheme, the final results will not be able to be 
!       specified in certain spectral bands, set num of bands to 1 for simplicity.

      integer, parameter :: NBDLW       = 1   

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

