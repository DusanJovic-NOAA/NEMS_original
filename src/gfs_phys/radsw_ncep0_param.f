!!!!!  ==========================================================  !!!!!
!!!!!            sw-ncep0 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the sw-ncep0 package includes these parts:                         !
!                                                                      !
!      'radsw_ncep0_param.f'                                           !
!      'radsw_ncep0_datatb.f'                                          !
!      'radsw_ncep0_main.f'                                            !
!                                                                      !
!   the 'radsw_ncep0_param.f' contains:                                !
!                                                                      !
!      'module_radsw_parameters'  -- band parameters set up            !
!      'module_radsw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radsw_ncep0_datatb.f' contains:                               !
!                                                                      !
!      'module_radsw_co2tab'      -- co2 lookup tables                 !
!      'module_radsw_cldprtb'     -- cloud property coefficients table !
!                                                                      !
!   the 'radsw_ncep0_main.f' contains:                                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main ncep0 sw radiation routine                 !
!      'rswinit'    -- initialization routine                          !
!                                                                      !
!   all the sw radiation subprograms become contained subprograms      !
!   in module 'module_radsw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_ncep0_param.f'                                           !
!      'radsw_ncep0_datatb.f'                                          !
!      'radsw_ncep0_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use sw modules     !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radsw_cntr_para      !
!........................................!
!
        implicit   none
!
        integer :: iaersw, ioxysw, ico2sw, ih2osw, ioznsw
        integer :: iswrate, ibndsw, iflagc

!
!  ---  set up control parameters for sw radiation
!
        parameter ( iswrate=2 )     !===> ... flag for heating rate unit
                                    ! =1: output in k/day
                        ! (default) ! =2: output in k/second
        parameter ( iaersw=1 )      !===> ... flag for aerosols
                                    ! =0: without aerosol effect
                        ! (default) ! =1: include aerosol effect

        parameter ( ioxysw=1 )      !===> ... flag for o2 absorption
                                    ! =0: without o2 absorption
                        ! (default) ! =1: include o2 absorption

        parameter ( ico2sw=1 )      !===> ... flag for co2 absorption
                                    ! =0: without co2 absorption
                        ! (default) ! =1: include co2 absorption

        parameter ( ih2osw=1 )      !===> ... flag for h2o absorption
                                    ! =0: without h2o absorption
                        ! (default) ! =1: include h2o absorption

        parameter ( ioznsw=1 )      !===> ... flag for o3 absorption
                                    ! =0: without o3 absorption
                        ! (default) ! =1: include o3 absorption

        parameter ( ibndsw=0 )      !===> ... band selections in nir spectrum
                        ! (default) ! =0: faster computation, use 1 nir band
                                    ! =1: better accuracy,    use 3 nir bands

        parameter ( iflagc=1 )      !===> ... flag for cloud property method
                                    ! =0: input cloud opt depth, fixed ssa, asy
                        ! (default) ! =1: input cwp/cip, t-adj coeff for all bands
                                    ! =2: input cwp/cip, chou (1999) coeff for uv
                                    !     and 3 nir bands, or gfs's tuned 1 nir
                                    !     band coeff.
                                    ! =3: input cwp/cip, chou (2002) coeff for uv
                                    !     and 3 nir bands, or averaged 1 nir
                                    !     band coeff.

!
!........................................!
      end module module_radsw_cntr_para  !
!========================================!



!========================================!
      module module_radsw_parameters     !
!........................................!

      use machine,                 only : kind_phys
      use module_radsw_cntr_para,  only : ibndsw

      implicit   none
!
      public
!
!  ---  define type construct for radiation fluxes at toa
!
      type :: topfsw_type
        real (kind=kind_phys) :: upfxc         ! total sky upward flux at toa
        real (kind=kind_phys) :: dnfxc         ! total sky downward flux at toa
        real (kind=kind_phys) :: upfx0         ! clear sky upward flux at toa
      end type
!
!  ---  define type construct for radiation fluxes at surface
!
      type :: sfcfsw_type
        real (kind=kind_phys) :: upfxc         ! total sky upward flux at sfc
        real (kind=kind_phys) :: dnfxc         ! total sky downward flux at sfc
        real (kind=kind_phys) :: upfx0         ! clear sky upward flux at sfc
        real (kind=kind_phys) :: dnfx0         ! clear sky downward flux at sfc
      end type
!
!  ---  define type construct for optional radiation flux profiles
!
      type :: profsw_type
        real (kind=kind_phys) :: upfxc         ! total sky level upward flux
        real (kind=kind_phys) :: dnfxc         ! total sky level downward flux
        real (kind=kind_phys) :: upfx0         ! clear sky level upward flux
        real (kind=kind_phys) :: dnfx0         ! clear sky level downward flux
      end type
!
!  ---  define type construct for optional component downward fluxes at surface
!
      type :: cmpfsw_type
        real (kind=kind_phys) :: uvbfc         ! total sky downward uv-b flux at sfc
        real (kind=kind_phys) :: uvbf0         ! clear sky downward uv-b flux at sfc
        real (kind=kind_phys) :: nirbm         ! sfc downward nir direct beam flux
        real (kind=kind_phys) :: nirdf         ! sfc downward nir diffused flux
        real (kind=kind_phys) :: visbm         ! sfc downward uv+vis direct beam flx
        real (kind=kind_phys) :: visdf         ! sfc downward uv+vis diffused flux
      end type
!
!  ---  parameter constants for sw band structures
!
      integer, parameter :: NK0 = 10    ! num of k-values in nir spectrum 1-band config
      integer, parameter :: NKB = 29    ! num of k-values in nir spectrum 3-band config
      integer, parameter :: NVB = 8     ! num of bands in uv+vis spectrum
      integer, parameter :: NRB = 4     ! num of bands in nir spectrum (1-b + 3-b)

      integer, parameter :: NBK0 = NK0 +NVB  ! 10 k (1 nir); 8 uv,vis
      integer, parameter :: NBK1 = NBK0+NKB  ! 10 k (1 nir); 8 uv+vis; 29 k (3 nir)
      integer, parameter :: NBD0 = 1   +NRB  ! 1 uv+vis + 1 single nir + 3 nir
      integer, parameter :: NBD1 = NVB +NRB  ! 8 uv+vis + 1 singlw nir + 3 nir
      integer, parameter :: NK0P = NK0 +1
      integer, parameter :: NBK0P= NBK0+1

      integer, parameter :: NSWSTR = ibndsw + 1
      integer, parameter :: NSWEND = ibndsw*3 + 1 + NVB
      integer, parameter :: NBDSW  = NSWEND - NSWSTR + 1

!
!  ---  band spectrum structures (start and end wavenumbers for each band in cm**-1)
!       ** note: band-3 actual has 2 sub intervals (40820-44440;35700-38460)
!          when ibndsw=0 (1 nir band), NBDSW=9, use wvnum bnad 1-9
!               ibndsw=1 (3 nir bands), NBDSW=11, use wvnum bnad 2-12

      real (kind=kind_phys) :: wvnum1(NBD1), wvnum2(NBD1)
      data wvnum1  /    1000., 44441., 35701., 38461., 33901., 32261.,  &
     &                 31251., 25001., 14281.,  8201.,  4401.,  1000. /
      data wvnum2  /   14280., 57140., 44440., 40820., 35700., 33900.,  &
     &                 32260., 31250., 25000., 14280.,  8200.,  4400. /

!
!........................................!
      end module module_radsw_parameters !
!========================================!
