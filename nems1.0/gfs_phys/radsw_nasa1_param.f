!!!!!  ==========================================================  !!!!!
!!!!!            sw-nasa1 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the sw-nasa1 package includes these parts:                         !
!                                                                      !
!      'radsw_nasa1_param.f'                                           !
!      'radsw_nasa1_datatb.f'                                          !
!      'radsw_nasa1_main.f'                                            !
!                                                                      !
!   the 'radsw_nasa1_param.f' contains:                                !
!                                                                      !
!      'module_radsw_parameters'  -- band parameters set up            !
!      'module_radsw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radsw_nasa1_datatb.f' contains:                               !
!                                                                      !
!      'module_radsw_co2tab'      -- co2 lookup tables                 !
!      'module_radsw_cldsca'      -- cloud scaling coefficients table  !
!                                                                      !
!   the 'radsw_nasa1_main.f' contains:                                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main nasa1 sw radiation routine                 !
!      'rswinit'    -- initialization routine                          !
!                                                                      !
!   all the sw radiation subprograms become contained subprograms      !
!   in module 'module_radsw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_nasa1_param.f'                                           !
!      'radsw_nasa1_datatb.f'                                          !
!      'radsw_nasa1_main.f'                                            !
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
        integer :: iswrate, iaersw, ioxysw, ico2sw, iflagc

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

        parameter ( iflagc=2 )      !===> ... flag for cloud scaling method
                                    ! =1: no scaling, clouds either 0 or 1
                        ! (default) ! =2: fractional clouds with scaling

!
!........................................!
      end module module_radsw_cntr_para  !
!========================================!



!========================================!
      module module_radsw_parameters     !
!........................................!

      use machine,                 only : kind_phys

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
        real (kind=kind_phys) :: visbm         ! sfc downward vis(+uv) direct beam flx
        real (kind=kind_phys) :: visdf         ! sfc downward vis(+uv) diffused flux
!       real (kind=kind_phys) :: uvbm          ! sfc downward uv direct beam flx
!       real (kind=kind_phys) :: uvdf          ! sfc downward uv diffused flux
      end type
!
!  ---  parameter constants for sw band structures
!
      integer, parameter :: NK0    = 10    ! num of k-values in each nir band
      integer, parameter :: NUVVIS = 8     ! num of bands in uv+vis spectrum
      integer, parameter :: NIRBND = 3     ! num of bands in nir spectrum
      integer, parameter :: NUVBS  = 4     ! staring band num for uv-b spectrum
      integer, parameter :: NUVBE  = 6     ! ending band num for uv-b spectrum

      integer, parameter :: NSPCLD = 3     ! num of species in cloud condensates

      integer, parameter :: NBDSW = NUVVIS+NIRBND
      integer, parameter :: NSWSTR= 1
      integer, parameter :: NSWEND= NBDSW
!
!  ---  band spectrum structures (boundaries of wavenumber in cm**-1)
!       ** note: band-2 actual has 2 sub intervals (40820-44440;35700-38460)

      real (kind=kind_phys) :: wvnum1(NBDSW), wvnum2(NBDSW)
      data wvnum1  /   44441., 35701., 38461., 33901., 32261., 31251.,  &
     &                 25001., 14281.,  8201.,  4401.,  1000.   /
      data wvnum2  /   57140., 44440., 40820., 35700., 33900., 32260.,  &
     &                 31250., 25000., 14280.,  8200.,  4400.   /

!
!........................................!
      end module module_radsw_parameters !
!========================================!
