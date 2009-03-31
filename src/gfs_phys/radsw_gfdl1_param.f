!!!!!  ==========================================================  !!!!!
!!!!!            sw-gfdl1 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the sw-gfdl1 package includes these parts:                         !
!                                                                      !
!      'radsw_gfdl1_param.f'                                           !
!      'radsw_gfdl1_datatb.f'                                          !
!      'radsw_gfdl1_main.f'                                            !
!                                                                      !
!   the 'radsw_gfdl1_param.f' contains:                                !
!                                                                      !
!      'module_radsw_parameters'  -- band parameters set up            !
!      'module_radsw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radsw_gfdl1_datatb.f' contains:                               !
!                                                                      !
!      'module_radsw_bandtbl'     -- band structure and data           !
!      'module_radsw_cldprtb'     -- cloud property coefficients table !
!                                                                      !
!   the 'radsw_gfdl1_main.f' contains the main module:                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main gfdl sw radiation routine                  !
!      'rswinit'    -- initialization routine                          !
!                                                                      !
!   all the sw radiation subprograms become contained subprograms      !
!   in module 'module_radsw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_gfdl1_param.f'                                           !
!      'radsw_gfdl1_datatb.f'                                          !
!      'radsw_gfdl1_main.f'                                            !
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
        integer :: iswrate, iaersw, NSOLWG

!
!  ---  set up control parameters for sw radiation
!
        parameter ( iswrate=2 )     !===> ... flag for heating rate unit
                                    ! =1: output in k/day
                        ! (default) ! =2: output in k/second

        parameter ( iaersw=1 )      !===> ... flag for aerosols
                                    ! =0: without aerosol effect
                        ! (default) ! =1: include aerosol effect

        parameter ( NSOLWG=2 )      !===> ... num of gaussian quadrature weights
                                    ! =1, 2, 4, or 8 only
                        ! (default) ! =2

!
!........................................!
      end module module_radsw_cntr_para  !
!========================================!



!========================================!
      module module_radsw_parameters     !
!........................................!
!
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
        real (kind=kind_phys) :: visbm         ! sfc downward uv+vis direct beam flx
        real (kind=kind_phys) :: visdf         ! sfc downward uv+vis diffused flux
      end type
!
!  ---  parameter constants for sw band structures
!
      integer, parameter :: NBANDS      = 25   ! number of frequency bands
      integer, parameter :: NFRQPTS     = 72   ! number of pseudo-monochromatic frequencies
      integer, parameter :: NH2OBANDS   = 14   ! number of h2o bands
      integer, parameter :: NSTREAMS    =  4   ! number of streams
      integer, parameter :: FIRSTRAYBAND=  9   ! first band number where includes raleigh
                                               ! scattering contribution
      integer, parameter :: NIRBANDS    = 10   ! number of bands in NIR (for sfc albedo
      integer, parameter :: NINTSOLAR   =151   ! number of wavenumber regions where the
                                               ! solar flux is constant
      integer, parameter :: NLIQCLDV    = 24   ! number of scat intervals for cloud drops
      integer, parameter :: NICECLDV    = 25   ! number of scat intervals for ice crystals
      integer, parameter :: NRAINCLDV   = 4    ! number of scat intervals for rain drops
      integer, parameter :: NSNOWCLDV   = 6    ! number of scat intercals for snow
      integer, parameter :: NUVBSTR     = 17   ! starting band num for uv-b spectrum
      integer, parameter :: NUVBEND     = 21   ! ending band num for uv-b spectrum

      integer, parameter :: TOT_WVNUMS  = 57600

      integer, parameter :: NBDSW       = NBANDS
      integer, parameter :: NSWSTR      = 1
      integer, parameter :: NSWEND      = NBANDS


!
!  ---  band spectrum structures (boundaries of wavenumber in cm**-1)
!
      real (kind=kind_phys) :: wvnum1(NBDSW), wvnum2(NBDSW)
      data wvnum1  /     1.,  2501.,  2901.,  3401.,  4201.,  4701.,    &
     &        5601.,  6201.,  8201., 11501., 14601., 16701., 20001.,    &
     &       22301., 24601., 27501., 30001., 31901., 33001., 33801.,    &
     &       34501., 35301., 36501., 40001., 43301.  /
      data wvnum2  /  2500.,  2900.,  3400.,  4200.,  4700.,  5600.,    &
     &        6200.,  8200., 11500., 14600., 16700., 20000., 22300.,    &
     &       24600., 27500., 30000., 31900., 33000., 33800., 34500.,    &
     &       35300., 36500., 40000., 43300., 57600.  /

!
!........................................!
      end module module_radsw_parameters !
!========================================!

