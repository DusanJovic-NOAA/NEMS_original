!!!!!  ==========================================================  !!!!!
!!!!!            lw-nasa1 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the lw-nasa1 package includes these parts:                         !
!                                                                      !
!      'radlw_nasa1_param.f'                                           !
!      'radlw_nasa1_datatb.f'                                          !
!      'radlw_nasa1_main.f'                                            !
!                                                                      !
!   the 'radlw_nasa1_param.f' contains:                                !
!                                                                      !
!      'module_radlw_parameters'  -- band parameters set up            !
!      'module_radlw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radlw_nasa1_datatb.f' contains:                               !
!                                                                      !
!      'module_radlw_tables'      -- pre-computed transmittance tables !
!                                                                      !
!   the 'radlw_nasa1_main.f' contains:                                 !
!                                                                      !
!      'module_radlw_main'        -- main lw radiation transfer        !
!not!  'module_radlw_opaer'       -- compute aerosol lw opt properties !
!                                                                      !
!   in the main module 'module_radlw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'lwrad'      -- main nasa1 lw radiation routine                 !
!      'rlwinit'    -- initialization routine                          !
!                                                                      !
!   all the lw radiation subprograms become contained subprograms      !
!   in module 'module_radlw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radlw_nasa1_param.f'                                           !
!      'radlw_nasa1_datatb.f'                                          !
!      'radlw_nasa1_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use lw modules     !
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
        integer :: ilwrate, iaerlw, itrace, iflagc, iovcst, itrans

!
!  ---  set up control parameters for lw radiation
!
        parameter ( ilwrate=2 )     !===> ... flag for heating rate unit
                                    ! =1: output in k/day
                        ! (default) ! =2: output in k/second
        parameter ( iaerlw=1 )      !===> ... flag for aerosols
                                    ! =0: without aerosol effect
                        ! (default) ! =1: aeros opt prop are calc for each spectral band
                                    ! =2: broad band aeros opt prop are used for all bands

        parameter ( itrace=1 )      !===> ... flag for trace gases: n2o, ch4, cfcs
                                    ! =0: without trace gases absorption
                        ! (default) ! =1: include trace gases absorption

        parameter ( iflagc=2 )      !===> ... flag for cloud optical property
                                    ! =1: input cloud optical depth
                        ! (default) ! =2: input cloud cwp/cip/crp...

        parameter ( iovcst=2 )      !===> ... flag for cloud fraction assuption
                                    ! =1: overcast assumption, clouds either 0 or 1
                        ! (default) ! =2: allowing for fractional clouds

        parameter ( itrans=2 )      !===> ... flag for transm funct for h2o,co2,o3
                                    ! =1: k-distr method, faster, less accurate
                        ! (default) ! =2: table-lookup, slower, more accurate

!
!........................................!
      end module module_radlw_cntr_para  !
!========================================!



!========================================!
      module module_radlw_parameters     !
!........................................!

      use machine,                 only : kind_phys

      implicit   none
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
        real (kind=kind_phys) :: upfxc         ! total sky level upward flux
        real (kind=kind_phys) :: dnfxc         ! total sky level downward flux
        real (kind=kind_phys) :: upfx0         ! clear sky level upward flux
        real (kind=kind_phys) :: dnfx0         ! clear sky level downward flux
      end type
!
!  ---  parameter constants for lw band structures
!
      integer, parameter :: NBW   = 9   ! num of water vapor absorption bands
      integer, parameter :: NK0   = 6   ! num of k-dist exponentials for h2o and co2
      integer, parameter :: NK1   = 4   ! num of k-dist exponentials for n2o and ch4
      integer, parameter :: NK2   = 3   ! num of k-dist exponentials for h2o cont.
      integer, parameter :: NB3   = 3   ! num of sub bands in band-3 spectrum

      integer, parameter :: NBDLW = 10  ! total num of lw spectral bands

!  ---  parameters defining the size of the pre-computed tables for
!       transmittance using table look-up.

      integer, parameter :: NX = 26     ! num of intervals in pressure
      integer, parameter :: NC = 30     ! num of intervals in co2 amount
      integer, parameter :: NO = 21     ! num of intervals in o3 amount
      integer, parameter :: NH = 31     ! num of intervals in h2o amount

!  ---  band spectrum structures (boundaries of wavenumber in cm**-1)
      real (kind=kind_phys) :: wvnlw1(NBDLW),  wvnlw2(NBDLW)
      data wvnlw1  /   1.00,  341.0,  541.0,  801.0,  981.0,            &
     &               1101.0, 1216.0, 1381.0, 1901.0,  540.0  /
      data wvnlw2  /  340.0,  540.0,  800.0,  980.0, 1100.0,            &
     &               1215.0, 1380.0, 1900.0, 3000.0,  620.0  /

!
!........................................!
      end module module_radlw_parameters !
!========================================!
