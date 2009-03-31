!!!!!  ==========================================================  !!!!!
!!!!!              rrtm2 radiation package description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the rrtm2 package includes these parts:                           !
!                                                                      !
!       'radlw_rrtm2_param.f'                                          !
!       'radlw_rrtm2_datatb.f'                                         !
!       'radlw_rrtm2_main.f'                                           !
!                                                                      !
!    the 'radlw_rrtm2_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_rrtm2_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_refprof'     -- reference t and p profiles (mls) !
!       'module_radlw_avplank'     -- plank flux data                  !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16     !
!                                     bands, where nn = 01-16          !
!                                                                      !
!    the 'radlw_rrtm2_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main rrtm2 lw radiation routine                 !
!       'rlwinit'   -- to initialize rrtm2 lw radiation                !
!                                                                      !
!    all the lw radiation subprograms become contained subprograms     !
!    in module 'module_radlw_rrtm' and many of them are not directly   !
!    accessable from places outside the module.                        !
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_rrtm2_param.f'                                          !
!       'radlw_rrtm2_datatb.f'                                         !
!       'radlw_rrtm2_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!                                                                      !
!                                                                      !
!    the original program declarations:                                !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! Copyright 2002, 2003, Atmospheric & Environmental Research, Inc.(AER)!
! This software may be used, copied, or redistributed as long as it is !
! not sold and this copyright notice is reproduced on each copy made.  !
! This model is provided as is without any express or implied warranties
!                      (http://www.rtweb.aer.com/)                     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                               rrtm                                   !
!                                                                      !
!                   rapid radiative transfer model                     !
!                                                                      !
!            atmospheric and environmental research, inc.              !
!                        840 memorial drive                            !
!                        cambridge, ma 02139                           !
!                                                                      !
!                           eli j. mlawer                              !
!                         steven j. taubman~                           !
!                         shepard a. clough                            !
!                                                                      !
!                         ~currently at gfdl                           !
!                                                                      !
!                       email:  mlawer@aer.com                         !
!                                                                      !
!        the authors wish to acknowledge the contributions of the      !
!        following people:  patrick d. brown, michael j. iacono,       !
!        ronald e. farren, luke chen, robert bergstrom.                !
!                                                                      !
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
        integer :: ilwrate, iaerlw, irgaslw, icfclw, iflagliq, iflagice
        integer :: numangs

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

        parameter ( irgaslw=1 )     !===> ... control flag for rare gases (ch4,n2o,o2, etc.)
                                    ! =0: do not include rare gases
                        !(default)  ! =1: include all rare gases

        parameter ( icfclw=1  )     !===> ... control flag for halocarbon (cfc) gases
                                    ! =0: do not include cfc gases
                        !(default)  ! =1: include all cfc gases

        parameter ( iflagliq=3 )    !===> ... liq-cloud optical properties contrl flag
                                    ! =0: input cloud opt depth, ignor iflagice setting
                                    ! =1: input cwp,cip, (ccm2 method) ignor iflagice setting
                                    ! =2: input cwp rew, ccm3 method for liquid clouds
                        !(default)  ! =3: input cwp rew, hu and stamnes(1993) method for liq cld

        parameter ( iflagice=1 )    !===> ... ice-cloud optical properties contrl flag
                                    !         only used when iflagliq .ge. 2, else is ignored
                                    ! =0: input cip rei, ccm3 method for ice clouds
                        !(default)  ! =1: input cip rei, ebert and curry(1997) for ice clouds
                                    ! =2: input cip rei, streamer (1996) for ice clouds

        parameter ( numangs=3 )     !===> ... number of angles for radiance calculation
                                    ! =0: diffusivity approx, cos(ang) = 1/1.66
                        !(default=3)! =n: n=1-4, cos of the angles have the standard x-axis
                                    !     values in first-moment gaussian quadrature

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
        real (kind=kind_phys) :: upfxc         ! level up flux for total sky
        real (kind=kind_phys) :: dnfxc         ! level dn flux for total sky
        real (kind=kind_phys) :: upfx0         ! level up flux for clear sky
        real (kind=kind_phys) :: dnfx0         ! level dn flux for clear sky
      end type
!
!  ---  parameter constants for lw band structures
!
      integer, parameter :: NBANDS = 16         ! num of total spectral bands
      integer, parameter :: NGMX   = 16         ! max num of g-points per band
      integer, parameter :: NTBL   = 10000      ! num of elem in trans table
      integer, parameter :: MAXGAS = 7          ! max num of absorbing gases
      integer, parameter :: MAXXSEC= 4          ! num of halocarbon gases
      integer, parameter :: MAXANG = 4          ! max num of angles
      integer, parameter :: NPLNK  = 181        ! dim for plank function table

      integer, parameter :: NBDLW  = NBANDS

!  ---  number of g-point in each band
      integer  :: NG01, NG02, NG03, NG04, NG05, NG06, NG07, NG08,       &
     &            NG09, NG10, NG11, NG12, NG13, NG14, NG15, NG16
      parameter (NG01=16, NG02=16, NG03=16, NG04=16, NG05=16, NG06=16,  &
     &           NG07=16, NG08=16, NG09=16, NG10=16, NG11=16, NG12=16,  &
     &           NG13=16, NG14=16, NG15=16, NG16=16)

!  ---  begining index of each band
      integer  :: NS01, NS02, NS03, NS04, NS05, NS06, NS07, NS08,       &
     &            NS09, NS10, NS11, NS12, NS13, NS14, NS15, NS16
      parameter (NS01=00, NS02=16, NS03=32, NS04=48, NS05=64, NS06=80,  &
     &           NS07=96, NS08=112, NS09=128, NS10=144, NS11=160,       &
     &           NS12=176, NS13=192, NS14=208, NS15=224, NS16=240)

!  ---  band spectrum structures (wavenumber in cm**-1)
      real (kind=kind_phys) :: wvnlw1(NBANDS), wvnlw2(NBANDS)
      data wvnlw1  /                                                    &
     &         10.,  351.,  501.,  631.,  701.,  821.,  981., 1081.,    &
     &       1181., 1391., 1481., 1801., 2081., 2251., 2381., 2601. /
      data wvnlw2  /                                                    &
     &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
     &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /


!........................................!
      end module module_radlw_parameters !
!========================================!
