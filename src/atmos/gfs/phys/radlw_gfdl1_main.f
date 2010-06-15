!!!!!  ==========================================================  !!!!!
!!!!!              gfdl1 radiation package description             !!!!!
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
!       'module_radlw_banddata'    -- data for each lw spectral band   !
!       'module_radlw_cntdata'     -- data for continum coeff          !
!       'module_radlw_levdata'     -- std press level data fof lbl mdl !
!                                                                      !
!    the 'radlw_gfdl1_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main gfdl1 lw radiation routine                 !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                     !
!            clouds,iovr,aerosols,sfemis,                              !
!            NPTS, NLAY, NLP1, iflip, lprnt,                           !
!          outputs:                                                    !
!            hlwc,topflx,sfcflx,                                       !
!!         optional outputs:                                           !
!            HLW0,HLWB,FLXPRF)                                         !
!                                                                      !
!       'rlwinit'   -- initialization routine                          !
!          inputs:                                                     !
!           ( icwp, me, NLAY )                                         !
!          outputs:                                                    !
!           (none)                                                     !
!                                                                      !
!    all the lw radiation subprograms become contained subprograms     !
!    in module 'module_radlw_main' and many of them are not directly   !
!    accessable from places outside the module.                        !
!                                                                      !
!                                                                      !
!    derived data type constructs used:                                !
!                                                                      !
!     1. radiation flux at toa: (from module 'module_radlw_parameters')!
!          topflw_type   -  derived data type for toa rad fluxes       !
!            upfxc              total sky upward flux at toa           !
!            upfx0              clear sky upward flux at toa           !
!                                                                      !
!     2. radiation flux at sfc: (from module 'module_radlw_parameters')!
!          sfcflw_type   -  derived data type for sfc rad fluxes       !
!            upfxc              total sky upward flux at sfc           !
!            dnfxc              total sky downward flux at sfc         !
!            dnfx0              clear sky downward flux at sfc         !
!                                                                      !
!     3. radiation flux profiles(from module 'module_radlw_parameters')!
!          proflw_type    -  derived data type for rad vertical prof   !
!            netfc              level net lw flux for total sky        !
!            netf0              level net lw flux for clear sky        !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'                                               !
!       'module physcons'                                              !
!       'module_iounitdef'                                             !
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_gfdl1_param.f'                                          !
!       'radlw_gfdl1_datatb.f'                                         !
!       'radlw_gfdl1_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!    the original program description:                                 !
!                                                                      !
!     author:                                                          !
!       m.d.schwarzkopf                                                !
!       p.o.box 308                                                    !
!       princeton, nj 08542                                            !
!       (609) 452-6521                                                 !
!                                                                      !
!     statement:                                                       !
!       this code is distributed to users with the explicit permission !
!     of the author.  no additional distribution should be made without!
!     the expressed permission of the author.  improvements may be     !
!     be undertaken, but all such changes should be brought to the     !
!     attention of the author.                                         !
!                                                                      !
!     revised: 9/15/99                                                 !
!     certified:  radiation version 1.0                                !
!                                                                      !
!     general limitation:                                              !
!       this radiation code is degigned for terrestrial applications   !
!     and is not suitable for differing atmospheric simulations.       !
!     specifically:                                                    !
!                                                                      !
!      1) the temperature must lie between 100k and 370k (without code !
!         modifications).                                              !
!      2) the maximum pressure should be less than ~1.1 atm, while the !
!         minimum pressure should be at least 10E-6 atm.               !
!      3) h2o, o3, and co2 amounts should be in terrestrial ranges.    !
!         for co2, this means a variation from ~0.5 to ~4 times        !
!         terrestrial value.                                           !
!                                                                      !
!     comments:                                                        !
!       this code is generally designed to be independent of the       !
!     vertical level structure of the user.  However, there are a      !
!     number of places which are level-dependent.                      !
!                                                                      !
!     references:                                                      !
!      1) schwarzkopf, m.d., and s.b. fels, "the simplified            !
!         exchange method revisited: an accurate, rapid method for     !
!         computation of infrared cooling rates and fluxes," journal   !
!         of geophysical research, 96 (1981), 9075-9096.               !
!      2) schwarzkopf, m.d., and s.b. fels, "improvements to the       !
!         algorithm for computing co2 transmissivities and cooling     !
!         rates," journal geophysical research, 90 (1985) 10541-10550. !
!      3) fels, s.b., "simple strategies for inclusion of voigt        !
!         effects in infrared cooling calculations," application       !
!         optics, 18 (1979), 2634-2637.                                !
!      4) fels, s.b., and m.d. schwarzkopf, "the simplified exchange   !
!         approximation: a new method for radiative transfer           !
!         calculations," journal atmospheric science, 32 (1975),       !
!         1475-1488.                                                   !
!                                                                      !
!                                                                      !
!                                                                      !
!    ncep modifications history log:                                   !
!       sep 2000,  ken campana  -- received original code from gfdl    !
!       jan 2001,  ken campana                                         !
!                  successful running with column data input           !
!       aug 2005,  yu-tai hou                                          !
!                  recoded to fit in ncep unified radiation package    !
!       mar 2007,  yu-tai hou                                          !
!                  add aerosol effect for lw radiation                 !
!       apr 2007,  yu-tai hou                                          !
!                  add spectral band heating as optional output        !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radlw_main           !
!........................................!
!
      use machine,          only : kind_phys, kind_io4
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_rgas, con_p0
      use module_iounitdef, only : NIOFRAD

      use module_radlw_cntr_para
      use module_radlw_parameters
      use module_radlw_banddata
      use module_radlw_cntdata
      use module_radlw_levdata
!     use radiation_diag_mod,     only: radiag   ! -- diag output prog
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGLW='GFDL-LW v99.1   Mar 2007'
      character(24), parameter :: VTAGLW='GFDL-LW v99.1   Apr 2007'

!  ---  constant parameters
      real (kind=kind_phys), parameter :: f_zero = 0.0

!     integer, parameter :: NIOFRAD = 16     ! unit num for input and output
                                             ! gas transm func tables
!  ---  convert mks to cgs unit system
      real (kind=kind_phys), parameter :: grav    = 1.0e+2 * con_g
      real (kind=kind_phys), parameter :: pstd    = 1.0e+1 * con_p0
      real (kind=kind_phys), parameter :: rgas    = 1.0e+7 * con_rgas

      real (kind=kind_phys), parameter :: diffac  = 1.66
      real (kind=kind_phys), parameter :: rh2oair = con_amw/con_amd
      real (kind=kind_phys), parameter :: secday  = 86400.0

!     the values of the molecular weights of f11 and f12 are derived
!   from elemental atomic weights adopted by the International Union of
!   Pure and Applied Chemistry in 1961. These values are also used in
!   the US Standard Atmosphere, 1976.

      real (kind=kind_phys), parameter ::  wtmf11  = 1.3736855E+02
      real (kind=kind_phys), parameter ::  wtmf12  = 1.2091395E+02
      real (kind=kind_phys), parameter ::  wtmf113 = 1.873765E+02
      real (kind=kind_phys), parameter ::  wtmf22  = 8.646892E+01

!  ---  module variables
!  (1)  the following arrays are co2 transm funcs, temp and press
!       derivatives for the 560-800 cm-1 band, at std temp t0.
!
!     for p(surface)=1013.25mb:
!       co251   =  transmition functions for t0 (standard profile)
!       cdt51   =  first temperature derivative of co251.
!       c2d51   =  second temperature derivative of co251.
!       co2m51  =  transm funcs for t0 for adjacent press lvls, no press
!                  quadrature.  used for nearby layer comps.
!       cdtm51  =  first temperature derivative of co2m51.
!       c2dm51  =  second temperature derivative of co2m51.
!     for p(surface)=810mb:
!       co258   =  transmition functions for t0 (standard profile)
!       cdt58   =  first temperature derivative of co258.
!       c2d58   =  second temperature derivative of co258.
!       co2m58  =  transm funcs for t0 for adjacent press lvls, no press
!                  quadrature.  used for nearby layer comps.
!       cdtm51  =  first temperature derivative of co2m51.
!       cdtm58  =  first temperature derivative of co2m58.
!       c2dm51  =  second temperature derivative of co2m51.
!       c2dm58  =  second temperature derivative of co2m58.

      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       co251, cdt51, c2d51, co258, cdt58, c2d58
      real (kind=kind_io4 ), allocatable, dimension(:)   ::             &
     &       co2m51, cdtm51, c2dm51, co2m58, cdtm58, c2dm58

!  (2)  the following arrays are co2 transm funcs and temp and press
!       derivatives for (NBCO215) narrow bands in the 15um co2 band.
!
!     for p(surface)=1013.25mb:
!        co215nbps1    = transm funcs for USSTD profile
!        co2dt15nbps1  = temperature derivative of co215nbps1.
!        co2d2t15nbps1 = second temperature derivative of co215nbps1.
!     for p(surface)=810mb:
!        co215nbps8    = transm funcs for USSTD profile
!        co2dt15nbps8  = temperature derivative of co215nbps8.
!        co2d2t15nbps8 = second temperature derivative of co215nbps8.

      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       co215nbps1, co215nbps8, co2dt15nbps1, co2dt15nbps8,        &
     &       co2d2t15nbps1, co2d2t15nbps8

!  (3)  the following arrays are ch4 and n2o transmission functions for
!       the 1200-1400 cm-1 band.
!
!     for p(surface)=1013.25mb:
!       ch451   =  ch4 transmission functions for t0 (standard profile)
!       ch4dt51 =  first temperature derivative of ch4 transm funcs
!       ch4d2t51=  second temperature derivative of ch4 transm funcs
!       n2o51   =  n2o transmission functions for t0 (standard profile)
!       n2odt51 =  first temperature derivative of n2o transm funcs
!       n2od2t51=  second temperature derivative of n2o transm funcs
!     for p(surface)=810mb:
!       ch458   =  ch4 transmission functions for t0 (standard profile)
!       ch4dt58 =  first temperature derivative of ch4 transm funcs
!       ch4d2t58=  second temperature derivative of ch4 transm funcs
!       n2o58   =  n2o transmission functions for t0 (standard profile)
!       n2odt58 =  first temperature derivative of n2o transm funcs
!       n2od2t58=  second temperature derivative of n2o transm funcs

      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       ch451, ch4dt51, ch4d2t51, ch458, ch4dt58, ch4d2t58
      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       n2o51, n2odt51, n2od2t51, n2o58, n2odt58, n2od2t58

!  (4)  the following arrays are n2o transm funcs for the 560-630 cm-1 band.
!
!     for p(surface)=1013.25mb:
!       n2o71   =  n2o transmission functions for t0 (standard profile)
!       n2odt71 =  first temperature derivative of n2o transm funcs
!       n2od2t71=  second temperature derivative of n2o transm funcs
!     for p(surface)=810mb:
!       n2o78   =  n2o transmission functions for t0 (standard profile)
!       n2odt78 =  first temperature derivative of n2o transm funcs
!       n2od2t78=  second temperature derivative of n2o transm funcs

      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       n2o71, n2odt71, n2od2t71, n2o78, n2odt78, n2od2t78

!  (5)  the following arrays are n2o transm funcs for 1070-1200 cm-1 band.
!
!     for p(surface)=1013.25mb:
!       n2o91   =  n2o transmission functions for t0 (standard profile)
!       n2odt91 =  first temperature derivative of n2o transm funcs
!       n2od2t91=  second temperature derivative of n2o transm funcs
!     for p(surface)=810mb:
!       n2o98   =  n2o transmission functions for t0 (standard profile)
!       n2odt98 =  first temperature derivative of n2o transm funcs
!       n2od2t98=  second temperature derivative of n2o transm funcs

      real (kind=kind_io4 ), allocatable, dimension(:,:) ::             &
     &       n2o91, n2odt91, n2od2t91, n2o98, n2odt98, n2od2t98

!  (6)  variables for optical_path, setup in 'rlwinit'
!   fvj     = foreign-broadened ckd 2.1 coeff (including all corrections),
!             averaged over 7 specified wide freq bands in the 560-1200
!             cm-1 range. The average is weighted by the freq of the
!             individual 10 cm-1 bands used in the averaging process.
!   fvjinw  = band-averaged foreign coeff (as in fvj) over the 900-990,
!             1070-1200 cm-1 range.
!   fvjwd   = band-averaged foreign coeff (as in fvj) over 560-800 cm-1 range.
!   svj     = self-broadened ckd 2.1 coeff (including all corrections),
!             averaged over 7 specified wide freq bands in 560-1200 cm-1
!             range. The average is weighted by the freq of the individual
!             10 cm-1 bands used in the averaging process.
!   svjinw  = band-averaged self coeff (as in svj) over 900-990,1070-1200 cm-1.
!   svjwd   = band-averaged self coeff (as in svj) over 560-800 cm-1 range.
!   radfnbd = the radiation function (radfn) averaged over each of the 7
!             freq bands: assumed to be altitude-independent
!   radfnbdinw= same as radfnbd, but for the 560-800 cm-1 range.
!   radfnbdwd = same as radfnbd, but for the 900-990,1070-1200 cm-1 range

      real (kind=kind_phys) :: fvj(7), svj(7), radfnbd(7), fvjinw,      &
     &       svjinw, radfnbdinw, fvjwd,  svjwd, radfnbdwd, vvj(40),     &
     &       ab15wd

!  (7)  the following are table variables, values set in 'table'

      type (tab1_type) :: tab1, tab1a, tab2, tab2a, tab3, tab3a, tab1w
      type (tab3_type) :: tabsr

      type (axis_type) :: temp_1 = axis_type(1, 100.0, 370.0, 10.0),    &
     &                    mass_1 = axis_type(1, -16.0,   1.9,  0.1)

!  (8)  the following are band boundaries, and corresponding cloud
!       properties indexes

      integer, dimension(NBANDS) :: band_no_str, band_no_end
      data band_no_str   /  1, 25, 41, 42, 43, 44, 46, 47 /
      data band_no_end   / 24, 40, 41, 42, 43, 45, 46, 47 /

      integer, dimension(NBLY-1) :: cld_indx_table
      data  cld_indx_table  / 40*1, 2, 2, 2, 3, 4, 5, 6 /

      integer, dimension(NLWCLDB) :: cld_band_no, cld_indx

      real (kind=kind_phys), dimension(NLWCLDB,NBFL) :: fulwwts

!  (9)  index limits variables and conversion constant set in 'rlwinit'
      real (kind=kind_phys) :: radcon

      integer :: IXPRNLTE, IXPRKMINH2O

!  (10) standard layer-mean temperature and pressure coeff set in 'ptz'
      real (kind=kind_phys), allocatable, dimension(:) :: stemp, gtemp

      logical, save :: lfirst=.true.

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

!  ---  procedure interfaces

      interface looktab
        module procedure  looktab_type1, looktab_type3
      end interface

      interface table_alloc
        module procedure    table1_alloc, table3_alloc
      end interface

!  ---  public accessable interfaces

      public   lwrad, rlwinit


! ================
      contains
! ================


!-----------------------------------
      subroutine lwrad                                                  &
!...................................
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,iovr,aerosols,sfemis,                               &
     &       NPTS, NLAY, NLP1, iflip, lprnt,                            &
!  ---  outputs:
     &       hlwc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HLW0,HLWB,FLXPRF                                           &
     &     )
!-----------------------------------------------------------------------
!
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!     plyr   (NPTS,NLAY)    - layer pressures (mb)                      !
!     plvl   (NPTS,NLP1)    - interface pressures (mb)                  !
!     tlyr   (NPTS,NLAY)    - layer temperature (k)                     !
!     tlvl   (NPTS,NLP1)    - interface temperatures (k)                !
!     qlyr   (NPTS,NLAY)    - layer h2o mixing ratio (gm/gm)*see inside !
!     olyr   (NPTS,NLAY)    - layer o3 mixing ratio (gm/gm) *see inside !
!     gasvmr (NPTS,NLAY,:)  - atmospheric gases amount:                 !
!                       (check module_radiation_gases for definition)   !
!       gasvmr(:,:,1)   -      co2 volume mixing ratio                  !
!       gasvmr(:,:,2)   -      n2o volume mixing ratio                  !
!       gasvmr(:,:,3)   -      ch4 volume mixing ratio                  !
!       gasvmr(:,:,4)   -      o2  volume mixing ratio     (not used)   !
!       gasvmr(:,:,5)   -      co  volume mixing ratio     (not used)   !
!       gasvmr(:,:,6)   -      cfc11 volume mixing ratio                !
!       gasvmr(:,:,7)   -      cfc12 volume mixing ratio                !
!       gasvmr(:,:,8)   -      cfc22 volume mixing ratio                !
!       gasvmr(:,:,9)   -      cfccl4 volume mixing ratio  (not used)   !
!       gasvmr(:,:,10)  -      cfc113 volume mixing ratio               !
!     clouds (NPTS,NLAY,:)  - cloud profiles                            !
!                       (check module_radiation_clouds for definition)  !
!                ---  for  iflgcld > 0  ---                             !
!        clouds(:,:,1)  -  layer total cloud fraction                   !
!        clouds(:,:,2)  -  layer cloud liq water path      (g/m**2)     !
!        clouds(:,:,3)  -  mean eff radius for liq cloud   (micron)     !
!        clouds(:,:,4)  -  layer cloud ice water path      (g/m**2)     !
!        clouds(:,:,5)  -  mean eff radius for ice cloud   (micron)     !
!        clouds(:,:,6)  -  layer rain drop water path      (g/m**2)     !
!        clouds(:,:,7)  -  mean eff radius for rain drop   (micron)     !
!        clouds(:,:,8)  -  layer snow flake water path     (g/m**2)     !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!        clouds(:,:,9)  -  mean eff radius for snow flake  (micron)     !
!                ---  for  iflgcld = 0  ---                             !
!        clouds(:,:,1)  -  layer total cloud fraction                   !
!        clouds(:,:,2)  -  layer cloud optical depth                    !
!        clouds(:,:,3)  -  layer cloud single scattering albedo         !
!        clouds(:,:,4)  -  layer cloud asymmetry factor                 !
!     iovr                  - control flag for cloud overlapping        !
!                             =0: random overlapping clouds             !
!                             =1: max/ran overlapping clouds            !
!     aerosols(NPTS,NLAY,NBDLW,:) - aerosol optical prop.               !
!                       (check module_radiation_aerosols for definition)!
!        (:,:,:,1)      - optical depth                                 !
!        (:,:,:,2)      - single scattering albedo                      !
!        (:,:,:,3)      - asymmetry parameter                           !
!     sfemis (NPTS)         - surface emissivity         **not used!!** !
!     NPTS                  - total number of horizontal points         !
!     NLAY,NLP1             - total number of vertical layers, levels   !
!     iflip                 - control flag for in/out vertical index    !
!                             =0: index from toa to surface             !
!                             =1: index from surface to toa             !
!     lprnt                 - cntl flag for diagnostic print out        !
!                                                                       !
!  control parameters in module "module_radlw_cntr_para":               !
!     ilwrate               - heating rate unit selections              !
!                             =1: output in k/day                       !
!                             =2: output in k/second                    !
!     iaerlw                - control flag for aerosols                 !
!                             =0: do not include aerosol effect         !
!                             >0: include aerosol effect                !
!     iflgcld               - cloud optical properties control flag     !
!                             =0: input cloud optical depth             !
!                             =1: input cld condensates:cwp,cip,etc...  !
!     ithkcld               - control flag for thick cloud adjustment   !
!                             =0: no peseudo-conv adj for maxi ovlp cld !
!                             =1: do peseudo-conv adj for maxi ovlp cld !
!     iflgco2               - control flag for co2 gas effect           !
!                             =1: logarithmic press interp co2 coeff    !
!                             =2: linear press interp co2 coeff         !
!     ico2tfs               - control flag for co2 transmission funcs   !
!                             =1: calc and write out co2 transm funcs   !
!                             =2: read in pre calc co2 transm funcs     !
!     iflgcfc               - control flag for cfc gases effect         !
!                             =0: do not include cfc gases              !
!                             =1: include 4 types of cfc effect         !
!     ifch4n2o              - control flag for ch4 and n2o gases        !
!                             =0: do not include ch4 and n2o gases      !
!                             =1: ch4,n2o gases without lbl temp interp !
!                             =2: ch4,n2o gases with lbl temp interp    !
!     ich4n2otf             - control flag for ch4 transmission func    !
!                             =1: calc and write out ch4 and n2o transm !
!                             =2: read in ch4 and n2o transm funcs      !
!                                                                       !
!  output variables:                                                    !
!     hlwc   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)   !
!     topflx (NPTS)         - radiation fluxes at top, component:       !
!                       (check module_radlw_paramters for definition)   !
!        upfxc                 total sky upward flux at top (w/m2)      !
!        upfx0                 clear sky upward flux at top (w/m2)      !
!     sfcflx (NPTS)         - radiation fluxes at sfc, component:       !
!                       (check module_radlw_paramters for definition)   !
!        upfxc                 total sky upward flux at sfc (w/m2)      !
!        dnfxc                 total sky downward flux at sfc (w/m2)    !
!        dnfx0                 clear sky downward flux at sfc (w/m2)    !
!                                                                       !
!! optional output variables:                                           !
!     hlwb(NPTS,NLAY,NBDLW) - spectral band total sky heating rates     !
!     hlw0   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)   !
!     flxprf (NPTS,NLP1)    - level radiative fluxes (w/m2), components !
!                       (check module_radlw_paramters for definition)   !
!        netfc                 total sky net lw flux                    !
!        netf0                 clear sky net lw flux                    !
!                                                                       !
!                                                                       !
!                                                                       !
!  subroutine lwrad is called by : grrad - ncep radiation driver        !
!                                                                       !
!  subroutines called by lwrad : ptz, gases_stdtf, cloudrad_driver,     !
!                                locate_in_table, lwrad1, (radiag)      !
!                                                                       !
!  *******************************************************************  !
!

      implicit none

!  ---  inputs:
      integer, intent(in)   :: NPTS, NLAY, NLP1, iovr, iflip

      logical,  intent(in) :: lprnt

      real (kind=kind_phys), dimension(:,:),  intent(in) :: plvl, tlvl, &
     &       plyr, tlyr, qlyr, olyr
      real (kind=kind_phys), dimension(:,:,:),intent(in) :: gasvmr,     &
     &       clouds
      real (kind=kind_phys), dimension(:,:,:,:),intent(in) :: aerosols
      real (kind=kind_phys), dimension(:),    intent(in) :: sfemis

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: hlwc

      type (topflw_type),    dimension(:),   intent(out) :: topflx
      type (sfcflw_type),    dimension(:),   intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys),dimension(:,:,:),optional,intent(out):: hlwb
      real (kind=kind_phys),dimension(:,:),optional,intent(out):: hlw0
      type (proflw_type),   dimension(:,:),optional,intent(out):: flxprf

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY,NLWCLDB) ::            &
     &       emmxolw, emrndlw
      real (kind=kind_phys), dimension(NPTS,NLP1,NLWCLDB) ::            &
     &       fluxn, fluxncf
      real (kind=kind_phys), dimension(NPTS,NLAY,NBLY)    :: exctsn
      real (kind=kind_phys), dimension(NPTS,NLP1,NBDLW)   :: tauaertot

      real (kind=kind_phys), dimension(NPTS,NLP1) :: temp, press, pflux,&
     &       tflux, deltaz, tpl1, tpl2, dte1, dte2, flxnet, flxnetcf
      real (kind=kind_phys), dimension(NPTS,NLAY) :: qo3, rh2o, clds,   &
     &       cwt, cic, rwt, ric, cra, csn, rra, rsn, cmxolw,  crndlw,   &
     &       ctau, pdflux, pdfinv, excts, cts, ctsco2, ctso3, heatra,   &
     &       heatracf
      real (kind=kind_phys), dimension(NPTS,NBLY) :: fctsg

      real (kind=kind_phys), dimension(NLP1) :: pd, pd8, plm, plm8
      real (kind=kind_phys), dimension(NPTS) :: flx1e1, flx1e1f, gxcts

      real (kind=kind_phys) :: pss, pss8, prnlte, sigl(NLAY), rco2, tmp
      real (kind=kind_phys) :: rch4, rn2o, cfc11, cfc12, cfc22, cfc113

      integer, dimension(NPTS,NLP1) :: ixoe1, ixoe2
      integer, dimension(NPTS)      :: nmxolw,  nrndlw
      integer :: i, j, k, k1, prkminh2o, idxflx

!
!===> ...  begin here
!

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

!  ---  the internal array is always from top to surface
!       gfdl lw use cgs unit, need a factor of 1.0e3 for pressure
!       h2o mixing ratio should greater than 2.0e-7, and temp is
!       between 100K and 370K, the limits of the tables.

      if (iflip == 0) then        ! input from toa to sfc

        do k = 1, NLAY
          do i = 1, NPTS
            press(i,k) = plyr(i,k) * 1.0e+3
!           press(i,k) = plyr(i,k)
            temp (i,k) = max( min(tlyr(i,k), 370.0), 100.0 )
!test use
!           rh2o (i,k) = max(2.0e-7, qlyr(i,k))                 ! input mass mixing ratio
!ncep model use
            rh2o (i,k) = max(2.0e-7, qlyr(i,k)/(1.0-qlyr(i,k))) ! input specific humidity
            qo3  (i,k) = max(0.0, olyr(i,k))                    ! input mass mixing ratio
!           qo3  (i,k) = max(0.0, olyr(i,k)*ramdo3)             ! input vol mixing ratio
          enddo
        enddo

        do i = 1, NPTS
          press(i,NLP1) = plvl(i,NLP1) * 1.0e+3
!         press(i,NLP1) = plvl(i,NLP1)
          temp (i,NLP1) = max( min(tlvl(i,NLP1), 370.0), 100.0 )
        enddo

!  ---  compute pressure, temperature, altitude arrays for layer
!       boundaries that are used in radiation routines and in interface
!       routines between radiation and other physics (clouds, aerosols)

        do k = 1, NLP1
          do i = 1, NPTS
            pflux(i,k) = plvl(i,k) * 1.0e+3
!           pflux(i,k) = plvl(i,k)
            tflux(i,k) = max( min(tlvl(i,k), 370.0), 100.0 )
          enddo
        enddo

        if ( iflgcld > 0 ) then    ! input cloud condensates
          do k = 1, NLAY
            do i = 1, NPTS
              clds(i,k) = clouds(i,k,1)
              cwt (i,k) = clouds(i,k,2)
              rwt (i,k) = clouds(i,k,3)
              cic (i,k) = clouds(i,k,4)
              ric (i,k) = clouds(i,k,5)
              cra (i,k) = clouds(i,k,6)
              rra (i,k) = clouds(i,k,7)
              csn (i,k) = clouds(i,k,8)
              rsn (i,k) = clouds(i,k,9)
            enddo
          enddo
        else                       ! input cloud optical depth
          do k = 2, NLAY
            do i = 1, NPTS
              clds(i,k) = clouds(i,k,1)
              ctau(i,k) = clouds(i,k,2)
            enddo
          enddo
        endif                    ! end if_iflgcld

        if (iaerlw > 0) then
          do j = 1, NBDLW
            tauaertot(:,1,j) = f_zero

            do k = 1, NLAY
              do i = 1, NPTS
                tmp = 1.66*aerosols(i,k,j,1)*(1.0-aerosols(i,k,j,2))
                tauaertot(i,k+1,j) = tauaertot(i,k,j) + tmp
              enddo
            enddo
          enddo
        else
          tauaertot(:,:,:) = f_zero
        endif

      else                        ! input data from sfc to toa

        do k = 1, NLAY
          k1 = NLP1 - k
          do i = 1, NPTS
            press(i,k) = plyr(i,k1) * 1.0e+3
!           press(i,k) = plyr(i,k1)
            temp (i,k) = max( min(tlyr(i,k1), 370.0), 100.0 )
!test use
!           rh2o (i,k) = max(2.0e-7, qlyr(i,k1))                  ! input mass mixing ratio
!ncep model use
            rh2o (i,k) = max(2.0e-7, qlyr(i,k1)/(1.0-qlyr(i,k1))) ! input specific humidity
            qo3  (i,k) = max(0.0, olyr(i,k1))                     ! input mass mixing ratio
!           qo3  (i,k) = max(0.0, olyr(i,k1)*ramdo3)              ! input vol mixing ratio
          enddo
        enddo

        do i = 1, NPTS
          press(i,NLP1) = plvl(i,1) * 1.0e+3
!         press(i,NLP1) = plvl(i,1)
          temp (i,NLP1) = max( min(tlvl(i,1), 370.0), 100.0 )
        enddo

        do k = 1, NLP1
          k1 = NLP1 - k + 1
          do i = 1, NPTS
            pflux(i,k) = plvl(i,k1) * 1.0e+3
!           pflux(i,k) = plvl(i,k1)
            tflux(i,k) = max( min(tlvl(i,k1), 370.0), 100.0 )
          enddo
        enddo

        if ( iflgcld > 0 ) then   ! input cloud condensates
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              clds(i,k) = clouds(i,k1,1)
              cwt (i,k) = clouds(i,k1,2)
              rwt (i,k) = clouds(i,k1,3)
              cic (i,k) = clouds(i,k1,4)
              ric (i,k) = clouds(i,k1,5)
              cra (i,k) = clouds(i,k1,6)
              rra (i,k) = clouds(i,k1,7)
              csn (i,k) = clouds(i,k1,8)
              rsn (i,k) = clouds(i,k1,9)
            enddo
          enddo
        else                      ! input cloud optical depth
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              clds(i,k) = clouds(i,k1,1)
              ctau(i,k) = clouds(i,k1,2)
            enddo
          enddo
        endif                    ! end if_iflgcld

        if (iaerlw > 0) then
          do j = 1, NBDLW
            tauaertot(:,1,j) = f_zero

            do k = 1, NLAY
              k1 = NLP1 - k
              do i = 1, NPTS
                tmp = 1.66*aerosols(i,k1,j,1)*(1.0-aerosols(i,k1,j,2))
                tauaertot(i,k+1,j) = tauaertot(i,k,j) + tmp
              enddo
            enddo
          enddo
        else
          tauaertot(:,:,:) = f_zero
        endif

      endif                       ! end if_iflip

!  ---  gases volumn mixing ratio are fixed global values

      k1 = NLP1 / 2
      rco2  = gasvmr(1,k1,1)
      rn2o  = gasvmr(1,k1,2)
      rch4  = gasvmr(1,k1,3)
      cfc11 = gasvmr(1,k1,6)
      cfc12 = gasvmr(1,k1,7)
      cfc22 = gasvmr(1,k1,8)
      cfc113= gasvmr(1,k1,10)

!  ---  define deltaz in meters

      deltaz(:,1) = 2.0e-2 * rgas * temp(:,1) / (grav * con_amd)
      do k = 2, NLAY
        deltaz(:,k) = alog( pflux(:,k+1) / pflux(:,k) )                 &
     &              * 1.0e-2 * rgas * temp(:,k) / (grav * con_amd)
      enddo

!  ---  convert cloud condensate paths (g/m**2) to concentrations (g/m**3)

      clds(:,1) = f_zero      !  ---  do not allow clouds at model layer 1

      do k = 1, NLAY
        do i = 1, NPTS
          cwt(i,k) = cwt(i,k) / deltaz(i,k)
          cic(i,k) = cic(i,k) / deltaz(i,k)
          cra(i,k) = cra(i,k) / deltaz(i,k)
          csn(i,k) = csn(i,k) / deltaz(i,k)
        enddo
      enddo

!  ---  compute appropriate vertical level structure and input standard
!       pressures and temperatures for that structure

      if ( lfirst ) then

!  ---  define "sigma interface from input pressure
        if ( press(1,NLP1) >= 800.0 ) then
          do k = 1, NLAY
            sigl(k) = press(1,k) / press(1,NLP1)
          enddo
        else
          k1 = (NPTS + 2) / 2
          do k = 1, NLAY
            sigl(k) = press(k1,k) / press(k1,NLP1)
          enddo
        endif

!  ---  define layer-mean pressures (the factor 1.e-2 change pa to mb)
!       from "sigp_std"

        if ( lprnt ) print *,' sigl =',sigl

        pss  = con_p0 * 1.0e-2
        pss8 = con_p0 * 1.0e-2 * 0.8

        do k = 1, NLAY
          pd (k) = sigl(k) * pss
          pd8(k) = sigl(k) * pss8
        enddo
        pd (NLP1) = pss
        pd8(NLP1) = pss8

!  ---  define layer interface as the mean of adjacent layer mean pressure
!       this is necessary co2 interpolation

        plm (1) = f_zero
        plm8(1) = f_zero
        do k = 1, NLAY
          plm (k+1) = 0.5 * (pd (k) + pd (k+1))
          plm8(k+1) = 0.5 * (pd8(k) + pd8(k+1))
        enddo
        plm (NLP1) = pss
        plm8(NLP1) = pss8

!  ---  convert pressure specification for bottom (flux) pressure
!       level for nlte calculation into an index (IXPRNLTE)
!       from "co2_source_init"

        prnlte = 0.1
        IXPRNLTE = 1
        do k = 2, NLAY
          if ((plm(k) - prnlte) < f_zero) then
            IXPRNLTE = k
          else
            exit
          endif
        enddo

!  ---  call ptz to compute standard temps and a pressure coeff(gtemp)
!       used in the radiation algorithm.

        call ptz                                                        &
!  ---  inputs:
     &     ( pd, plm, NLAY, NLP1 )
!  ---  outputs: (to module variables)

!  ---  convert pressure specification for top (flux) pressure level
!       for nearby layer calculation into an index (IXPRKMINH2O)
!
!       note: min value of IXPRKMINH2O is 1 . (but plm(1) is zero,
!       so min value of KSRAD is at least 2). if all levels used
!       for radiative calcs are at pressures less than 28 mb, then
!       ignore nearby layer effects, so IXPRKMINH2O is set to NLP1.

        prkminh2o = 28.0         ! press(mb) above which h2o-co2 overlap
                                 ! affects nearby layer transmissivities
        if ( plm(1) >= prkminh2o ) then
          IXPRKMINH2O = 1
        else if ( plm(NLAY) < prkminh2o ) then
          IXPRKMINH2O = NLP1
        else
          do k = 2, NLAY
            if ( plm(k)-prkminh2o >= f_zero ) then
              IXPRKMINH2O = k
              exit
            endif
          enddo
        endif

        if ( ich4n2otf == 1 .or. ico2tfs == 1 ) then

          call gases_stdtf                                              &
!  ---  inputs:
     &     ( rch4, rn2o, rco2, pd, plm, pd8, plm8, NLAY, NLP1 )
!  ---  outputs: ( none, to module variables )

        endif

        if ( ifch4n2o > 0 ) then
          idxflx = 7
        else
          idxflx = 6
        endif

        lfirst = .false.
      endif    ! end if_lfirst_block

!  ---  call cloudrad_driver to obtain cloud radiative properties from
!       specified or computed microphysical cloud properties

      call cloudrad_driver                                              &
!  ---  inputs:
     &     ( clds,cwt,rwt,cic,ric,cra,rra,csn,rsn,deltaz,ctau,          &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emmxolw,emrndlw, cmxolw,crndlw, nmxolw,nrndlw              &
     &     )

!  ---  compute difference between flux level pressures (pdflux)

      do k = 1, NLAY
        do i = 1, NPTS
          pdflux(i,k) = pflux(i,k+1) - pflux(i,k)
          pdfinv(i,k) = 1.0 / pdflux(i,k)
        enddo
      enddo

!  ---  compute mean temperature in the "nearby layer" between a flux
!       level and the first data level below the flux level (tpl1) or the
!       first data level above the flux level (tpl2)

      do i = 1, NPTS
        tpl1(i,1     ) = temp(i,NLAY)
        tpl1(i,2:NLAY) = tflux(i,2:NLAY)
        tpl1(i,NLP1  ) = 0.5 * (tflux(i,NLP1) + temp(i,NLAY))
        tpl2(i,2:NLAY) = tflux(i,2:NLAY)
        tpl2(i,NLP1  ) = 0.5 * (tflux(i,NLAY) + temp(i,NLAY))
      enddo

      call locate_in_table                                              &
!  ---  inputs:
     &     ( temp_1, temp, 1, NLP1,                                     &
!  ---  outputs:
     &       dte1, ixoe1                                                &
     &     )

      call locate_in_table                                              &
!  ---  inputs:
     &     ( temp_1, tflux, 1, NLP1,                                    &
!  ---  outputs:
     &       dte2, ixoe2                                                &
     &     )

      do k = 1, NLAY
        do i = 1, NPTS
          ixoe2(i,k) = ixoe2(i,k+1)
          dte2 (i,k) = dte2 (i,k+1)
        enddo
      enddo

      do i = 1, NPTS
        ixoe2(i,NLP1) = ixoe1(i,NLAY)
        dte2 (i,NLP1) = dte1 (i,NLAY)
      enddo

!  ---  compute longwave radiation.

      call lwrad1                                                       &
!  ---  inputs:
     &     ( emmxolw,emrndlw,cmxolw,crndlw,nmxolw,nrndlw,               &
     &       temp,press,tflux,pflux,pdflux,pdfinv,                      &
     &       tpl1,tpl2,dte1,dte2,ixoe1,ixoe2,tauaertot,                 &
     &       rco2,qo3,rh2o,cfc11,cfc12,cfc22,cfc113,deltaz,             &
     &       NPTS, NLAY, NLP1, IDXFLX,                                  &
!  ---  outputs:
     &       flx1e1, flx1e1f, heatra, heatracf, flxnet, flxnetcf,       &
     &       fluxn, fluxncf, exctsn, fctsg, excts, gxcts,               &
     &       cts, ctsco2, ctso3                                         &
     &     )

!  ---  call radiag to compute radiation diagnostics at desired points

!     call radiag                                                       &
!  ---  inputs:
!    &     ( NBLY, NBLW, IOFFSET, radcon,                               &
!    &       NPTS, NLAY, NLP1, ifch4n2o, IDXFLX,                        &
!    &       bdlocm, bdhicm, iband, bdbond,                             &
!    &       pdfinv, pflux, temp, press, rh2o, qo3,                     &
!    &       cfc11, cfc12, cfc22, cfc113, rch4, rn2o, rco2,             &
!    &       flx1e1, flx1e1f, heatra, heatracf, flxnet, flxnetcf,       &
!    &       fluxn, fluxncf, exctsn, fctsg, excts, gxcts,               &
!    &       cts, ctsco2, ctso3                                         &
!  ---  outputs: ( none )
!    &     )

!  ---  output total-sky and clear-sky fluxes and heating rates.
!       convert flux to w/m**2

      do i = 1, NPTS
        topflx(i)%upfxc = flxnet  (i,1) * 1.0e-3
        topflx(i)%upfx0 = flxnetcf(i,1) * 1.0e-3

        sfcflx(i)%upfxc = 5.673e-8*tflux(i,NLP1)**4
        sfcflx(i)%dnfxc = sfcflx(i)%upfxc - flxnet  (i,NLP1)*1.0e-3
        sfcflx(i)%dnfx0 = sfcflx(i)%upfxc - flxnetcf(i,NLP1)*1.0e-3
      enddo

      if (iflip == 0) then        ! output from toa to sfc

        do k = 1, NLAY
          do i = 1, NPTS
            hlwc(i,k) = heatra(i,k)
          enddo
        enddo

!! ---  optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
          do i = 1, NPTS
            hlw0(i,k) = heatracf(i,k)
          enddo
          enddo
        endif

!! ---  optional spectral heating rate
!!      *** currently only one broad band
        if ( lhlwb ) then
          do k = 1, NLAY
          do i = 1, NPTS
            hlwb(i,k,1) = heatra(i,k)
          enddo
          enddo
        endif

!! ---  optional fluxes
        if ( lflxprf ) then
          do k = 1, NLP1
          do i = 1, NPTS
            flxprf(i,k)%netfc = flxnet  (i,k) * 1.0e-3
            flxprf(i,k)%netf0 = flxnetcf(i,k) * 1.0e-3
          enddo
          enddo
        endif

      else                        ! output from sfc to toa

        do k = 1, NLAY
          k1 = NLP1 - k
          do i = 1, NPTS
            hlwc(i,k) = heatra(i,k1)
          enddo
        enddo

!! ---  optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              hlw0(i,k) = heatracf(i,k1)
            enddo
          enddo
        endif

!! ---  optional spectral heating rate
!!      *** currently only one broad band
        if ( lhlwb ) then
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              hlwb(i,k,1) = heatra(i,k1)
            enddo
          enddo
        endif

!! ---  optional fluxes
        if ( lflxprf ) then
          do k = 1, NLP1
            k1 = NLP1 - k + 1
            do i = 1, NPTS
              flxprf(i,k)%netfc = flxnet  (i,k1) * 1.0e-3
              flxprf(i,k)%netf0 = flxnetcf(i,k1) * 1.0e-3
            enddo
          enddo
        endif

      endif                       ! if_iflip

!
      return
!...................................
      end subroutine lwrad
!-----------------------------------


!-----------------------------------
      subroutine rlwinit                                                &
!...................................
!  ---  inputs:
     &     ( icwp, me, NLAY )
!  ---  outputs:( none )

! --------------------------------------------------------------------- !
!                                                                       !
!   subroutine rlwinit                                                  !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
!                                                                       !
!  inputs:                                                              !
!    icwp     -  flag of cloud schemes used by model                    !
!                =0: diagnostic scheme gives cloud tau, omiga, and g    !
!                =1: prognostic scheme gives cloud liq/ice path, etc.   !
!    me       - print control for parallel process                      !
!    NLAY     - number of vertical layers                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!                                                                       !
!  subroutines called by rlwinit :  table_alloc, lwtable                !
!                                                                       !
!                                                                       !
! --------------------------------------------------------------------- !
!
      implicit none

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: v1sh2o_296 = 5.0
      real (kind=kind_phys), parameter :: dvsh2o_296 = 10.0

!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: ( none )

!  ---  locals:

!  ---  parameters for Fu and Liou lw snow water parameterization
!       wavenumber ranges  for Fu-Liou ice crystal parameterizations
!       these apply to the ice crystal (El), cloud rain (Furainlw)
!       and cloud snow (Fusnowlw) parameterizations.
!       note: the cloud liquid drop parameterization (Cliqlw) is
!       frequency-independent,
!
!       endfubands = high wavenumber boundary of wavenumber bands used
!             in Fu-Liou parameterization. the order of increasing
!             wavenumber, thus bands 1-12 correspond to Fu's
!             bands 18-7.
!       iendfubands = index of model 10 cm-1 bands corresponding to
!              the value of endfubands. computed in the code.
!
!     data endfubands/   280, 400, 540, 670, 800, 980,                  &
!                       1100,1250,1400,1700,1900,2200 /

      integer, dimension(NBFL) :: iendfubands
      data iendfubands /  28,  40,  54,  67,  80,  98,                  &
     &                   110, 125, 140, 170, 190, 220 /

      integer, dimension(NLWCLDB+1) :: istartcldband, iendcldband,      &
     &                                 nivl1lwicecld, nivl2lwicecld

      real (kind=kind_phys), dimension(NLWCLDB+1,NBFL):: planckivlicecld
      real (kind=kind_phys), dimension(NLWCLDB)       :: planckcldband

!  ---  note: the cloud properties for wavelengths beyond 1400
!       wavenumbers are included in the results for the first band,
!       ie, that band actually is 0-560, 1400-2200 cm-1. thus the
!       indices include an extra band.

      data istartcldband /  1,  57,  81,  91, 100, 108, 121, 141 /
      data iendcldband   / 56,  80,  90,  99, 107, 120, 140, 220 /

      integer :: i, j, n, nb, NLP1

      real (kind=kind_phys) :: del, sum, a1, a2, b1, b2
      real (kind=kind_phys) :: src1nb(NBLW), vjtab(300)

!
!====> ... begin here
!
      NLP1 = NLAY + 1

      if ( me == 0 ) then
        print *,' Using GFDL Longwave Radiation, Version: ', VTAGLW

        if ( iaerlw > 0 ) then
          print *,'   --- Using input aerosol parameters for LW'
        else
          print *,'   --- Aerosol effect is NOT included in LW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        endif

        if ( iflgco2 == 1 ) then
          print *,'   --- Use logarithmic press interpolation ',        &
     &            'for co2 coefficients'
        else if ( iflgco2 == 2 ) then
          print *,'   --- Use Rlinear press interpolation ',            &
     &            'for co2 coefficients'
        endif

        if ( ico2tfs == 1 ) then
          print *,'   --- Compute and write out co2 transmission ',     &
     &            'function tables'
        else if ( ico2tfs == 2 ) then
          print *,'   --- Read in pre_computed co2 transmission ',      &
     &            'function tables'
        endif

        if ( ifch4n2o == 0 ) then
          print *,'   --- CH4 and N2O effects are NOT included in LW'
        else
          if ( ifch4n2o == 1 ) then
            print *,'   --- Include CH4, N2O effect in LW, but ',       &
     &              'without LBL temp interpolation'
          else if (ifch4n2o == 2) then
            print *,'   --- Include CH4, N2O effect in LW with ',       &
     &              'LBL temp interpolation'
          endif

          if ( ich4n2otf == 1 ) then
            print *,'       Compute and write out ch4 and n2o ',        &
     &              'transmission function tables'
          else if ( ich4n2otf == 2 ) then
            print *,'       Read in pre-computed ch4 and n2o ',         &
     &              'transmission function tables'
          endif
        endif

        if ( iflgcfc == 1 ) then
          print *,'   --- Include CFC gases absorptions in LW'
        else
          print *,'   --- CFC gases effect is NOT included in LW'
        endif

        print *,'   --- Use ckd2.1 h2o continuum method'

        print *,'   --- No up/down, only net flux profiles produced'
      endif

      if ( ilwrate == 1 ) then
!       radcon = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        radcon = con_g * 864.0 / con_cp            !   (in k/day)
      else
        radcon = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!  ---  allocate module arrays

      if ( .not. allocated(stemp) ) then
        allocate ( stemp(NLAY+1), gtemp(NLAY+1) )
      endif

      if ( .not. allocated(co251) ) then
        allocate (co251(NLP1,NLP1), cdt51(NLP1,NLP1), c2d51(NLP1,NLP1), &
     &            co258(NLP1,NLP1), cdt58(NLP1,NLP1), c2d58(NLP1,NLP1) )

        allocate (cdtm51(NLAY), co2m51(NLAY), c2dm51(NLAY),             &
     &            cdtm58(NLAY), co2m58(NLAY), c2dm58(NLAY) )

        allocate ( co215nbps1   (NLP1,NBCO215),                         &
     &             co215nbps8   (NLP1,NBCO215),                         &
     &             co2dt15nbps1 (NLP1,NBCO215),                         &
     &             co2dt15nbps8 (NLP1,NBCO215),                         &
     &             co2d2t15nbps1(NLP1,NBCO215),                         &
     &             co2d2t15nbps8(NLP1,NBCO215) )
      endif

      if ( ifch4n2o == 1 ) then
       if ( .not. allocated(ch451) ) then
        allocate ( ch451   (NLP1,NLP1), ch458   (NLP1,NLP1),            &
     &             ch4dt51 (NLP1,NLP1), ch4dt58 (NLP1,NLP1),            &
     &             ch4d2t51(NLP1,NLP1), ch4d2t58(NLP1,NLP1) )

        allocate ( n2o51   (NLP1,NLP1), n2o58   (NLP1,NLP1),            &
     &             n2odt51 (NLP1,NLP1), n2odt58 (NLP1,NLP1),            &
     &             n2od2t51(NLP1,NLP1), n2od2t58(NLP1,NLP1) )

        allocate ( n2o71   (NLP1,NLP1), n2o78   (NLP1,NLP1),            &
     &             n2odt71 (NLP1,NLP1), n2odt78 (NLP1,NLP1),            &
     &             n2od2t71(NLP1,NLP1), n2od2t78(NLP1,NLP1) )

        allocate ( n2o91   (NLP1,NLP1), n2o98   (NLP1,NLP1),            &
     &             n2odt91 (NLP1,NLP1), n2odt98 (NLP1,NLP1),            &
     &             n2od2t91(NLP1,NLP1), n2od2t98(NLP1,NLP1) )
       endif
      endif

!  ---  input gases transfer coeff at std pressures and temperatures

      if ( ifch4n2o > 0 .and. ich4n2otf == 2 ) then
          open(NIOFRAD, file='stdch4n2otfs', form='unformatted',        &
     &          access='sequential')
          read(NIOFRAD) ch451,ch458, ch4dt51,ch4dt58, ch4d2t51,ch4d2t58
          read(NIOFRAD) n2o51,n2o58, n2odt51,n2odt58, n2od2t51,n2od2t58
          read(NIOFRAD) n2o71,n2o78, n2odt71,n2odt78, n2od2t71,n2od2t78
          read(NIOFRAD) n2o91,n2o98, n2odt91,n2odt98, n2od2t91,n2od2t98
          close(NIOFRAD)
      endif

      if ( ico2tfs == 2 ) then
          open(NIOFRAD, file='stdco2tfs', form='unformatted',           &
     &           access='sequential')
          read(NIOFRAD) co215nbps1, co215nbps8, co2dt15nbps1,           &
     &                  co2dt15nbps8, co2d2t15nbps1, co2d2t15nbps8
          read(NIOFRAD) co251,co258, cdt51,cdt58, c2d51,c2d58,          &
     &                  co2m51,co2m58, cdtm51,cdtm58, c2dm51,c2dm58
          close(NIOFRAD)
      endif

!  ---  allocate coeff table spaces

      call table_alloc                                                  &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab1                                                       &
     &     )

      call table_alloc                                                  &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab2                                                       &
     &     )

      call table_alloc                                                  &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab3                                                       &
     &     )

      call table_alloc                                                  &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab1w                                                      &
     &     )

      call table_alloc                                                  &
!  ---  inputs:
     &     ( NTTABH2O, NBLY,                                            &
!  ---  in/outputs:
     &       tabsr                                                      &
     &     )

      if ( ifch4n2o > 0 ) then
        call table_alloc                                                &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab1a                                                      &
     &     )

        call table_alloc                                                &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab2a                                                      &
     &     )

        call table_alloc                                                &
!  ---  inputs:
     &     ( NTTABH2O, NUTABH2O,                                        &
!  ---  in/outputs:
     &       tab3a                                                      &
     &     )
      endif

!  ---  define the cloud band index to be used in the longwave 
!       transmission calculations for each cloud band.
!       from "longwave_clouds_init" and "longwave_driver_init"

      do i = 1, NLWCLDB
        cld_band_no(i) = i
        cld_indx(i) = min( i, NLWCLDB )
      end do

!  ---  compute a*b for computational frequency bands for the 15 um
!       region, as 1 band (ab15wd)
!       from "optical_path_init" and "longwave_tables_init"

      ab15wd = awide_c * bwide_c

!  ---  vvj are the frequencies for calculation of h2o coefficients.
!       by assumption, no other frequencies are used.
!       from "optical_ckd2p1_init"

      do i = 1, NBCNTN
        vvj(i) = v1sh2o_296 + dvsh2o_296*float(i+IOFFH2O-1)
      enddo

!  ---  define tables for lw radiation
!       from "longwave_tables_init"

      call lwtable
!  ---  inputs: ( none )
!  ---  outputs:( none )

!  ---  compute band-averaged coeffs for microphysics-to-radiation
!       interface for cloud species in infrared freq ranges. the
!       actual ext coeffts (emiss and other coeffs) are calc as time-
!       dependent quantities.
!       at present, the species included are:
!         1) cloud drops; 2) ice crystal; 3) cloud rain; 4) cloud snow.
!       from "cloudrad_package_init"

!  ---  calculation for ice crystals (Fu-Liou)
!       compute weighting function. according to Fu and Liou, this
!       should be the Planck function at -40C.

      do n = 1, NBLW
        del  = 10.0
        a1   = 5.0 + (n - 1)*del
        a2   = 3.7412e-5 * a1**3
        b1   = 1.4387 * a1 / 233.15
        b2   = exp( b1 )
        src1nb(n) = del * a2 / (b2 - 1.0)
      enddo

!  ---  compute summed weighting function over NLWCLDB cloud bands

      planckcldband(:) = f_zero

      do n = 1, NLWCLDB
        do i = istartcldband(n), iendcldband(n)
          planckcldband(n) = planckcldband(n) + src1nb(i)
        enddo
      enddo

!  ---  add contribution of 1400-2200 cm-1 region to first band

      do i = istartcldband(NLWCLDB+1), iendcldband(NLWCLDB+1)
        planckcldband(1) = planckcldband(1) + src1nb(i)
      enddo

      n  = 1
      nb = 1

      sum = f_zero
      planckivlicecld(:,:) = f_zero
      nivl1lwicecld(1) = 1

      lab_do_j : do j = 1, NBLW
        sum = sum + src1nb(j)
        if ( j == iendfubands(n) ) then
          planckivlicecld(nb,n) = sum
          sum = f_zero
        end if

        if ( j == iendcldband(nb) ) then
          if ( j <> iendfubands(n) ) then
            planckivlicecld(nb,n) = sum
            sum = f_zero
          end if

          nivl2lwicecld(nb) = n
          nb = nb + 1

          if ( nb <= NLWCLDB+1 ) then
            if ( j == iendfubands(n) ) then
              nivl1lwicecld(nb) = n + 1
            else
              nivl1lwicecld(nb) = n
            end if
          end if
        end if

        if ( j == iendfubands(n) ) n = n + 1
        if ( j >= iendcldband(NLWCLDB+1) ) then
          exit lab_do_j
        endif
      enddo lab_do_j

!  ---  compute planck-weighted band weights for Fu lw microphysics

      fulwwts(:,:) = f_zero

      do n = 1, NLWCLDB
        do i = nivl1lwicecld(n), nivl2lwicecld(n)
          fulwwts(n,i) = planckivlicecld(n,i) / planckcldband(n)
        enddo
      enddo

!  ---  add band (NLWCLDB+1) to band 1 weights

      do i = nivl1lwicecld(NLWCLDB+1), nivl2lwicecld(NLWCLDB+1)
        fulwwts(1,i) = planckivlicecld(NLWCLDB+1,i) / planckcldband(1)
      enddo

!  --- ...  from optical_ckd2p1_init

!  ---  vjtab is the frequency points at which tabulations occurred.

        do j = 1, 300
          vjtab(j) = 5.0 + 10.0*(j-1)
        enddo

!  ---  compute h2o coeff averaged over the broad bands used in the 560
!       -1200 cm-1 range.
!
!     1) 560-630            2) 630-700 (assuming 3 bands in 15um complex)
!     3) 700-800            4) 560-800 (1 band for entire complex)
!     5) 800-900            6) 900-990
!     7) 990-1070           8) 1070-1200
!     9) 800-900,1070-1200   (until this band is broken into 2)
!
!       we assume that, for best accuracy:
!     the quantity required is <svj> and <fvj) where angle brackets are
!     averages over frequency, s and f are self- and foreign coeff-
!     icients, including corrections, and vj is frequency (from vjtab).
!     notations for special bands attempt similarity with that
!     previously used in the radiation code.
!        we also assume that one value may be used (at all altitudes)
!     for the radiation correction term radfn, in each frequency band.
!     the values used below result from experimentation.

      svj    = f_zero
      fvj    = f_zero
      svjwd  = f_zero
      fvjwd  = f_zero
      svjinw = f_zero
      fvjinw = f_zero

!  ---  560-630 band:

      do j = 57, 63
        svj(1) = svj(1) + vjtab(j)*ssh2o_296(j)*sfac (j)/7.0
        fvj(1) = fvj(1) + vjtab(j)*sfh2o(j)    *fscal(j)/7.0
      enddo

      radfnbd(1) = 0.90

!  ---  630-700 band:

      do j = 64, 70
        svj(2) = svj(2) + vjtab(j)*ssh2o_296(j)*sfac (j)/7.0
        fvj(2) = fvj(2) + vjtab(j)*sfh2o(j)    *fscal(j)/7.0
      enddo

      radfnbd(2) = 0.92

!  ---  700-800 band:

      do j = 71, 80
        svj(3) = svj(3) + vjtab(j)*ssh2o_296(j)*sfac (j)/10.0
        fvj(3) = fvj(3) + vjtab(j)*sfh2o(j)    *fscal(j)/10.0
      enddo

      radfnbd(3) = 0.95

!  ---  800-900 band:

      do j = 81, 90
        svj(4) = svj(4) + vjtab(j)*ssh2o_296(j)*sfac (j)/10.0
        fvj(4) = fvj(4) + vjtab(j)*sfh2o(j)    *fscal(j)/10.0
      enddo

      radfnbd(4) = 0.97

!  ---  900-990 band:

      do j = 91, 99
        svj(5) = svj(5) + vjtab(j)*ssh2o_296(j)*sfac (j)/9.0
        fvj(5) = fvj(5) + vjtab(j)*sfh2o(j)    *fscal(j)/9.0
      enddo

      radfnbd(5) = 0.98

!  ---  990-1070 band:

      do j = 100, 107
        svj(6) = svj(6) + vjtab(j)*ssh2o_296(j)*sfac (j)/8.0
        fvj(6) = fvj(6) + vjtab(j)*sfh2o(j)    *fscal(j)/8.0
      enddo

      radfnbd(6) = 0.99

!  ---  1070-1200 band:

      do j = 108, 120
        svj(7) = svj(7) + vjtab(j)*ssh2o_296(j)*sfac (j)/13.0
        fvj(7) = fvj(7) + vjtab(j)*sfh2o(j)    *fscal(j)/13.0
      enddo

      radfnbd(7) = 0.992

!  ---  560-800 combined band:

      do j = 57, 80
        svjwd = svjwd + vjtab(j)*ssh2o_296(j)*sfac (j)/24.0
        fvjwd = fvjwd + vjtab(j)*sfh2o(j)    *fscal(j)/24.0
      enddo

      radfnbdwd = 0.92

!  ---  800-990,1070-1200 combined band:

      do j = 81, 99
        svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac (j)/22.0
        fvjinw = fvjinw + vjtab(j)*sfh2o(j)    *fscal(j)/32.0
      enddo

      do j = 108, 120
        svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac (j)/32.0
        fvjinw = fvjinw + vjtab(j)*sfh2o(j)    *fscal(j)/32.0
      enddo

      radfnbdinw = 0.98


      return
!...................................
      end subroutine rlwinit
!-----------------------------------



!  =========================================
!  *****   cloudrad_package section    *****
!  =========================================


!-----------------------------------
      subroutine cloudrad_driver                                        &
!...................................
!  ---  inputs:
     &     ( ccover,conc_drop,size_drop,conc_ice,size_ice,              &
     &       conc_rain,size_rain,conc_snow,size_snow,deltaz,ctau,       &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emmxolw,emrndlw, cmxolw,crndlw, nmxolw, nrndlw             &
     &     )

!---------------------------------------------------------------------
!
      implicit none

!  ---  inputs:
      integer,              intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: ccover,      &
     &       conc_drop, size_drop, conc_ice, size_ice, conc_rain,       &
     &       size_rain, conc_snow, size_snow, deltaz, ctau 

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       emmxolw, emrndlw
      real (kind=kind_phys), dimension(:,:),   intent(out) ::           &
     &       cmxolw, crndlw

      integer, dimension(:),  intent(out) :: nmxolw, nrndlw

!  ---  locals:
      integer :: i, k, n

!
!===> ...  begin here
!
!  ---  for the predicted clouds, initialize the cloud and cloud index
!       fields, assuming that clouds are not present.

      cmxolw(:,:) = f_zero
      crndlw(:,:) = f_zero
      emmxolw(:,:,:) = f_zero
      emrndlw(:,:,:) = f_zero
      nmxolw (:) = 0
      nrndlw (:) = 0

!  ---  define the number of clouds, their tops and bottoms and the
!       cloud amount present in the grid box.

      do i = 1, NPTS
        n = 0

        do k = NLAY, 2, -1

!  ---  search out the distinct clouds in a column. allow for thick clouds
!       when cloud is present at two adjacent levels.

          if ( n == 0 ) then
            if ( ccover(i,k) <> f_zero ) then
              if ( ccover(i,k-1) == f_zero ) then
                nrndlw(i) = nrndlw(i) + 1
                crndlw(i,k) = ccover(i,k)
              else
                n = 1
                cmxolw(i,k) = ccover(i,k)
              endif
            endif
          else
            if ( ccover(i,k) <> f_zero ) then
              cmxolw(i,k) = ccover(i,k)
              if ( k == 2 ) then
                nmxolw(i) = nmxolw(i) + 1
              endif

              if ( ccover(i,k-1) == f_zero ) then
                nmxolw(i) = nmxolw(i) + 1
                n = 0
              endif

            endif
          endif   ! end if_n_block

        end do   ! end do_k_loop
      end do   ! end do_i_loop

!  ---  start cloud_driver part

!     if ( lprnt ) then
!       print *,'size_drop,size_ice,size_rain,conc_rain,conc_snow=',    &
!    &           size_drop(1,1),size_ice(1,1),size_rain(1,1),           &
!    &           conc_rain(1,1),conc_snow(1,1)
!     endif

      if ( iflgcld > 0 ) then   ! input cloud condensates

        call cloud_lwpar
!  ---  inputs:  (from parent variables)
!  ---  outputs: (to   parent variables)

      else                      ! input cloud optical depth

        do n = 1, NLWCLDB
          emmxolw(:,:,n)  = 1.0 - exp( -diffac*ctau(:,:) )
        enddo
        emrndlw  = emmxolw

      endif   ! end if_iflgcld_block


! ================
      contains
! ================

!-----------------------------------
      subroutine cloud_lwpar
!...................................
!  ---  inputs: use parent variables
!  ---  outputs: to parent variables

!----------------------------------------------------------------------
!     determine the infrared cloud emissivities for specified wavenumber
!     bands from parameterizations for absorption coefficients due to
!     cloud drops, cloud ice crystals, rain and snow. conceptually one
!     could have separate concentrations and sizes for "thin" or randomly
!     overlapped and for maximally overlapped clouds. for now, there is
!     one concentration and size, therefore the two emissivities are set
!     equal.
!
!  intent in:
!    size_drop = the cloud drop effective diameter in microns
!    size_ice  = the ice crystal effective size in microns
!    size_rain = the rain drop effective diameter in microns
!    conc_drop = the cloud drop liquid water concentration in g/m**3
!    conc_ice = the ice water concentation in g/m**3
!    conc_rain = the rain drop water concentration in g/m**3
!    conc_snow = the snow concentration in g/m**3
!    deltaz    = the thickness of pressure layers (in m)
!
!  intent out:
!    emmxolw   = the infrared cloud emissivity for the (NLWCLDB)
!                bands. used for maximally overlapped clouds.
!    emrndlw   = the infrared cloud emissivity for the (NLWCLDB)
!                bands. used for randomly overlapped clouds.
!
!  intent local:
!    cldextbndrainlw = the specified values of the extinction
!                      coefficient for rain water in kilometers**(-1)
!                      over wavenumber bands used by the radiation code
!    cldssalbbndrainlw = the specified values of the single-
!                        scattering albedo for rain water
!                        over wavenumber bands used by the radiation code
!    cldasymmbndrainlw = the specified values of the asymmetry
!                        factor for rain water
!                        over wavenumber bands used by the radiation code
!    cldextbndsnowlw = the specified values of the extinction
!                      coefficient for snow water in kilometers**(-1)
!                      over wavenumber bands used by the radiation code
!    cldssalbbndsnowlw = the specified values of the single-
!                        scattering albedo for snow water
!                        over wavenumber bands used by the radiation code
!    cldasymmbndsnowlw = the specified values of the asymmetry
!                        factor for snow water
!                        over wavenumber bands used by the radiation code
!    cldextbndicelw = the specified values of the extinction
!                     coefficient for ice particles in kilometers**(-1)
!                     over wavenumber bands used by the radiation code
!    cldssalbbndicelw = the specified values of the single-
!                       scattering albedo for ice particles
!                       over wavenumber bands used by the radiation code
!    cldasymmbndicelw = the specified values of the asymmetry
!                       factor for ice particles
!                       over wavenumber bands used by the radiation code
!    cldextbnddroplw = the specified values of the extinction
!                      coefficient for cloud drops in kilometers**(-1)
!                      over wavenumber bands used by the radiation code
!----------------------------------------------------------------------

      implicit none

!  ---  inputs: (use parent variables)

!  ---  outputs: (to parent variables)

!  ---  locals:
      real (kind=kind_phys), dimension (NPTS,NLAY,NLWCLDB) ::           &
     &       cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw,     &
     &       cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw,     &
     &       cldextbndicelw,  cldssalbbndicelw,  cldasymmbndicelw,      &
     &       cldextbnddroplw, abscoeff

      integer       :: k, i, n

!
!===> ...  begin here
!
!  ---  compute extinction coefficient, single scattering coefficient
!       asymmetry parameter for rain

      call furainlw                                                     &
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw )

!  ---  compute extinction coefficient, single scattering coefficient
!       asymmetry parameter for snow

      call fusnowlw                                                     &
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw )

!  ---  compute extinction coefficient for cloud drops

      call cliqlw                                                       &
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbnddroplw )

!  ---  compute extinction coefficient, single scattering coefficient
!       asymmetry parameter for cloud ice crystals

      call el                                                           &
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw )

!  ---  compute absorption coefficient (in km-1) for each species
!       as the product of the extinction coefficient and (1 - single
!       scattering albedo).  the total absorption coefficient is the
!       sum of the species absorption coefficients. over a single
!       frequency band (see goody and yung, eq. 6.72) the emissivity
!       in a layer is defined as: 1 - T(f)  where T(f) is the flux
!       transmissivity, which may be computed as:
!         EXP(-(1.66)*(abs. coeff)*(layer thickness))
!       the factor 1.66 is the diffusivity factor (diffac).

      abscoeff = cldextbnddroplw                                        &
     &         + cldextbndicelw  * (1.0 - cldssalbbndicelw )            &
     &         + cldextbndsnowlw * (1.0 - cldssalbbndsnowlw)            &
     &         + cldextbndrainlw * (1.0 - cldssalbbndrainlw)

!  ---  1.0E-3 is conversion factor from (m) to (km).

      do n = 1, NLWCLDB
        emmxolw(:,:,n)  = 1.0                                           &
     &              - exp( -diffac*abscoeff(:,:,n)*deltaz(:,:)*1.0e-3 )
      enddo
      emrndlw  = emmxolw
!

      return
!...................................
      end subroutine cloud_lwpar
!----------------------------------- contained in cloudrad_driver


!-----------------------------------
      subroutine furainlw                                               &
!...................................
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw )
!...................................

!----------------------------------------------------------------------
!      Calculates absorption coefficient for cloud rain water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for rain water with radii between 60 um and 1.8 mm.
!      See also notes from Q. Fu (4 Sept 98)
!
!      note: the size of the rain water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!----------------------------------------------------------------------
!
!  inputs:
!    conc_rain = the rain drop water concentration in grams / meter**3
!
!  outputs:
!    cldextbndrainlw   = the specified values of the extinction
!                        coefficient for rain water in kilometers**(-1)
!                        over wavenumber bands used by the radiation code
!    cldssalbbndrainlw = the specified values of the single-scattering
!                        albedo for rain water over wavenumber bands
!                        used by the radiation code
!    cldasymmbndrainlw = the specified values of the asymmetry factor
!                        for rain water over wavenumber bands used by
!                        the radiation code
!
!  locals:
!    boundaries for spectral bands(cm**(-1)) (7-18):
!       2200,1900,1700,1400,1250,1100,980,800,670,540,400,280,0
!
!    brn  = empirical coefficients for extinction coefficient
!           parameterization (km**-1)
!    wrnf = empirical coefficients for single scattering albedo
!           parameterization
!    grn  = empirical coefficients for asymmetry parameter
!           parameterization
!    rwc0 = rain water content (g/m**3) used to obtain above
!           empirical coefficients.
!
!    note: since model wavenumber bands are in order of increasing
!    wavenumber, the Fu coefficients have been reversed (thus are
!    in order  bands 18-7)
!
!    cldextivlrain   = the specified spectral values of the extinction
!                      coefficient for rain water in kilometers**(-1)
!    cldssalbivlrain = the specified spectral values of the single-
!                      scattering albedo for rain water
!    cldasymmivlrain = the specified spectral values of the asymmetry
!                      factor for rain water
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs: (use in-scope variables)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY,NBFL) ::               &
     &       cldextivlrain, cldssalbivlrain, cldasymmivlrain

      real (kind=kind_phys), dimension(NBFL) :: brn, wrnf, grn
      data brn  / 1.6765,  1.6149,  1.5993,  1.5862,  1.5741,  1.5647,  &
     &            1.5642,  1.5600,  1.5559,  1.5512,  1.5478,  1.5454 /
      data wrnf / .55218,  .55334,  .55488,  .55169,  .53859,  .51904,  &
     &            .52321,  .52716,  .52969,  .53192,  .52884,  .53233 /
      data grn  / .90729,  .92990,  .93266,  .94218,  .96374,  .98584,  &
     &            .98156,  .97745,  .97467,  .97216,  .97663,  .97226 /

      real (kind=kind_phys) :: rwc0
      data rwc0 / 0.5 /

      integer :: i, k, n, ni

!
!===> ...  begin here
!
!  ---  Calculate extinction coefficient (km**(-1)) over wavenumber
!       bands of the Fu-Liou parameterization (not the radiation
!       code wavenumber bands) as of 4/8/99, the asymmetry parameter
!       is not used in the infrared code. therefore the calculations
!       are commented out.

      do n = 1, NBFL
        cldextivlrain  (:,:,n) = brn(n) * conc_rain(:,:) / rwc0
        cldssalbivlrain(:,:,n) = wrnf(n)
!       cldasymmivlrain(:,:,n) = grn(n)
      enddo

!  ---  use band weighting factors (computed in initialization module)
!       to derive band quantities for these quantities

      cldextbndrainlw   = f_zero
      cldssalbbndrainlw = f_zero
!     cldasymmbndrainlw = f_zero

      do n = 1, NLWCLDB
        do ni = 1, NBFL
          cldextbndrainlw  (:,:,n) = cldextbndrainlw  (:,:,n)           &
     &                       + cldextivlrain  (:,:,ni)*fulwwts(n,ni)
          cldssalbbndrainlw(:,:,n) = cldssalbbndrainlw(:,:,n)           &
     &                       + cldssalbivlrain(:,:,ni)*fulwwts(n,ni)
!         cldasymmbndrainlw(:,:,n) = cldasymmbndrainlw(:,:,n)           &
!    &                       + cldasymmivlrain(:,:,ni)*fulwwts(n,ni)
        enddo
      enddo
!

      return
!...................................
      end subroutine furainlw
!----------------------------------- contained in cloudrad_driver


!-----------------------------------
      subroutine fusnowlw                                               &
!...................................
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw )

!----------------------------------------------------------------------
!      Calculates absorption coefficient for cloud snow water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for snow water with radii between 60 um and 5.0 mm.
!      See also notes from Q. Fu (4 Sept 98)
!      note: the size of the snow water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!-----------------------------------------------------------------------
!
!  inputs:
!    conc_snow = the snow drop water concentration in grams / meter**3
!
!  outputs:
!    cldextbndsnowlw   = the specified values of the extinction
!                        coefficient for snow water in kilometers**(-1)
!                        over wavenumber bands used by the radiation code
!    cldssalbbndsnowlw = the specified values of the single-scattering
!                        albedo for snow water over wavenumber bands
!                        used by the radiation code
!    cldasymmbndsnowlw = the specified values of the asymmetry factor
!                        for snow water over wavenumber bands used by
!                        the radiation code
!
!  locals:
!    boundaries for spectral bands(cm**(-1)) (7-18):
!       2200,1900,1700,1400,1250,1100,980,800,670,540,400,280,0
!
!    brn  = empirical coefficients for extinction coefficient
!           parameterization (km**-1)
!    wrnf = empirical coefficients for single scattering albedo
!           parameterization
!    grn  = empirical coefficients for asymmetry parameter
!           parameterization
!    swc0 = snow water content (g/m**3) used to obtain above
!           empirical coefficients.
!
!    note: since model wavenumber bands are in order of increasing
!    wavenumber, the Fu coefficients have been reversed (thus are
!    in order  bands 18-7)
!
!    cldextivlsnow   = the specified spectral values of the extinction
!                      coefficient for snow water in kilometers**(-1)
!    cldssalbivlsnow = the specified spectral values of the single-
!                      scattering albedo for snow water
!    cldasymmivlsnow = the specified spectral values of the asymmetry
!                      factor for snow water
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs: (use in-scope variables)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY,NBFL)     ::           &
     &       cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow

      real (kind=kind_phys), dimension(NBFL) :: brn, wrnf, grn
      data brn  / .87477,  .85421,  .84825,  .84418,  .84286,  .84143,  &
     &            .84097,  .84058,  .84029,  .83995,  .83979,  .83967 /
      data wrnf / .55474,  .53160,  .54307,  .55258,  .54914,  .52342,  &
     &            .52446,  .52959,  .53180,  .53182,  .53017,  .53296 /
      data grn  / .93183,  .97097,  .95539,  .94213,  .94673,  .98396,  &
     &            .98274,  .97626,  .97327,  .97330,  .97559,  .97173 /

      real (kind=kind_phys) :: swc0
      data swc0 / 0.5 /

      integer :: i, k, n, ni

!
!===> ...   begin here
!
!  ---  Calculate extinction coefficient (km**(-1)) over wavenumber
!       bands of the Fu-Liou parameterization (not the radiation
!       code wavenumber bands) as of 4/8/99, the asymmetry parameter
!       is not used in the infrared code. therefore the calculations
!       are commented out.

      do n = 1, NBFL
        cldextivlsnow  (:,:,n) = brn (n) * conc_snow(:,:) / swc0
        cldssalbivlsnow(:,:,n) = wrnf(n)
      enddo

!  ---  use band weighting factors (computed in initialization module)
!       to derive band quantities for these quantities

      cldextbndsnowlw   = f_zero
      cldssalbbndsnowlw = f_zero

      do n = 1, NLWCLDB
        do ni = 1, NBFL
          cldextbndsnowlw  (:,:,n) = cldextbndsnowlw  (:,:,n)           &
     &                        + cldextivlsnow(:,:,ni)*fulwwts(n,ni)
          cldssalbbndsnowlw(:,:,n) = cldssalbbndsnowlw(:,:,n)           &
     &                        + cldssalbivlsnow(:,:,ni)*fulwwts(n,ni)
        enddo
      enddo
!

      return
!...................................
      end subroutine fusnowlw
!----------------------------------- contained in cloudrad_driver


!-----------------------------------
      subroutine cliqlw                                                 &
!...................................
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbnddroplw )

!----------------------------------------------------------------------
!     Calculates longwave absorption optical depth for liquid.
!     Follows Held et al. (J. Atmos. Sci., 1993).
!
!     Leo Donner, GFDL, 1 Feb 1999
!-----------------------------------------------------------------------
!
!  inputs:
!     conc_drop = the cloud drop concentration in grams / meter**3
!
!  outputs:
!    cldextbnddroplw = the specified values of the extinction coefficient
!                      for cloud drops in kilometers**(-1) over wavenumber
!                      bands used by the radiation code
!
!  locals:
!    alpha = frequency-independent parameter (in m**2/g) for absorption
!            due to cloud drops in the infrared. this value is given in
!            held et al, JAS, 1993.
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs: (use in-scope variables)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cldextbnddroplw

!  ---  locals:
      real (kind=kind_phys) :: alpha
      data alpha / 0.1 /

      integer :: n

!
!===> ...  begin here
!
      do n = 1, NLWCLDB
        cldextbnddroplw(:,:,n) = 1.0e3 * alpha * conc_drop(:,:)
      enddo
!

      return
!...................................
      end subroutine cliqlw
!----------------------------------- contained in cloudrad_driver


!-----------------------------------
      subroutine el                                                     &
!...................................
!  ---  inputs: (use in-scope variables)
!  ---  outputs:
     &     ( cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw )

!----------------------------------------------------------------------
!     Calculates total optical depth and scattering optical depth
!     for infrared radiation using Fu and Liou (1993,
!     JAS). To be used for crystal effective sizes from 20 to 130 um.
!     limits changed to 18.6 to 130.2 um on 2 April 1999 to
!     match shortwave limits.
!
!     Leo Donner, GFDL, 29 Aug 98
!-----------------------------------------------------------------------
!
!  intent in:
!    conc_ice = the ice crystal concentation in grams / meter**3
!    size_ice = the ice crystal effective size in microns
!
!  intent out:
!    cldextbndicelw   = the specified values of the extinction coefficient
!                       for ice particles in kilometers**(-1) over
!                       wavenumber bands used by the radiation code
!    cldssalbbndicelw = the specified values of the single-scattering
!                       albedo for ice particles over wavenumber bands
!                       used by the radiation code
!    cldasymmbndicelw = the specified values of the asymmetry factor for
!                       ice particles over wavenumber bands used by the
!                       radiation code
!
!  locals:
!    cldextivlice   = the specified spectral values of the extinction
!                     coefficient for ice particles in kilometers**(-1)
!    cldssalbivlice = the specified spectral values of the single-
!                     scattering albedo for ice particles
!    cldasymmivlice = the specified spectral values of the asymmetry
!                     factor for ice particles
!
!    a0,a1,a2 = empirical coefficients for extinction coefficient
!               parameterization
!    b0,...b3 = empirical coefficients for single scattering albedo
!               parameterization
!
!    note: since model wavenumber bands are in order of increasing
!    wavenumber, the Fu coefficients have been reversed (thus are
!    in order  bands 18-7)
!
!    NBFL = number of frequency bands in the lw snow water para-
!    meterization.corresponds to bands 7-18 in Fu and Liou (Table 2).
!    NBA, NBB, NBC = number of terms in parameterization for series
!    expansions (not counting the 0th power term) for ai, bi, ci
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs: (use in-scope variables)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY,NBFL)     ::           &
     &       cldextivlice, cldssalbivlice, cldasymmivlice

      real (kind=kind_phys), dimension(NBFL) :: a0,a1,a2, b0,b1,b2,b3

      data a0 /                                                         &
     &    -7.752e-3,-1.741e-2,-1.704e-2,-1.151e-2,-1.026e-2,-8.294e-3,  &
     &    -1.153e-2,-9.609e-3,-9.061e-3,-8.441e-3,-8.088e-3,-7.770e-3 /
      data a1 /  4.624,   5.541,   4.830,   4.182,   4.105,   3.925,    &
     &           4.109,   3.768,   3.741,   3.715,   3.717,   3.734 /
      data a2 /-42.010, -58.420,  16.270,  31.130,  16.360,   1.315,    &
     &          17.320,  34.110,  26.480,  19.480,  17.170,  11.850 /

      data b0 /  0.8079,  0.3964,  0.1028,  0.3254,  0.5207,  0.5631,   &
     &           0.2307,  0.2037,  0.3105,  0.3908,  0.3014,  0.1996 /
      data b1 / -0.7004e-2, -0.3155e-2,  0.5019e-2,  0.3434e-2,         &
     &          -0.9778e-3, -0.1434e-2,  0.3830e-2,  0.4247e-2,         &
     &           0.2603e-2,  0.1272e-2,  0.2639e-2,  0.3780e-2 /
      data b2 /  0.5209e-4,  0.6417e-4, -0.2024e-4, -0.3081e-4,         &
     &           0.3725e-5,  0.6298e-5, -0.1616e-4, -0.1810e-4,         &
     &          -0.1139e-4, -0.5564e-5, -0.1116e-4, -0.1491e-4 /
      data b3 / -0.1425e-6, -0.2979e-6,  0.0000e+0,  0.9143e-7,         &
     &           0.0000e+0,  0.0000e+0,  0.0000e+0,  0.0000e+0,         &
     &           0.0000e+0,  0.0000e+0,  0.0000e+0,  0.0000e+0 /

      integer :: n, ni

!
!===> ...  begin here
!
!  ---  Calculate extinction coefficient (km**(-1)) over wavenumber
!       bands of the Fu-Liou parameterization (not the radiation
!       code wavenumber bands)

      do n = 1, NBFL
        cldextivlice(:,:,n) = 1.0e+3 * conc_ice(:,:)                    &
     &          * (a0(n) + a1(n)/size_ice(:,:) + a2(n)/size_ice(:,:)**2)
      enddo

!  ---  Calculate single-scattering albedo and asymmetry parameter.
!       these are dimensionless. as of 4/8/99, the asymmetry parameter
!       is not used in the infrared code. therefore the calculations
!       are commented out.

      do n = 1, NBFL
        cldssalbivlice(:,:,n) = 1.0 - ( b0(n) + b1(n)*size_ice(:,:)     &
     &           + b2(n)*size_ice(:,:)**2 + b3(n)*size_ice(:,:)**3 )
      enddo

!  ---  use band weighting factors (computed in initialization module)
!       to derive band quantities for these quantities

      cldextbndicelw   = f_zero
      cldssalbbndicelw = f_zero

      do n = 1, NLWCLDB
        do ni = 1, NBFL
          cldextbndicelw  (:,:,n) = cldextbndicelw(:,:,n)               &
     &                       + cldextivlice(:,:,ni)*fulwwts(n,ni)
          cldssalbbndicelw(:,:,n) = cldssalbbndicelw(:,:,n)             &
     &                       + cldssalbivlice(:,:,ni)*fulwwts(n,ni)
        enddo
      enddo
!

      return
!...................................
      end subroutine el
!----------------------------------- contained in cloudrad_driver
!
!...................................
      end subroutine cloudrad_driver
!-----------------------------------



!  =========================================
!  *****    longwave_driver section    *****
!  =========================================


!-----------------------------------
      subroutine lwrad1                                                 &
!...................................
!  ---  inputs:
     &     ( emmxolw,emrndlw,cmxolw,crndlw,nmxolw,nrndlw,               &
     &       temp,press,tflux,pflux,pdflux,pdfinv,                      &
     &       tpl1,tpl2,dte1,dte2,ixoe1,ixoe2,tauaertot,                 &
     &       rco2,qo3,rh2o,cfc11,cfc12,cfc22,cfc113,deltaz,             &
     &       NPTS, NLAY, NLP1, IDXFLX,                                  &
!  ---  outputs: 
     &       flx1e1, flx1e1f, heatra, heatracf, flxnet, flxnetcf,       &
     &       fluxn, fluxncf, exctsn, fctsg, excts, gxcts,               &
     &       cts, ctsco2, ctso3                                         &
     &     )

!---------------------------------------------------------------------
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1991), 9075-9096.
!
!     (2) schwarzkopf, m. d., and s. b. fels, "improvements to the
!         algorithm for computing co2 transmissivities and cooling
!         rates," journal geophysical research, 90 (1985) 10541-10550.
!
!     (3) fels, s.b., "simple strategies for inclusion of voigt
!         effects in infrared cooling calculations," application
!         optics, 18 (1979), 2634-2637.
!
!     (4) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 10/7/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
!  intent in:
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from 1 to NLAY.
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from 1 to NLAY.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from 1 to NLAY. (default is one).
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from 1 to NLAY. (default is one).
!     kmxolw =  maximum number of maximally overlapped longwave clouds.
!     krndlw =  maximum number of randomly overlapped longwave clouds.
!     nmxolw  =  number of maximally overlapped longwave clouds
!                at each grid point.
!     nrndlw  =  number of maximally overlapped longwave clouds
!                at each grid point.
!     press  =  pressure at data levels of model.
!     qo3    =  mass mixing ratio of o3 at model data levels.
!     rh2o   =  mass mixing ratio of h2o at model data levels.
!     temp   =  temperature at data levels of model.
!
!  intent out:
!     flxnet  =  net longwave flux at model flux levels (including the
!                ground and the top of the atmpsphere).
!     heatra  =  heating rate at data levels.
!     pflux   =  pressure at flux levels of model.
!     tflux   =  temperature assigned to model flux levels.
!
!  intent local:
!     atmden   =  atmospheric density, in gm/cm**2, for each of the
!                 NLAY layers.
!     co21c    =  transmission function for the 560-800 cm-1 band
!                 from levels k through NLP1 to level k. used for flux
!                 at level k arising from exchange with levels k
!                 through NLP1. includes co2 (from tables), h2o (after
!                 multiplication with over).
!     co21diag =  transmission function for co2 only in the 560-800
!                 cm-1 band, from levels k to k, where k ranges from
!                 1 to NLP1.
!     co21r    =  transmission function for the 560-800 cm-1 band
!                 from level k to levels k+1 through NLP1. used for
!                 flux at levels k+1 through NLP1 arising from
!                 exchange with level k. includes co2 (from tables),
!                 h2o (after multiplication with over).
!     co2spnb  =  transmission functions from level 1 to levels 2
!                 through NLP1, for the (NBCO215) narrow bands
!                 comprising the 15um range (560-800 cm-1). used for
!                 exact CTS computations. includes co2 and possibly
!                 n2o (in the first such band).
!       to3cnt =  transmission function for the 990-1070 cm-1 band
!                 including o3(9.6 um) + h2o continuum (no lines)
!                 and possibly cfcs.
!       sorc15 =  planck function for 560-800 cm-1 bands (sum over
!                 bands 9 and 10).
!       sorc   =  planck function, at model temperatures, for all
!                 bands;  used in cool-to-space calculations.
!       tmp1   =  temporary array, used for computational purposes
!                 in various places. should have no consequences
!                 beyond the immediate module wherein it is defined.
!    tmp2,tmp3 =  temporary arrays, similar to tmp1
!         vv   =  layer-mean pressure in atmospheres. due to quadra-
!                 ture considerations, this does not equal the pressure
!                 at the data level (press).
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, IDXFLX

      real (kind=kind_phys), dimension(:,:,:), intent(in) ::            &
     &       emmxolw, emrndlw, tauaertot
      real (kind=kind_phys), dimension(:,:),   intent(in) ::            &
     &       cmxolw, crndlw, press, temp, pflux, tflux, pdflux, pdfinv, &
     &       tpl1, tpl2, dte1, dte2, qo3, rh2o, deltaz

      real (kind=kind_phys),               intent(in) :: rco2, cfc11,   &
     &       cfc12, cfc22, cfc113

      integer, dimension(:),  intent(in) :: nmxolw,  nrndlw
      integer, dimension(:,:),intent(in) :: ixoe1, ixoe2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: heatra,     &
     &       heatracf, flxnet,flxnetcf, fctsg, excts,cts, ctsco2,ctso3
      real (kind=kind_phys), dimension(:),   intent(out) :: flx1e1,     &
     &       flx1e1f, gxcts

      real (kind=kind_phys), dimension(:,:,:), intent(out) :: fluxn,    &
     &       fluxncf, exctsn

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1,NLWCLDB) ::            &
     &       clddiag, cldtf
      real (kind=kind_phys), dimension(NPTS,NLP1,NBCO215) :: co2spnb
      real (kind=kind_phys), dimension(NPTS,NLP1,NBANDS)  :: sorc

      real (kind=kind_phys), dimension(NPTS,NLP1)         :: cnttaub1,  &
     &       cnttaub2, cnttaub3, co21diag, co21c, co21r, dsorc15,       &
     &       dsorc93, dsorcb1, dsorcb2, dsorcb3, t4, dt4, emiss,        &
     &       heatem, overod, sorc15, to3cnt, tmp1, tmp2, soe2, soe3,    &
     &       soe4, soe5, e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx,         &
     &       contodb1, contodb2, contodb3, emissb, emisdg, to3dg,       &
     &       flx, tdav, tstdav, tlsqu, tmpdiff, n2o9c, tn2o17,          &
     &       empl1, empl2, totvo2, avephi, totch2obdwd, totphi,         &
     &       cntval, toto3, tphio3
      real (kind=kind_phys), dimension(NPTS,NLP1) :: heatemcf, flxcf,   &
     &       taero8, taero8kp, emissf, emissbf, emisdgf, tch4n2oe,      &
     &       tch4n2o_diag, e1cts1f, e1cts2f, e1flxf, empl1f, empl2f,    &
     &       vrpfh2o, avephif, tphfh2o, tcfc8, totf11, totf12,          &
     &       totf113, totf22

      real (kind=kind_phys), dimension(NPTS,NLAY) :: cts_sum, cts_sumcf,&
     &       co2nbl, cts_tem, var1, var2, var3, var4, xch2obdwd, rh2os, &
     &       rfrgn, tfac, wk

      real (kind=kind_phys), dimension(NPTS)   :: s1a, flxge1, a1, a2,  &
     &       emx1, emx2, flx1e1cf, flxge1cf, gxctscf, flx1e1fcf,        &
     &       flxge1fcf, flxge1f, emx1f, emx2f
      real (kind=kind_phys), dimension(NPTS,2) :: emspec, emspecf
      real (kind=kind_phys), dimension(NPTS,NLP1,3) :: contdg
      real (kind=kind_phys), dimension(NPTS,NLAY,7) :: xch2obd

      real (kind=kind_phys), dimension(NPTS,NLAY,NLWCLDB) ::            &
     &       taucld_mxolw, taucld_rndlw, taunbl_mxolw

      integer :: n, k, kp, m, j, kmxolw

!
!===> ...  begin here
!
!  ---  find the maximum number of maximally overlapped clouds for
!       longwave radiation (thickcld).

      kmxolw = maxval( nmxolw )

!  ---  initialization flux arrays

      fluxn  (:,:,:) = f_zero
      fluxncf(:,:,:) = f_zero

      exctsn(:,:,:) = f_zero
      fctsg (:,:)   = f_zero

      flxnet  (:,:) = f_zero
      flxnetcf(:,:) = f_zero

!  ---  setup optical path parameters

      call optical_path_setup
!  ---  inputs: ( from in-scope variables )
!  ---  outputs:( to  in-scope variables )

!  ---  call co2coef to compute some co2 temperature and pressure
!       interpolation quantities.

      call co2coef                                                      &
!  ---  inputs:
     &     ( press, temp, tflux, pdflux,                                &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       tdav, tstdav, tlsqu, tmpdiff, a1, a2                       &
     &     )

!  ---  call cldtau to compute cloud layer transmission functions for
!       all layers.

      call cldtau                                                       &
!  ---  inputs:
     &     ( cmxolw, crndlw, emmxolw, emrndlw, NPTS, NLAY, NLWCLDB)
!  ---  outputs: ( none )

!  ---  call transfn to compute temperature-corrected co2 transmission
!        functions (co2spnb and co2nbl). these fields remain in 
!        co2_tf_mod, from where they are accessed as needed.

      call transfn                                                      &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co2nbl, co2spnb                                            &
     &     )

!  ---  compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!       mission functions and n2o 560-670 cm-1 transmission functions
!       appropriate for level 1.

      call transcolrow                                                  &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       1, 1, 1, NLP1, 2, NLP1,                                    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, tch4n2oe, n2o9c, tn2o17                      &
     &     )

!  ---  save co2 and ch4n2o diagonal elements for use in nearby layer
!       co2 transmission function calculations.

      co21diag    (:,1)      = co21c   (:,1)
      tch4n2o_diag(:,1:NLP1) = tch4n2oe(:,1:NLP1)

!  ---  go into optical_path_mod to obtain the optical path functions
!       needed for use from level 1. the 15um band transmission
!       functions between levels 1 and 2 are stored in overod and
!       co2nbl; they will not be overwritten, as long as calculations
!       are made for pressure levels increasing from 1.

      call optical_trans_funct_from_KS                                  &
!  ---  inputs:
     &     ( tn2o17, NPTS, NLAY, NLP1,                                  &
!  ---  outputs:
     &       to3cnt, overod, cnttaub1, cnttaub2, cnttaub3               &
     &     )

!  ---  obtain combined transmission functions for 560-800 cm-1 band,
!       pertaining to level 1.

      co21r (:,2:NLP1) = overod(:,1:NLAY)*co21r(:,2:NLP1)
      co21c (:,2:NLP1) = overod(:,1:NLAY)*co21c(:,2:NLP1)

!  ---  compute cloud transmission functions between level 1 and all
!       other levels.

      call cloud                                                        &
!  ---  inputs:
     &     ( 1, cmxolw, crndlw, NPTS, NLP1, NLWCLDB,                    &
!  ---  outputs:
     &       cldtf                                                      &
     &     )

!  ---  save "nearby layer" cloud transmission function for later use.

      do n = 1,NLWCLDB
        clddiag(:,1,n) = cldtf(:,1,n)
      enddo

!  ---  calculate the source function

      call co2_source_calc                                              &
!  ---  inputs:
     &     ( press, rco2, pdfinv, dte1, ixoe1,                          &
     &       tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       sorc15, soe2, soe3, soe4, soe5, sorc                       &
     &     )

!  ---  define "source function" appropriate for emissivity calculations
!       (temp**4), source functions for selected ranges including more
!       than 1 frequency band (sorc15 for 15 um co2 band) and
!       differences in source functions (deltab) over pressure layers.
!
!       note: the values of sorc, sorc15, sorcwin, and derivatives
!       depend on the no. of freq. bands!

      t4(:,1:NLP1) = temp (:,1:NLP1)**4

      dsorc93 (:,2:NLP1) = soe2  (:,2:NLP1) - soe2  (:,1:NLAY)
      dsorcb1 (:,2:NLP1) = soe3  (:,2:NLP1) - soe3  (:,1:NLAY)
      dsorcb2 (:,2:NLP1) = soe4  (:,2:NLP1) - soe4  (:,1:NLAY)
      dsorcb3 (:,2:NLP1) = soe5  (:,2:NLP1) - soe5  (:,1:NLAY)
      dsorc15 (:,2:NLP1) = sorc15(:,2:NLP1) - sorc15(:,1:NLAY)
      dt4     (:,2:NLP1) = t4    (:,2:NLP1) -  t4   (:,1:NLAY)

!  ---  obtain cool-to-space flux at the top by integration of heating
!       rates and using cool-to-space flux at the bottom (current value
!       of gxcts).  note that the pressure quantities and conversion
!       factors have not been included either in excts or in gxcts.
!       these cancel out, thus reducing computations.

      call cool_to_space_exact                                          &
!  ---  inputs:
     &     ( cldtf, press, temp, to3cnt, pflux, pdfinv,                 &
     &       dte1, ixoe1, co2spnb, n2o9c, tn2o17, sorc,                 &
     &       tfac, rh2os, rfrgn, wk, xch2obd, totvo2, var1, var2,       &
     &       totf11, totf12, totf113, totf22, tauaertot,                &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       gxcts, gxctscf, cts_sum, cts_sumcf,                        &
     &       exctsn, fctsg, excts                                       &
     &     )

!  ---  compute the emissivity fluxes for k=1.

      if ( ifch4n2o > 0 ) then

        call e1e290                                                     &
!  ---  inputs:
     &     ( tab1, tab2, tab1a, tab2a, tab1w,                           &
     &       dte1, ixoe1, dte2, ixoe2, mass_1, avephi,                  &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss               &
!  ---  optional in/out:
     &,      avephif, e1cts1f, e1cts2f, e1flxf, emissf                  &
     &     )

      else

        call e1e290                                                     &
!  ---  inputs:
     &     ( tab1, tab2, tab1a, tab2a, tab1w,                           &
     &       dte1, ixoe1, dte2, ixoe2, mass_1, avephi,                  &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss               &
     &     )

      endif

!  ---  add the effects of other radiative gases on these flux arrays.
!       the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!       thus, tch4e and tn2oe are obtained directly from the transmission
!       functions.

      if ( ifch4n2o > 0 ) then
        e1flxf (:,2:NLP1) = e1flxf (:,2:NLP1) * tch4n2oe(:,2:NLP1)
        e1cts1f(:,1:NLP1) = e1cts1f(:,1:NLP1) * tch4n2oe(:,1:NLP1)
        e1cts2f(:,1:NLAY) = e1cts2f(:,1:NLAY) * tch4n2oe(:,2:NLP1)
        emissf (:,1:NLAY) = emissf (:,1:NLAY) * tch4n2oe(:,2:NLP1)

!  ---  add cfc transmissivities if species which absorb in this 
!       frequency range ( presently index 8) are present.

        if ( iflgcfc > 0 ) then
          tcfc8(:,1:NLP1) = 1.0                                         &
     &                     - CFC_basic%strf113(8)*totf113(:,1:NLP1)     &
     &                     - CFC_basic%strf22 (8)*totf22 (:,1:NLP1)

          e1flxf (:,2:NLP1) = e1flxf (:,2:NLP1) * tcfc8(:,2:NLP1)
          e1cts1f(:,1:NLP1) = e1cts1f(:,1:NLP1) * tcfc8(:,1:NLP1)
          e1cts2f(:,1:NLAY) = e1cts2f(:,1:NLAY) * tcfc8(:,2:NLP1)
          emissf (:,1:NLAY) = emissf (:,1:NLAY) * tcfc8(:,2:NLP1)
        endif

!  ---  compute aerosol transmission function for 1200-1400 cm-1 region

        if ( iaerlw > 0 ) then
!         call get_totaerooptdep(8, totaer_tmp)
!         taero8(:,1:NLP1) = EXP( -1.0*totaer_tmp(:,1:NLP1) )

          if ( NBDLW == 1 ) then
            taero8(:,1:NLP1) = EXP( -1.0*tauaertot(:,1:NLP1,1) )
          else                                       ! *** to be developed !!!
            taero8(:,1:NLP1) = EXP( -1.0*tauaertot(:,1:NLP1,1) )
          endif

          e1flxf (:,2:NLP1) = e1flxf (:,2:NLP1) * taero8(:,2:NLP1)
          e1cts1f(:,1:NLP1) = e1cts1f(:,1:NLP1) * taero8(:,1:NLP1)
          e1cts2f(:,1:NLAY) = e1cts2f(:,1:NLAY) * taero8(:,2:NLP1)
          emissf (:,1:NLAY) = emissf (:,1:NLAY) * taero8(:,2:NLP1)
        endif
      endif

!  ---  the following is a rewrite of the original code largely to
!       eliminate three-dimensional arrays.  the code works on the
!       following principles.  let k be a fixed flux level and kp
!       be a varying flux level, then
!
!       flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=1,NLP1.
!
!       if we assume a symmetrical array tau(k,kp)=tau(kp,k), we can
!       break down the calculations for k=1,NLP1 as follows:
!
!       flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=k+1,NLP1           (1)
!       flux(kp) =   (deltab(k )*tau(kp,k)) for kp=k+1,NLP1.          (2)
!
!       plus deltab(k)*tau(k,k) for all k.
!
!       if we compute a one-dimensional array tauod(kp) for
!       kp=k+1,NLP1, equations (1) and (2) become:
!
!       tauod(kp) = tau(kp,k)                                         (3)
!
!       flux (k ) = sum(deltab(kp)*tauod(kp)) for kp=k+1,NLP1         (4)
!       flux (kp) =    (deltab(k )*tauod(kp)) for kp=k+1,NLP1         (5)
!
!       where tau(k,k) and nearby layer terms are handled separately.
!
!       compute fluxes at level k = 1
!       compute the terms for flux at levels 2 to NLP1 from level 1.
!       compute terms for flux at level 1 from level 1.
!       compute the terms for flux at level 1 due to levels KP from 2
!       to NLP1.

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( t4, e1flx, 1, dt4, emiss, 0, cldtf, cld_indx(1),           &
     &       cld_band_no(1), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( sorc15, co21r, 1, dsorc15, co21c, 1, cldtf, cld_indx(2),   &
     &       cld_band_no(2), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( soe3, cnttaub1,0, dsorcb1, cnttaub1,0, cldtf, cld_indx(3), &
     &       cld_band_no(3), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( soe4, cnttaub2,0, dsorcb2, cnttaub2,0, cldtf, cld_indx(4), &
     &       cld_band_no(4), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( soe2, to3cnt, 0, dsorc93, to3cnt, 0, cldtf, cld_indx(5),   &
     &       cld_band_no(5), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_ks                                           &
!  ---  inputs:
     &     ( soe5, cnttaub3,0, dsorcb3, cnttaub3,0, cldtf, cld_indx(6), &
     &       cld_band_no(6), NPTS, NLP1 )
!  ---  outputs: ( none )

      if (ifch4n2o > 0) then

        call longwave_fluxes_ks                                         &
!  ---  inputs:
     &     ( t4, e1flxf(:,:),1, dt4, emissf(:,:),0, cldtf, cld_indx(7), &
     &       cld_band_no(7), NPTS, NLP1 )
!  ---  outputs: ( none )

      endif

!  ---  compute approximate cool-to-space heating rates for 1 wide band
!       in the 15um  range (560-800 cm-1) (ctsco2) and for 1 band in
!       the 9.6 um band (ctso3).

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 1, sorc15, co21r, 0, cldtf(:,:,cld_indx(2)),               &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       ctsco2                                                     &
     &     )

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 2, soe2, to3cnt, 1, cldtf(:,:,cld_indx(5)),                &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       ctso3                                                      &
     &     )

!  ---  compute the emissivity cool-to-space heating rates for the
!       frequency ranges: 160-560, 800-990, and 1070-1200 cm-1. (the
!       latter 2 are combined in calculations).

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 3, t4, e1ctw2, 1, cldtf(:,:,cld_indx(1)),                  &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       cts_tem                                                    &
!  ---  optional arguments:
     &,      e1ctw1, 0                                                  &
     &     )

      cts(:,:) = cts_tem(:,:)

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 4, soe3, cnttaub1, 1, cldtf(:,:,cld_indx(3)),              &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       cts_tem                                                    &
     &     )

      cts(:,:) = cts(:,:) + cts_tem(:,:)

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 5, soe4, cnttaub2, 1, cldtf(:,:,cld_indx(4)),              &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       cts_tem                                                    &
     &     )

      cts(:,:) = cts(:,:) + cts_tem(:,:)

      call cool_to_space_approx                                         &
!  ---  inputs:
     &     ( 6, soe5, cnttaub3, 1, cldtf(:,:,cld_indx(6)),              &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       cts_tem                                                    &
     &     )

      cts(:,:) = cts(:,:) + cts_tem(:,:)

!  ---  obtain the flux at the top of the atmosphere in the 0-160,
!       1200-2200 cm-1 frequency ranges, where heating rates and fluxes
!       are derived from h2o emissivity calculations (flx1e1) by:
!       1) obtaining the surface flux (flxge1); 2) summing the
!       emissivity flux divergence for these ranges (tmp1) over all
!       pressure layers.
!  ifdef ch4n2o
!       if the 1200-1400 cm-1 range is computed separately, flux calcu-
!       lations are done separately in this range, then combined with
!       those from the other frequency range.

      s1a   (:) = t4(:,NLP1)*(e1cts1(:,NLP1) - e1ctw1(:,NLP1))
      flxge1(:) = s1a(:)*cldtf(:,NLP1,1)
      tmp1(:,1:NLAY) = t4(:,1:NLAY)*(e1cts1(:,1:NLAY)-e1ctw1(:,1:NLAY))
      tmp2(:,1:NLAY) = t4(:,1:NLAY)*(e1cts2(:,1:NLAY)-e1ctw2(:,1:NLAY))
      flx1e1(:) = flxge1(:)

      do k = 1, NLAY
        flx1e1(:) = flx1e1(:) + tmp1(:,k)*cldtf(:,k,1)                  &
     &            - tmp2(:,k)*cldtf(:,k+1,1)
      enddo

      flxge1cf(:) = s1a(:)
      flx1e1cf(:) = flxge1cf(:)

      do k=1,NLAY
        flx1e1cf(:) = flx1e1cf(:) + tmp1(:,k) - tmp2(:,k)
      enddo

      if ( ifch4n2o > 0 ) then
        s1a(:) = t4(:,NLP1)*e1cts1f(:,NLP1)
        flxge1f(:) = s1a(:)*cldtf(:,NLP1,cld_indx(7))
        flx1e1f(:) = flxge1f(:)

        do k = 1, NLAY
          tmp1(:,k) = t4(:,k)*e1cts1f(:,k)
          tmp2(:,k) = t4(:,k)*e1cts2f(:,k)
          flx1e1f(:) = flx1e1f(:) + tmp1(:,k)*cldtf(:,k,cld_indx(7))    &
     &               - tmp2(:,k)*cldtf(:,k+1,cld_indx(7))
        enddo

        flx1e1(:) = flx1e1(:) + flx1e1f(:)

        flxge1fcf(:) = s1a(:)
        flx1e1fcf(:) = s1a(:)

        do k = 1, NLAY
          flx1e1fcf(:) = flx1e1fcf(:) + tmp1(:,k) - tmp2(:,k)
        enddo

        flx1e1cf(:) = flx1e1cf(:) + flx1e1fcf(:)
      endif

!  ---  perform flux calculations for the flux levels 2 to NLAY-1.
!       calculations for flux levels NLAY and NLP1 are done separately,
!       as all calculations are special cases or nearby layers.

      do 100 k = 2, NLAY-1

!  ---  compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!       mission functions and n2o 560-670 cm-1 transmission functions
!       appropriate for level k.

        call transcolrow                                                &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       k, k, k, NLP1, k+1, NLP1,                                  &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, tch4n2oe, n2o9c, tn2o17                      &
     &     )

!  ---  save co2 diagonal element and ch4n2o tf term for use in nearby
!       layer co2 transmission function calculations.

        co21diag(:,k) = co21c(:,k)
        tch4n2o_diag(:,k+1:NLP1) = tch4n2oe(:,k+1:NLP1)

!  ---  the 15 um band transmission functions between levels k and k+1
!       are stored in overod and co2nbl; they will not be overwritten,
!       as long as calculations are made for pressure levels increasing
!       from k.
!
        call optical_trans_funct_k_down                                 &
!  ---  inputs:
     &     ( k, cnttaub1, cnttaub2, cnttaub3, tn2o17,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       to3cnt, overod, contodb1, contodb2, contodb3               &
     &     )

!  ---  add the overod factor to then co21 transmission functions.

        co21c(:,k+1:NLP1) = overod(:,k:NLAY)*co21c(:,k+1:NLP1)
        co21r(:,k+1:NLP1) = overod(:,k:NLAY)*co21r(:,k+1:NLP1)

!  ---  compute cloud transmission functions between level k and all
!       other levels greater or equal to k.

        call cloud                                                      &
!  ---  inputs:
     &     ( k, cmxolw, crndlw, NPTS, NLP1, NLWCLDB,                    &
!  ---  outputs:
     &       cldtf                                                      &
     &     )

!  ---  save "nearby layer" cloud transmission function for later use.

        do n = 1,NLWCLDB
          clddiag(:,k,n) = cldtf(:,k,n)
        enddo

!  ---  compute the exchange terms in the flux equation (except the
!       nearby layer (k,k) terms, done later).

        if ( ifch4n2o > 0) then

          call e290                                                     &
!  ---  inputs:
     &     ( k, tab2, tab2a, dte2, ixoe2, mass_1, avephi,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emiss, emissb                                              &
!  --- optional in/out:
     &,      avephif, emissf, emissbf                                   &
     &     )

        else

          call e290                                                     &
!  ---  inputs:
     &     ( k, tab2, tab2a, dte2, ixoe2, mass_1, avephi,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emiss, emissb                                              &
     &     )

        endif

! ---  add the effects of other radiative gases on these flux arrays.
!      the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!      thus, tch4e and tn2oe are obtained directly from the transmission
!      functions.

        if ( ifch4n2o > 0 ) then

          emissbf(:,k:NLAY) = emissbf(:,k:NLAY)*tch4n2oe(:,k+1:NLP1)
          emissf(:,k:NLAY) = emissf(:,k:NLAY)*tch4n2oe(:,k+1:NLP1)

!  ---  add cfc transmissivities if species which absorb in this
!       frequency range are present.

          if ( iflgcfc > 0 ) then
            do kp = k, NLAY
              tcfc8(:,kp+1) = 1.0                                       &
     &          - CFC_basic%strf113(8)*(totf113(:,kp+1) - totf113(:,k)) &
     &          - CFC_basic%strf22 (8)*(totf22 (:,kp+1) - totf22 (:,k))
            enddo

            emissbf(:,k:NLAY) = emissbf(:,k:NLAY)*tcfc8(:,k+1:NLP1)
            emissf(:,k:NLAY) = emissf(:,k:NLAY)*tcfc8(:,k+1:NLP1)
          endif

!  ---  compute aerosol transmission function for 1200-1400 cm-1 region
!       (as quotient of 2 exponentials), taero8kp(k) contains the
!       (k+1,k) transmissivities for all k in the 1200-1400 cm-1
!       frequency range.

          if ( iaerlw > 0 ) then
            do kp = k+1, NLP1
              taero8kp(:,kp) = taero8(:,kp)/taero8(:,k)
            enddo

            emissbf(:,k:NLAY) = emissbf(:,k:NLAY)*taero8kp(:,k+1:NLP1)
            emissf(:,k:NLAY) = emissf(:,k:NLAY)*taero8kp(:,k+1:NLP1)
          endif
        endif

!  ---  compute the terms for flux at levels k+1 to NLP1 from level k.
!       compute the terms for flux at level k due to levels
!       kp from k+1 to NLP1.

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dt4, emissb, emiss, 0, cldtf(:,:,cld_indx(1)),          &
     &       cld_band_no(1), NPTS, NLP1 )
!  ---  outputs: ( none )

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dsorc15, co21r, co21c, 1, cldtf(:,:,cld_indx(2)),       &
     &       cld_band_no(2), NPTS, NLP1 )
!  ---  outputs: ( none )

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dsorcb1, contodb1, contodb1,0, cldtf(:,:,cld_indx(3)),  &
     &       cld_band_no(3), NPTS, NLP1 )
!  ---  outputs: ( none )

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dsorcb2, contodb2, contodb2,0, cldtf(:,:,cld_indx(4)),  &
     &       cld_band_no(4), NPTS, NLP1 )
!  ---  outputs: ( none )

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dsorc93, to3cnt, to3cnt, 0, cldtf(:,:,cld_indx(5)),     &
     &       cld_band_no(5), NPTS, NLP1 )
!  ---  outputs: ( none )

        call longwave_fluxes_k_down                                     &
!  ---  inputs:
     &     ( k, dsorcb3, contodb3, contodb3,0, cldtf(:,:,cld_indx(6)),  &
     &       cld_band_no(6), NPTS, NLP1 )
!  ---  outputs: ( none )

        if ( ifch4n2o > 0 ) then

          call longwave_fluxes_k_down                                   &
!  ---  inputs:
     &     ( k, dt4,emissbf(:,:),emissf(:,:),0, cldtf(:,:,cld_indx(7)), &
     &       cld_band_no(7), NPTS, NLP1 )
!  ---  outputs: ( none )

        endif

100   continue

!  ---  compute remaining flux terms. these include:
!       1) the (k,k) terms, for pressure levels k from 2 to NLAY-1
!          (the 1,1 term was handled earlier);
!       2) terms for pressure level NLAY. these include the (NLAY,NLAY) term,
!          computed as in (1), and the (NLAY,NLP1) and (NLP1,NLAY) terms,
!          computed somewhat differently from the similar terms at
!          higher levels, owing to the proximity to the surface layer
!          NLP1;
!       3) the term for pressure level NLP1 (the (NLP1,NLP1 term).
!
!       compute k=NLAY case.  since the kp loop is length one, many
!       simplifications occur.  the co2 quantities and the emissivity
!       quantities) are computed in the nbl section. therefore, we want
!       to compute over, to3cnt, and contod; according to our notation
!       over(:,NLAY), to3cnt(:,NLAY), and contod(:,NLAY).  the boundary
!       layer and nearby layer corrections to the transmission functions
!       are obtained above.  the following ratios are used in various nbl
!       nbl calculations.  the remaining calculations are for:
!
!       1) the (k,k) terms, k=2,NLAY-1;
!       2) the (NLAY,NLAY) term;
!       3) the (NLAY,NLP1) term;
!       4) the (NLP1,NLAY) term;
!       5) the (NLP1,NLP1) term.
!
!       each is uniquely handled.  different flux terms are computed
!       differently the fourth section obtains water transmission
!       functions used in q(approximate) calculations and also makes nbl
!       corrections:
!
!       1) emiss (:,:) is the transmission function matrix obtained
!          using E2spec;
!       2) "nearby layer" corrections (emiss(i,i)) are obtained
!          using E3v88;
!       3) special values at the surface (emiss(NLAY,NLP1),
!          emiss(NLP1,NLAY), emiss(NLP1,NLP1)) are calculated.
!
!       compute temperature and/or scaled amount indices and residuals
!       for nearby layer and special layer lookup tables.
!
!          calculation for special cases (NLAY,NLP1) and (NLP1,NLAY)
!
!       compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!       mission functions and n2o 560-670 cm-1 transmission functions
!       appropriate for level NLAY.

      call transcolrow                                                  &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       NLAY,NLAY,NLAY,NLP1,NLP1,NLP1,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, tch4n2oe, n2o9c, tn2o17                      &
     &     )

!  ---  get optical path terms for NLAY

      call optical_trans_funct_KE                                       &
!  ---  inputs:
     &     ( cnttaub1, cnttaub2, cnttaub3, tn2o17,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       to3cnt, overod, contodb1, contodb2, contodb3               &
     &     )

!  ---  compute cloud transmission functions between level NLAY and NLAY
!       and NLP1

      call cloud                                                        &
!  ---  inputs:
     &     ( NLAY, cmxolw, crndlw, NPTS, NLP1, NLWCLDB,                 &
!  ---  outputs:
     &       cldtf                                                      &
     &     )

!  ---  save "nearby layer" cloud transmission function for later use.
!       save ch4n2o tf term

      do n = 1,NLWCLDB
        clddiag(:,NLAY,n) = cldtf(:,NLAY,n)
      enddo

      tch4n2o_diag(:,NLP1) = tch4n2oe(:,NLP1)

!  ---  call enear to calculate emissivity arrays

      if ( ifch4n2o > 0) then

        call enear                                                      &
!  ---  inputs:
     &     ( tab2, tab3, tab2a, tab3a, tpl1, tpl2,                      &
     &       dte1, dte2, ixoe1, ixoe2, temp_1, mass_1,                  &
     &       empl1, empl2, emx1, emx2, var2,                            &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emisdg, emspec                                             &
!  ---  optional in/out:
     &,      empl1f, empl2f, emx1f, emx2f, vrpfh2o,                     &
     &       emisdgf, emspecf                                           &
     &     )

      else

        call enear                                                      &
!  ---  inputs:
     &     ( tab2, tab3, tab2a, tab3a, tpl1, tpl2,                      &
     &       dte1, dte2, ixoe1, ixoe2, temp_1, mass_1,                  &
     &       empl1, empl2, emx1, emx2, var2,                            &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emisdg, emspec                                             &
     &     )

      endif

!  ---  add the effects of other radiative gases on these flux arrays.
!       the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!       thus, tch4e and tn2oe are obtained directly from the transmission
!       functions.

      if ( ifch4n2o > 0 ) then
        do k = 1, 2
          emspecf(:,K) = emspecf(:,K)*tch4n2oe(:,NLP1)
        enddo

        emisdgf(:,2:NLP1) = emisdgf(:,2:NLP1)*tch4n2o_diag(:,2:NLP1)

!  ---  add cfc transmissivities if species which absorb in this
!       frequency range are present.

        if ( iflgcfc > 0 ) then
          tcfc8(:,NLP1) = 1.0                                           &
     &      - CFC_basic%strf113(8)*(totf113(:,NLP1) - totf113(:,NLAY))  &
     &      - CFC_basic%strf22 (8)*(totf22 (:,NLP1) - totf22 (:,NLAY))

          do k = 1, 2
            emspecf(:,k) = emspecf(:,k)*tcfc8(:,NLP1)
          enddo

          emisdgf(:,2:NLP1) = emisdgf(:,2:NLP1) * tcfc8(:,2:NLP1)
        endif
      endif

!  ---  compute nearby layer transmission functions for 15 um band, 
!       continuum bands, and 9.3 um band in subroutine Nearbylyrtf.
!       transmission functions for the special cases (NLAY,NLP1) and
!       (NLP1,NLAY) are also computed for the 15 um band.

      call trans_nearby                                                 &
!  ---  inputs:
     &     ( press, pflux, overod, pdflux, pdfinv, co2nbl,              &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, co21diag                                     &
     &     )

!  ---  obtain fluxes for the two terms (NLAY,NLP1) and (NLP1,NLAY), 
!       both using the same cloud transmission functions (from layer NLAY)

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dt4, emspec(:,2), emspec(:,1),                             &
     &       cldtf(:,:,cld_indx(1)), cld_band_no(1), NPTS )
!  ---  outputs: ( none )

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dsorc15, co21c(:,NLP1), co21r(:,NLP1),                     &
     &       cldtf(:,:,cld_indx(2)), cld_band_no(2), NPTS )
!  ---  outputs: ( none )

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dsorcb1, contodb1(:,NLAY), contodb1(:,NLAY),               &
     &       cldtf(:,:,cld_indx(3)), cld_band_no(3), NPTS )
!  ---  outputs: ( none )

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dsorcb2, contodb2(:,NLAY), contodb2(:,NLAY),               &
     &       cldtf(:,:,cld_indx(4)), cld_band_no(4), NPTS )
!  ---  outputs: ( none )

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dsorc93, to3cnt(:,NLAY), to3cnt(:,NLAY),                   &
     &       cldtf(:,:,cld_indx(5)), cld_band_no(5), NPTS )
!  ---  outputs: ( none )

      call longwave_fluxes_KE_KEp1                                      &
!  ---  inputs:
     &     ( dsorcb3, contodb3(:,NLAY), contodb3(:,NLAY),               &
     &       cldtf(:,:,cld_indx(6)), cld_band_no(6), NPTS )
!  ---  outputs: ( none )

      if (ifch4n2o > 0) then

        call longwave_fluxes_KE_KEp1                                    &
!  ---  inputs:
     &     ( dt4, emspecf(:,2), emspecf(:,1),                           &
     &       cldtf(:,:,cld_indx(7)), cld_band_no(7), NPTS )
!  ---  outputs: ( none )

      endif

!  ---  obtain optical path transmission functions for diagonal terms

      call optical_trans_funct_diag
!  ---  inputs: ( use in-scope variables )
!  ---  outputs:( use in-scope variables )

!  ---  compute cloud transmission functions between level NLP1 and NLP1

      call cloud                                                        &
!  ---  inputs:
     &     ( NLP1, cmxolw, crndlw, NPTS, NLP1, NLWCLDB,                 &
!  ---  outputs:
     &       cldtf                                                      &
     &     )

!  ---  save "nearby layer" cloud transmission function for later use.

      do n = 1,NLWCLDB
        clddiag(:,NLP1,n) = cldtf(:,NLP1,n)
      enddo

!  ---  obtain fluxes for the diagonal terms at all levels.

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dt4, emisdg, clddiag(:,:,cld_indx(1)),                     &
     &       cld_band_no(1), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dsorc15, co21diag, clddiag(:,:,cld_indx(2)),               &
     &       cld_band_no(2), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dsorcb1, contdg(:,:,1), clddiag(:,:,cld_indx(3)),          &
     &       cld_band_no(3), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dsorcb2, contdg(:,:,2), clddiag(:,:,cld_indx(4)),          &
     &       cld_band_no(4), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dsorc93, to3dg, clddiag(:,:,cld_indx(5)),                  &
     &       cld_band_no(5), NPTS, NLP1 )
!  ---  outputs: ( none )

      call longwave_fluxes_diag                                         &
!  ---  inputs:
     &     ( dsorcb3, contdg(:,:,3), clddiag(:,:,cld_indx(6)),          &
     &       cld_band_no(6), NPTS, NLP1 )
!  ---  outputs: ( none )

      if ( ifch4n2o > 0 ) then
        call longwave_fluxes_diag                                       &
!  ---  inputs:
     &     ( dt4, emisdgf(:,:), clddiag(:,:,cld_indx(7)),               &
     &       cld_band_no(7), NPTS, NLP1 )
!  ---  outputs: ( none )
      endif

!  ---  sum up fluxes over bands

      call longwave_fluxes_sum                                          &
!  ---  inputs:
     &     ( NPTS,                                                      &
!  ---  outputs:
     &       flx, flxcf                                                 &
     &     )

!  ---  compute emissivity heating rates.

      heatem  (:,1:NLAY) = radcon*(flx  (:,2:NLP1)-flx  (:,1:NLAY))     &
     &                   * pdfinv(:,1:NLAY)
      heatemcf(:,1:NLAY) = radcon*(flxcf(:,2:NLP1)-flxcf(:,1:NLAY))     &
     &                   * pdfinv(:,1:NLAY)

!  ---  compute total heating rates.

      heatra  (:,1:NLAY) = heatem  (:,1:NLAY) + cts_sum  (:,1:NLAY)
      heatracf(:,1:NLAY) = heatemcf(:,1:NLAY) + cts_sumcf(:,1:NLAY)

!  ---  compute the flux at each flux level using the flux at the
!       top (flx1e1 + gxcts) and the integral of the heating rates

      flxnet  (:,1) = flx1e1  (:) + gxcts  (:)
      flxnetcf(:,1) = flx1e1cf(:) + gxctscf(:)

      tmp1(:,1:NLAY) = heatra  (:,1:NLAY)*pdflux(:,1:NLAY)/radcon
      tmp2(:,1:NLAY) = heatracf(:,1:NLAY)*pdflux(:,1:NLAY)/radcon

      do k = 2, NLP1
        flxnet  (:,k) = flxnet  (:,k-1) + tmp1(:,k-1)
        flxnetcf(:,k) = flxnetcf(:,k-1) + tmp2(:,k-1)
      enddo

!  ---  call Thickcld to perform "pseudo-convective adjustment" for
!       maximally overlapped clouds, if desired.

      if ( ithkcld == 1 ) then

        call thickcld                                                   &
!  ---  inputs:
     &     ( pflux, cmxolw, emmxolw, kmxolw, pdfinv, NPTS, NLAY,        &
!  ---  outputs:
     &       flxnet, heatra                                             &
     &     )

      endif


! ================
      contains
! ================


!  =========================================
!  *****    longwave_clouds section    *****
!  =========================================


!-----------------------------------
      subroutine cldtau                                                 &
!...................................
!  ---  inputs:
     &     ( cmxolw, crndlw, emmxolw, emrndlw, NPTS, NLAY, NLWCLDB )
!  ---  outputs: ( none )

!--------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLWCLDB

      real (kind=kind_phys), dimension(:,:),   intent(in) ::            &
     &       cmxolw, crndlw
      real (kind=kind_phys), dimension(:,:,:), intent(in) ::            &
     &       emmxolw, emrndlw

!  ---  outputs: ( none )

!  ---  locals:
         integer :: i, k, n
!
!===> ...  begin here
!

!  --- define max overlap layer transmission function over layers 1,NLAY

      do n = 1, NLWCLDB
        do k = 1, NLAY
          taucld_mxolw(:,k,n) = 1.0 - emmxolw(:,k,n)
        enddo
      enddo

!  --- define "weighted random cloud" layer transmission function
!      over layers 1,NLAY

      do n = 1, NLWCLDB
        do k = 1, NLAY
          do i = 1, NPTS
            if ( crndlw(i,k) > f_zero ) then
              taucld_rndlw(i,k,n) = (crndlw(i,k)/(1.0 - cmxolw(i,k)))   &
     &                            * (1.0 - emrndlw(i,k,n)) + 1.0        &
     &                            - crndlw(i,k)/(1.0 - cmxolw(i,k))
            else
              taucld_rndlw(i,k,n) = 1.0
            endif
          enddo
        enddo
      enddo

!  ---  define "nearby layer" cloud transmission function for max
!       overlapped clouds (if emissivity not equal to one)

      do n = 1, NLWCLDB
        do k=1, NLAY
          taunbl_mxolw(:,k,n) = f_zero
        enddo
      enddo

!
      return
!...................................
      end subroutine cldtau
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine cloud                                                  &
!...................................
!  ---  inputs:
     &     ( kl, cmxolw, crndlw, NPTS, NLP1, NLWCLDB,                   &
!  ---  outputs: ( none )
     &       cldtf                                                      &
     &     )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in)  :: kl, NPTS, NLP1, NLWCLDB

      real (kind=kind_phys), dimension(:,:), intent(in) :: cmxolw,crndlw

!  ---  outputs: ( none )
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: cldtf

!  ---  locals:
        real (kind=kind_phys), dimension(NPTS,NLP1) :: cldtfmo, cldtfrd
        integer   ::   n, i, kp

!
!===> ...  begin here
!

!  ---  the definition of "within a max overlapped cloud" is:
!       at pressure level k (separating layers k and (k-1)), the max
!       overlap cloud amounts for layers k and (k-1) must be 1) nonzero
!       (else no such cloud) and 2) equal (else one such cloud ends at
!       level k and another begins). Another way to define this is: if
!       considering the transmission across layer kp (between levels
!       kp and (kp+1)) the max overlap cloud amounts for layers kp and
!       (kp-1) must be nonzero and equal.

      do n = 1, NLWCLDB

!  ---  cloud "nearby layer" transmission functions

        cldtfmo(:,kl) = f_zero
        cldtfrd(:,kl) = 1.0

!  ---  if level kl is within a maximum overlapped cloud, the cloud
!       "nearby layer" transmission function may be non-unity. Exception:
!       at levels 1, NLP1 the function must be unity.

        if ( kl>1 .and. kl<NLP1 ) then
          do i = 1, NPTS
            if ( cmxolw(i,kl-1) /= f_zero .and.                         &
     &           cmxolw(i,kl) == cmxolw(i,kl-1) ) then
              cldtfmo(i,kl) = cmxolw(i,kl) * taunbl_mxolw(i,kl,n)
              cldtfrd(i,kl) = 1.0 - cmxolw(i,kl)
            endif
          enddo
        endif

        cldtf(:,kl,n) = cldtfmo(:,kl) + cldtfrd(:,kl)

!  ---  cloud transmission functions between level kl and higher
!       levels ( when kl le NLP1)

        if ( kl < NLP1 ) then

          cldtfmo(:,kl) = f_zero
          cldtfrd(:,kl) = 1.0

!  --- for first layer below  level kl, assume flux at level kl
!      is unity and is apportioned between (cmxolw) max. overlap cld,
!      (crndlw) rnd overlap cld, and remainder as clear sky.

          cldtfmo(:,kl+1) = cmxolw(:,kl)*taucld_mxolw(:,kl,n)
          cldtfrd(:,kl+1) = (1.0 - cmxolw(:,kl))*taucld_rndlw(:,kl,n)
          cldtf(:,kl+1,n) = cldtfmo(:,kl+1) + cldtfrd(:,kl+1)

          do kp = kl+2, NLP1

!  ---  if layers above and below level (kp-1) have no max overlap cloud,
!       or their amounts differ (ie, either top of new max overlap cld or
!       no max overlap cld at all), then apportion total "flux" (or,
!       cloud tf (cldtf)) between any max overlap cloud in layer(kp-1),
!       any rnd overlap cloud and clear sky.

            do i = 1, NPTS
              if ( cmxolw(i,kp-2) == f_zero .or.                        &
     &             cmxolw(i,kp-2) /= cmxolw(i,kp-1) ) then
                cldtfmo(i,kp) = cldtf(i,kp-1,n) * cmxolw(i,kp-1)        &
     &                        * taucld_mxolw(i,kp-1,n)
                cldtfrd(i,kp) = cldtf(i,kp-1,n) * (1.0-cmxolw(i,kp-1))  &
     &                        * taucld_rndlw(i,kp-1,n)
                cldtf(i,kp,n) = cldtfmo(i,kp) + cldtfrd(i,kp)

!  ---  if layer above level (kp-1) has a max overlap cloud, and layer
!       layer below level (kp-1) also does (ie, within max overlap cld)
!       obtain separate cloud tfs for max overlap cloud and for
!       remainder (which may contain a random overlap cloud).

              else
                cldtfmo(i,kp) = cldtfmo(i,kp-1)*taucld_mxolw(i,kp-1,n)
                cldtfrd(i,kp) = cldtfrd(i,kp-1)*taucld_rndlw(i,kp-1,n)
                cldtf(i,kp,n) = cldtfmo(i,kp) + cldtfrd(i,kp)
              endif
            enddo

          enddo
        endif
      enddo

!
      return
!...................................
      end subroutine cloud
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine thickcld                                               &
!  ---  inputs:
     &     ( pflux, cmxolw, emmxolw, kmxolw, pdfinv, NPTS, NLAY,        &
!  ---  outputs:
     &       flxnet, heatra                                             &
     &     )

!------------------------------------------------------------------
!     this module recomputes cloud fluxes in "thick" clouds assuming
!     that df/dp is constant. the effect is to reduce top-of-cloud
!     cooling rates, thus performing a "pseudo-convective adjustment"
!     by heating (in a relative sense) the cloud top.
!
!     NOTE: this module cannot handle a frequency-dependent emissivity.
!     Therefore, it assumes that emissivity quantities (emmxolw) are
!     from frequency band 1 (normally unity).
!---------------------------------------------------------------------
!     input variables:
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from 1 to NLAY.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from 1 to NLAY. (default is one).
!     kmxolw  =  maximum number of maximally overlapped longwave clouds.
!     pdfinv  =  inverse of pressure difference between flux levels.
!     pflux   =  pressure at flux levels of model.
!
!----------------------------------------------------------------------
!     output variables:
!     flxnet  =  net flux at flux levels of model (also an input
!                variable)
!     heatra  =  heating rate at data levels. (also an input variable)
!
!----------------------------------------------------------------------


      implicit none

!  ---  inputs:
      integer, intent(in) :: kmxolw, NPTS, NLAY

      real (kind=kind_phys), dimension(:,:),   intent(in) :: pflux,     &
     &       cmxolw, pdfinv
      real (kind=kind_phys), dimension(:,:,:), intent(in) :: emmxolw

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),intent(out) :: flxnet,heatra

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: tmp1
      real (kind=kind_phys), dimension(NPTS)      :: delptc, fbtm,      &
     &       ftop, pbtm, ptop

      integer, dimension(NPTS,NLAY) :: ktopmxo, kbtmmxo
      integer, dimension(NPTS)      :: itopmxo, ibtmmxo

      integer :: i, k, kc, kc1, kc2
!
!===> ...  begin here
!

!  ---  determine levels at which max overlap clouds start and stop

      itopmxo(:) = 0
      ibtmmxo(:) = 0
      ktopmxo(:,:) = 0
      kbtmmxo(:,:) = 0

!  ---  max overlap cloud in first layer (not likely)

      do i = 1, NPTS
        if ( cmxolw(i,1) > f_zero ) then
          itopmxo(i) = itopmxo(i) + 1
          ktopmxo(i,itopmxo(i)) = 1
        endif
      enddo

!  ---  k-level for which top of max overlap cloud is defined

      do k = 2, NLAY
        do i = 1, NPTS
          if ( cmxolw(i,k) > f_zero .and.                               &
     &         cmxolw(i,k-1) /= cmxolw(i,k) ) then
            itopmxo(i) = itopmxo(i) + 1
            ktopmxo(i,itopmxo(i)) = k
          endif
        enddo
      enddo

!  ---  k-level for which bottom of max overlap cloud is defined

      do k = 1, NLAY-1
        do i = 1, NPTS
          if ( cmxolw(i,k) > f_zero .and.                               &
     &         cmxolw(i,k+1) /= cmxolw(i,k) ) then
            ibtmmxo(i) = ibtmmxo(i) + 1
            kbtmmxo(i,ibtmmxo(i)) = k+1
          endif
        enddo
      enddo

!  ---  bottom of max overlap cloud in NLAY'th level

      do i = 1, NPTS
        if ( cmxolw(i,NLAY) > f_zero ) then
          ibtmmxo(i) = ibtmmxo(i) + 1
          kbtmmxo(i,ibtmmxo(i)) = NLP1
        endif
      enddo

      if ( kmxolw /= 0 ) then

!  ---  obtain the pressures and fluxes of the top and bottom of the cloud.

        do kc = 1, kmxolw
          do i = 1, NPTS
            kc1 = ktopmxo(i,kc)
            kc2 = kbtmmxo(i,kc)

            if ( kc2 > kc1 ) then
              ptop(i) = pflux (i,kc1)
              pbtm(i) = pflux (i,kc2)
              ftop(i) = flxnet(i,kc1)
              fbtm(i) = flxnet(i,kc2)

!  ---  compute the "flux derivative" df/dp delptc.

              delptc(i) = (ftop(i) - fbtm(i)) / (ptop(i) - pbtm(i))

!  ---  compute the total flux change from the top of the cloud.

              do k = kc1+1, kc2-1
                tmp1(i,k) = (pflux(i,k) - ptop(i))*delptc(i) + ftop(i)
                flxnet(i,k) = flxnet(i,k)                               &
     &                      * (1.0 - cmxolw(i,k)*emmxolw(i,k,1))        &
     &                      + tmp1(i,k)*cmxolw(i,k)*emmxolw(i,k,1)
              enddo
            endif
          enddo
        enddo

      endif

!  ---  recompute the heating rates based on the revised fluxes.

      heatra(:,1:NLAY) = radcon*(flxnet(:,2:NLAY+1)                     &
     &                 - flxnet(:,1:NLAY))*pdfinv(:,1:NLAY)

!
      return
!...................................
      end subroutine thickcld
!----------------------------------- contained in lwrad1


!  =========================================
!  *****    longwave_fluxes section    *****
!  =========================================


!-----------------------------------
      subroutine longwave_fluxes_ks                                     &
!...................................
!  ---  inputs:
     &     ( source, trans, iof, source2, trans2, iof2,                 &
     &       cld_trans, cld_ind, m, NPTS, NLP1 )
!  ---  outputs: ( none )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: iof, iof2, m, cld_ind, NPTS, NLP1

      real (kind=kind_phys), dimension(:,:),   intent(in) :: source,    &
     &       source2, trans2, trans
      real (kind=kind_phys), dimension(:,:,:), intent(in) :: cld_trans

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: flux_tmp,flux_tmp2
      integer :: k

!
!===> ...  begin here
!
      do k = 2, NLP1
        flux_tmp (:,k) = source (:,1)*trans (:,k-1+iof)
        flux_tmp2(:,k) = source2(:,k)*trans2(:,k-1+iof2)
      enddo

      if ( m==1 .or. m>=7 ) then
        fluxn(:,1,m) = fluxn(:,1,m) + source(:,1)* trans(:,1)
      else
        fluxn(:,1,m) = fluxn(:,1,m) + source(:,1)
      endif

      do k = 2, NLP1
        fluxn(:,k,m) = fluxn(:,k,m)                                     &
     &               + flux_tmp(:,k) * cld_trans(:,k,cld_ind)
        fluxn(:,1,m) = fluxn(:,1,m)                                     &
     &               + flux_tmp2(:,k)* cld_trans(:,k,cld_ind)
      enddo

      if ( m==1 .or. m>=7 ) then
        fluxncf(:,1,m) =  source(:,1)*trans(:,1)
      else
        fluxncf(:,1,m) =  source(:,1)
      endif

      do k = 2, NLP1
        fluxncf(:,k,m) = fluxncf(:,k,m) + flux_tmp(:,k)
        fluxncf(:,1,m) = fluxncf(:,1,m) + flux_tmp2(:,k)
      enddo

!
      return
!...................................
      end subroutine longwave_fluxes_ks
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine longwave_fluxes_k_down                                 &
!...................................
!  ---  inputs:
     &     ( klevel, source, trans, trans2, iof, cld_trans, m,          &
     &       NPTS, NLP1 )
!  ---  outputs: ( none )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: iof, klevel, m, NPTS, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: cld_trans,   &
     &       source, trans2, trans

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: flux_tmp, flux_tmp2

      integer :: kp

!
!===> ...  begin here
!
      do kp = klevel+1, NLP1
        flux_tmp (:,kp) = source(:,klevel)*trans (:,kp-1+iof)
        flux_tmp2(:,kp) = source(:,kp)    *trans2(:,kp-1+iof)
      enddo

      do kp = klevel+1, NLP1
        fluxn(:,kp,m) = fluxn(:,kp,m) + flux_tmp(:,kp)*cld_trans(:,kp)
        fluxn(:,klevel,m) = fluxn(:,klevel,m)                           &
     &                    + flux_tmp2(:,kp)*cld_trans(:,kp)
      enddo

      do kp = klevel+1, NLP1
        fluxncf(:,kp,    m) = fluxncf(:,kp,    m) + flux_tmp (:,kp)
        fluxncf(:,klevel,m) = fluxncf(:,klevel,m) + flux_tmp2(:,kp)
      enddo

!
      return
!...................................
      end subroutine longwave_fluxes_k_down
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine longwave_fluxes_KE_KEp1                                &
!...................................
!  ---  inputs:
     &     ( source, trans, trans2, ctrans, m, NPTS )
!  ---  outputs: ( none )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: m, NPTS

      real (kind=kind_phys), dimension(:),   intent(in) :: trans2,trans
      real (kind=kind_phys), dimension(:,:), intent(in) :: source,ctrans

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS) :: flux_tmp, flux_tmp2

!
!===> ...  begin here
!
      flux_tmp (:) = source(:,NLP1) * trans (:)
      flux_tmp2(:) = source(:,NLAY) * trans2(:)

      fluxn(:,NLAY,m) = fluxn(:,NLAY,m) + flux_tmp(:) *ctrans(:,NLP1)
      fluxn(:,NLP1,m) = fluxn(:,NLP1,m) + flux_tmp2(:)*ctrans(:,NLP1)

      fluxncf(:,NLAY,m) = fluxncf(:,NLAY,m) + flux_tmp (:)
      fluxncf(:,NLP1,m) = fluxncf(:,NLP1,m) + flux_tmp2(:)

!
      return
!...................................
      end subroutine longwave_fluxes_KE_KEp1
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine longwave_fluxes_diag                                   &
!...................................
!  ---  inputs:
     &     ( source, trans, cld_trans, m, NPTS, NLP1 )
!  ---  outputs: ( none )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: m, NPTS, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: source,      &
     &       trans, cld_trans

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: flux_tmp

      integer :: k
!
!===> ...  begin here
!
      do k = 2, NLP1
        flux_tmp(:,k) = source(:,k) * trans(:,k)
      enddo

      do k = 2, NLP1
        fluxn(:,k,m) = fluxn(:,k,m) + flux_tmp(:,k)*cld_trans(:,k)
      enddo

      do k = 2, NLP1 
        fluxncf(:,k,m) = fluxncf(:,k,m) + flux_tmp(:,k)
      enddo

!
      return
!...................................
      end subroutine longwave_fluxes_diag
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine longwave_fluxes_sum                                    &
!...................................
!  ---  inputs:
     &     ( NPTS,                                                      &
!  ---  outputs:
     &       flux, fluxcf                                               &
     &     )

!--------------------------------------------------------------------
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: flux,fluxcf

! ---  locals:
      integer ::  m, j

!
!===> ... begin here
!
      flux   = f_zero
      fluxcf = f_zero
      do m = 1, IDXFLX
        flux  (:,:) = flux  (:,:) + fluxn  (:,:,m)
        fluxcf(:,:) = fluxcf(:,:) + fluxncf(:,:,m)
      enddo

!
      return
!...................................
      end subroutine longwave_fluxes_sum
!----------------------------------- contained in lwrad1


!  =========================================
!  *****     optical path section      *****
!  =========================================


!-----------------------------------
      subroutine optical_path_setup
!...................................
!  ---  inputs: ( from parent variables )
!  ---  outputs:( to   parent variables )

! --------------------------------------------------------------------- !
!   input variables:                                                    !
!     qo3   (NPTS,NLAY)  - o3 mass mixing ratio at model data levels    !
!     cfc11,cfc12,cfc22,cfc113                                          !
!                        - volume mixing ratio of cfc gases             !
!     temp  (NPTS,NLAY)  - temperature at data levels of model          !
!     press (NPTS,NLAY)  - pressure at data levels of model             !
!     pflux (NPTS,NLP1)  - pressure at flux levels of model             !
!     pdflux(NPTS,NLAY)  - pressure at                                  !
!     tpl1  (NPTS,NLP1)  - temperature at "upper" nearby layers of model!
!     tpl2  (NPTS,NLP1)  - temperature at "lower" nearby layers of model!
!                                                                       !
!   output variables:                                                   !
!     empl1 (NPTS,NLP1)  - h2o press scaled opt path between k lev/lay  !
!     empl2 (NPTS,NLP1)  - h2o press scaled opt path between k,k+1 lv/la!
!     emx1  (NPTS)       - h2o press scaled opt path between NLAY lv/la !
!     emx2  (NPTS)       - h2o press scaled opt path between NLP1,NLAY  !
!     var1  (NPTS,NLAY)  - h2o optical path in model layers             !
!     var2  (NPTS,NLAY)  - press-weighted h2o opt path in model layers  !
!     var3  (NPTS,NLAY)  - o3 optical path in model layers              !
!     var4  (NPTS,NLAY)  - press-weighted o3 opt path in model layers   !
!     totch2obdwd(:,:P)  - self & forgn continuum opt-path from top to k!
!     xch2obdwd  (:,:)   - self & forgn continuum opt-path for layer k  !
!     xch2obd    (:,:,:) - self & forgn continuum opt-path for layer k  !
!     totphi(NPTS,NLP1)  - total opt-path for h2o                       !
!     cntval(NPTS,NLP1)  - h2o cont path for 800-990/1070-1200 cm-1 band!
!     toto3 (NPTS,NLP1)  - total o3 opt-path from top to k level        !
!     tphio3(NPTS,NLP1)  - total press-wghted o3 opt-path from top to k !
!     totvo2(NPTS,NLP1)  - total h2o continuum path from top to k level !
!     rh2os (NPTS,NLAY)  -                                              !
!     rfrgn (NPTS,NLAY)  -                                              !
!     tfac  (NPTS,NLAY)  -                                              !
!     wk    (NPTS,NLAY)  -                                              !
!                                                                       !
!   optional output variables (ifch4n2o>0 in the 1200-1400 range)       !
!     empl1f(NPTS,NLP1)  -                                              !
!     empl2 (NPTS,NLP1)  -                                              !
!     emx1  (NPTS)       -                                              !
!     emx2  (NPTS)       -                                              !
!     vrpfh2o(NPTS,NLP1) -                                              !
!     tphfh2o(NPTS,NLP1) -                                              !
!                                                                       !
!   local variables:                                                    !
!     atmden(NPTS,NLP1)  - atmospheric density, in gm/cm**2             !
!     vv    (NPTS,NLP1)  - layer-mean pressure in atmospheres           !
! --------------------------------------------------------------------- !

      implicit none

!  ---  inputs: ( from parent variables )

!  ---  outputs:( to   parent variables )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: atmden, vv
      real (kind=kind_phys) :: f11, f12, f22, f113

      integer :: k
!
!===> ...  begin here
!

!  ---  define the layer-mean pressure in atmospheres (vv) and the layer
!       density (atmden).

      do k = 1, NLAY
        atmden(:,k) = pdflux(:,k) / grav
        vv(:,k)     = 0.5 * (pflux(:,k+1) + pflux(:,k)) / pstd
      enddo

!  ---  compute optical paths.

      if ( ifch4n2o > 0 ) then

        call optical_h2o                                                &
!  ---  inputs:
     &     ( press, pflux, rh2o, temp, tpl1, tpl2, atmden, vv,          &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       empl1, empl2, emx1, emx2, totphi, var1, var2               &
!  ---  optional in/out:
     &,      empl1f, empl2f, emx1f, emx2f, vrpfh2o, tphfh2o             &
     &     )

      else

        call optical_h2o                                                &
!  ---  inputs:
     &     ( press, pflux, rh2o, temp, tpl1, tpl2, atmden, vv,          &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       empl1, empl2, emx1, emx2, totphi, var1, var2               &
     &     )

      endif

!  ---  call Optical_ckd2p1 to determine self- and foreign-broadened h2o
!       continuum paths, for the given temperature, pressure and mixing
!       ratio, over the predetermined frequency range.

      call optical_path_ckd2p1                                          &
!  ---  inputs:
     &     ( press, rh2o, temp, atmden,                                 &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       xch2obd, xch2obdwd, totch2obdwd, rh2os, rfrgn, tfac, wk    &
     &     )

      call optical_o3                                                   &
!  ---  inputs:
     &     ( atmden, qo3, vv,                                           &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       toto3, tphio3, var3, var4                                  &
     &     )

      if ( iflgcfc == 1 ) then
!  ---  change vmr to mmr
        f11 = cfc11 * (wtmf11 / con_amd)
        f12 = cfc12 * (wtmf12 / con_amd)
        f22 = cfc22 * (wtmf22 / con_amd)
        f113= cfc113* (wtmf113/ con_amd)

        totf11 (:,1) = f_zero
        totf12 (:,1) = f_zero
        totf113(:,1) = f_zero
        totf22 (:,1) = f_zero

        do k = 2, NLP1
          totf11 (:,k) = totf11 (:,k-1) + atmden(:,k-1)*f11 * 2.0
          totf12 (:,k) = totf12 (:,k-1) + atmden(:,k-1)*f12 * 2.0
          totf113(:,k) = totf113(:,k-1) + atmden(:,k-1)*f113* 2.0
          totf22 (:,k) = totf22 (:,k-1) + atmden(:,k-1)*f22 * 2.0
        enddo
      endif

!
      return
!...................................
      end subroutine optical_path_setup
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_h2o                                            &
!...................................
!  ---  inputs:
     &     ( press, pflux, rh2o, temp, tpl1, tpl2, atmden, vv,          &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       empl1, empl2, emx1, emx2, totphi, var1, var2               &
!  ---  optional in/out:
     &,      empl1f, empl2f, emx1f, emx2f, vrpfh2o, tphfh2o             &
     &     )

! -------------------------------------------------------------------- !
!                                                                      !
!     Optical_h2o computes optical paths for h2o.                      !
!                                                                      !
!     author: m. d. schwarzkopf                                        !
!     revised: 1/1/93                                                  !
!     certified:  radiation version 1.0                                !
!                                                                      !
! -------------------------------------------------------------------- !
!   inputs from parent program:                                        !
!     press   -  pressure at data levels of model.                     !
!     pflux   -  pressure at flux levels of model.                     !
!     rh2o    -  mass mixing ratio of h2o at model data levels         !
!     temp    -  temperature at data levels of model.                  !
!     tpl1    -  temperature at "upper" nearby layers of model.        !
!     tpl2    -  temperature at "lower" nearby layers of model.        !
!                                                                      !
!   outputs to parent program:                                         !
!     empl1   -  h2o press scaled optical path between flux            !
!                level k and nearest data level k?                     !
!     empl2   -  h2o pressure scaled optical path between flux         !
!                level k and nearest data level k+1?                   !
!     emx1    -  h2o pressure scaled optical path between flux level   !
!                NLAY and data level NLAY.                             !
!     emx2    -  h2o pressure scaled optical path between flux level   !
!                NLP1 and data level NLAY.                             !
!     totphi  -  summed pressure weighted h2o optical path from top    !
!                of atmosphere to flux level k.                        !
!     var1    -  h2o optical path in model layers.                     !
!     var2    -  pressure-weighted h2o optical path in model layers.   !
!                                                                      !
!   optional output variables (ifch4n2o>0 in the 1200-1400 range)      !
!     empl1f  -                                                        !
!     empl2   -                                                        !
!     emx1    -                                                        !
!     emx2    -                                                        !
!     vrpfh2o -                                                        !
!     tphfh2o -                                                        !
!                                                                      !
!   local variables:                                                   !
!     qh2o    -  h2o mass mixing ratio, multiplied by the diffusivity  !
!                factor diffac.                                        !
!                                                                      !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: press, pflux,&
     &       rh2o, temp, tpl1, tpl2, atmden, vv

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out):: empl1, empl2,&
     &       var1, var2, totphi
      real (kind=kind_phys), dimension(:),   intent(out):: emx1, emx2

!  ---  optional outputs:
      real (kind=kind_phys), dimension(:,:), optional, intent(out) ::   &
     &       empl1f, empl2f, vrpfh2o, tphfh2o
      real (kind=kind_phys), dimension(:),   optional, intent(out) ::   &
     &       emx1f, emx2f

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: qh2o, tdif, tdif2
      real (kind=kind_phys), dimension(2)         :: csfah2o

      integer :: m, k
!
!===> ...  begin here
!

!  ---  compute optical paths for h2o, using the diffusivity approximation
!       1.66 for the angular integration.  obtain unweighted values var1,
!       and weighted values var2. the quantities 0.0003 (.0003) appearing
!       in the var2 expressions are the approximate voigt corrections for
!       h2o.  vv is the layer-mean pressure (in atmosphere), which is not
!       the same as the level pressure press.

      qh2o(:,:) = rh2o(:,:) * diffac
      var1(:,:) = atmden(:,:) * qh2o(:,:)
      var2(:,:) = var1(:,:) * (vv(:,:) + 3.0e-4)

!  ---  compute summed optical paths for h2o.

      totphi(:,1) = f_zero
      do k = 2, NLP1
        totphi(:,k) = totphi(:,k-1) + var2(:,k-1)
      enddo

!  ---  emx1 is the additional pressure-scaled mass from press(NLAY) to
!       pflux(NLAY).  it is used in nearby layer and emiss calculations.
!       emx2 is the additional pressure-scaled mass from press(NLAY) to
!       pflux(NLP1).  it is used in calculations between flux levels k
!       and NLP1.

      emx1(:) = qh2o(:,NLAY)*press(:,NLAY)                              &
     &        * (press(:,NLAY) - pflux(:,NLAY)) / (grav*pstd)
      emx2(:) = qh2o(:,NLAY)*press(:,NLAY)                              &
     &        * (pflux(:,NLP1) - press(:,NLAY)) / (grav*pstd)

!  ---  empl is the pressure scaled mass from pflux(k) to press(k) or to
!       press(k+1).

      empl1(:,1) = var2(:,NLAY)
      empl1(:,2:NLP1) = qh2o(:,1:NLAY)*pflux(:,2:NLP1)                  &
     &      * (pflux(:,2:NLP1) - press(:,1:NLAY)) / (grav*pstd)
      empl2(:,2:NLAY) = qh2o(:,2:NLAY)*pflux(:,2:NLAY)                  &
     &      * (press(:,2:NLAY) - pflux(:,2:NLAY)) / (grav*pstd)
      empl2(:,NLP1)   = empl2(:,NLAY)

      if ( ifch4n2o > 0 ) then
        csfah2o(1) = CN_basic%csf1h2o(1)
        csfah2o(2) = CN_basic%csf1h2o(2)

!  ---  compute h2o optical paths for use in the 1200-1400 cm-1 range if
!       temperature dependence of line intensities is accounted for.

        tdif(:,:) = temp(:,:) - 2.5e+2

        vrpfh2o(:,1:NLAY) = var2(:,1:NLAY)                              &
     &        * exp( csfah2o(1)*(tdif(:,1:NLAY) )                       &
     &        + csfah2o(2)*( tdif(:,1:NLAY))**2 )

        tphfh2o(:,1) = f_zero
        do k = 2, NLP1
          tphfh2o(:,k) = tphfh2o(:,k-1) + vrpfh2o(:,k-1)
        enddo

        tdif2(:,2:NLP1) = tpl2(:,2:NLP1)-2.5e+2
        tdif (:,2:NLP1) = tpl1(:,2:NLP1)-2.5e+2

!  ---  compute this additional mass, for use in the 1200-1400 cm-1 range,
!       if temperature dependence of line intensities is accounted for.

        emx1f(:) = emx1(:) * exp( csfah2o(1)*tdif2(:,NLP1)              &
     &        + csfah2o(2)*tdif2(:,NLP1)**2 )
        emx2f(:) = emx2(:) * exp( csfah2o(1)*tdif(:,NLP1)               &
     &        + csfah2o(2)*tdif(:,NLP1)**2 )

!  ---  compute this additional mass, for use in the 1200-1400 cm-1 range,
!       if temperature dependence of line intensities is accounted for.

        empl1f(:,2:NLP1) = empl1(:,2:NLP1)                              &
     &        * exp( csfah2o(1)*tdif(:,2:NLP1)                          &
     &        + csfah2o(2)*tdif(:,2:NLP1)**2 )
        empl2f(:,2:NLAY)  = empl2(:,2:NLAY)                             &
     &        * exp( csfah2o(1)*tdif2(:,2:NLAY)                         &
     &        + csfah2o(2)*tdif2(:,2:NLAY)**2 )
        empl1f(:,1) = vrpfh2o(:,NLAY)
        empl2f(:,NLP1) = empl2f(:,NLAY)
      endif

!
      return
!...................................
      end subroutine optical_h2o
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_path_ckd2p1                                    &
!...................................
!  ---  inputs:
     &     ( press, rh2o, temp, atmden,                                 &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       xch2obd, xch2obdwd, totch2obdwd, rh2os, rfrgn, tfac, wk    &
     &     )

! -------------------------------------------------------------------- !
!     subroutine Optical_ckd2p1 computes h2o continuum optical paths   !
!     (self + foreign) over the frequency range specified by IOFFH2O   !
!     and nptch2o using the ckd2.1 algorithm, modified for the gcm     !
!     parameterization. (this routine is previously called contnm.F)   !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: press, rh2o, &
     &       temp, atmden

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: xch2obd
      real (kind=kind_phys), dimension(:,:),   intent(out) :: rh2os,    &
     &       rfrgn, tfac, wk, totch2obdwd, xch2obdwd

!  ---  locals:
      real(kind=kind_phys), parameter ::  t0 = 296.0

      real (kind=kind_phys), dimension(NPTS,NLAY) :: xch2obdinw, rvh2o, &
     &       rhoave, tmpexp
      real (kind=kind_phys), dimension(NPTS,NLP1) :: totch2obdinw

      integer :: n, k, nu
!
!===> ...  begin here
!

!  ---  define the volume mixing ratio of h2o

      rvh2o(:,:) = rh2o(:,:) * con_amd / con_amw

!  ---  define input arguments to optical_ckd2p1

      wk    (:,:) = rvh2o(:,:) * con_avgd / con_amd * atmden(:,:)
      rhoave(:,:) = (press(:,:)/pstd) * (t0/temp(:,:))
      rh2os (:,:) = rvh2o(:,:) * rhoave(:,:)
      rfrgn (:,:) = rhoave(:,:) - rh2os(:,:)
      tfac  (:,:) = temp(:,:) - t0

!  ---  compute self-broadened temperature-dependent continuum coefficient
!       using the single coefficient -.020 for all frequencies in the 560-
!       1200 cm-1 range. experiments with the mid-latitude summer profile
!       show errors of < .01 W/m**2 (in the net broadband flux, 0-2200 cm-1)
!       using this value. this value is used instead of tmpfctrs at each
!       frequency band.

      tmpexp(:,:) = exp( -.020*tfac(:,:) )

!  ---  compute h2o self- and foreign- broadened continuum optical path
!       for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!       the top of the atmosphere through layer k (totch2obd, totch2obdinw,
!       totch2obdwd).

      do nu = 1, 7
        do k = 1, NLAY
          xch2obd(:,k,nu) = wk(:,k)*1.0e-20                             &
     &        * ( svj(nu)*rh2os(:,k)*tmpexp(:,k) + fvj(nu)*rfrgn(:,k) ) &
     &        * radfnbd(nu)
        enddo
      enddo

      do k = 1, NLAY
        xch2obdinw(:,k) = wk(:,k)*1.0e-20                               &
     &        * ( svjinw*rh2os(:,k)*tmpexp(:,k) + fvjinw*rfrgn(:,k) )   &
     &        * radfnbdinw
        xch2obdwd (:,k) = wk(:,k)*1.0e-20                               &
     &        * ( svjwd*rh2os(:,k)*tmpexp(:,k) + fvjwd*rfrgn(:,k) )     &
     &        * radfnbdwd
      enddo

      totch2obdinw(:,1) = f_zero
      totch2obdwd (:,1) = f_zero
      do k = 2, NLP1
        totch2obdinw(:,k) = totch2obdinw(:,k-1) + xch2obdinw(:,k-1)
        totch2obdwd (:,k) = totch2obdwd (:,k-1) + xch2obdwd (:,k-1)
      enddo

!
      return
!...................................
      end subroutine optical_path_ckd2p1
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_o3                                             &
!...................................
!  ---  inputs:
     &     ( atmden, qo3, vv,                                           &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       toto3, tphio3, var3, var4                                  &
     &     )

! -------------------------------------------------------------------- !
!                                                                      !
!     optical_o3 computes optical paths for o3.                        !
!                                                                      !
!     author: m. d. schwarzkopf                                        !
!     revised: 1/1/93                                                  !
!     certified:  radiation version 1.0                                !
!                                                                      !
! -------------------------------------------------------------------- !
!   inputs from parent program:                                        !
!     qo3     - mass mixing ratio of o3 at model data levels.          !
!     atmden  -                                                        !
!     vv      -                                                        !
!                                                                      !
!   outputs to parent program:                                         !
!     toto3   - summed o3 optical path from the top of atmosphere to   !
!               flux level k.                                          !
!     tphio3  - summed pressure weighted o3 optical path from top of   !
!               atmosphere to flux level k.                            !
!     var3    - o3 optical path in model layers.                       !
!     var4    - pressure-weighted o3 optical path in model layers.     !
!                                                                      !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:),intent(in) :: qo3,atmden,vv

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),intent(out) :: toto3,       &
     &       tphio3, var3, var4

!  ---  locals:
      integer :: k
!
!===> ...  begin here
!

!  ---  compute optical paths for o3, using the diffusivity approximation
!       1.66 for the angular integration.  obtain unweighted values var3
!       and weighted values  var4. the quantities  0.003 (.003) appearing
!       in the var4 expression are the approximate voigt corrections for o3.

      var3(:,:) = atmden(:,:) * qo3(:,:) * diffac
      var4(:,:) = var3(:,:) * (vv(:,:) + 3.0e-3)

!  ---  compute summed optical paths for o3.

      toto3 (:,1) = f_zero
      tphio3(:,1) = f_zero
      do k = 2, NLP1
        toto3 (:,k) = toto3 (:,k-1) + var3(:,k-1)
        tphio3(:,k) = tphio3(:,k-1) + var4(:,k-1)
      enddo

!
      return
!...................................
      end subroutine optical_o3
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_trans_funct_from_KS                            &
!...................................
!  ---  inputs:
     &     ( tn2o17, NPTS, NLAY, NLP1,                                  &
!  ---  outputs:
     &       to3cnt, overod, cnttaub1, cnttaub2, cnttaub3               &
     &     )

! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      integer :: NPTS, NLAY, NLP1

      real (kind=kind_phys), intent(in) :: tn2o17(:,:)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: to3cnt,     &
     &       overod, cnttaub1, cnttaub2, cnttaub3

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: totch2o_tmp, tmp1, &
     &       tmp2, tmp3
      real (kind=kind_phys), dimension(NPTS,NLAY) :: cfc_tf

      integer :: k, m

!
!===> ...  begin here
!

!  ---  compute transmission functions in 990-1070 cm-1 range, including
!       ozone and h2o continuum, from level 1 to all other levels.

      tmp1(:,1:NLAY) = bo3rnd(2)*tphio3(:,2:NLP1)/toto3(:,2:NLP1)
      tmp2(:,1:NLAY) = 0.5 * ( tmp1(:,1:NLAY)                           &
     &      * ( sqrt(1.0 + (4.0*ao3rnd(2)*toto3(:,2:NLP1))              &
     &      / tmp1(:,1:NLAY) ) - 1.0 ) )

      totch2o_tmp(:,1) = f_zero
      do k = 2, NLP1
        totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,6)
      enddo

      tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + diffac*totch2o_tmp(:,2:NLP1)

      if ( iaerlw > 0 ) then
!       call get_totaerooptdep(6, totaer_tmp)
        if ( NBDLW == 1 ) then
          tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        else                                           ! to be developed !!
          tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        endif
      endif

      to3cnt(:,1:NLAY) = exp( -1.0*tmp2(:,1:NLAY) )

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in to3cnt.

      if ( iflgcfc == 1 ) then

        cfc_tf(:,1:NLAY) = 1.0                                          &
     &                   - CFC_basic%strf113(6)*totf113(:,2:NLP1)       &
     &                   - CFC_basic%strf22 (6)*totf22 (:,2:NLP1)       &
     &                   - CFC_basic%strf11 (6)*totf11 (:,2:NLP1)       &
     &                   - CFC_basic%strf12 (6)*totf12 (:,2:NLP1)

        to3cnt(:,1:NLAY) = to3cnt(:,1:NLAY)* cfc_tf(:,1:NLAY)
      endif

!  ---  compute transmission function in the 560-800 cm-1 range evaluate
!       optical depth contributions add contributions from h2o(lines)
!       and h2o(continuum), either Roberts or CKD2.1

      tmp1(:,1:NLAY) = sqrt( ab15wd*totphi(:,2:NLP1) )

      tmp1(:,1:NLAY) = tmp1(:,1:NLAY) + diffac*totch2obdwd(:,2:NLP1)

!  ---  add contribution from longwave aerosols (if desired).

      if ( iaerlw > 0 ) then
!       tmp1(:,1:NLAY) = tmp1(:,1:NLAY) + totaerooptdep_15(:,2:NLP1)
        if ( NBDLW == 1 ) then
          tmp1(:,1:NLAY) = tmp1(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        else                                     ! to be developed !!
          tmp1(:,1:NLAY) = tmp1(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        endif
      endif

!  ---  compute transmission function due to these contributions. the
!       effects of co2, n2o  and  cfc's (not exponentials) are added
!       later.

      overod(:,1:NLAY) = exp( -1.0*tmp1(:,1:NLAY) )

!  ---  add contribution from the 17 um n2o band (if desired). the
!       expression with tn2o17 retains the 560-630 cm-1 equivalent
!       widths in evaluating 560-800 cm-1 transmissivities.

      if ( ifch4n2o > 0 ) then
        if     ( NBCO215 == 2 ) then
          overod(:,1:NLAY) = overod(:,1:NLAY)                           &
     &         * ( 130.0/240.0 + 110.0/240.0*tn2o17(:,2:NLP1) )
        elseif ( NBCO215 == 3 ) then
          overod(:,1:NLAY) = overod(:,1:NLAY)                           &
     &         * ( 170.0/240.0 + 70.0/240.0*tn2o17(:,2:NLP1) )
        endif
      endif

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in overod .

      if ( iflgcfc == 1 ) then

        cfc_tf(:,1:NLAY) = 1.0                                          &
     &                   - CFC_basic%sf11315*totf113(:,2:NLP1)          &
     &                   - CFC_basic%sf2215 *totf22 (:,2:NLP1)          &
     &                   - CFC_basic%sf1115 *totf11 (:,2:NLP1)          &
     &                   - CFC_basic%sf1215 *totf12 (:,2:NLP1)

        overod(:,1:NLAY) = overod(:,1:NLAY)*cfc_tf(:,1:NLAY)
      endif

!  ---  compute continuum band transmission functions from level 1 to
!       other levels (cnttau). the continuum transmission function from
!       level k to kp (contod) equals cnttau for k=1, so is not evaluated
!       here. for all other levels k, contod is obtained by division of
!       relevant values of cnttau.

      totch2o_tmp(:,1) = f_zero
      do k = 2, NLP1
        totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,4)
      enddo
      tmp1(:,1:NLAY) = diffac*totch2o_tmp(:,2:NLP1)

      totch2o_tmp(:,1) = f_zero
      do k = 2, NLP1
        totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,5)
      enddo
      tmp2(:,1:NLAY) = diffac*totch2o_tmp(:,2:NLP1)

      totch2o_tmp(:,1) = f_zero
      do k = 2, NLP1
        totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,7)
      enddo
      tmp3(:,1:NLAY) = diffac*totch2o_tmp(:,2:NLP1)

      if ( iaerlw > 0 ) then
!       call get_totaerooptdep(4, totaer_tmp)
!       tmp1(:,1:NLAY) =  tmp1(:,1:NLAY) + totaer_tmp(:,2:NLP1)
!       call get_totaerooptdep(5, totaer_tmp)
!       tmp2(:,1:NLAY) =  tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
!       call get_totaerooptdep(7, totaer_tmp)
!       tmp3(:,1:NLAY) =  tmp3(:,1:NLAY) + totaer_tmp(:,2:NLP1)
        if ( NBDLW == 1 ) then
          tmp1(:,1:NLAY) =  tmp1(:,1:NLAY) + tauaertot(:,2:NLP1,1)
          tmp2(:,1:NLAY) =  tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
          tmp3(:,1:NLAY) =  tmp3(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        else                                ! to be developed !!
          tmp1(:,1:NLAY) =  tmp1(:,1:NLAY) + tauaertot(:,2:NLP1,1)
          tmp2(:,1:NLAY) =  tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
          tmp3(:,1:NLAY) =  tmp3(:,1:NLAY) + tauaertot(:,2:NLP1,1)
        endif
      endif

      cnttaub1(:,1:NLAY) = exp( -1.0*tmp1(:,1:NLAY) )
      cnttaub2(:,1:NLAY) = exp( -1.0*tmp2(:,1:NLAY) )
      cnttaub3(:,1:NLAY) = exp( -1.0*tmp3(:,1:NLAY) )

!  ---  if cfcs are included, add transmission functions for f11, f12,
!       f113, and f22.

      if ( iflgcfc == 1 ) then

        cfc_tf(:,1:NLAY) = 1.0                                          &
     &                   - CFC_basic%strf113(4)*totf113(:,2:NLP1)       &
     &                   - CFC_basic%strf22 (4)*totf22 (:,2:NLP1)       &
     &                   - CFC_basic%strf11 (4)*totf11 (:,2:NLP1)       &
     &                   - CFC_basic%strf12 (4)*totf12 (:,2:NLP1)
        cnttaub1(:,1:NLAY) = cnttaub1(:,1:NLAY)*cfc_tf(:,1:NLAY)

        cfc_tf(:,1:NLAY) = 1.0                                          &
     &                   - CFC_basic%strf113(5)*totf113(:,2:NLP1)       &
     &                   - CFC_basic%strf22 (5)*totf22 (:,2:NLP1)       &
     &                   - CFC_basic%strf11 (5)*totf11 (:,2:NLP1)       &
     &                   - CFC_basic%strf12 (5)*totf12 (:,2:NLP1)
        cnttaub2(:,1:NLAY) = cnttaub2(:,1:NLAY)*cfc_tf(:,1:NLAY)

        cfc_tf(:,1:NLAY) = 1.0                                          &
     &                   - CFC_basic%strf113(7)*totf113(:,2:NLP1)       &
     &                   - CFC_basic%strf22 (7)*totf22 (:,2:NLP1)       &
     &                   - CFC_basic%strf11 (7)*totf11 (:,2:NLP1)       &
     &                   - CFC_basic%strf12 (7)*totf12 (:,2:NLP1)
        cnttaub3(:,1:NLAY) = cnttaub3(:,1:NLAY)*cfc_tf(:,1:NLAY)

      endif

!  ---  evaluate h2o (mbar*phibar) between level 1 and other levels.

      avephi(:,1:NLAY) = totphi(:,2:NLP1)

!  ---  the evaluation of emiss over the layer between data level (1)
!       and flux level (NLP1) is done by averaging E2 functions
!       referring to the top and bottom of the layer. a special value
!       of (mbar*phibar) is required; it is stored in the (otherwise
!       vacant) NLP1'th position of avephi.

      avephi(:,NLP1) = avephi(:,NLAY-1) + emx1(:)

!  ---  if h2o lines in the 1200-1400 range are assumed to have a 
!       temperature dependent intensity, similar evaluation for 
!       (mbar*phibar) is performed, with a special value for the
!       lowest layer

      if ( ifch4n2o > 0 ) then
        avephif(:,1:NLAY) = tphfh2o(:,2:NLP1)
        avephif(:,NLP1)  = avephif(:,NLAY-1) + emx1f(:)
      endif

!
      return
!...................................
      end subroutine optical_trans_funct_from_KS
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_trans_funct_k_down                             &
!...................................
!  ---  inputs:
     &     ( k, cnttaub1, cnttaub2, cnttaub3, tn2o17,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       to3cnt, overod, contodb1, contodb2, contodb3               &
     &     )

! -------------------------------------------------------------------- !
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, k

      real (kind=kind_phys), dimension(:,:), intent(in) :: cnttaub1,    &
     &       cnttaub2, cnttaub3, tn2o17

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: to3cnt,     &
     &       overod, contodb1, contodb2, contodb3

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: avmo3, avpho3,     &
     &       tmp1, tmp2, avckdwd, avckdo3, totch2o_tmp, avaero3
      real (kind=kind_phys), dimension(NPTS,NLAY) :: cfc_tf

      integer :: kp, m
!
!===> ...  begin here
!

      totch2o_tmp(:,1) = f_zero
      do kp = 2, NLP1
        totch2o_tmp(:,kp) = totch2o_tmp(:,kp-1)+xch2obd(:,kp-1,6)
      enddo

      if ( iaerlw > 0 ) then
!       call get_totaerooptdep(6, totaer_tmp)
        if ( NBDLW == 1 ) then
          do kp = 1, NLP1-k
            avaero3(:,kp+k-1) = tauaertot(:,kp+k,1) - tauaertot(:,k,1)
          enddo
        else                                ! to be developed !!
          do kp = 1, NLP1-k
            avaero3(:,kp+k-1) = tauaertot(:,kp+k,1) - tauaertot(:,k,1)
          enddo
        endif
      endif

       do kp = 1, NLP1-k
         avmo3 (:,kp+k-1) = toto3 (:,kp+k) - toto3 (:,k)
         avpho3(:,kp+k-1) = tphio3(:,kp+k) - tphio3(:,k)

         avckdwd(:,kp+k-1) = totch2obdwd(:,kp+k) - totch2obdwd(:,k)
         avckdo3(:,kp+k-1) = totch2o_tmp(:,kp+k) - totch2o_tmp(:,k)
       enddo

       do kp = 1, NLP1-k
         contodb1(:,kp+k-1) = cnttaub1(:,kp+k-1)/cnttaub1(:,k-1)
         contodb2(:,kp+k-1) = cnttaub2(:,kp+k-1)/cnttaub2(:,k-1)
         contodb3(:,kp+k-1) = cnttaub3(:,kp+k-1)/cnttaub3(:,k-1)
         avephi  (:,kp+k-1) = totphi(:,kp+k) - totphi(:,k)
       enddo
       avephi(:,NLP1) = avephi(:,NLAY-1) + emx1(:)

!  ---  if h2o lines in the 1200-1400 range are assumed to have a temp-
!       erature dependent intensity, similar evaluation for (mbar*phibar)
!       is performed, with a special value for the lowest layer

       if ( ifch4n2o > 0 ) then
         do kp = 1, NLP1-k
           avephif(:,kp+k-1) = tphfh2o(:,kp+k) - tphfh2o(:,k)
         enddo

         avephif(:,NLP1) = avephif(:,NLAY-1) + emx1f(:)
       endif

!  ---  compute transmission function in the 560-800 cm-1 range evaluate
!       optical depth contributions, add contributions from h2o(lines) and
!       h2o(continuum) (either Roberts or CKD2.1)

       tmp1(:,k:NLAY) = sqrt(ab15wd*avephi(:,k:NLAY))

       tmp1(:,k:NLAY) = tmp1(:,k:NLAY) + diffac*avckdwd(:,k:NLAY)

!  ---  add contribution from longwave aerosols (if desired).

       if ( iaerlw > 0 ) then
!       do kp = k, NLAY
!         tmp1(:,kp) = tmp1(:,kp) + ( totaerooptdep_15(:,kp+1)          &
!                    - totaerooptdep_15(:,k) )
!       enddo
        if ( NBDLW == 1 ) then
          do kp = k, NLAY
            tmp1(:,kp) = tmp1(:,kp) + ( tauaertot(:,kp+1,1)             &
     &                 - tauaertot(:,k,1) )
          enddo
        else                                ! to be developed !!
          do kp = k, NLAY
            tmp1(:,kp) = tmp1(:,kp) + ( tauaertot(:,kp+1,1)             &
     &                 - tauaertot(:,k,1) )
          enddo
        endif
       endif

!  ---  compute transmission function due to these contributions. the
!       effects of co2, n2o  and  cfc's (not exponentials) are added later.

      overod(:,k:NLAY) = exp( -1.0*tmp1(:,k:NLAY) )

!  ---  add contribution from the 17 um n2o band (if desired). the expression
!       with tn2o17 retains the 560-630 cm-1 equivalent widths in evaluating
!       560-800 cm-1 transmissivities.

      if ( ifch4n2o > 0 ) then
        if ( NBCO215 == 2 ) then
          overod(:,k:NLAY) = overod(:,k:NLAY)                           &
     &         * ( 130.0/240.0 + 110.0/240.0*tn2o17(:,k+1:NLP1) )
        elseif (NBCO215 == 3 ) then
          overod(:,k:NLAY) = overod(:,k:NLAY)                           &
     &         * ( 170.0/240.0 + 70.0/240.0*tn2o17(:,k+1:NLP1) )
        endif
      endif

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in overod .

      if ( iflgcfc == 1 ) then
        do kp = k, NLAY
          cfc_tf(:,kp) = 1.0                                            &
     &            - CFC_basic%sf11315*(totf113(:,kp+1) - totf113(:,k))  &
     &            - CFC_basic%sf2215 *(totf22 (:,kp+1) - totf22 (:,k))  &
     &            - CFC_basic%sf1115 *(totf11 (:,kp+1) - totf11 (:,k))  &
     &            - CFC_basic%sf1215 *(totf12 (:,kp+1) - totf12 (:,k))
        enddo

        overod(:,k:NLAY) = overod(:,k:NLAY)*cfc_tf(:,k:NLAY)
      endif

!  ---  compute transmission functions in 990-1070 cm-1 range, including
!       ozone and h2o continuum, from level k to all other levels.

      tmp1(:,k:NLAY) = bo3rnd(2)*avpho3(:,k:NLAY)/avmo3(:,k:NLAY)
      tmp2(:,k:NLAY) = 0.5 * ( tmp1(:,k:NLAY)                           &
     &       * ( sqrt( 1.0 + ( 4.0*ao3rnd(2)*avmo3(:,k:NLAY) )          &
     &       / tmp1(:,k:NLAY) ) - 1.0 ) )

      tmp2(:,k:NLAY) = tmp2(:,k:NLAY) + diffac*avckdo3(:,k:NLAY)

      if ( iaerlw > 0 ) then
        tmp2(:,k:NLAY) = tmp2(:,k:NLAY) + avaero3(:,k:NLAY)
      endif

      to3cnt(:,k:NLAY) = exp( -1.0*tmp2(:,k:NLAY) )

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in to3cnt.

      if ( iflgcfc == 1 ) then
        do kp = k, NLAY
          cfc_tf(:,kp) = 1.0                                            &
     &         - CFC_basic%strf113(6)*(totf113(:,kp+1) - totf113(:,k))  &
     &         - CFC_basic%strf22 (6)*(totf22 (:,kp+1) - totf22 (:,k))  &
     &         - CFC_basic%strf11 (6)*(totf11 (:,kp+1) - totf11 (:,k))  &
     &         - CFC_basic%strf12 (6)*(totf12 (:,kp+1) - totf12 (:,k))
        enddo

        to3cnt(:,k:NLAY) = to3cnt(:,k:NLAY)*cfc_tf(:,k:NLAY)
      endif

!
      return
!...................................
      end subroutine optical_trans_funct_k_down
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_trans_funct_KE                                 &
!...................................
!  ---  inputs:
     &     ( cnttaub1, cnttaub2, cnttaub3, tn2o17,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       to3cnt, overod, contodb1, contodb2, contodb3               &
     &     )

! --------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: cnttaub1,    &
     &       cnttaub2, cnttaub3, tn2o17

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: to3cnt,     &
     &       overod, contodb1, contodb2, contodb3

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: tmp1, tmp2
      real (kind=kind_phys), dimension(NPTS,NLAY) :: cfc_tf

      integer :: k
!
!===> ...  begin here
!

!  ---  compute transmission function in the 560-800 cm-1 range
!       add optical depth contributions from h2o(lines) and h2o(continuum).
!       h2o(continuum) contributions are either Roberts or CKD2.1

      tmp1(:,NLAY) = sqrt( ab15wd*var2(:,NLAY) )
      tmp1(:,NLAY) = tmp1(:,NLAY) + diffac*xch2obdwd(:,NLAY)

!  ---  add contribution from longwave aerosols (if desired).

      if ( iaerlw > 0 ) then
!       tmp1(:,NLAY) = tmp1(:,NLAY) + aerooptdep_KE_15(:)
        if (NBDLW == 1 ) then
          tmp1(:,NLAY) = tmp1(:,NLAY)                                   &
     &                 + (tauaertot(:,NLP1,1) - tauaertot(:,NLAY,1))
        else                               ! to be developed !!
          tmp1(:,NLAY) = tmp1(:,NLAY)                                   &
     &                 + (tauaertot(:,NLP1,1) - tauaertot(:,NLAY,1))
        endif
      endif

!  ---  compute transmission function due to these contributions. the
!       effects of co2, n2o  and  cfc's (not exponentials) are added later.

      overod(:,NLAY) = exp( -1.0*tmp1(:,NLAY) )

!  ---  add contribution from the 17 um n2o band.  the expression with
!       tn2o17 retains the 560-630 cm-1 equivalent widths in evaluating
!       560-800 cm-1 transmissivities.

      if ( ifch4n2o > 0 ) then
        if ( NBCO215 == 2 ) then
          overod(:,NLAY) = overod(:,NLAY)                               &
     &       * ( 130.0/240.0 + 110.0/240.0*tn2o17(:,NLP1) )
        elseif ( NBCO215 == 3 ) then
          overod(:,NLAY) = overod(:,NLAY)                               &
     &       * ( 170.0/240.0 + 70.0/240.0*tn2o17(:,NLP1) )
        endif
      endif

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in overod .

      if ( iflgcfc == 1 ) then
        cfc_tf(:,NLAY) = 1.0                                            &
     &        - CFC_basic%sf11315*(totf113(:,NLP1) - totf113(:,NLAY))   &
     &        - CFC_basic%sf2215 *(totf22 (:,NLP1) - totf22 (:,NLAY))   &
     &        - CFC_basic%sf1115 *(totf11 (:,NLP1) - totf11 (:,NLAY))   &
     &        - CFC_basic%sf1215 *(totf12 (:,NLP1) - totf12 (:,NLAY))

        overod(:,NLAY) = overod(:,NLAY)*cfc_tf(:,NLAY)
      endif

!  ---  compute transmission functions in 990-1070 cm-1 range, including
!       ozone and h2o continuum, from level 1 to all other levels.

      tmp1(:,NLAY) = bo3rnd(2)*var4(:,NLAY)/var3(:,NLAY)
      tmp2(:,NLAY) = 0.5 * ( tmp1(:,NLAY)                               &
     &       * ( sqrt( 1.0 + ( 4.0*ao3rnd(2)*var3(:,NLAY) )             &
     &       / tmp1(:,NLAY) ) - 1.0 ) )

      tmp2(:,NLAY) = tmp2(:,NLAY) + diffac*xch2obd(:,NLAY,6)

      to3cnt(:,NLAY) = exp( -1.0*tmp2(:,NLAY) )

!  ---  if cfcs are included, also include the transmission functions for
!       f11, f12, f113, and f22 in overod and to3cnt.

      if ( iflgcfc == 1 ) then
        cfc_tf(:,NLAY) = 1.0                                            &
     &      - CFC_basic%strf113(6)*(totf113(:,NLP1) - totf113(:,NLAY))  &
     &      - CFC_basic%strf22 (6)*(totf22 (:,NLP1) - totf22 (:,NLAY))  &
     &      - CFC_basic%strf11 (6)*(totf11 (:,NLP1) - totf11 (:,NLAY))  &
     &      - CFC_basic%strf12 (6)*(totf12 (:,NLP1) - totf12 (:,NLAY))

        to3cnt(:,NLAY) = to3cnt(:,NLAY)*cfc_tf(:,NLAY)
      endif

      contodb1(:,NLAY) = cnttaub1(:,NLAY) / cnttaub1(:,NLAY-1)
      contodb2(:,NLAY) = cnttaub2(:,NLAY) / cnttaub2(:,NLAY-1)
      contodb3(:,NLAY) = cnttaub3(:,NLAY) / cnttaub3(:,NLAY-1)

!
      return
!...................................
      end subroutine optical_trans_funct_KE
!----------------------------------- contained in lwrad1


!-----------------------------------
      subroutine optical_trans_funct_diag
!...................................
!  ---  inputs:  ( use in-scope parent program variables )
!  ---  outputs: ( use in-scope parent program variables )

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: ca, cb,            &
     &       csuba, csubb, ctmp2, ctmp3, delpr1, delpr2

!
!===> ...  begin here
!
      delpr1(:,2:NLAY) = pdfinv(:,2:NLAY)*(press(:,2:NLAY)              &
     &                 - pflux(:,2:NLAY))
      delpr2(:,2:NLP1) = pdfinv(:,1:NLAY)*(pflux(:,2:NLP1)              &
     &                  - press(:,1:NLAY))

!  ---  compute nearby-layer transmissivities for the o3 band and for the
!       one-band continuum band.  the sf function is used. the method is
!       the same as described for co2 in reference(4).
!  ---  compute sf2.  continuum band 1

      csuba(:,2:NLAY) = diffac*xch2obd(:,2:NLAY,4)*delpr1(:,2:NLAY)
      csubb(:,2:NLAY) = diffac*xch2obd(:,1:NLAY-1,4)*delpr2(:,2:NLAY)

      ca(:,2:NLAY) = csuba(:,2:NLAY) * ( -0.5 + csuba(:,2:NLAY)         &
     &     * (0.166666 - csuba(:,2:NLAY)*0.416666e-1) )
      cb(:,2:NLAY) = csubb(:,2:NLAY) * ( -0.5 + csubb(:,2:NLAY)         &
     &     * (0.166666 - csubb(:,2:NLAY)*0.416666e-1) )

      contdg(:,NLP1,1)  = 1.0 + cb(:,NLAY)
      contdg(:,2:NLAY,1) = 1.0 + 0.5*( ca(:,2:NLAY) + cb(:,2:NLAY) )

!  ---  continuum band 2

      csuba(:,2:NLAY) = diffac*xch2obd(:,2:NLAY,5)*delpr1(:,2:NLAY)
      csubb(:,2:NLAY) = diffac*xch2obd(:,1:NLAY-1,5)*delpr2(:,2:NLAY)

      ca(:,2:NLAY) = csuba(:,2:NLAY) * ( -0.5 + csuba(:,2:NLAY)         &
     &     * (0.166666 - csuba(:,2:NLAY)*0.416666e-1) )
      cb(:,2:NLAY) = csubb(:,2:NLAY) * ( -0.5 + csubb(:,2:NLAY)         &
     &     * (0.166666 - csubb(:,2:NLAY)*0.416666e-1) )

      contdg(:,NLP1,2)  = 1.0 + cb(:,NLAY)
      contdg(:,2:NLAY,2) = 1.0 + 0.5*( ca(:,2:NLAY) + cb(:,2:NLAY) )

!  ---  continuum band 3

      csuba(:,2:NLAY) = diffac*xch2obd(:,2:NLAY,7)*delpr1(:,2:NLAY)
      csubb(:,2:NLAY) = diffac*xch2obd(:,1:NLAY-1,7)*delpr2(:,2:NLAY)

      ca(:,2:NLAY) = csuba(:,2:NLAY) * ( -0.5 + csuba(:,2:NLAY)         &
     &     * (0.166666 - csuba(:,2:NLAY)*0.416666e-1) )
      cb(:,2:NLAY) = csubb(:,2:NLAY) * ( -0.5 + csubb(:,2:NLAY)         &
     &     * (0.166666 - csubb(:,2:NLAY)*0.416666e-1) )

      contdg(:,NLP1,3)  = 1.0 + cb(:,NLAY)
      contdg(:,2:NLAY,3) = 1.0 + 0.5*( ca(:,2:NLAY) + cb(:,2:NLAY) )

!  ---  ozone band

      csuba(:,2:NLAY) = diffac*xch2obd(:,2:NLAY,6)*delpr1(:,2:NLAY)
      csubb(:,2:NLAY) = diffac*xch2obd(:,1:NLAY-1,6)*delpr2(:,2:NLAY)

      ca(:,2:NLAY) = csuba(:,2:NLAY) * ( -0.5 + csuba(:,2:NLAY)         &
     &     * (0.166666 - csuba(:,2:NLAY)*0.416666e-1) )
      cb(:,2:NLAY) = csubb(:,2:NLAY) * ( -0.5 + csubb(:,2:NLAY)         &
     &     * (0.166666 - csubb(:,2:NLAY)*0.416666e-1) )

      to3dg(:,NLP1)  = 1.0 + cb(:,NLAY)
      to3dg(:,2:NLAY) = 1.0 + 0.5*( ca(:,2:NLAY) + cb(:,2:NLAY) )

!
      return
!...................................
      end subroutine optical_trans_funct_diag
!----------------------------------- contained in lwrad1
!
!...................................
      end subroutine lwrad1
!-----------------------------------



!-----------------------------------
      subroutine get_totch2o                                            &
!...................................
!  ---  inputs:
     &     ( n, vvj, dte1, ixoe1, tfac, rh2os, rfrgn, wk,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       totch2o                                                    &
     &     )

! -------------------------------------------------------------------- !
!                                                                      !
!                                                                      !
!                                                                      !
! -------------------------------------------------------------------- !

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, n, ixoe1(:,:)

      real (kind=kind_phys), dimension(:,:), intent(in) :: dte1, tfac,  &
     &       rh2os, rfrgn, wk
      real (kind=kind_phys), dimension(:),   intent(in) :: vvj

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: totch2o

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: radf, sh2o, tmpexp
      real (kind=kind_phys) :: fh2o0, sh2o0

      integer :: k, nu

      type (tab3_type) :: radfunc
!
!===> ...  begin here
!

      call table_alloc                                                  &
!  ---  inputs:
     &     ( 40, 300,                                                   &
!  ---  in/outputs:
     &       radfunc                                                    &
     &     )

      radfunc%vae(:,:) = func_vae(:,:)
      radfunc%td (:,:) = func_td (:,:)

!  ---  compute self-broadened temperature-dependent continuum coefficient
!       using the single coefficient -.013 for all frequencies in the 160-
!       560 cm-1 range. experiments with the mid-latitude summer profile
!       show errors of < .01 W/m**2 (in the net broadband flux, 0-2200 cm-1)
!       using this value. this value is used instead of tmpfctrs at each
!       frequency band.

      tmpexp(:,:) = exp( -.013*tfac(:,:) )

!  ---  compute source function for frequency bands (IOFFH2O+1 to IOFFH2O
!       +nptch2o) at layer temperatures using table lookup. note that ixoe1
!       can be used for temp index, and dte1 for deltat, as the table 
!       extent for radf is the same as for the e1 tables of the model.

      nu = n

      call looktab                                                      &
!  ---  inputs:
     &     ( radfunc, ixoe1, dte1, 1, NLAY, nu+IOFFH2O,                 &
!  ---  outputs:
     &       radf                                                       &
     &     )

      sh2o0 = ssh2o_296(nu+IOFFH2O) * sfac(nu+IOFFH2O)

      do k = 1, NLAY
        sh2o(:,k) = sh2o0 * tmpexp(:,k)
      enddo

!  ---  compute h2o self- and foreign- broadened continuum optical path,
!       summed from the top of the atmosphere through layer k.

      fh2o0 = sfh2o(nu+IOFFH2O) * fscal(nu+IOFFH2O)

      totch2o(:,1) = f_zero
      do k = 2, NLP1
        totch2o(:,k) = wk(:,k-1)*1.0e-20                                &
     &               * (sh2o(:,k-1)*rh2os(:,k-1) + fh2o0*rfrgn(:,k-1))  &
     &               * vvj(nu)*radf(:,k-1) + totch2o(:,k-1)
      enddo

!
      return
!...................................
      end subroutine get_totch2o
!-----------------------------------



!  =========================================
!  *****      co2_source section       *****
!  =========================================


!-----------------------------------
      subroutine co2_source_calc                                        &
!...................................
!  ---  inputs:
     &     ( press, rco2, pdfinv, dte1, ixoe1,                          &
     &       tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       soe1, soe2, soe3, soe4, soe5, sorc                         &
     &     )

!----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: press,       &
     &       pdfinv, dte1, tdav, tstdav, tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2
      real (kind=kind_phys), intent(in) :: rco2

      integer, dimension(:,:), intent(in) :: ixoe1

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: sorc
      real (kind=kind_phys), dimension(:,:),   intent(out) :: soe1,     &
     &       soe2, soe3, soe4, soe5

!  ---  locals:
      integer :: n, n1

!
!===> ...  begin here
!
!  ---   compute source function for frequency bands (9+IOFFSET to
!        NBLY-1) at layer temperatures using table lookup.

      do n = 9+IOFFSET, NBLY-1
        n1 = n - 8 - IOFFSET

        call looktab                                                    &
!  ---  inputs:
     &     ( tabsr, ixoe1, dte1, 1, NLP1, n,                            &
!  ---  outputs:
     &       sorc(:,:,n1)                                               &
     &     )

      enddo
 
!  ---  compute the nlte source function for co2.

      call nlte                                                         &
!  ---  inputs:
     &     ( press,rco2,pdfinv,tdav,tstdav,tlsqu,tmpdiff,a1,a2,         &
!  ---  in/outputs:
     &       sorc                                                       &
     &     )

!  ---  pass the values needed in longwave_driver back there. sorc will
!       also be used in cool_to_space_exact.

      soe1(:,:)  = sorc(:,:,1) + sorc(:,:,2) + sorc(:,:,3)
      soe2(:,:)  = sorc(:,:,6)
      soe3(:,:)  = sorc(:,:,4)
      soe4(:,:)  = sorc(:,:,5)
      soe5(:,:)  = sorc(:,:,7)

! ================
      contains
! ================

!-----------------------------------
      subroutine nlte                                                   &
!...................................
!  ---  inputs:
     &     ( press,rrvco2,pdfinv,tdav,tstdav,tlsqu,tmpdiff,a1,a2,       &
!  ---  in/outputs:
     &       sorc                                                       &
     &     )

!-----------------------------------------------------------------------
!
!     nlte is the present formulation of an nlte calculation of the
!     source function in the 15 um region (two bands).
!
!     the essential theory is:  phi = C*j
                    !             j = b + E*phi
!     where:  C = Curtis matrix
!             E = NLTE contribution (diagonal matrix)
!           phi = heating rate vector
!             b = LTE source function vector
!             j = NLTE source function vector
!
!             j = b (by assumption) for pressure layers > ixnltr
!             j = b (by assumption) for pressure layers > IXPRNLTE
!      E is obtained using a formulation devised by Fels (denoted
!      Ri in his notes).
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93       certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
!
!  inputs:
!
!  outputs:
!
!  locals:
!     degeneracy factor = 0.5
!     fnlte  = NLTE contribution: (E in above notes)
!     phifx  = fixed portion of PHI (contributions from
!              layers > ixnltr, where j(k) = b(k))
!              layers > IXPRNLTE, where j(k) = b(k))
!     phivar = varying portion of PHI (contributions
!              from layers <= IXPRNLTE).
!              from layers <= ixnltr).
!
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      real (kind=kind_phys), dimension(:,:), intent(in) :: press,       &
     &       pdfinv, tdav, tstdav, tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2
      real (kind=kind_phys), intent(in) :: rrvco2

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(inout) :: sorc

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,IXPRNLTE,NBCO215) :: fnlte
      real (kind=kind_phys), dimension(NPTS,NLAY,IXPRNLTE)    :: cmtrx

      real (kind=kind_phys), dimension(NPTS,IXPRNLTE) :: ag, az, bdenom,&
     &       cdiag, tcoll, phifx, phivar
      real (kind=kind_phys), dimension(NBCO215) :: c1b7, c2b7, cent, del

      real (kind=kind_phys) :: degen = 0.5
      integer :: n, n1, k, inb, kp

!
!===> ...  begin here
!
!  ---  elements of the source function for bands in the 15 um range

      do n = 1, NBCO215
        n1 = n + 8 + IOFFSET
        cent(n) = 0.5 * ( bdlocm(n1) + bdhicm(n1) )
        del (n) = bdhicm(n1) - bdlocm(n1)
        c1b7(n) = 3.7412e-5 * cent(n)*cent(n)*cent(n) * del(n)
        c2b7(n) = 1.4387 * cent(n)
      enddo

!  ---  compute curtis matrix for both frequency bands.

      call co2curt                                                      &
!  ---  inputs:
     &     ( pdfinv, tdav, tstdav, tlsqu, tmpdiff, a1, a2,              &
!  ---  outputs:
     &       cmtrx                                                      &
     &     )

      do k = 1, IXPRNLTE
        cdiag(:,k) = cmtrx(:,k,k)
      enddo

!  ---  collisional relaxation time (see fels notes for "tcoll")

      do k = 1, IXPRNLTE
        tcoll(:,k) = degen*1.5e-5*press(:,NLP1) / (secday*press(:,k))
      enddo

!  ---  compute NLTE contribution for eack band at each pressure level
!       <= IXPRNLTE. fnlte = zero by assumption at other levels.

      do n = 1, NBCO215
        fnlte(:,1:IXPRNLTE,n) = 3.5 * tcoll(:,1:IXPRNLTE)               &
     &                        * c1b7(n) / (rrvco2*c2b7(n))
      enddo

!  ---  begin computations for (NBCO215) bands in 15um range.

      do inb = 1, NBCO215

        bdenom(:,1:IXPRNLTE) = 1.0 / (1.0 - fnlte(:,1:IXPRNLTE,inb)     &
     &                       * cdiag(:,:))
        phifx (:,:) = f_zero

        do k = 1, IXPRNLTE
          do kp = IXPRNLTE+1, NLAY
            phifx(:,k) = phifx(:,k) + cmtrx(:,kp,k)*sorc(:,kp,inb)
          enddo
        enddo

        az(:,1:IXPRNLTE) = sorc(:,1:IXPRNLTE,inb)                       &
     &                   + fnlte(:,1:IXPRNLTE,inb)*phifx(:,1:IXPRNLTE)

!  ---  first iteration. (J(k) = B(k)) as initial guess)

        phivar(:,1:IXPRNLTE) = f_zero
        do k = 1, IXPRNLTE
          do kp = 1, IXPRNLTE
            phivar(:,k) = phivar(:,k) + cmtrx(:,kp,k)*sorc(:,kp,inb)
          enddo
        enddo

        ag(:,1:IXPRNLTE) = fnlte(:,1:IXPRNLTE,inb)                      &
     &                   * (phivar(:,1:IXPRNLTE) - cdiag(:,1:IXPRNLTE)  &
     &                   * sorc(:,1:IXPRNLTE,inb))

        sorc(:,1:IXPRNLTE,inb) = bdenom(:,1:IXPRNLTE)                   &
     &                         * (az(:,1:IXPRNLTE) + ag(:,1:IXPRNLTE))

!  ---  second iteration.  (J(k) = result of first iteration as guess)

        phivar(:,1:IXPRNLTE) = f_zero
        do k = 1, IXPRNLTE
          do kp = 1, IXPRNLTE
            phivar(:,k) = phivar(:,k) + cmtrx(:,kp,k)*sorc(:,kp,inb)
          enddo
        enddo

        ag(:,1:IXPRNLTE) = fnlte(:,1:IXPRNLTE,inb)                       &
     &                   * (phivar(:,1:IXPRNLTE) - cdiag(:,1:IXPRNLTE)   &
     &                   * sorc(:,1:IXPRNLTE,inb))

        sorc(:,1:IXPRNLTE,inb) = bdenom(:,1:IXPRNLTE)                    &
     &                         * (az(:,1:IXPRNLTE) + ag(:,1:IXPRNLTE))

      enddo   ! end do_inb_loop

!
!...................................
      end subroutine nlte
!----------------------------------- contained in co2_source_calc


!-----------------------------------
      subroutine co2curt                                                &
!...................................
!  ---  inputs:
     &     ( pdfinv, tdav, tstdav, tlsqu, tmpdiff, a1, a2,              &
!  ---  outputs:
     &       cmtrx                                                      &
     &     )

!----------------------------------------------------------------------
!
!     co2curt computes Curtis matrix elements derived from co2
!     transmission functions.
!
!     author: m. d. schwarzkopf
!     revised: 8/18/94       certified:  radiation version 1.0
!
!----------------------------------------------------------------------
!
!  inputs:
!
!  outputs:
!     cmtrx  = cutris matrix.
!
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      real (kind=kind_phys), dimension(:,:), intent(in) :: pdfinv,      &
     &       tdav, tstdav, tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: cmtrx

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: co2row, co2rowp

      integer :: k, krow, kp

!
!===> ...  begin here
!
!  ---  compute co2 transmission functions.

      co2row (:,:) = 1.0
      co2rowp(:,:) = 1.0

!  ---  compute curtis matrix for rows from 1 to ixprnlte

      do k = 1, IXPRNLTE
        krow = k

        call transcol                                                   &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       1, krow, 1, NLP1, NPTS, NLAY, NLP1,                        &
!  ---  outputs:
     &       co2row                                                     &
     &     )

        call transcol                                                   &
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       1, krow+1, 1, NLP1, NPTS, NLAY, NLP1,                      &
!  ---  outputs:
     &       co2rowp                                                    &
     &     )

        do kp = 1, NLAY-1
          cmtrx(:,kp,k) = radcon * pdfinv(:,k)                          &
     &                  * (co2rowp(:,kp) - co2rowp(:,kp+1)              &
     &                  -  co2row(:,kp) + co2row(:,kp+1))
        enddo

        cmtrx(:,NLAY,k) = radcon * pdfinv(:,k)                          &
     &                  * (co2rowp(:,NLAY) - co2row(:,NLAY))
      enddo

!
!...................................
      end subroutine co2curt
!----------------------------------- contained in co2_source_calc
!
!...................................
      end subroutine co2_source_calc
!-----------------------------------



!  =========================================
!  *****        gas_tf section         *****
!  =========================================


!-----------------------------------
      subroutine ptz                                                    &
!...................................
!  ---  inputs:
!    &     ( pd, plm, NLAY, NLP1,                                       &
     &     ( pd, plm, NLAY, NLP1 )
!  ---  outputs: (to module variables)

!---------------------------------------------------------------------
!
!     this program calculates temperatures at up to 200 user specified
!     pressures. it makes use of an analytical function which 
!     approximates  the us standard atm(1976). this is calculated in
!     function 'antemp'.  the form of the analytical function was
!     suggested to me (S.B. Fels) in 1971 by richard s. lindzen.
!
!---------------------------------------------------------------------
!  inputs:
!     pd    - press (mb) for layer boundaries. (also known as flx levs).
!     plm   - press at midpoint of layer (avg of adjacent pd values)
!
!  outputs: (to module variables)
!     stemp - std layer-mean temp with index 1 at the top. appl for nq=5
!     gtemp - pressure coeff
!
!  locals:
!     plmcgs- plm in cgs units. needed for gtemp calc.
!     prsint- same as pd, but with indices reversed (index 1 at sfc).
!     press - press used for quadratures (4 pts.) indices as in prsint.
!     tempquad- temperature at quad. pts. index 1 is at the sfc.
!     tmpint- temperature at quad. pts, saved over both height and
!             quadrature index. values come from tempquad.
!     tmpout- as in temp, but with indices reversed (index 1 at top
!     altquad- height (km) generated by antemp. lowest index = surface.
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, NLP1

      real (kind=kind_phys), dimension(:), intent(in) :: pd, plm

!  ---  outputs: (none)

!  ---  locals:
      real (kind=kind_phys), parameter :: delzap = 0.5

      real (kind=kind_phys), dimension(NLP1) :: pres1, altquad,         &
     &       tempquad, prsint, plmcgs, tmpint

      real (kind=kind_phys) :: dlogp, znint, dz, ht, rk1, rk2, rk3, rk4
      integer :: k, nint1, m

!
!===> ... begin here
!
!  ---  the gtemp code below assumes plmcgs in cgs units

      do k = 1, NLP1
        prsint(k) = plm(NLAY+2-k)
        plmcgs(k) = pd (k)*1.0e+3
      enddo

      do k = 1, NLAY
        gtemp(k) = plmcgs(k)**0.2                                       &
     &           * (1.0 + plmcgs(k)/30000.0)**0.8 / 1013250.0
      enddo
      gtemp(NLP1) = 1.0

      altquad(1) = f_zero
      tempquad(1) = antemp( f_zero )

!  ---  obtain layer-mean quantities by quadrature over the layers.
!       the calc is made  to find the temperature at the layer-mean
!       (plm). calcs are done (oddly!) 1 quad. interval at a time,
!       with each going from the sfc upward to the top layer.

      pres1(1) = prsint(1)
      do k = 2, NLP1
        pres1(k) = pd(NLAY+2-k)
      enddo

!  ---  pres1 is the pressure at the quadrature point; alt and temp
!       are the heights and pressures for each such quad. pt. these
!       are saved as tmpint and a.

      do k = 1, NLAY

!  ---  establish comp levs between user levs at intervals of approx.
!       'delzap' km. special care is needed for the topmost layer,
!       which usually goes to zero pressure.

        dlogp = 7.0 * ALOG( pres1(k)/pres1(k+1) )
        nint1 = dlogp / delzap
        nint1 = nint1 + 1
        znint = nint1

!  ---  the conversion factor is used to convert dz from cm (using the
!       model's values for rgas and grav) to km (as in this program)

        dz = 1.0e-5 * rgas * dlogp / (7.0 * con_amd * grav * znint)
        ht = altquad(k)

!  ---  calculate height at next user level by means of runge-kutta
!       integration.

        do m = 1, nint1
          rk1 = antemp(ht)         * dz
          rk2 = antemp(ht+0.5*rk1) * dz
          rk3 = antemp(ht+0.5*rk2) * dz
          rk4 = antemp(ht+rk3)     * dz
          ht  = ht + 0.16666667 * (rk1+rk2+rk2+rk3+rk3+rk4)
        enddo

        altquad(k+1) = ht
        tempquad(k+1)= antemp(ht)
      enddo

!  ---  save temperature (tmpint) at quad. pts for layer-mean 
!       evaluations by simpsons rule.

      do k = 1, NLP1
        tmpint(k) = tempquad(k)
      enddo

!  ---  stemp is layer-mean temp with index 1 at the top. appl for nq=5

      do k = 1, NLP1
        stemp(k) = tmpint(NLAY+2-k)
      enddo


! ================
      contains
! ================

!-----------------------------------
      real (kind=kind_phys) function antemp                             &
!...................................
!  ---  inputs:
     &     ( z )

!----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      real (kind=kind_phys), intent(in) :: z

!  ---  output: (to antemp)

!  ---  locals:
      real (kind=kind_phys) :: zb(10), delta(10), c(11)
      data zb    /  11.0,    20.0,    32.0,    47.0,    51.0,           &
     &              71.0,    84.8520, 90.0,    91.0,    92.0 /

      data delta/   0.3,     1.0,     1.0,     1.0,     1.0,            &
     &              1.0,     1.0,     1.0,     1.0,     1.0  /

      data c    /  -6.5,     0.0,     1.0,     2.80,    0.0,            &
     &             -2.80,   -2.00,    0.0,     0.0,     0.0,    0.0 /

      real (kind=kind_phys), parameter :: tstar = 288.15

      real (kind=kind_phys) :: temp, expo, x, y, zlog, expp, faclog

      integer :: n, nlast

!
!===> ...  begin here
!
      temp = tstar + c(1) * z
      nlast = 10

      do n = 1, nlast
        expo = (z - zb(n)) / delta(n)
        if (abs(expo) <= 60.0) then
          x = exp( expo )
          y = x + 1.0 / x
          zlog = alog( y )
        else
          zlog = abs( expo )
        endif

        expp = zb(n) / delta(n)
        if (abs(expp) <= 60.0) then
          x = exp( expp )
          y = x + 1.0 / x
          faclog = alog( y )
        else
          faclog = abs( expp )
        endif

        temp = temp + (c(n+1) - c(n)) * 0.5                             &
     &       * ( z + delta(n)*(zlog - faclog) )
      enddo

      antemp = temp
!
!...................................
      end function antemp
!----------------------------------- contained in ptz
!
!...................................
      end subroutine ptz
!-----------------------------------



!-----------------------------------
      subroutine co2coef                                                &
!...................................
!  ---  inputs:
     &     ( press, temp, tflux, pdflux,                                &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       tdav, tstdav, tlsqu, tmpdiff, a1, a2                       &
     &     )

!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: press, temp, &
     &       tflux, pdflux

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: tdav,       &
     &       tstdav, tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(out) :: a1, a2

!  ---  define 4 coeffs (formerly in Id3). b0, b1, b2, b3 are coeffs used
!       to correct for the use of 250k in the planck function used in
!       evaluating planck-weighted co2 transm functions. (see reference(1).)

      real (kind=kind_phys), parameter :: b0 = -0.51926410e-4
      real (kind=kind_phys), parameter :: b1 = -0.18113332e-3
      real (kind=kind_phys), parameter :: b2 = -0.10680132e-5
      real (kind=kind_phys), parameter :: b3 = -0.67303519e-7

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: tdif
      real (kind=kind_phys) :: palog8, alogps8

      integer :: i, k

!
!===> ...  begin here
!
!  ---  compute temperature difference between model profile and
!       standard profile

      do k = 1, NLP1
        tmpdiff(:,k) = temp(:,k) - stemp(k)
      enddo

!  ---  compute weighted temperature difference (tdav) and pressure
!       integrals (tstdav) from level 1 to level NLAY. the quotient
!       will be used to obtain  the difference (dift) between the
!       model temperature profile and the standard profile.

      tstdav(:,1) = f_zero
      tdav  (:,1) = f_zero
      do k = 1, NLAY
        tstdav(:,k+1) = tstdav(:,k) + gtemp(k)*pdflux(:,k)
        tdav  (:,k+1) = tdav  (:,k) + gtemp(k)*pdflux(:,k)*tmpdiff(:,k)
      enddo

!  ---  a logarithmic interpolation is presently assumed, with the 2
!       2nd pressure profile having pressures 0.8* the first, thus
!       accounting for the 0.8 and 0.2 factors. The denominator, which
!       is (log(pstd)-log(0.8*pstd)) is a constant (-log(0.8)) so the
!       expression can be replaced by the quantity palog8.

      if ( iflgco2 == 1 ) then
        alogps8 = alog( pstd*0.8 )
        palog8  = -alog( 0.8 )

        a1(:) = (alog( press(:,NLP1) ) - alogps8) / palog8
        a2(:) = 1.0 - a1(:)

!  ---  evaluate coefficients for co2 pressure interpolation (a1, a2).
!       a linear interpolation is presently assumed, with the 2nd
!       pressure profile having pressures 0.8* the first, thus
!       accounting for the 0.8 and 0.2 factors.

      else
        a1(:) = (press(:,NLP1) - pstd*0.8) / (pstd*0.2)
        a2(:) = (pstd - press(:,NLP1)) / (pstd*0.2)
      endif

!  ---  compute temperature coefficient based on tflux. see fels and
!       schwarzkopf (1981) for details.

      tdif(:,:) = tflux(:,:) - 2.5e+2
      do k = 1, NLP1
        do i = 1, NPTS
          if ( tflux(i,k) <= 2.5e+2 ) then
            tlsqu(i,k) = b0 + tdif(i,k) * (b1 + tdif(i,k)               &
     &                 * (b2 + b3*tdif(i,k)) )
          else
            tlsqu(i,k) = b0
          endif
        enddo
      enddo

!
      return
!...................................
      end subroutine co2coef
!-----------------------------------



!-----------------------------------
      subroutine transfn                                                &
!...................................
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co2nbl, co2spnb                                            &
     &     )

!----------------------------------------------------------------------
!
!     transfn computes the temperature-corrected co2 nearby layer
!     transmission functions.
!
!     author: m. d. schwarzkopf
!     revised: 1/1/93        certified:  radiation version 1.0
!
!--------------------------------------------------------------------
!
!  inputs:
!
!  outputs:
!     co2nbl =  co2 transmission functions (not pressure-integrated)
!               for adjacent levels, over the 560-800 cm-1 range.
!     co2spnb = co2 transmission functions between a flux level and
!               space, for each of (NBCO215) frequency bands in
!               the 15 um range. used for cool-to-space calculations.
!
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: tdav,        &
     &       tstdav, tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: co2spnb
      real (kind=kind_phys), dimension(:,:),   intent(out) :: co2nbl

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1,NBCO215) :: co215nb,   &
     &       co2dt15nb, co2d2t15nb
      real (kind=kind_phys), dimension(NPTS,NLP1) :: dift
      real (kind=kind_phys), dimension(NPTS,NLAY) :: co2m2d,co2md,co2mr

      integer :: inb, k

!
!===> ...  begin here
!
!  ---  perform co2 pressure interpolation on all inputted transmission
!       functions and temperature derivatives successively computing
!       co2r, dco2dt, and d2cdt2.

      do inb = 1, NBCO215
        do k = 1, NLP1
          co215nb   (:,k,inb) = a1(:)*co215nbps1(k,inb)                 &
     &                        + a2(:)*co215nbps8(k,inb)
          co2dt15nb (:,k,inb) = 1.0e-2 * (a1(:)*co2dt15nbps1(k,inb)     &
     &                        +           a2(:)*co2dt15nbps8(k,inb))
          co2d2t15nb(:,k,inb) = 1.0e-3 * (a1(:)*co2d2t15nbps1(k,inb)    &
     &                        +           a2(:)*co2d2t15nbps8(k,inb))
        enddo
      enddo

      do k = 1, NLAY
        co2mr (:,k) = a1(:)*co2m51(k) + a2(:)*co2m58(k)
        co2md (:,k) = 1.0e-2 * (a1(:)*cdtm51(k) + a2(:)*cdtm58(k))
        co2m2d(:,k) = 1.0e-3 * (a1(:)*c2dm51(k) + a2(:)*c2dm58(k))
      enddo

!  ---  perform the temperature interpolation for these transmissivities

      dift(:,2:NLP1) = tdav(:,2:NLP1) / tstdav(:,2:NLP1)

      do inb = 1, NBCO215
        co2spnb(:,1,inb) = 1.0
        co2spnb(:,2:NLP1,inb) = co215nb(:,2:NLP1,inb)                   &
     &           + dift(:,2:NLP1) * ( co2dt15nb(:,2:NLP1,inb)           &
     &           + 0.5*dift(:,2:NLP1)*co2d2t15nb(:,2:NLP1,inb) )

        do k = 1, NLP1
          co2spnb(:,k,inb) = co2spnb(:,k,inb) * (1.0 - tlsqu(:,1))      &
     &                     + tlsqu(:,1)
        enddo
      enddo

!  ---  compute special nearby layer transmission functions for combined
!       band in 15 um range. the transmissivities are not layer-averaged.

      co2nbl(:,1:NLAY) = co2mr(:,1:NLAY)                                &
     &                 + tmpdiff(:,1:NLAY)*( co2md(:,1:NLAY)            &
     &                 + 0.5*tmpdiff(:,1:NLAY)*co2m2d(:,1:NLAY) )

      co2nbl(:,1:NLAY) = co2nbl(:,1:NLAY)                               &
     &                 * (1.0 - tlsqu(:,1:NLAY)) + tlsqu(:,1:NLAY)

!
      return
!...................................
      end subroutine transfn
!-----------------------------------



!-----------------------------------
      subroutine transcol                                               &
!...................................
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       kcol, krow, kcols, kcole, NPTS, NLAY, NLP1,                &
!  ---  outputs:
     &       co21c                                                      &
     &     )

!---------------------------------------------------------------------
!
!     transcol computes temperature-corrected co2 transmission
!     functions at a particular (krow).
!
!     author: c. l. kerr
!     revised: 11/11/93       certified:  radiation version 1.0
!
!---------------------------------------------------------------------
!
!  inputs:
!
!  outputs:
!     co21c  = column of transmission functions.
!
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, krow, kcol, kcols, kcole

      real (kind=kind_phys), dimension(:,:), intent(in) :: tdav, tstdav,&
     &       tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: co21c

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: co2r, dift, d2cdt2,&
     &       dco2dt

      integer :: k, kp

!
!===> ...  begin here
!
      co21c(:,:) = 1.0

      do kp = kcols, kcole
        if ( kp <> krow ) then
          dift(:,kp) = (tdav  (:,kp) - tdav  (:,krow))                   &
     &               / (tstdav(:,kp) - tstdav(:,krow))
        elseif (krow <> 1) then
          dift(:,kp) = 0.5 * (tmpdiff(:,kp) + tmpdiff(:,kp-1))
        else
          dift(:,kp) = f_zero
        endif
      end do

!  ---  obtain transmission functions used for the flux at a fixed level
!       (krow). ie, tf's  from varying flux levels (kp) to (krow)

!  ---  pressure interpolation

      do kp = kcols, kcole
        co2r  (:,kp) = a1(:)*co251(kp,krow) + a2(:)*co258(kp,krow)
        dco2dt(:,kp) = 1.0e-2 * ( a1(:)*cdt51(kp,krow)                  &
     &               +            a2(:)*cdt58(kp,krow) )
        d2cdt2(:,kp) = 1.0e-3 * ( a1(:)*c2d51(kp,krow)                  &
     &               +            a2(:)*c2d58(kp,krow) )
      enddo

!  ---  temperature interpolation

      do kp = kcols, kcole
        co21c (:,kp) = co2r(:,kp) + dift(:,kp)*(dco2dt(:,kp)            &
     &               + 0.5*dift(:,kp)*d2cdt2(:,kp))
      enddo

!  ---  correction for finite width of co2 bands (Eqs. 7a-7c, Ref. (2))

      do kp = kcols, kcole
        co21c(:,kp) = co21c(:,kp)*(1.0 - tlsqu(:,kp)) + tlsqu(:,kp)
      enddo

!
      return
!...................................
      end subroutine transcol
!-----------------------------------



!-----------------------------------
      subroutine transcolrow                                            &
!...................................
!  ---  inputs:
     &     ( tdav, tstdav, tlsqu, tmpdiff, a1, a2,                      &
     &       kcol, krow, kcols, kcole, krows, krowe,                    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, tch4n2oe, n2o9c, tn2o17                      &
     &     )

!-----------------------------------------------------------------------
!
!     transcolrow computes the temperature-corrected co2 transmission
!     functions for a particular (krow) (varying column index) and for
!     a particular (kcol) (varying row index).
!
!     transcolrow also computes the pressure-interpolated ch4 and n2o
!     transmission functions for a particular (krow) and for a parti-
!     cular (kcol). By assumption, no correction for finite bandwidth
!     is performed.
!
!     author: c. l. kerr
!     revised: 11/11/93       certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
!
!  inputs:
!
!  outputs:
!     co21c  = column of transmission functions (fixed krow).
!     co21r  = column of transmission functions (fixed kcol).
!     tch4n2oe =
!     n2o9c  = column of n2o 9 um transmission functions (fixed krow).
!     tn2o17 =
!
!  locals:
!     ch41c  = column of ch4 transmission functions (fixed krow).
!     n2o1c  = column of n2o transmission functions (fixed krow).
!     n2o17c = column of n2o 17 um transmission functions (fixed krow).
!
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1,                          &
     &                       kcol, krow, kcols, kcole, krows, krowe

      real (kind=kind_phys), dimension(:,:), intent(in) :: tdav, tstdav,&
     &       tlsqu, tmpdiff
      real (kind=kind_phys), dimension(:),   intent(in) :: a1, a2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: co21c,      &
     &       co21r, tch4n2oe, n2o9c, tn2o17

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: co2p, dift, d2cdt2,&
     &       dco2dt, ch41c, n2o1c, n2o17c, ch4p, d2ch4dt2, dch4dt,      &
     &       d2n2odt2, dn2odt, d2n2o17dt2, dn2o17dt, d2n2o9dt2,         &
     &       dn2o9dt, n2op, n2o17p, n2o9p

      integer :: kp
!
!===> ...  begin here
!

!  ---  initialization

      co21c(:,:) = 1.0
      co21r(:,:) = 1.0
      ch41c(:,:) = 1.0
      n2o1c(:,:) = 1.0
      n2o17c(:,:)= 1.0
      n2o9c(:,:) = 1.0

!  ---  temperature difference averaged between levels k and kp

      do kp = kcols, kcole
        if ( kp <> krow ) then
          dift(:,kp) = (tdav  (:,kp) - tdav  (:,krow))                  &
     &               / (tstdav(:,kp) - tstdav(:,krow))
        elseif (krow <> 1) then
          dift(:,kp) = 0.5 * (tmpdiff(:,kp) + tmpdiff(:,kp-1))
        else
          dift(:,kp) = f_zero
        endif
      end do

!  ---  obtain transmission functions used for the flux at a fixed level
!       (krow). ie, tf's  from varying flux levels (kp) to (krow)

!  ---  pressure interpolation

      do kp = kcols, kcole
        co2p  (:,kp) = a1(:)*co251(kp,krow) + a2(:)*co258(kp,krow)
        dco2dt(:,kp) = 1.0e-2 * ( a1(:)*cdt51(kp,krow)                  &
     &               +            a2(:)*cdt58(kp,krow) )
        d2cdt2(:,kp) = 1.0e-3 * ( a1(:)*c2d51(kp,krow)                  &
     &               +            a2(:)*c2d58(kp,krow) )

        if ( ifch4n2o == 1 ) then

          ch41c(:,kp) = a1(:)*ch451(kp,krow) + a2(:)*ch458(kp,krow)
          n2o1c(:,kp) = a1(:)*n2o51(kp,krow) + a2(:)*n2o58(kp,krow)
          n2o17c(:,kp)= a1(:)*n2o71(kp,krow) + a2(:)*n2o78(kp,krow)
          n2o9c(:,kp) = a1(:)*n2o91(kp,krow) + a2(:)*n2o98(kp,krow)

        else if ( ifch4n2o == 2 ) then

          ch4p (:,kp) = a1(:)*ch451(kp,krow) + a2(:)*ch458(kp,krow)
          n2op (:,kp) = a1(:)*n2o51(kp,krow) + a2(:)*n2o58(kp,krow)
          n2o17p(:,kp)= a1(:)*n2o71(kp,krow) + a2(:)*n2o78(kp,krow)
          n2o9p(:,kp) = a1(:)*n2o91(kp,krow) + a2(:)*n2o98(kp,krow)

          dch4dt(:,kp)    = 1.0e-2 * ( a1(:)*ch4dt51(kp,krow)           &
     &                    +            a2(:)*ch4dt58(kp,krow) )
          dn2odt(:,kp)    = 1.0e-2 * ( a1(:)*n2odt51(kp,krow)           &
     &                    +            a2(:)*n2odt58(kp,krow) )
          dn2o17dt(:,kp)  = 1.0e-2 * ( a1(:)*n2odt71(kp,krow)           &
     &                    +            a2(:)*n2odt78(kp,krow) )
          dn2o9dt(:,kp)   = 1.0e-2 * ( a1(:)*n2odt91(kp,krow)           &
     &                    +            a2(:)*n2odt98(kp,krow) )
          d2ch4dt2(:,kp)  = 1.0e-3 * ( a1(:)*ch4d2t51(kp,krow)          &
     &                    +            a2(:)*ch4d2t58(kp,krow) )
          d2n2odt2(:,kp)  = 1.0e-3 * ( a1(:)*n2od2t51(kp,krow)          &
     &                    +            a2(:)*n2od2t58(kp,krow) )
          d2n2o17dt2(:,kp)= 1.0e-3 * ( a1(:)*n2od2t71(kp,krow)          &
     &                    +            a2(:)*n2od2t78(kp,krow) )
          d2n2o9dt2(:,kp) = 1.0e-3 * ( a1(:)*n2od2t91(kp,krow)          &
     &                    +            a2(:)*n2od2t98(kp,krow) )
        endif
      enddo

!  ---  temperature interpolation

      do kp = kcols, kcole
        co21c (:,kp) = co2p(:,kp) + dift(:,kp)*( dco2dt(:,kp)           &
     &               + 0.5*dift(:,kp)*d2cdt2(:,kp) )

        if ( ifch4n2o == 2 ) then
          ch41c (:,kp) = ch4p(:,kp) + dift(:,kp)*( dch4dt(:,kp)         &
     &                 + 0.5*dift(:,kp)*d2ch4dt2(:,kp) )
          n2o1c (:,kp) = n2op(:,kp) + dift(:,kp)*( dn2odt(:,kp)         &
     &                 + 0.5*dift(:,kp)*d2n2odt2(:,kp) )
          n2o17c(:,kp) = n2o17p(:,kp) + dift(:,kp)*( dn2o17dt(:,kp)     &
     &                 + 0.5*dift(:,kp)*d2n2o17dt2(:,kp) )
          n2o9c (:,kp) = n2o9p(:,kp) + dift(:,kp)*( dn2o9dt(:,kp)       &
     &                 + 0.5*dift(:,kp)*d2n2o9dt2(:,kp) )
        endif
      enddo

!  ---  correction for finite width of co2 bands (Eqs. 7a-7c, Ref. (2))

      do kp = kcols, kcole
        co21c(:,kp) = co21c(:,kp)*(1.0 - tlsqu(:,kp)) + tlsqu(:,kp)
      enddo

!  ---  obtain transmission functions used for the flux for varying
!       levels (krow) from a fixed level (kcol). ie, tf's  from a fixed
!       flux level (kcol) to varying levels (krow).

!  ---  temperature difference averaged between levels k and kp. This
!       computation is made unless krow = kcol, and range (krows,krowe)
!       is entirely within (kcols,kcole), in which case the dift
!       computed for column tfs is applicable to row tfs.

      if ( kcol <> krow .or. krows < kcols .or. krowe > kcole ) then
        do kp = krows, krowe
          if ( kp <> krow ) then
            dift(:,kp) = (tdav  (:,kp) - tdav  (:,krow))                &
     &                 / (tstdav(:,kp) - tstdav(:,krow))
          else if ( krow <> 1 ) then
            dift(:,kp) = 0.5 * (tmpdiff(:,kp) + tmpdiff(:,kp-1))
          else
            dift(:,kp) = f_zero
          endif
        enddo
      endif

!  ---  pressure interpolation

      do kp = krows, krowe
        co2p  (:,kp) = a1(:)*co251(kcol,kp) + a2(:)*co258(kcol,kp)
        dco2dt(:,kp) = 1.0e-2 * ( a1(:)*cdt51(kcol,kp)                  &
     &               +            a2(:)*cdt58(kcol,kp) )
        d2cdt2(:,kp) = 1.0e-3 * ( a1(:)*c2d51(kcol,kp)                  &
     &               +            a2(:)*c2d58(kcol,kp) )

        if ( ifch4n2o == 2 ) then
          ch4p      (:,kp) = a1(:)*ch451(kcol,kp) + a2(:)*ch458(kcol,kp)
          n2op      (:,kp) = a1(:)*n2o51(kcol,kp) + a2(:)*n2o58(kcol,kp)
          n2o17p    (:,kp) = a1(:)*n2o71(kcol,kp) + a2(:)*n2o78(kcol,kp)
          n2o9p     (:,kp) = a1(:)*n2o91(kcol,kp) + a2(:)*n2o98(kcol,kp)

          dch4dt    (:,kp) = 1.0e-2 * ( a1(:)*ch4dt51(kcol,kp)          &
     &                     +            a2(:)*ch4dt58(kcol,kp) )
          dn2odt    (:,kp) = 1.0e-2 * ( a1(:)*n2odt51(kcol,kp)          &
     &                     +            a2(:)*n2odt58(kcol,kp) )
          dn2o17dt  (:,kp) = 1.0e-2 * ( a1(:)*n2odt71(kcol,kp)          &
     &                     +            a2(:)*n2odt78(kcol,kp) )
          dn2o9dt   (:,kp) = 1.0e-2 * ( a1(:)*n2odt91(kcol,kp)          &
     &                     +            a2(:)*n2odt98(kcol,kp) )
          d2ch4dt2  (:,kp) = 1.0e-3 * ( a1(:)*ch4d2t51(kcol,kp)         &
     &                     +            a2(:)*ch4d2t58(kcol,kp) )
          d2n2odt2  (:,kp) = 1.0e-3 * ( a1(:)*n2od2t51(kcol,kp)         &
     &                     +            a2(:)*n2od2t58(kcol,kp) )
          d2n2o17dt2(:,kp) = 1.0e-3 * ( a1(:)*n2od2t71(kcol,kp)         &
     &                     +            a2(:)*n2od2t78(kcol,kp) )
          d2n2o9dt2 (:,kp) = 1.0e-3 * ( a1(:)*n2od2t91(kcol,kp)         &
     &                     +            a2(:)*n2od2t98(kcol,kp) )
        endif
      enddo

!  ---  temperature interpolation

      do kp = krows, krowe
        co21r(:,kp) = co2p(:,kp) + dift(:,kp) * ( dco2dt(:,kp)          &
     &              + 0.5*dift(:,kp)*d2cdt2(:,kp) )
      enddo

!  ---  correction for finite width of co2 bands (Eqs. 7a-7c, Ref. (2))

      do kp = krows, krowe
        co21r(:,kp) = co21r(:,kp)*(1.0 - tlsqu(:,kcol)) + tlsqu(:,kcol)
      enddo

!  ---  tn2o17 results are for 560-630 cm-1 band. (if NBCO215=3)

      tch4n2oe(:,:) = ch41c(:,:) * n2o1c(:,:)
      tn2o17  (:,:) = n2o17c(:,:)

!
      return
!...................................
      end subroutine transcolrow
!-----------------------------------



!-----------------------------------
      subroutine trans_nearby                                           &
!...................................
!  ---  inputs:
     &     ( press, pflux, overod, pdflux, pdfinv, co2nbl,              &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       co21c, co21r, co21diag                                     &
     &     )

!-------------------------------------------------------------------
!
!     compute "nearby  layer" transmission functions at level k
!     ( tau(p(k),p(k))) in the frequency band at 15 um. include all
!     gases (co2, h2o, h2o cont) used in computing fluxes in this band.
!     the algorithm assumes that at pressures (p') near the pressure
!     at level k (p(k)), the transmission function may be written as:
!
!              tau(p',p(k)) = EXP(-alpha*SQRT(p'-p(k)))
!
!     with alpha determined by the boundary condition(s) tau(p(k+1),p(k))
!     and tau(p(k-1),p(k)) = the values from "normal" calculations.
!     An integration is performed over the "nearby" pressure layer to
!     obtain tau(p(k),p(k)). the computation is not done for levels
!     from 1 to KMINH2O-1 (if different), where it is assumed that the
!     h2o transmissivities are near unity, and that the precomputed co2
!     transmissivities may be used.
!
!     two "special case" transmissivities, viz., tau(p(NLAY),p(NLP1))
!     and tau(p(NLP1),p(NLAY)) are also evaluated using the above
!     assumptions and an integration.
!
!-------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: press, pflux,&
     &       overod, pdflux, pdfinv, co2nbl

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: co21c,      &
     &       co21r, co21diag

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: alpa, alpb, ca, cb,&
     &       rlog, delpr1, delpr2

      integer :: k, km, kmp1

!
!===> ...  begin here
!
      km  = max( IXPRKMINH2O-1, 1 )
      kmp1= max( IXPRKMINH2O-1, 2 )

      rlog  (:,km  :NLAY) = log( co2nbl(:,km:NLAY)*overod(:,km:NLAY) )
      delpr1(:,kmp1:NLAY) = pdfinv(:,kmp1:NLAY)                         &
     &                    * ( press(:,kmp1:NLAY) - pflux(:,kmp1:NLAY) )
      alpb  (:,kmp1:NLAY) = -sqrt( delpr1(:,kmp1:NLAY) )                &
     &                    * rlog(:,kmp1:NLAY)
      delpr2(:,kmp1:NLP1) = pdfinv(:,kmp1-1:NLAY)                       &
     &                    * ( pflux(:,kmp1:NLP1)-press(:,kmp1-1:NLAY) )
      alpa  (:,km  :NLAY) = -sqrt( delpr2(:,km+1:NLP1) )                &
     &                    * rlog(:,km:NLAY)
      alpa  (:,NLP1)      = -rlog(:,NLAY)
      alpb  (:,NLP1)      = -rlog(:,NLAY) * sqrt( pdfinv(:,NLAY)        &
     &                    * ( pflux(:,NLP1) - press(:,NLAY-1) ))

      ca(:,km  :NLP1) = alpa(:,km:NLP1) * ( -0.66667 + alpa(:,km:NLP1)  &
     &                * ( 0.25 + alpa(:,km:NLP1)*(-0.066667) ))
      cb(:,kmp1:NLP1) = alpb(:,kmp1:NLP1)*( -0.66667 + alpb(:,kmp1:NLP1)&
     &                * ( 0.25 + alpb(:,kmp1:NLP1) * (-0.066667) ))

      do k = IXPRKMINH2O, NLAY
        co21diag(:,k) = 1.0 + 0.5*(cb(:,k) + ca(:,k-1))
      enddo

      co21diag(:,NLP1) = 1.0 + ca(:,NLAY)

      co21c(:,NLP1) = 1.0 + ( pdflux(:,NLAY)*ca(:,NLP1)                 &
     &               - (press(:,NLAY) - pflux(:,NLAY))*cb(:,NLAY) )     &
     &               / (pflux(:,NLP1) - press(:,NLAY))

      co21r(:,NLP1) = 1.0                                               &
     &               + ((pflux(:,NLP1) - press(:,NLAY-1))*cb(:,NLP1)    &
     &               -  (pflux(:,NLP1) - press(:,NLAY)  )*ca(:,NLAY ))  &
     &               / (press(:,NLAY) - press(:,NLAY-1))

!
      return
!...................................
      end subroutine trans_nearby
!-----------------------------------



!  =========================================
!  *****     cool_to_space section     *****
!  =========================================


!-----------------------------------
      subroutine cool_to_space_exact                                    &
!...................................
!  ---  inputs:
     &     ( cldtf, press, temp, to3cnt, pflux, pdfinv,                 &
     &       dte1, ixoe1, co2spnb, n2o9c, tn2o17, sorc,                 &
     &       tfac, rh2os, rfrgn, wk, xch2obd, totvo2, var1, var2,       &
     &       totf11, totf12, totf113, totf22, tauaertot,                &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       gxcts, gxctscf, cts_sum, cts_sumcf,                        &
     &       exctsn, fctsg, excts                                       &
     &     )

!---------------------------------------------------------------------
!
!     cool_to_space calculates the cool-to-space cooling rate for
!     a band n.
!
!     author: m. d. schwarzkopf
!     revised: 7/21/94
!     certified:  radiation version 2.0
!
!-----------------------------------------------------------------------
!
!   intent in:
!     cldtf    = cloud transmission function between levels k and
!                  level 1.
!     co2spnb  =  co2 transmission functions between levels k and
!                  level 1 for the freq bands in the 15 um range.
!                 the first band may include n2o 17um transmissivities.
!     pdfinv   =  inverse of pressure difference between flux levels.
!     pflux    =  pressure at flux levels of model.
!     press    =  pressure at data levels of model.
!     sorc     =   band-integrated Planck function, for each combined
!                band in the 160-1200 cm-1 region.
!     tfac     =
!     rh2os    =
!     rfrgn    =
!     wk       =
!     xch2obd  = 
!     temp     =  temperature at data levels of model.
!     to3cnt   =   transmission functions between levels k and
!                   level 1 for the 990-1070 cm-1 range.
!     totvo2   =  summed h2o continuum path from top of atmosphere to
!                   flux level k.
!     var1     =  h2o optical path in model layers.
!     var2     =  pressure-weighted h2o optical path in model layers.
!     vvj      =  frequencies for calculation of h2o coefficients
!     apcm, bpcm   =   capphi coefficients for NBLY bands.
!     atpcm, btpcm =   cappsi coefficients for NBLY bands.
!     acomb        =   random "a" parameter for NBLY bands.
!     bcomb        =   random "b" parameter for NBLY bands.
!
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension (:,:,:), intent(in) :: cldtf,    &
     &       co2spnb, sorc, xch2obd, tauaertot
      real (kind=kind_phys), dimension (:,:),   intent(in) :: press,    &
     &       temp, to3cnt, pflux, pdfinv, dte1, n2o9c, tn2o17, tfac,    &
     &       rh2os, rfrgn, wk, totvo2, var1, var2, totf11, totf12,      &
     &       totf113, totf22

      integer, dimension(:,:), intent(in) :: ixoe1

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: exctsn
      real (kind=kind_phys), dimension(:,:),   intent(out) :: cts_sum,  &
     &       cts_sumcf, fctsg, excts
      real (kind=kind_phys), dimension(:), intent(out) :: gxcts, gxctscf

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1,NBCO215) :: co2spnb_d
      real (kind=kind_phys), dimension(NPTS,NLP1) :: sorc_tmp, ctmp,    &
     &       totch2o_tmp
      real (kind=kind_phys), dimension(NPTS,NLAY) :: exctscf, tt, x, y, &
     &       topm, topphi, phitmp, psitmp, ag, agg, f, ff, fac1, fac2,  &
     &       tmp1, tmp2, cfc_tf
      real (kind=kind_phys), dimension(NPTS)      :: gxctsbd, gxctsbdcf,&
     &       pfac1, pfac2
      real (kind=kind_phys), dimension(NPTS,NBLY) :: fctsgcf
      real (kind=kind_phys), dimension(NPTS,NLAY,NBLY) :: exctsncf

      integer :: n, k, n1

!
!===> ...  begin here
!
!  ---  initialize quantities.

      excts(:,:) = f_zero
      gxcts(:)   = f_zero

      exctscf(:,:) = f_zero
      gxctscf(:)   = f_zero

!  ---  retrieve co2spnb, tn2o17 and n2o9c from co2_tf_mod.

      co2spnb_d(:,:,1)   = co2spnb(:,:,1)*tn2o17(:,:)
      co2spnb_d(:,:,2:3) = co2spnb(:,:,2:3)

!   ---  compute temperature quantities.

      x(:,:) = temp(:,:) - 2.5e2
      y(:,:) = x(:,:) * x(:,:)
      ctmp(:,1) = 1.0

      lab_do_n : do n = 1, NBLY-1

!  ---  obtain temperature correction capphi, cappsi, then multiply
!       by optical path var1, var2 to compute temperature-corrected
!       optical path and mean pressure for a layer: phitmp, psitmp.

        f  (:,:) = 0.44194e-1 * (apcm (n)*x(:,:) + bpcm (n)*y(:,:))
        ff (:,:) = 0.44194e-1 * (atpcm(n)*x(:,:) + btpcm(n)*y(:,:))
        ag (:,:) = (1.418191e0 + f (:,:)) * f (:,:) + 1.0e0
        agg(:,:) = (1.418191e0 + ff(:,:)) * ff(:,:) + 1.0e0

        phitmp(:,:) = var1(:,:) * ((((ag (:,:)*ag (:,:))**2)**2)**2)
        psitmp(:,:) = var2(:,:) * ((((agg(:,:)*agg(:,:))**2)**2)**2)

!  ---  obtain optical path and mean pressure from the top of the
!       atmosphere to the level k.

        topm  (:,1) = phitmp(:,1)
        topphi(:,1) = psitmp(:,1)

        do k = 2, NLAY
          topm  (:,k) = topm  (:,k-1) + phitmp(:,k)
          topphi(:,k) = topphi(:,k-1) + psitmp(:,k)
        enddo

!  ---  tt is the cloud-free cool-to-space transmission function.

        fac1(:,:) = acomb(n)  * topm(:,:)
        fac2(:,:) = fac1(:,:) * topm(:,:) / (bcomb(n)*topphi(:,:))
        tmp1(:,:) = fac1(:,:) / SQRT(1.0e0 + fac2(:,:))

!  ---  if use ckd2p1, calculation for combined bands 1-24.
!       else,          calculation for combined bands 1-4.

        if ( n >= band_no_str(1) .and. n <= band_no_end(1) ) then

          call get_totch2o                                              &
!  ---  inputs:
     &     ( n, vvj, dte1, ixoe1, tfac, rh2os, rfrgn, wk,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       totch2o_tmp                                                &
     &     )

          tt(:,1:NLAY) = exp( -1.0 * (tmp1(:,1:NLAY) + diffac           &
     &                 * totch2o_tmp(:,2:NLP1)) )

!  ---  if use ckd2p1, calculation for combined bands 25-40.
!       else,          calculation for combined bands 5-8.

        elseif (n >= band_no_str(2) .and. n <= band_no_end(2)) then

          call get_totch2o                                              &
!  ---  inputs:
     &     ( n, vvj, dte1, ixoe1, tfac, rh2os, rfrgn, wk,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       totch2o_tmp                                                &
     &     )

          tt(:,1:NLAY) = exp( -1.0 * (tmp1(:,1:NLAY) + diffac           &
     &                 * totch2o_tmp(:,2:NLP1)) )

!  ---  if use ckd2p1, calculation for combined band 41. (first co2 band)
!       else,          calculation for combined band 9.  (first co2 band)

        elseif ( n >= band_no_str(3) .and. n <= band_no_end(3) ) then

          totch2o_tmp(:,1) = f_zero
          do k = 2, NLP1
            totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,n-40)
          enddo

          tmp2(:,1:NLAY) = tmp1(:,1:NLAY) + diffac                      &
     &                   * totch2o_tmp(:,2:NLP1)

          if ( iaerlw > 0 ) then
!           call get_totaerooptdep(1, totaer_tmp)
!           tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
            if ( NBDLW == 1 ) then
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            else                                    ! to be developed !!
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            endif
          endif

          tt(:,1:NLAY) = exp( -1.0e0 * tmp2(:,1:NLAY) )                 &
     &                 * co2spnb_d(:,2:NLP1,1)

!  ---  if use ckd2p1, calculation for combined band 42. (second co2 band)
!       else,          calculation for combined band 10. (second co2 band)

        elseif ( n >= band_no_str(4) .and. n <= band_no_end(4) ) then

          totch2o_tmp(:,1) = f_zero
          do k = 2, NLP1
            totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,n-40)
          enddo

          tmp2(:,1:NLAY) = tmp1(:,1:NLAY) + diffac                      &
     &                   * totch2o_tmp(:,2:NLP1)

          if ( iaerlw > 0 ) then
!           call get_totaerooptdep(2, totaer_tmp)
!           tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
            if ( NBDLW == 1 ) then
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            else                                  ! to be developed !!
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            endif
          endif

          tt(:,1:NLAY) = exp( -1.0 * tmp2(:,1:NLAY) )                   &
     &                 * co2spnb_d(:,2:NLP1,2)

!  ---  if use ckd2p1, calculation for combined band 43. (third co2 band)
!       else,          calculation for combined band 11. (third co2 band)

        elseif ( n >= band_no_str(5) .and. n <= band_no_end(5) ) then

          totch2o_tmp(:,1) = f_zero
          do k = 2, NLP1
            totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,n-40)
          enddo

          tmp2(:,1:NLAY) = tmp1(:,1:NLAY) + diffac                      &
     &                   * totch2o_tmp(:,2:NLP1)

          if ( iaerlw > 0 ) then
!           call get_totaerooptdep(3, totaer_tmp)
!           tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
            if ( NBDLW == 1 ) then
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            else                                  ! to be developed !!
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            endif
          endif

          tt(:,1:NLAY) = exp( -1.0 * tmp2(:,1:NLAY) )                   &
     &                 * co2spnb_d(:,2:NLP1,3)

!  ---  if use ckd2p1, calculation for combined bands 44-45.
!       else,          calculation for combined bands 12-13.

        elseif ( n >= band_no_str(6) .and. n <= band_no_end(6) ) then

          totch2o_tmp(:,1) = f_zero
          do k = 2, NLP1
            totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,n-40)
          enddo

          tmp2(:,1:NLAY) = tmp1(:,1:NLAY) + diffac                      &
     &                   * totch2o_tmp(:,2:NLP1)

          if ( iaerlw > 0 ) then
!           call get_totaerooptdep(n-8-IOFFSET, totaer_tmp)
!           tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
            if ( NBDLW == 1 ) then
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            else                                  ! to be developed !!
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            endif
          endif

          tt(:,1:NLAY) = exp( -1.0 * tmp2(:,1:NLAY) )

!  ---  if use ckd2p1, calculation for combined band 46.
!       else,          calculation for combined band 14.

        elseif ( n >= band_no_str(7) .and. n <= band_no_end(7) ) then

          tt(:,1:NLAY) = exp(-1.0*tmp1(:,1:NLAY)) * to3cnt(:,1:NLAY)

!  ---  if use ckd2p1, calculation for combined band 47.
!       else,          calculation for combined band 15.

        elseif (n >= band_no_str(8) .and. n <= band_no_end(8)) then

          totch2o_tmp(:,1) = f_zero
          do k = 2, NLP1
            totch2o_tmp(:,k) = totch2o_tmp(:,k-1)+xch2obd(:,k-1,n-40)
          enddo

          tmp2(:,1:NLAY) = tmp1(:,1:NLAY) + diffac                      &
     &                   * totch2o_tmp(:,2:NLP1)

          if ( iaerlw > 0 ) then
!           call get_totaerooptdep(7, totaer_tmp)
!           tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + totaer_tmp(:,2:NLP1)
            if ( NBDLW == 1 ) then
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            else                                  ! to be developed !!
              tmp2(:,1:NLAY) = tmp2(:,1:NLAY) + tauaertot(:,2:NLP1,1)
            endif
          endif

          tt(:,1:NLAY) = exp( -1.0*tmp2(:,1:NLAY)) * n2o9c(:,2:NLP1)

        endif

!  ---  calculate or retrieve the source function for the current band.

        if ( n <= 8+IOFFSET ) then

          call looktab                                                  &
!  ---  inputs:
     &     ( tabsr, ixoe1, dte1, 1, NLP1, n,                            &
!  ---  outputs:
     &       sorc_tmp                                                   &
     &     )

        else
          n1 = n - (8 + IOFFSET)
          sorc_tmp(:,:) = sorc(:,:,n1)
        endif

!  ---  retrieve the cfc effect if cfcs are activated.

        if ( iflgcfc == 1 .and. n >= 9+IOFFSET ) then
          n1 = n - 8 - IOFFSET
          cfc_tf(:,1:NLAY) = 1.0                                        &
     &                     - CFC_basic%strf113(n1)*totf113(:,2:NLP1)    &
     &                     - CFC_basic%strf22 (n1)*totf22 (:,2:NLP1)    &
     &                     - CFC_basic%strf11 (n1)*totf11 (:,2:NLP1)    &
     &                     - CFC_basic%strf12 (n1)*totf12 (:,2:NLP1)

          tt(:,1:NLAY) = tt(:,1:NLAY)*cfc_tf(:,1:NLAY)
        endif

!  ---  define some near-surface pressure functions that are needed

        pfac1(:) = 0.5 * pdfinv(:,NLAY)                                 &
     &           * (pflux(:,NLP1)-press(:,NLAY)) * tt(:,NLAY-1)
        pfac2(:) = 0.5 * pdfinv(:,NLAY)                                 &
     &           * (pflux(:,NLP1)+press(:,NLAY)-2.0*pflux(:,NLAY))      &
     &           * tt(:,NLAY)

!  ---  calculate the ground fluxes (?)

        gxctsbdcf(:) = tt(:,NLAY) * sorc_tmp(:,NLAY)                    &
     &               + (pfac1(:) + pfac2(:))                            &
     &               * (sorc_tmp(:,NLP1) - sorc_tmp(:,NLAY))

        gxctsbd(:) = cldtf(:,NLP1,cld_indx_table(n+32-IOFFSET))         &
     &             * gxctsbdcf(:)
        gxctscf(:) = gxctscf(:) + gxctsbdcf(:)
        gxcts(:) = gxcts(:) + gxctsbd(:)

        fctsg(:,n) = gxctsbd(:)
        fctsgcf(:,n) = gxctsbdcf(:)

!  ---  include the effect of the cloud transmission function

        ctmp(:,2:NLP1) = tt(:,1:NLAY)                                   &
     &                  * cldtf(:,2:NLP1,cld_indx_table(n+32-IOFFSET))

!  ---  if diagnostics is on, save each band's contribution separately.
!       exctsn is the cool-to-space heating rate for each frequency
!       band. fctsg is the "exact" surface flux for each frequency
!       band in the 160-1200 cm-1 range.

        exctsn(:,1:NLAY,n) = sorc_tmp(:,1:NLAY)                         &
     &                     * (ctmp(:,2:NLP1) - ctmp(:,1:NLAY))

        exctsncf(:,1,n) =  sorc_tmp(:,1) * (tt(:,1) - 1.0e0)
        exctsncf(:,2:NLAY,n) = sorc_tmp(:,2:NLAY)                       &
     &                       * (tt(:,2:NLAY) - tt(:,1:NLAY-1))

!  ---  excts is the cool-to-space cooling rate accumulated over
!       frequency bands.

        excts(:,1:NLAY) = excts(:,1:NLAY) + sorc_tmp(:,1:NLAY)          &
     &                  * (ctmp(:,2:NLP1) - ctmp(:,1:NLAY))

        exctscf(:,1) = exctscf(:,1)                                     &
     &               + sorc_tmp(:,1) * (tt(:,1) - 1.0e0)
        exctscf(:,2:NLAY) = exctscf(:,2:NLAY) + sorc_tmp(:,2:NLAY)      &
     &                    * (tt(:,2:NLAY) - tt(:,1:NLAY-1))

      enddo  lab_do_n

!  ---  gxcts is the "exact" surface flux accumulated over the
!       frequency bands in the 160-1200 cm-1 range. obtain cool-to-
!       space flux at the top by integration of heating rates and
!       using cool-to-space flux at the bottom (current value of
!       gxcts). note that the pressure quantities and conversion
!       factors have not been included either in excts or in gxcts.
!       these cancel out, thus reducing computations.

      do k = 1, NLAY
        gxcts(:) = gxcts(:) - excts(:,k)
      enddo

      do k=1,NLAY
        gxctscf(:) = gxctscf(:) - exctscf(:,k)
      enddo

!  ---  now scale the cooling rate excts by including the pressure
!       factor pdfinv and the conversion factor radcon.

      do n = 1, NBLY-1
        exctsn  (:,1:NLAY,n) = exctsn(:,1:NLAY,n) * radcon              &
     &                       * pdfinv(:,1:NLAY)
        exctsncf(:,1:NLAY,n) = exctsncf(:,1:NLAY,n) * radcon            &
     &                       * pdfinv(:,1:NLAY)
      enddo

      excts  (:,1:NLAY) = excts(:,1:NLAY)*radcon*pdfinv(:,1:NLAY)
      exctscf(:,1:NLAY) = exctscf(:,1:NLAY)*radcon*pdfinv(:,1:NLAY)

!  ---  this is the end of the exact cool-to-space computations; at this
!       point excts has its appropriate value.

!  ---  save the heating rates to be later sent to longwave_driver_mod.

      cts_sum  (:,:) = excts(:,:)
      cts_sumcf(:,:) = exctscf(:,:)

!
      return
!...................................
      end subroutine cool_to_space_exact
!-----------------------------------



!-----------------------------------
      subroutine cool_to_space_approx                                   &
!...................................
!  ---  inputs:
     &     ( index, source, trans, iof, cld_trans,                      &
     &       pdfinv, NPTS, NLAY, NLP1,                                  &
!  ---  in/outputs:
     &       cts_sum, cts_sumcf,                                        &
!  ---  outputs:
     &       cts_out                                                    &
!  ---  optional arguments:
     &,      trans2, iof2                                               &
     &     )

!---------------------------------------------------------------------
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, iof, index
      real (kind=kind_phys), dimension(:,:), intent(in) :: source,      &
     &       trans, cld_trans, pdfinv

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: cts_sum,  &
     &       cts_sumcf

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out)   :: cts_out

!  ---  optional arguments:
      integer,               intent(in), optional :: iof2
      real (kind=kind_phys), intent(in), optional :: trans2(:,:)

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: cts_outcf
      integer :: j

!
!===> ...  begin here
!

      if ( present(trans2) ) then
        cts_out(:,1) = radcon * pdfinv(:,1) * source(:,1)               &
     &               * ( trans(:,2-iof)*cld_trans(:,2)                  &
     &               -   trans2(:,1   )*cld_trans(:,1) )

        cts_out(:,2:NLAY) = radcon*pdfinv(:,2:NLAY)*source(:,2:NLAY)    &
     &             * ( trans(:,3-iof:NLP1-iof)*cld_trans(:,3:NLP1)      &
     &             -   trans2(:,2-iof2 :NLAY-iof2)*cld_trans(:,2:NLAY) )
      else
        cts_out(:,1) = radcon * pdfinv(:,1) * source(:,1)               &
     &               * ( trans(:,2-iof)*cld_trans(:,2) - 1.0 )

        cts_out(:,2:NLAY) = radcon*pdfinv(:,2:NLAY)*source(:,2:NLAY)    &
     &             * ( trans(:,3-iof:NLP1-iof)*cld_trans(:,3:NLP1)      &
     &             -   trans(:,2-iof :NLAY-iof)*cld_trans(:,2:NLAY) )
      endif


      if ( present(trans2) ) then
        cts_outcf(:,1) = radcon * pdfinv(:,1) * source(:,1)             &
     &                 * ( trans(:,2-iof) - trans2(:,1-iof2) )

        cts_outcf(:,2:NLAY) = radcon*pdfinv(:,2:NLAY)*source(:,2:NLAY)  &
     &                 * ( trans(:,3-iof:NLP1-iof)                      &
     &                 -   trans2(:,2-iof2:NLAY-iof2) )
      else
        cts_outcf(:,1) = radcon * pdfinv(:,1) * source(:,1)             &
     &                 * ( trans(:,2-iof) -  1.0 )

        cts_outcf(:,2:NLAY) = radcon*pdfinv(:,2:NLAY)*source(:,2:NLAY)  &
     &                 * ( trans(:,3-iof:NLP1-iof)                      &
     &                 -   trans(:,2-iof:NLAY-iof) )
      endif

!  ---  cts_sum is the sum of the values from cool_to_space_exact and
!       the values defined here in cool_to_space_approx. it will be
!       used by longwave_driver_mod.

      cts_sum  (:,:) = cts_sum  (:,:) - cts_out  (:,:)
      cts_sumcf(:,:) = cts_sumcf(:,:) - cts_outcf(:,:)

!
      return
!...................................
      end subroutine cool_to_space_approx
!-----------------------------------



!-----------------------------------
      subroutine e1e290                                                 &
!  ---  inputs:
     &     ( tab1, tab2, tab1a, tab2a, tab1w,                           &
     &       dte1, ixoe1, dte2, ixoe2, mass_1, avephi,                  &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss               &
!  --- optional in/out:
     &,      avephif, e1cts1f, e1cts2f, e1flxf, emissf                  &
     &     )

!-----------------------------------------------------------------------
!
!     E1e290 computes two different quantities.
!
!     1) emissivities used to compute the exchange terms for flux at the
!     top of the atmosphere (level 1). (the top layer, isothermal by
!     assumption, does not contribute to photon exchanges with other
!     layers). these terms are obtained using precomputed e2 functions
!     (see ref. (2)).
!
!     2) emissivities used to obtain the cool-to-space heating rates
!     for all pressure layers. these are obtained using precomputed
!     e1 functions (see ref. (2)).
!
!     the frequency ranges for the e2 calculations are 0-560 and 1200-
!     2200 cm-1. the CTS calculations also require calculations in the
!     160-560 cm-1 range. (see refs. (1) and (2)).
!ifdef ch4n2o
!
!     if ch4 and n2o are included, the frequency range for emissivities
!     is 1400-2200 cm-1, with separate emissivity quantities for the
!     1200-1400 cm-1 range.
!endif ch4n2o
!
!     the reason for combining these calculations is that both use
!     the same scaled h2o amount (avephi) as input, thus reducing
!     some calculation time for obtaining index quantities.
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      type (tab1_type), intent(in) :: tab1, tab1a, tab1w, tab2, tab2a
      type (axis_type), intent(in) :: mass_1

      real (kind=kind_phys), dimension(:,:), intent(in) :: dte1, dte2,  &
     &       avephi
      integer, dimension(:,:), intent(in) :: ixoe1, ixoe2

      real (kind=kind_phys), dimension(:,:), optional, intent(in) ::    &
     &       avephif

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: e1cts1,     &
     &       e1cts2, e1ctw1, e1ctw2, e1flx, emiss

      real (kind=kind_phys), dimension(:,:), optional, intent(out) ::   &
     &       emissf, e1flxf, e1cts1f, e1cts2f

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1):: avphilog,dt1,du,dup
      integer,               dimension(NPTS,NLP1):: ixo1, iyo, iyop

      integer :: k, m

!
!===> ...  begin here
!

!  ---  obtain the "exchange" emissivities as a function of temperature
!       (fxo) and water amount (fyo). the temperature indices have
!       been obtained in longwave_setup_mod.

      avphilog(:,1:NLP1) = LOG10(avephi(:,1:NLP1))

      call locate_in_table                                              &
!  ---  inputs:
     &     ( mass_1, avphilog, 1, NLP1,                                 &
!  ---  outputs:
     &       du, iyo                                                    &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab2, ixoe2, iyo, dte2, du, 1, NLP1,                       &
!  ---  outputs:
     &       emiss                                                      &
     &     )

!  ---  the special case emiss(:,NLAY) for layer NLAY is obtained by
!       averaging the values for NLAY and NLAY+1.

      emiss(:,NLAY) = 0.5 * (emiss(:,NLAY) + emiss(:,NLP1))

!  ---  perform calculations for the e1 function. the terms involving
!       top layer du are not known.  we use index two to represent
!       index one in previous calculations.

      iyop(:,1)      = 1
      iyop(:,2:NLP1) = iyo(:,1:NLAY)
      dup (:,1)      = f_zero
      dup (:,2:NLP1) = du (:,1:NLAY)

      do k = 1, NLP1
        ixo1(:,k) = ixoe1(:,1)
        dt1 (:,k) = dte1 (:,1)
      enddo

!  ---  e1flx(:,1) equals e1cts1(:,1).

      call looktab                                                      &
!  ---  inputs:
     &     ( tab1, ixoe1, iyop, dte1, dup, 1, NLP1,                     &
!  ---  outputs:
     &       e1cts1                                                     &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab1, ixoe1, iyo, dte1, du, 1, NLAY,                       &
!  ---  outputs:
     &       e1cts2                                                     &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab1, ixo1, iyop, dt1, dup, 1, NLP1,                       &
!  ---  outputs:
     &       e1flx                                                      &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab1w, ixoe1, iyop, dte1, dup, 1, NLP1,                    &
!  ---  outputs:
     &       e1ctw1                                                     &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab1w, ixoe1, iyo, dte1, du, 1, NLAY,                      &
!  ---  outputs:
     &       e1ctw2                                                     &
     &     )

!  ---  calculations with ch4 and n2o require NBTRGE separate emissivity
!       bands for h2o.

      if ( ifch4n2o > 0 ) then
        avphilog(:,1:NLP1) = LOG10(avephif(:,1:NLP1))

        call locate_in_table                                            &
!  ---  inputs:
     &     ( mass_1, avphilog, 1, NLP1,                                 &
!  ---  outputs:
     &       du, iyo                                                    &
     &     )

        iyop(:,1)      = 1
        iyop(:,2:NLP1) = iyo(:,1:NLAY)
        dup (:,1)      = f_zero
        dup (:,2:NLP1) = du (:,1:NLAY)

        call looktab                                                    &
!  ---  inputs:
     &     ( tab2a, ixoe2, iyo, dte2, du, 1, NLP1,                      &
!  ---  outputs:
     &       emissf                                                     &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab1a, ixoe1, iyop, dte1, dup, 1, NLP1,                    &
!  ---  outputs:
     &       e1cts1f                                                    &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab1a, ixoe1, iyo, dte1, du, 1, NLAY,                      &
!  ---  outputs:
     &       e1cts2f                                                    &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab1a, ixo1, iyop, dt1, dup, 1, NLP1,                      &
!  ---  outputs:
     &       e1flxf                                                     &
     &     )

!  ---  the special case emissf(:,NLAY) for layer NLAY is obtained by
!       averaging the values for NLAY and NLAY+1.

        emissf(:,NLAY) = 0.5e0 * (emissf(:,NLAY) +emissf(:,NLP1))
      endif

!
      return
!...................................
      end  subroutine e1e290
!-----------------------------------



!  =========================================
!  *****    longwave_tables section    *****
!  =========================================


!-----------------------------------
      subroutine lwtable
!...................................
!  ---  inputs: ( none )
!  ---  outputs:( none )

!---------------------------------------------------------------------
!
!     table computes table entries used in longwave radiation.
!
!     author: m. d. schwarzkopf
!     revised: 1/1/93
!     certified:  radiation version 1.0
!
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:  through module variables

!  ---  outputs: through module variables

!  ---  locals:
      real (kind=kind_phys), dimension(NBLW) :: alfanb, arotnb,         &
     &       centnb, delnb
      real (kind=kind_phys), dimension(30)   :: cnusb, dnusb

      real (kind=kind_phys), dimension(NTTABH2O,NBLW) :: dbdtnb, src1nb
      real (kind=kind_phys), dimension(NTTABH2O,NBLX) :: srcwd

      real (kind=kind_phys), dimension(NTTABH2O, NUTABH2O) :: sumdbe,   &
     &       sum, sum3, sumwde, suma, sumdbea, sum3a

      real (kind=kind_phys), dimension(NTTABH2O) :: ddsc, fortcu, r1,   &
     &       r1wd, r2, s2, sc, srcs, sum4, sum4wd, sum6, sum7, sum8,    &
     &       t3, tfour, x, x1, xtemv, r1a, r2a, s2a, t3a, sum4a,        &
     &       sum6a, sum7a, sum8a
      real (kind=kind_phys), dimension(NUTABH2O) :: expo, fac, x2,      &
     &       zmass, zroot

      integer  :: n, m, itab, jtab, nsubds, nsb,  iter
      real (kind=kind_phys) :: zmassincr, cent, del, bdlo, bdhi, anu,   &
     &       c1, freq_cutoff, bdlah2o, bdhah2o

!
!===>  ...  start here
!
!  ---  compute a*b and sqrt(a*b) for all 10 cm-1 frequency bands.

      do n = 1, NBLW
        centnb(n) = 0.5 * (bdbond(n) + bdbond(n+1))
        delnb (n) = bdbond(n+1) - bdbond(n)
        alfanb(n) = brndm(n)*arndm(n)
        arotnb(n) = sqrt( alfanb(n) )
      enddo

      if ( ifch4n2o > 0 ) then
!  ---  select, from random band coeff in 1200-1400 cm-1 range,
!       those appropriate for NBTRGE h2o bands.

        bdlah2o = CN_basic%bdl1h2o
        bdhah2o = CN_basic%bdh1h2o

!  ---  define critical frequency (cutoff for wide band ?? )

        freq_cutoff = 1400.
      else
        freq_cutoff = 1200.
      endif

!  ---  begin table computations here.  compute temperatures and 
!       masses for table entries.
!
!       note: the dimensioning and initialization of xtemv and other
!       arrays with dimension of NTTABH2O=28 imply a restriction of
!       model temperatures from 100k to 370k.
!
!       the dimensioning of zmass, zroot and other arrays with
!       dimension of NUTABH2O=180 imply a restriction of model h2o
!       amounts such that optical paths are between 10**-16 and 10**2,
!       in cgs units.

      zmass(1)  = 10.0 ** mass_1%min_val
      zmassincr = 10.0 ** mass_1%tab_inc

!  ---  the definition of zmassincr as 10**0.1 is slightly different
!       from all previous versions, in which it is 1.258925411.
!       This produces slightly different answers (fluxes differ by
!       1.0e-6 W/m2).

      do jtab = 2, NUTABH2O
        zmass(jtab) = zmass(jtab-1) * zmassincr
      enddo

      do jtab = 1, NUTABH2O
        zroot(jtab) = sqrt( zmass(jtab) )
      enddo

      do itab = 1, NTTABH2O
        xtemv (itab) = temp_1%min_val + temp_1%tab_inc*(itab-1)
        tfour (itab) = xtemv(itab) ** 4
        fortcu(itab) = 4.0 * xtemv(itab)**3
      enddo

!  ---  the computation of source, dsrce is needed only for the
!       combined wide band case.  to obtain them,  the source must
!       be computed for each of the NBLX wide bands srcwd then
!       combined using iband into source.

      do n = 1, NBLY
        do itab = 1, NTTABH2O
          tabsr%vae(itab,n) = f_zero
        enddo
      enddo

      do n = 1, NBLX
        do itab = 1, NTTABH2O
          srcwd(itab,n) = f_zero
        enddo
      enddo
!
!  ---  begin frequency loop.
!

      do n = 1, NBLX

!  ---  the 160-560 cm-1 region

        if ( n <= 40 ) then
          cent = centnb(n+16)
          del  = delnb (n+16)
          bdlo = bdbond(n+16)
          bdhi = bdbond(n+16+1)

!  ---  the 560-1200 cm-1 region, and the 2270-2380 cm-1 region

        else
          cent = 0.5 * (bdlocm(n-32+IOFFSET) + bdhicm(n-32+IOFFSET))
          del  = bdhicm(n-32+IOFFSET) - bdlocm(n-32+IOFFSET)
          bdlo = bdlocm(n-32+IOFFSET)
          bdhi = bdhicm(n-32+IOFFSET)
        endif

!  ---  for purposes of accuracy, all evaluations of planck functions
!       are made on 10 cm-1 intervals, then summed into the NBLX wide
!       bands.  the last subband may be narrower than 10 cm-1.

        nsubds = (del - 1.0e-3) / 10 + 1
        do nsb = 1, nsubds
          if ( nsb <> nsubds ) then
            cnusb(nsb) = 10.0 * (nsb - 1) + bdlo + 5.0
            dnusb(nsb) = 10.0
          else
            cnusb(nsb) = 0.5 * (10.0*(nsb - 1) + bdlo + bdhi)
            dnusb(nsb) = bdhi -   (10.0*(nsb - 1) + bdlo)
          endif
          c1 = 3.7412e-5 * cnusb(nsb)**3

!  ---  begin temperature loop.

          do itab = 1, NTTABH2O
            x    (itab)   = 1.4387 * cnusb(nsb) / xtemv(itab)
            x1   (itab)   = exp( x(itab) )
            srcs (itab)   = c1 / (x1(itab) - 1.0)
            srcwd(itab,n) = srcwd(itab,n) + srcs(itab)*dnusb(nsb)
          enddo

        enddo
      enddo

!  ---  the following loops create the combined wide band quantities
!       for robert scheme: the first 40 bands map to bands 1 to 8 in
!       for ckd2.1 scheme: the first 40 bands map to bands 1 to 40 in

      do n = 1, 40
        do itab = 1, NTTABH2O
          tabsr%vae(itab,iband(n)) = tabsr%vae(itab,iband(n)) +         &
     &                               srcwd(itab,n)
        enddo
      enddo

      do n = 9+IOFFSET, NBLY
        do itab = 1, NTTABH2O
          tabsr%vae(itab,n) = srcwd(itab,n+32-IOFFSET)
        enddo
      enddo

      do n = 1, NBLY
        do itab = 1, NTTABH2O-1
          tabsr%td(itab,n) = 0.1*(tabsr%vae(itab+1,n)-tabsr%vae(itab,n))
        enddo
      enddo

!  ---  first compute planck functions src1nb and derivatives dbdtnb
!       for use in table evaluations.  these are different from
!       source, dsrce because different frequency points are used in
!       evaluation, the frequency ranges are different, and the
!       derivative algorithm is different.

      do n = 1, NBLW
        cent = centnb(n)
        del  = delnb (n)

!  ---  note: at present, the iter loop is only used for iter=2.
!       the loop structure is kept so that in the future, we may
!       use a quadrature scheme for the planck function evaluation,
!       rather than use the mid-band frequency.

        do iter = 2, 2
          anu = cent + 0.5 * (iter - 2) * del
          c1  = 3.7412e-5 * anu**3

!  ---  temperature loop.

          do itab = 1, NTTABH2O
            x  (itab)      = 1.4387 * anu / xtemv(itab)
            x1 (itab)      = exp( x(itab) )
            sc (itab)      = c1 / ((x1(itab) - 1.0) + 1.0e-20)
            sc (itab)      = c1 / ( x1(itab) - 1.0)
            ddsc(itab)     = sc(itab) / (x1(itab)-1.0)                  &
     &                     * x1(itab) * x(itab) / xtemv(itab)
            src1nb(itab,n) = del * sc(itab)
            dbdtnb(itab,n) = del * ddsc(itab)
          enddo

        enddo
      enddo

!  ---  next compute r1, r2, s2, and t3 coefficients used for e3
!       function when the optical path is less than 10**-4.  in this
!       case, we assume a different dependence on zmass.  also obtain
!       r1wd, which is r1 summed over the 160-560 cm-1 range.

      do itab = 1, NTTABH2O
        sum4  (itab) = f_zero
        sum6  (itab) = f_zero
        sum7  (itab) = f_zero
        sum8  (itab) = f_zero
        sum4wd(itab) = f_zero
      enddo

      if ( ifch4n2o > 0 ) then
        sum4a (:)    = f_zero
        sum6a (:)    = f_zero
        sum7a (:)    = f_zero
        sum8a (:)    = f_zero
      endif

      do n = 1, NBLW
        cent = centnb(n)

!  ---  if include ch4-n2o: sum for freq ranges of 0-560, 1400-2200 cm-1
!       otherwise         : sum for freq ranges of 0-560, 1200-2200 cm-1

        if ( cent<5.6e+2 .or. cent>freq_cutoff .and. cent<=2.2e+3 ) then
          do itab = 1, NTTABH2O
            sum4(itab) = sum4(itab) + src1nb(itab,n)
            sum6(itab) = sum6(itab) + dbdtnb(itab,n)
            sum7(itab) = sum7(itab) + dbdtnb(itab,n)*arotnb(n)
            sum8(itab) = sum8(itab) + dbdtnb(itab,n)*alfanb(n)
          enddo
        endif

        if ( ifch4n2o > 0 ) then

!  ---  perform summations for frequency range of 1200-1400 cm-1
!       for sum4a, sum6a, sum7a, and sum8a. the computation depends
!       on the value of NBTRGE.

          if ( cent > 1.2e+3 .and. cent <= 1.4e+3 ) then
            if ( cent > bdlah2o .and. cent <= bdhah2o ) then
              sum4a(:) = sum4a(:) + src1nb(:,n)
              sum6a(:) = sum6a(:) + dbdtnb(:,n)
              sum7a(:) = sum7a(:) + dbdtnb(:,n)*arotnb(n)
              sum8a(:) = sum8a(:) + dbdtnb(:,n)*alfanb(n)
            endif
          endif
        endif

!  ---  perform summations over 160-560 cm-1 frequency range for e1
!       calculations sum4wd.

        if ( cent > 1.6e+2 .and. cent < 5.6e+2 ) then
          do itab = 1, NTTABH2O
            sum4wd(itab) = sum4wd(itab) + src1nb(itab,n)
          enddo
        endif

      enddo

      do itab = 1, NTTABH2O
        r1(itab)   = sum4(itab) / tfour (itab)
        r2(itab)   = sum6(itab) / fortcu(itab)
        s2(itab)   = sum7(itab) / fortcu(itab)
        t3(itab)   = sum8(itab) / fortcu(itab)
        r1wd(itab) = sum4wd(itab)/tfour(itab)
      enddo

      do jtab = 1, NUTABH2O
        do itab = 1, NTTABH2O
          sum   (itab,jtab) = f_zero
          sumdbe(itab,jtab) = f_zero
          sum3  (itab,jtab) = f_zero
          sumwde(itab,jtab) = f_zero
        enddo
      enddo

      if ( ifch4n2o > 0 ) then
        r1a(:)   = sum4a(:) / tfour(:)
        r2a(:)   = sum6a(:) / fortcu(:)
        s2a(:)   = sum7a(:) / fortcu(:)
        t3a(:)   = sum8a(:) / fortcu(:)

        suma   (:,:) = f_zero
        sumdbea(:,:) = f_zero
        sum3a  (:,:) = f_zero
      endif

!  ---  frequency loop begins.

      do n = 1, NBLW
        cent = centnb(n)

!  ---  perform calculations for frequency ranges of 0-560,
!       if include ch4-n2o: 1400-2200 cm-1.
!       else              : 1200-2200 cm-1.

        if ( cent<5.6e+2 .or. cent>freq_cutoff .and. cent<=2.2e+3 ) then
          do jtab = 1, NUTABH2O
            x2  (jtab) = arotnb(n) * zroot(jtab)
            expo(jtab) = exp( - x2(jtab) )
          enddo

          do jtab = 121, NUTABH2O
            fac(jtab) = (1.0 - (1.0 + x2(jtab)) * expo(jtab))           &
     &                / (alfanb(n) * zmass(jtab))
          enddo

          do jtab = 1, NUTABH2O
            do itab = 1, NTTABH2O
              sum   (itab,jtab) = sum   (itab,jtab)                     &
     &                          + src1nb(itab,n)*expo(jtab)
              sumdbe(itab,jtab) = sumdbe(itab,jtab)                     &
     &                          + dbdtnb(itab,n)*expo(jtab)
            enddo
          enddo

          do jtab = 121, NUTABH2O
            do itab = 1, NTTABH2O
              sum3  (itab,jtab) = sum3  (itab,jtab)                     &
     &                          + dbdtnb(itab,n)*fac(jtab)
            enddo
          enddo
        endif

!  ---  perform calculations over the frequency range 1200-1400 cm-1.
!       the calculations depend on the value of NBTRGE.

        if ( ifch4n2o > 0 ) then
          if ( cent > 1.2e+3 .and. cent <= 1.4e+3 ) then
            if ( cent > bdlah2o .and. cent <= bdhah2o ) then
              x2  (:) = arotnb(n) * zroot(:)
              expo(:) = exp( - x2(:) )

              do jtab = 121, NUTABH2O
                fac(jtab) = (1.0 - (1.0 + x2(jtab))*expo(jtab))         &
     &                    / (alfanb(n)*zmass(jtab))
              enddo

              do jtab = 1, NUTABH2O
                suma   (:,jtab) = suma   (:,jtab)                       &
     &                          + src1nb(:,n)*expo(jtab)
                sumdbea(:,jtab) = sumdbea(:,jtab)                       &
     &                          + dbdtnb(:,n)*expo(jtab)
              enddo

              do jtab = 121, NUTABH2O
                sum3a(:,jtab)   = sum3a(:,jtab)                         &
     &                          + dbdtnb(:,n)*fac(jtab)
              enddo
            endif
          endif
        endif

!  ---  compute sum over 160-560 cm-1 range for use in e1 calculations
!       sumwde.

        if ( cent > 1.6e+2 .and. cent < 5.6e+2 ) then
          do jtab = 1, NUTABH2O
            do itab = 1, NTTABH2O
              sumwde(itab,jtab) = sumwde(itab,jtab)                     &
     &                          + src1nb(itab,n)*expo(jtab)
            enddo
          enddo
        endif

      enddo

!  ---  frequency loop ends

      do jtab = 1, NUTABH2O
        do itab = 1, NTTABH2O
          tab1%vae(itab,jtab) = sum   (itab,jtab) / tfour (itab)
          tab2%vae(itab,jtab) = sumdbe(itab,jtab) / fortcu(itab)
        enddo
      enddo

      do jtab = 121, NUTABH2O
        do itab = 1, NTTABH2O
          tab3%vae(itab,jtab) = sum3  (itab,jtab) / fortcu(itab)
        enddo
      enddo

      do jtab = 1, 2
        do itab = 1, NTTABH2O
          tab1%vae(itab,jtab) = r1(itab)
        enddo
      enddo

      do jtab = 1, 120
        do itab = 1, NTTABH2O
          tab3%vae(itab,jtab) = r2(itab)/2.0 - s2(itab)*zroot(jtab)/3.0 &
     &                        + t3(itab)*zmass(jtab)/8.0
        enddo
      enddo

!  ---  compute e1 tables for 160-560 cm-1 bands.

      do jtab = 1, NUTABH2O
        do itab = 1, NTTABH2O
          tab1w%vae(itab,jtab) = sumwde(itab,jtab) / tfour(itab)
        enddo
      enddo

      do jtab = 1, 2
        do itab = 1, NTTABH2O
          tab1w%vae(itab,jtab) = r1wd(itab)
        enddo
      enddo

!  ---  initialize all derivative table entries.

      do jtab = 1, NUTABH2O
        do itab = 1, NTTABH2O
          tab1%td (itab,jtab) = f_zero
          tab1w%td(itab,jtab) = f_zero
          tab2%td (itab,jtab) = f_zero
          tab3%td (itab,jtab) = f_zero
          tab1%md (itab,jtab) = f_zero
          tab1w%md(itab,jtab) = f_zero
          tab2%md (itab,jtab) = f_zero
          tab3%md (itab,jtab) = f_zero
          tab1%cd (itab,jtab) = f_zero
          tab1w%cd(itab,jtab) = f_zero
          tab2%cd (itab,jtab) = f_zero
          tab3%cd (itab,jtab) = f_zero
        enddo
      enddo

!  ---  compute table entries for temperature derivatives.

      do jtab = 1, NUTABH2O
        do itab = 1, NTTABH2O-1
          tab1%td (itab,jtab) = (tab1%vae(itab+1,jtab)                  &
     &                        -  tab1%vae(itab,jtab)) / temp_1%tab_inc

          tab1w%td(itab,jtab) = (tab1w%vae(itab+1,jtab)                 &
     &                        -  tab1w%vae(itab,jtab)) / temp_1%tab_inc

          tab2%td (itab,jtab) = (tab2%vae(itab+1,jtab)                  &
     &                        -  tab2%vae(itab,jtab)) / temp_1%tab_inc

          tab3%td (itab,jtab) = (tab3%vae(itab+1,jtab)                  &
     &                        -  tab3%vae(itab,jtab)) / temp_1%tab_inc
        enddo
      enddo

!  ---  compute table entries for mass derivatives.

      do jtab = 1, NUTABH2O-1
        do itab = 1, NTTABH2O
          tab1%md (itab,jtab) = (tab1%vae(itab,jtab+1)                  &
     &                        -  tab1%vae(itab,jtab)) / mass_1%tab_inc

          tab1w%md(itab,jtab) = (tab1w%vae(itab,jtab+1)                 &
     &                        -  tab1w%vae(itab,jtab)) / mass_1%tab_inc

          tab2%md (itab,jtab) = (tab2%vae(itab,jtab+1)                  &
     &                        -  tab2%vae(itab,jtab)) / mass_1%tab_inc

          tab3%md (itab,jtab) = (tab3%vae(itab,jtab+1)                  &
     &                        -  tab3%vae(itab,jtab)) / mass_1%tab_inc
        enddo
      enddo

!  ---  compute table entries for cross derivatives.

      do jtab = 1, NUTABH2O-1
        do itab = 1, NTTABH2O-1
          tab1%cd (itab,jtab) =                                         &
     &        (tab1%vae(itab+1,jtab+1) - tab1%vae(itab+1,jtab)          &
     &       - tab1%vae(itab  ,jtab+1) + tab1%vae(itab  ,jtab))         &
     &       / (temp_1%tab_inc*mass_1%tab_inc)

          tab1w%cd(itab,jtab) =                                         &
     &        (tab1w%vae(itab+1,jtab+1) - tab1w%vae(itab+1,jtab)        &
     &       - tab1w%vae(itab  ,jtab+1) + tab1w%vae(itab  ,jtab))       &
     &       / (temp_1%tab_inc*mass_1%tab_inc)

          tab3%cd (itab,jtab) =                                         &
     &        (tab3%vae(itab+1,jtab+1) - tab3%vae(itab+1,jtab)          &
     &       - tab3%vae(itab  ,jtab+1) + tab3%vae(itab  ,jtab))         &
     &       / (temp_1%tab_inc*mass_1%tab_inc)
        enddo
      enddo

      if ( ifch4n2o > 0 ) then
        do jtab = 1, NUTABH2O
          tab1a%vae(:,jtab) = suma(:,jtab)/tfour(:)
          tab2a%vae(:,jtab) = sumdbea(:,jtab)/fortcu(:)
        enddo

        do jtab = 121, NUTABH2O
          tab3a%vae(:,jtab) = sum3a(:,jtab)/fortcu(:)
        enddo

        do jtab = 1, 2
          tab1a%vae(:,jtab) = r1a(:)
        enddo

        do jtab = 1, 120
          tab3a%vae(:,jtab) = r2a(:)/2.0 - s2a(:)*zroot(jtab)/3.0       &
     &                      + t3a(:)*zmass(jtab)/8.0
        enddo

        tab1a%td(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab2a%td(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab3a%td(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab1a%md(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab2a%md(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab3a%md(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab1a%cd(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab2a%cd(1:NTTABH2O,1:NUTABH2O) = f_zero
        tab3a%cd(1:NTTABH2O,1:NUTABH2O) = f_zero

        tab1a%td(1:NTTABH2O-1,1:NUTABH2O) =                             &
     &           (tab1a%vae(2:NTTABH2O,1:NUTABH2O)                      &
     &         -  tab1a%vae(1:NTTABH2O-1,1:NUTABH2O)) / temp_1%tab_inc
        tab2a%td(1:NTTABH2O-1,1:NUTABH2O) =                             &
     &           (tab2a%vae(2:NTTABH2O,1:NUTABH2O)                      &
     &         -  tab2a%vae(1:NTTABH2O-1,1:NUTABH2O)) / temp_1%tab_inc
        tab3a%td(1:NTTABH2O-1,1:NUTABH2O) =                             &
     &           (tab3a%vae(2:NTTABH2O,1:NUTABH2O)                      &
     &         -  tab3a%vae(1:NTTABH2O-1,1:NUTABH2O)) / temp_1%tab_inc
        tab1a%md(1:NTTABH2O,1:NUTABH2O-1) =                             &
     &           (tab1a%vae(1:NTTABH2O,2:NUTABH2O)                      &
     &         -  tab1a%vae(1:NTTABH2O,1:NUTABH2O-1)) / mass_1%tab_inc
        tab2a%md(1:NTTABH2O,1:NUTABH2O-1) =                             &
     &           (tab2a%vae(1:NTTABH2O,2:NUTABH2O)                      &
     &         -  tab2a%vae(1:NTTABH2O,1:NUTABH2O-1)) / mass_1%tab_inc
        tab3a%md(1:NTTABH2O,1:NUTABH2O-1) =                             &
     &           (tab3a%vae(1:NTTABH2O,2:NUTABH2O)                      &
     &         -  tab3a%vae(1:NTTABH2O,1:NUTABH2O-1)) / mass_1%tab_inc

        tab1a%cd(1:NTTABH2O-1,1:NUTABH2O-1) =                           &
     &           (tab1a%vae(2:NTTABH2O,2:NUTABH2O)                      &
     &         -  tab1a%vae(2:NTTABH2O,1:NUTABH2O-1)                    &
     &         -  tab1a%vae(1:NTTABH2O-1,2:NUTABH2O)                    &
     &         +  tab1a%vae(1:NTTABH2O-1,1:NUTABH2O-1))                 &
     &         / (temp_1%tab_inc*mass_1%tab_inc)
        tab3a%cd(1:NTTABH2O-1,1:NUTABH2O-1) =                           &
     &           (tab3a%vae(2:NTTABH2O,2:NUTABH2O)                      &
     &         -  tab3a%vae(2:NTTABH2O,1:NUTABH2O-1)                    &
     &         -  tab3a%vae(1:NTTABH2O-1,2:NUTABH2O)                    &
     &         +  tab3a%vae(1:NTTABH2O-1,1:NUTABH2O-1))                 &
     &         / (temp_1%tab_inc*mass_1%tab_inc)
      endif

!
      return
!...................................
      end subroutine lwtable
!-----------------------------------



!-----------------------------------
      subroutine e290                                                   &
!...................................
!  ---  inputs:
     &     ( k, tab2, tab2a, dte2, ixoe2, mass_1, avephi,               &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emiss, emissb                                              &
!  --- optional in/out:
     &,      avephif, emissf, emissbf                                   &
     &     )


!-----------------------------------------------------------------------
!
!     e290 computes the exchange terms in the flux equation for longwave
!     radiation for all terms except the exchange with the top of the
!     atmosphere.  the method is a table lookup on a pre-computed e2
!     function (defined in reference (2)).  calculation are done in the
!     frequency range: 0-560, 1200-2200 cm-1 for q(approximate).
!     motivation for these calculations is in references (1) and (2).
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
      implicit none

!  ---  inputs:
      integer, intent(in) ::  NPTS, NLAY, NLP1, k

      type (tab1_type), intent(in) :: tab2, tab2a
      type (axis_type), intent(in) :: mass_1

      real (kind=kind_phys), dimension(:,:), intent(in) :: dte2, avephi
      integer, dimension(:,:), intent(in) :: ixoe2

      real (kind=kind_phys), dimension(:,:), optional, intent(in) ::    &
     &       avephif

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) ::emiss,emissb

      real (kind=kind_phys), dimension(:,:), optional, intent(out) ::   &
     &       emissf, emissbf

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: avphilog, dtk, du
      integer,               dimension(NPTS,NLP1) :: ixok, iyo

      integer                              :: kp, m

!
!===> ...  begin here
!
!  ---  obtain the "exchange" emissivities as a function of temperature
!       (fxo) and water amount (avephi). the temperature indices have
!       been obtained in Lwrad1. calculations are for flux level k, with
!       kp = k+1 to NLAY+1. the case k=1 is excluded (done in E1e290).
!       calculations are also made for flux levels k to NLAY, for
!       contributions from flux level k. in this case, the temperature
!       index (ixok) represents tflux(:,:,k-1); the water index (iyo)
!       has the same values as in the case with varying kp.

      do kp = k, NLAY
        ixok(:,kp) = ixoe2(:,k-1)
        dtk (:,kp) = dte2 (:,k-1)
      enddo

      avphilog(:,k:NLP1) = LOG10(avephi(:,k:NLP1))

      call locate_in_table                                              &
!  ---  inputs:
     &     ( mass_1, avphilog, k, NLP1,                                 &
!  ---  outputs:
     &       du, iyo                                                    &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab2, ixoe2, iyo, dte2, du, k, NLP1,                       &
!  ---  outputs:
     &       emiss                                                      &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab2, ixok, iyo, dtk, du, k, NLAY,                         &
!  ---  outputs:
     &       emissb                                                     &
     &     )

!  ---  the special case emiss(:,NLAY) for layer NLAY is obtained by
!       averaging the values for NLAY and NLAY+1. note that emiss(:,NLAY+1)
!       is not useful after this point.

      emiss(:,NLAY) = 0.5 * (emiss(:,NLAY) +emiss(:,NLP1))

!  ---  calculations with ch4 and n2o require NBTRGE separate emissivity
!       bands for h2o. reqults are in emissf (flux level k) and
!       emissbf (other levels).

      if ( ifch4n2o > 0 ) then
        avphilog(:,k:NLP1) = LOG10(avephif(:,k:NLP1))

        call locate_in_table                                            &
!  ---  inputs:
     &     ( mass_1, avphilog, k, NLP1,                                 &
!  ---  outputs:
     &       du, iyo                                                    &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab2a, ixoe2, iyo, dte2, du, k, NLP1,                      &
!  ---  outputs:
     &       emissf                                                     &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab2a, ixok, iyo, dtk, du, k, NLAY,                        &
!  ---  outputs:
     &       emissbf                                                    &
     &     )

!  ---  the special case emissf(:,NLAY) for layer NLAY is obtained by
!       averaging the values for NLAY and NLAY+1. note that emissf(:,NLAY+1)
!       is not useful after this point.

        emissf(:,NLAY) = 0.5e0 * (emissf(:,NLAY) +emissf(:,NLP1))
      endif

!
      return
!...................................
      end subroutine e290
!-----------------------------------



!-----------------------------------
      subroutine enear                                                  &
!...................................
!  ---  inputs:
     &     ( tab2, tab3, tab2a, tab3a, tpl1, tpl2,                      &
     &       dte1, dte2, ixoe1, ixoe2, temp_1, mass_1,                  &
     &       empl1, empl2, emx1, emx2, var2,                            &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       emisdg, emspec                                             &
!  ---  optional in/out:
     &,      empl1f, empl2f, emx1f, emx2f, vrpfh2o,                     &
     &       emisdgf, emspecf                                           &
     &     )

!--------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      type (tab1_type), intent(in) :: tab2, tab2a, tab3, tab3a
      type (axis_type), intent(in) :: temp_1, mass_1

      real (kind=kind_phys), dimension(:,:), intent(in) :: tpl1, tpl2,  &
     &       dte1, dte2, empl1, empl2, var2
      real (kind=kind_phys), dimension(:),   intent(in) :: emx1, emx2

      integer, dimension(:,:), intent(in) :: ixoe1, ixoe2

      real (kind=kind_phys), dimension(:,:), optional, intent(in) ::    &
     &       empl1f, empl2f, vrpfh2o
      real (kind=kind_phys), dimension(:),   optional, intent(in) ::    &
     &       emx1f, emx2f

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) ::emisdg,emspec

      real (kind=kind_phys), dimension(:,:), optional, intent(out) ::   &
     &       emisdgf,emspecf

!  ---  locals:

      real (kind=kind_phys), dimension(NPTS,NLP1) :: dxsp, ylog, dysp,  &
     &       emiss, emd1, emd2, emd2f, emd1f, emissf

      integer,               dimension(NPTS,NLP1) :: ixsp, iysp
      integer :: m

!
!===> ...  begin here
!

      ixsp(:,NLAY) = ixoe2(:,NLAY-1)
      ixsp(:,NLP1) = ixoe1(:,NLAY-1)
      dxsp(:,NLAY) = dte2 (:,NLAY-1)
      dxsp(:,NLP1) = dte1 (:,NLAY-1)

      ylog(:,NLAY) = ALOG10(var2(:,NLAY))
      ylog(:,NLP1) = ALOG10(var2(:,NLAY) + empl1(:,NLAY))

      call locate_in_table                                              &
!  ---  inputs:
     &     ( mass_1, ylog, NLAY, NLP1,                                  &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

!  ---  compute exchange terms in the flux equation for two terms used
!       for nearby layer computations.

      call looktab                                                      &
!  ---  inputs:
     &     ( tab2, ixsp, iysp, dxsp, dysp, NLAY, NLP1,                  &
!  ---  outputs:
     &       emiss                                                      &
     &     )

!  ---  obtain index values of h2o pressure-scaled mass for each band
!       in the 1200-1400 range.

      if ( ifch4n2o > 0 ) then
        ylog(:,NLAY) = ALOG10(vrpfh2o(:,NLAY))
        ylog(:,NLP1) = ALOG10(vrpfh2o(:,NLAY) + empl1f(:,NLAY))

        call locate_in_table                                            &
!  ---  inputs:
     &     ( mass_1, ylog, NLAY, NLP1,                                  &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

!  ---  compute exchange terms in the flux equation for two terms used
!       for nearby layer computations.

        call looktab                                                    &
!  ---  inputs:
     &     ( tab2a, ixsp, iysp, dxsp, dysp, NLAY, NLP1,                 &
!  ---  outputs:
     &       emissf                                                     &
     &     )

      endif

!  ---  compute nearby layer transmissivities for h2o.

      call locate_in_table                                              &
!  ---  inputs:
     &     ( temp_1, tpl1, 1, NLP1,                                     &
!  ---  outputs:
     &       dxsp, ixsp                                                 &
     &     )

      ylog(:,:) = ALOG10(empl1(:,:))

      call locate_in_table                                              &
!  ---  inputs:
     &     ( mass_1, ylog, 1, NLP1,                                     &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab3, ixsp, iysp, dxsp, dysp, 1, NLP1,                     &
!  ---  outputs:
     &       emd1                                                       &
     &     )

!  ---  obtain index values of h2o pressure-scaled mass for each band
!       in the 1200-1400 range.

      if ( ifch4n2o > 0 ) then
        ylog(:,1:NLP1) = ALOG10(empl1f(:,1:NLP1))

        call locate_in_table                                            &
!  ---  inputs:
     &     ( mass_1, ylog, 1, NLP1,                                     &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab3a, ixsp, iysp, dxsp, dysp, 1, NLP1,                    &
!  ---  outputs:
     &       emd1f                                                      &
     &     )

      endif

      call locate_in_table                                              &
!  ---  inputs:
     &     ( temp_1, tpl2, 2, NLP1,                                     &
!  ---  outputs:
     &       dxsp, ixsp                                                 &
     &     )

      ylog(:,2:NLP1) = ALOG10(empl2(:,2:NLP1))

      call locate_in_table                                              &
!  ---  inputs:
     &     ( mass_1, ylog, 2, NLP1,                                     &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

      call looktab                                                      &
!  ---  inputs:
     &     ( tab3, ixsp, iysp, dxsp, dysp, 2, NLP1,                     &
!  ---  outputs:
     &       emd2                                                       &
     &     )

!  ---  obtain index values of h2o pressure-scaled mass for each band
!       in the 1200-1400 range.

      if ( ifch4n2o > 0 ) then
        ylog(:,2:NLP1) = ALOG10(empl2f(:,2:NLP1))

        call locate_in_table                                            &
!  ---  inputs:
     &     ( mass_1, ylog, 2, NLP1,                                     &
!  ---  outputs:
     &       dysp, iysp                                                 &
     &     )

        call looktab                                                    &
!  ---  inputs:
     &     ( tab3a, ixsp, iysp, dxsp, dysp, 2, NLP1,                    &
!  ---  outputs:
     &       emd2f                                                      &
     &     )

      endif

!  ---  compute nearby layer and special-case transmissivities for
!       emissivity using methods for h2o given in reference (4).

      emisdg(:,2:NLAY) = emd2(:,2:NLAY) + emd1(:,2:NLAY)
      emisdg(:,NLP1)  = 2.0e0 * emd1(:,NLP1)
      emspec(:,1) = (emd1(:,1)*empl1(:,1)-emd1(:,NLP1)*empl1(:,NLP1))   &
     &            / emx1(:) + 0.25e0 * (emiss(:,NLAY) + emiss(:,NLP1))

      emspec(:,2) = 2.0e0 * (emd1(:,1)*empl1(:,1)                       &
     &            - emd2(:,NLP1)*empl2(:,NLP1)) / emx2(:)

      if ( ifch4n2o > 0 ) then
        emisdgf(:,2:NLAY) = emd2f(:,2:NLAY) + emd1f(:,2:NLAY)
        emisdgf(:,NLP1)  = 2.0e0 * emd1f(:,NLP1)
        emspecf(:,1) = (emd1f(:,1)*empl1f(:,1)                          &
     &               - emd1f(:,NLP1)*empl1f(:,NLP1)) / emx1f(:)         &
     &               + 0.25e0*(emissf(:,NLAY) + emissf(:,NLP1))
        emspecf(:,2) = 2.0e0 * (emd1f(:,1)*empl1f(:,1)                  &
     &               - emd2f(:,NLP1)*empl2f(:,NLP1)) / emx2f(:)
      endif

!
      return
!...................................
      end subroutine enear
!-----------------------------------



!  =========================================
!  *****    lw_gases_stdtf section     *****
!  =========================================


!-----------------------------------
      subroutine gases_stdtf                                            &
!...................................
!  ---  inputs:
     &     ( rch4, rn2o, rco2, pd, plm, pd8, plm8, NLAY, NLP1 )
!  ---  outputs: ( none )

!---------------------------------------------------------------------
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, NLP1

      real (kind=kind_phys), intent(in) :: rch4, rn2o, rco2
      real (kind=kind_phys), dimension(:), intent(in) :: pd, plm,       &
     &       pd8, plm8

!  ---  local constants:
      integer, parameter :: nfreq_bands_sea_ch4 = 1
      integer, parameter :: nfreq_bands_sea_n2o = 3
      integer, parameter :: nfreq_bands_sea_co2 = 5

!  ---  locals:
      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS,3) ::    &
     &       trns_std_hi_nf, trns_std_lo_nf

      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS)   ::    &
     &       pressint_hiv_std_pt1, pressint_lov_std_pt1,                &
     &       pressint_hiv_std_pt2, pressint_lov_std_pt2
      integer,            dimension(NSTDCO2LVLS,NSTDCO2LVLS)   ::       &
     &         indx_pressint_hiv_std_pt1, indx_pressint_lov_std_pt1,    &
     &         indx_pressint_hiv_std_pt2, indx_pressint_lov_std_pt2

      real (kind=kind_phys), dimension(NLP1,NLP1,3) ::                  &
     &       trns_interp_lyr_ps_nf, trns_interp_lyr_ps8_nf,             &
     &       trns_interp_lvl_ps_nf, trns_interp_lvl_ps8_nf

      real (kind=kind_phys) :: temp1_lvl(NLP1), temp2_lvl(NLP1,NLP1)

      real (kind=kind_phys) :: rch4_vmr, rn2o_vmr, rco2_vmr

      integer :: k

!
!===> ...  begin here
!

!  ---  comp ch4 and n2o transmission function tables and write out

      if ( ifch4n2o > 0 .and. ich4n2otf == 1 ) then

        rch4_vmr = rch4 * 1.0e+9
        rn2o_vmr = rn2o * 1.0e+9

        call Ch4_lblinterp

        call N2o_lblinterp

        open (NIOFRAD, file='stdch4n2otfs', form='unformatted',         &
     &         access='sequential')
        write(NIOFRAD) ch451,ch458, ch4dt51,ch4dt58, ch4d2t51,ch4d2t58
        write(NIOFRAD) n2o51,n2o58, n2odt51,n2odt58, n2od2t51,n2od2t58
        write(NIOFRAD) n2o71,n2o78, n2odt71,n2odt78, n2od2t71,n2od2t78
        write(NIOFRAD) n2o91,n2o98, n2odt91,n2odt98, n2od2t91,n2od2t98
        close(NIOFRAD)

      endif   ! end if_ifch4n2o_block

!  ---  comp co2 transmission function table and write out

      if ( ico2tfs == 1 ) then

        rco2_vmr = rco2 * 1.0e+6

        call Co2_lblinterp

        open (NIOFRAD,  file='stdco2tfs', form='unformatted',           &
     &         access='sequential')
        write(NIOFRAD) co215nbps1, co215nbps8, co2dt15nbps1,            &
     &                 co2dt15nbps8, co2d2t15nbps1, co2d2t15nbps8
        write(NIOFRAD) co251, co258, cdt51, cdt58, c2d51, c2d58,        &
     &                 co2m51, co2m58, cdtm51, cdtm58, c2dm51, c2dm58
        close(NIOFRAD)

      endif   ! end if_iflgco2_block


! ================
      contains
! ================

!-----------------------------------
      subroutine ch4_lblinterp
!...................................

!  -------------------------------------------------------------------
!
!     this module is 1) a standalone program for a ch4 interpolation
!     to user-specified pressure levels and ch4 concentration;
!     2) an interface between a GCM and ch4 interpolation
!
!     input files:
!
!         1)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (ch4_vmr).
!         2)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (ch4_vmr). may not be
!                  used, depending on value of (ch4_vmr).
!
!     output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these
!                  files are created if ifdef (writeinterpch4) is
!                  on. otherwise, it is assumed the data is fed
!                 directly back to the parent model.
!  -------------------------------------------------------------------
!

      implicit none

!  ---  local constants:
      real (kind=kind_phys), dimension(5) ::  ch4_std_vmr
      data  ch4_std_vmr / 700.0, 1250.0, 1750.0, 2250.0, 2800.0 /

      integer, dimension(nfreq_bands_sea_ch4) ::  ntbnd_ch4
      data  ntbnd_ch4(:)  / 3 /

      logical, dimension(nfreq_bands_sea_ch4) :: do_lyrcalc_nf,         &
     &      do_lvlcalc_nf, do_ctscalc_nf
      data  do_lyrcalc_nf(:) / .true.  /
      data  do_lvlcalc_nf(:) / .false. /
      data  do_ctscalc_nf(:) / .false. /

      character(len=3), parameter :: gas_type = 'ch4'

!  ---  local variables:
      logical ::  callrctrns, do_lyrcalc, do_lvlcalc, do_ctscalc

      real (kind=kind_phys), dimension(NLP1,NLP1) :: ch4p10_lyr,        &
     &       ch4p8_lyr, d2ch4t10_lyr, d2ch4t8_lyr, dch4dt10_lyr,        &
     &       dch4dt8_lyr

      real (kind=kind_phys) :: std_lo, std_hi
      integer :: k, kp, n, nf, nt, nstd_lo, nstd_hi

!
!===>  ...  begin here
!

!  ---  call gasint for ch4 interpolations.
!
!     using the value of the ch4 volume mixing ratio (ch4_vmr) and
!     the available standard ch4 mixing ratios (with lbl
!     transmission functions) obtain the two standard mixing ratios
!     which bracket (ch4_vmr). if (as in a fixed ch4 experiment) the
!     difference between (ch4_vmr) and the higher of the standard
!     mixing ratios is less than a tolerance (taken as 0.1 ppbv)
!     we will assume that that standard ch4 transmissivity applies
!     to (ch4_vmr) without interpolation. otherwise, interpolation
!     to (ch4_vmr) will be performed, in rctrns.F

      callrctrns = .false.

      if ( rch4_vmr<ch4_std_vmr(1) .or. rch4_vmr>ch4_std_vmr(5) ) then

        print *, ' Error! ch4 volume mixing ratio is out of range'
        print *, 'rch4_vmr =',rch4_vmr,' ppbv. the acceptable range',   &
     &           ' is ',ch4_std_vmr(1),' to ', ch4_std_vmr(5),' ppbv.'
        stop

      elseif ( rch4_vmr == ch4_std_vmr(1) ) then

        std_lo = ch4_std_vmr(1)
        std_hi = ch4_std_vmr(1)
        nstd_lo = 1
        nstd_hi = 1

      else 

        lab_do_n : do n = 1, 4
          if ( rch4_vmr > ch4_std_vmr(n) .and.                          &
     &         rch4_vmr <= ch4_std_vmr(n+1) ) then

            std_lo = ch4_std_vmr(n)
            std_hi = ch4_std_vmr(n+1)
            nstd_lo = n
            nstd_hi = n+1

            if ( abs(rch4_vmr-std_hi) > 1.0e-1 ) callrctrns=.true.
            exit lab_do_n
          endif
        enddo  lab_do_n

      endif   ! end if_rch4_vmr_block

!  --- ...  loop on frequency bands. in the 1996 SEA formulation, 
!           there are 1 frequency ranges for lbl ch4 transmissions:
!       nf = 1:  lbl transmissions over 1200-1400 cm-1

      do nf = 1, nfreq_bands_sea_ch4

        do_lyrcalc = do_lyrcalc_nf(nf)
        do_lvlcalc = do_lvlcalc_nf(nf)
        do_ctscalc = do_ctscalc_nf(nf)

!  ---  read in ch4 transmissivities, data is for all temperature profiles
!       required for the freq band (at 1 or 2 appropriate concentrations).
!       the number of temperature profiles for band (nf) is ntbnd(nf).
!       in 1996 SEA formulation, ntbnd_ch4=3 (USSTD,1976; USSTD,1976 +- 25).

        call read_lbltfs                                                &
!  ---  inputs:
     &     ( gas_type, callrctrns, nstd_lo, nstd_hi, nf, ntbnd_ch4,     &
!  ---  outputs:
     &       trns_std_hi_nf,trns_std_lo_nf                              &
     &     )

!  ---  loop on temperature structures needed for each frequency band.

        do nt = 1, ntbnd_ch4(nf)

!         print *,' nt,trns_std_hi_nf =',nt,trns_std_hi_nf(1,1,nt)

          call gasint                                                   &
!  ---  inputs:
     &     ( pd, plm, pd8, plm8, gas_type, rch4_vmr, std_lo, std_hi,    &
     &       callrctrns, do_lvlcalc, do_ctscalc, do_lyrcalc,            &
     &       nf, nt, NLAY, NLP1                                         &
     &     )

        enddo   ! end do_nt_loop

!  ---  perform final processing for each frequency band.

        call gasins                                                     &
!  ---  inputs:
     &     ( gas_type, do_lvlcalc, do_ctscalc, do_lyrcalc,              &
     &       nf, ntbnd_ch4(nf), NLP1, NLP1,                             &
!  ---  outputs:
     &       temp2_lvl, temp1_lvl, dch4dt10_lyr,                        &
     &       temp2_lvl, temp1_lvl, ch4p10_lyr,                          &
     &       temp2_lvl, temp1_lvl, d2ch4t10_lyr,                        &
     &       temp2_lvl, temp1_lvl, dch4dt8_lyr,                         &
     &       temp2_lvl, temp1_lvl, ch4p8_lyr,                           &
     &       temp2_lvl, temp1_lvl, d2ch4t8_lyr                          &
     &     )

!  ---  define arrays for the SEA module- the SEA model nomenclature
!       has been used here

        do k = 1, NLP1
          do kp = 1, NLP1
            ch4dt51 (kp,k) = dch4dt10_lyr(kp,k)
            ch451   (kp,k) = ch4p10_lyr  (kp,k)
            ch4d2t51(kp,k) = d2ch4t10_lyr(kp,k)
            ch4dt58 (kp,k) = dch4dt8_lyr (kp,k)
            ch458   (kp,k) = ch4p8_lyr   (kp,k)
            ch4d2t58(kp,k) = d2ch4t8_lyr (kp,k)
          enddo
        enddo

      enddo   ! end do_nf_loop
!
!...................................
      end subroutine ch4_lblinterp
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine n2o_lblinterp
!...................................

!---------------------------------------------------------------------
!
!     this module is 1) a standalone program for a n2o interpolation
!     to user-specified pressure levels and n2o concentration;
!     2) an interface between a GCM and n2o interpolation
!
!     input files:
!
!         1)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (n2o_vmr).
!         2)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (n2o_vmr). may not be
!                  used, depending on value of (n2o_vmr).
!
!     output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these
!                  files are created if ifdef (writeinterpn2o) is
!                  on. otherwise, it is assumed the data is fed
!                 directly back to the parent model.
!  -------------------------------------------------------------------
!

      implicit none

!  ---  local constants:
      real (kind=kind_phys), dimension(4) :: n2o_std_vmr
      data  n2o_std_vmr / 275.0, 310.0, 340.0, 375.0 /

      integer, dimension(nfreq_bands_sea_n2o) ::  ntbnd_n2o
      data  ntbnd_n2o(:)  / 3, 3, 3 /

      logical, dimension(nfreq_bands_sea_n2o) :: do_lyrcalc_nf,         &
     &      do_lvlcalc_nf, do_ctscalc_nf
      data  do_lyrcalc_nf(:) / .true.,  .true.,  .true.  /
      data  do_lvlcalc_nf(:) / .false., .false., .false. /
      data  do_ctscalc_nf(:) / .false., .false., .false. /

      character(len=3), parameter :: gas_type = 'n2o'

!  ---  local variables:
      logical ::  callrctrns, do_lyrcalc, do_lvlcalc, do_ctscalc

      real (kind=kind_phys), dimension(NLP1,NLP1) :: n2op10_lyr,        &
     &       n2op8_lyr, d2n2ot10_lyr, d2n2ot8_lyr, dn2odt10_lyr,        &
     &       dn2odt8_lyr

      real (kind=kind_phys) :: std_lo, std_hi
      integer :: k, kp, n, nf, nt, nstd_lo, nstd_hi

!
!===>  ...  begin here
!

!  ---  call gasint for n2o interpolations.
!
!     using the value of the n2o volume mixing ratio (n2o_vmr) and
!     the available standard n2o mixing ratios (with lbl
!     transmission functions) obtain the two standard mixing ratios
!     which bracket (n2o_vmr). if (as in a fixed n2o experiment) the
!     difference between (n2o_vmr) and the higher of the standard
!     mixing ratios is less than a tolerance (taken as 0.1 ppbv)
!     we will assume that that standard n2o transmissivity applies
!     to (n2o_vmr) without interpolation. otherwise, interpolation
!     to (n2o_vmr) will be performed, in rctrns.F

      callrctrns = .false.

      if ( rn2o_vmr<n2o_std_vmr(1) .or. rn2o_vmr>n2o_std_vmr(4) ) then

        print *, ' Error! n2o volume mixing ratio is out of range'
        print *, 'rn2o_vmr =',rn2o_vmr,' ppbv. the acceptable range',   &
     &           ' is ',n2o_std_vmr(1),' to ', n2o_std_vmr(4),' ppbv.'
        stop

      elseif ( rn2o_vmr == n2o_std_vmr(1) ) then

        std_lo = n2o_std_vmr(1)
        std_hi = n2o_std_vmr(1)
        nstd_lo = 1
        nstd_hi = 1

      else

        lab_do_n : do n = 1, 3
          if ( rn2o_vmr > n2o_std_vmr(n) .and.                          &
     &         rn2o_vmr <= n2o_std_vmr(n+1) ) then

            std_lo = n2o_std_vmr(n)
            std_hi = n2o_std_vmr(n+1)
            nstd_lo = n
            nstd_hi = n+1

            if ( abs(rn2o_vmr-std_hi) > 1.0e-1 ) callrctrns=.true.
            exit lab_do_n
          endif
        enddo  lab_do_n

      endif   ! end if_rn2o_vmr_block

!  --- ...  loop on frequency bands. in the 1996 SEA formulation, there
!           are 3 frequency ranges for lbl n2o transmissions:
!       nf = 1:  lbl transmissions over 1200-1400 cm-1
!       nf = 2:  lbl transmissions over 1070-1200 cm-1
!       nf = 3:  lbl transmissions over 560-630 cm-1

      do nf = 1, nfreq_bands_sea_n2o

        do_lyrcalc = do_lyrcalc_nf(nf)
        do_lvlcalc = do_lvlcalc_nf(nf)
        do_ctscalc = do_ctscalc_nf(nf)

!  ---  read in n2o transmissivities, data is read for all temperature
!       profiles required for the freq band (at 1 or 2 appropriate
!       concentrations). the number of temp profiles is ntbnd(nf).
!       in 1996 SEA formulation, ntbnd_n2o=3(USSTD,1976; USSTD,1976 +- 25).

        call read_lbltfs                                                &
!  ---  inputs:
     &     ( gas_type, callrctrns, nstd_lo, nstd_hi, nf, ntbnd_n2o,     &
!  ---  outputs:
     &       trns_std_hi_nf, trns_std_lo_nf                             &
     &     )

!  ---  loop on temperature structures needed for each frequency band.

        do nt = 1, ntbnd_n2o(nf)

          call gasint                                                   &
!  ---  inputs:
     &     ( pd, plm, pd8, plm8, gas_type, rn2o_vmr, std_lo, std_hi,    &
     &       callrctrns, do_lvlcalc, do_ctscalc, do_lyrcalc,            &
     &       nf, nt, NLAY, NLP1                                         &
     &     )

        enddo   ! end do_nt_loop

!  ---  perform final processing for each frequency band.

        call gasins                                                     &
!  ---  inputs:
     &     ( gas_type, do_lvlcalc, do_ctscalc, do_lyrcalc,              &
     &       nf, ntbnd_n2o(nf), NLP1, NLP1,                             &
!  ---  outputs:
     &       temp2_lvl, temp1_lvl, dn2odt10_lyr,                        &
     &       temp2_lvl, temp1_lvl, n2op10_lyr,                          &
     &       temp2_lvl, temp1_lvl, d2n2ot10_lyr,                        &
     &       temp2_lvl, temp1_lvl, dn2odt8_lyr,                         &
     &       temp2_lvl, temp1_lvl, n2op8_lyr,                           &
     &       temp2_lvl, temp1_lvl, d2n2ot8_lyr                          &
     &     )

!  ---  define arrays for the SEA module- the SEA model nomenclature
!       has been used here

        if ( nf == 1 ) then

          do k = 1, NLP1
            do kp = 1, NLP1
              n2odt51 (kp,k) = dn2odt10_lyr(kp,k)
              n2o51   (kp,k) = n2op10_lyr  (kp,k)
              n2od2t51(kp,k) = d2n2ot10_lyr(kp,k)
              n2odt58 (kp,k) = dn2odt8_lyr (kp,k)
              n2o58   (kp,k) = n2op8_lyr   (kp,k)
              n2od2t58(kp,k) = d2n2ot8_lyr (kp,k)
            enddo
          enddo

        elseif ( nf == 2 ) then

          do k = 1, NLP1
            do kp = 1, NLP1
              n2odt91 (kp,k) = dn2odt10_lyr(kp,k)
              n2o91   (kp,k) = n2op10_lyr  (kp,k)
              n2od2t91(kp,k) = d2n2ot10_lyr(kp,k)
              n2odt98 (kp,k) = dn2odt8_lyr (kp,k)
              n2o98   (kp,k) = n2op8_lyr   (kp,k)
              n2od2t98(kp,k) = d2n2ot8_lyr (kp,k)
            enddo
          enddo

        elseif ( nf == 3 ) then

          do k = 1, NLP1
            do kp = 1, NLP1
              n2odt71 (kp,k) = dn2odt10_lyr(kp,k)
              n2o71   (kp,k) = n2op10_lyr  (kp,k)
              n2od2t71(kp,k) = d2n2ot10_lyr(kp,k)
              n2odt78 (kp,k) = dn2odt8_lyr (kp,k)
              n2o78   (kp,k) = n2op8_lyr   (kp,k)
              n2od2t78(kp,k) = d2n2ot8_lyr (kp,k)
            enddo
          enddo

        endif   ! end if_nf_block

      enddo   ! end do_nf_loop
!
!...................................
      end subroutine n2o_lblinterp
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine co2_lblinterp
!...................................

!---------------------------------------------------------------------
!
!     this module is 1) a standalone program for a co2 interpolation
!     to user-specified pressure levels and co2 concentration;
!     2) an interface between a GCM and co2 interpolation
!
!     input files:
!
!         1)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (co2_vmr).
!         2)     : gas transmission function at higher of 2
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (co2_vmr). may not be
!                  used, depending on value of (co2_vmr).
!
!     output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these
!                  files are created if ifdef (writeinterpco2) is
!                  on. otherwise, it is assumed the data is fed
!                 directly back to the parent model.
!---------------------------------------------------------------------
!

      implicit none

!  ---  local constants:
      real (kind=kind_phys), dimension(9) :: co2_std_vmr
      data  co2_std_vmr / 165.0, 300.0, 330.0, 348.0, 356.0,            &
     &                    360.0, 600.0, 660.0, 1320.0 /

      integer, dimension(nfreq_bands_sea_co2) ::  ntbnd_co2
      data  ntbnd_co2(:)  / 3, 3, 3, 3, 1 /

      logical, dimension(nfreq_bands_sea_co2) :: do_lyrcalc_nf,         &
     &      do_lvlcalc_nf, do_ctscalc_nf
      data  do_lyrcalc_nf(:) / .true., .false.,.false.,.false.,.true.  /
      data  do_lvlcalc_nf(:) / .true., .true., .true., .true., .true.  /
      data  do_ctscalc_nf(:) / .false.,.true., .true., .true., .false. /

      character(len=3), parameter :: gas_type = 'co2'

!  ---  local variables:
      logical ::  callrctrns, do_lyrcalc, do_lvlcalc, do_ctscalc

      real (kind=kind_phys), dimension(NLP1,NLP1) :: co2p10_lyr,        &
     &       co2p8_lyr, d2ct10_lyr, d2ct8_lyr, dcdt10_lyr, dcdt8_lyr,   &
     &       co2p10_lvl, co2p8_lvl, d2ct10_lvl, d2ct8_lvl, dcdt10_lvl,  &
     &       dcdt8_lvl

      real (kind=kind_phys), dimension(NLP1)      :: co2p10_cts,        &
     &       co2p8_cts, d2ct10_cts, d2ct8_cts, dcdt10_cts, dcdt8_cts,   &
     &       temp1_lvl, temp2_lvl

      real (kind=kind_phys) :: std_lo, std_hi
      integer :: k, kp, n, nf, nt, nstd_lo, nstd_hi

!
!===>  ...  begin here
!

!  ---  call gasint for co2 interpolations.
!
!     using the value of the co2 volume mixing ratio (co2_vmr) and
!     the available standard co2 mixing ratios (with lbl
!     transmission functions) obtain the two standard mixing ratios
!     which bracket (co2_vmr). if (as in a fixed co2 experiment) the
!     difference between (co2_vmr) and the higher of the standard
!     mixing ratios is less than a tolerance (taken as .0001 ppmv)
!     we will assume that that standard co2 transmissivity applies
!     to (co2_vmr) without interpolation. otherwise, interpolation
!     to (co2_vmr) will be performed, in rctrns.F

      callrctrns = .false.

      if ( rco2_vmr<co2_std_vmr(1) .or. rco2_vmr>co2_std_vmr(9) ) then

        print *, ' Error! co2 volume mixing ratio is out of range'
        print *, 'rco2_vmr = ',rco2_vmr,' ppmv. the acceptable range',  &
     &           ' is ',co2_std_vmr(1),' to ', co2_std_vmr(9),' ppmv.'
        stop

      elseif ( rco2_vmr == co2_std_vmr(1) ) then

        std_lo = co2_std_vmr(1)
        std_hi = co2_std_vmr(1)
        nstd_lo = 1
        nstd_hi = 1

      else

        lab_do_n : do n = 1, 8
          if ( rco2_vmr > co2_std_vmr(n) .and.                          &
     &         rco2_vmr <= co2_std_vmr(n+1) ) then

            std_lo = co2_std_vmr(n)
            std_hi = co2_std_vmr(n+1)
            nstd_lo = n
            nstd_hi = n+1

            if ( abs(rco2_vmr-std_hi) > 1.0e-4 ) callrctrns=.true.
            exit lab_do_n
          endif
        enddo  lab_do_n

      endif

!  --- ... loop on frequency bands. in the 1996 SEA formulation, there
!          are 5 frequency ranges for lbl co2 transmissions:
!       nf = 1:  lbl transmissions over 490-850 cm-1
!       nf = 2:  lbl transmissions over 490-630 cm-1
!       nf = 3:  lbl transmissions over 630-700 cm-1
!       nf = 4:  lbl transmissions over 700-800 cm-1
!       nf = 5:  lbl transmissions over 2270-2380 cm-1

      do nf = 1, nfreq_bands_sea_co2

        do_lyrcalc = do_lyrcalc_nf(nf)
        do_lvlcalc = do_lvlcalc_nf(nf)
        do_ctscalc = do_ctscalc_nf(nf)

!  ---  read in co2 transmissivities, data is read for all temperature
!       profiles required for the freq band (at 1 or 2 appropriate
!       concentrations). the number of temp profiles is ntbnd(nf).
!       in 1996 SEA formulation, ntband_co2=3(USSTD,1976;USSTD,1976 +- 25)
!       except for the 4.3 um band (nf=2) where the number is one.

        call read_lbltfs                                                &
!  ---  inputs:
     &     ( gas_type,callrctrns,nstd_lo,nstd_hi,nf,ntbnd_co2,          &
!  ---  outputs:
     &       trns_std_hi_nf, trns_std_lo_nf                             &
     &     )

!  ---  loop on temperature structures needed for each frequency band.

        do nt = 1, ntbnd_co2(nf)

          call gasint                                                   &
!  ---  inputs:
     &     ( pd, plm, pd8, plm8, gas_type, rco2_vmr, std_lo, std_hi,    &
     &       callrctrns, do_lvlcalc,do_ctscalc,do_lyrcalc,              &
     &       nf, nt, NLAY, NLP1                                         &
     &     )

        enddo   ! end_do_nt_loop

!  ---  perform final processing for each frequency band.

        call gasins                                                     &
!  ---  inputs:
     &     ( gas_type, do_lvlcalc, do_ctscalc, do_lyrcalc,              &
     &       nf, ntbnd_co2(nf), NLP1, NLP1,                             &
!  ---  outputs:
     &       dcdt10_lvl, dcdt10_cts, dcdt10_lyr,                        &
     &       co2p10_lvl, co2p10_cts, co2p10_lyr,                        &
     &       d2ct10_lvl, d2ct10_cts, d2ct10_lyr,                        &
     &       dcdt8_lvl,  dcdt8_cts,  dcdt8_lyr,                         &
     &       co2p8_lvl,  co2p8_cts,  co2p8_lyr,                         &
     &       d2ct8_lvl,  d2ct8_cts,  d2ct8_lyr                          &
     &     )

!  ---  define arrays for the SEA module- the SEA model nomenclature
!       has been used here

        if ( nf == 1 ) then

          do k = 1, NLP1
            do kp = 1, NLP1
              cdt51(kp,k) = dcdt10_lyr(kp,k)
              co251(kp,k) = co2p10_lyr(kp,k)
              c2d51(kp,k) = d2ct10_lyr(kp,k)
              cdt58(kp,k) = dcdt8_lyr (kp,k)
              co258(kp,k) = co2p8_lyr (kp,k)
              c2d58(kp,k) = d2ct8_lyr (kp,k)
            enddo
          enddo

          do k=1,NLAY
            cdtm51(k) = dcdt10_lvl(k,k+1)
            co2m51(k) = co2p10_lvl(k,k+1)
            c2dm51(k) = d2ct10_lvl(k,k+1)
            cdtm58(k) = dcdt8_lvl (k,k+1)
            co2m58(k) = co2p8_lvl (k,k+1)
            c2dm58(k) = d2ct8_lvl (k,k+1)
          enddo

        elseif ( nf >= 2 .and. nf <= 4 ) then

          do kp = 1, NLP1
            co2dt15nbps1 (kp,nf-1) = dcdt10_cts(kp)
            co215nbps1   (kp,nf-1) = co2p10_cts(kp)
            co2d2t15nbps1(kp,nf-1) = d2ct10_cts(kp)
            co2dt15nbps8 (kp,nf-1) = dcdt8_cts (kp)
            co215nbps8   (kp,nf-1) = co2p8_cts (kp)
            co2d2t15nbps8(kp,nf-1) = d2ct8_cts (kp)
          enddo

        endif   ! end if_nf_block

      enddo   ! end do_nf_loop
!
!...................................
      end subroutine co2_lblinterp
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine gasint                                                 &
!...................................
!  ---  inputs:
     &     ( pd, plm, pd8, plm8, gas_type, gas_vmr, std_lo, std_hi,     &
     &       callrctrns, do_lvlcalc, do_ctscalc, do_lyrcalc,            &
     &       nf, nt, NLAY, NLP1                                         &
!  ---  outputs: (to parent variables)
     &     )

!---------------------------------------------------------------------
!
!        inputs:
!
!           trns_std_lo : array of gas transmission functions at the
!                         lower of two standard gas concentrations
!                         at a given temperature profile.
!                         used if interpolation to the actual gas
!                         mixing ratio is required (callrctrns = true).
!                         dimensions: (NSTDCO2LVLS,NSTDCO2LVLS)
!           trns_std_hi : array of gas transmission functions at the
!                         higher of two standard gas concentrations
!                         at a given temperature profile.
!                         dimensions: (NSTDCO2LVLS,NSTDCO2LVLS)
!              gas_vmr  : actual gas concentration (in ppmv)
!               std_lo  : gas concentration (ppmv) of lower of
!                         two standard concentrations.
!               std_hi  : gas concentration (ppmv) of higher of
!                         two standard concentrations.
!           callrctrns  : call rctrns.f if true.
! pressint_hiv_std_pt1  : allocated array used for rctrns hi pressure
! pressint_lov_std_pt1  : allocated array used for rctrns low pressure
! pressint_hiv_std_pt2  : allocated array used for rctrns hi pressure
! pressint_lov_std_pt2  : allocated array used for rctrns low pressure
! indx_pressint_hiv_std_pt1  : allocated index array used in rctrns
! indx_pressint_lov_std_pt1  : allocated index array used in rctrns
! indx_pressint_hiv_std_pt2  : allocated index array used in rctrns
! indx_pressint_lov_std_pt2  : allocated index array used in rctrns
!           do_lvlcalc  : compute level gas transmissivities if true.
!           do_lyrcalc  : compute layer gas transmissivities if true.
!        do_ctscalc  : compute cts level gas transmissivities if true
!                   nf  : frequency band index
!                   nt  : temperature index (for the freq band)
!        ndimkp, ndimk  : extents of dimensions for output interp
!                         transmissivity arrays.
!              pd, plm  : see description below. note that the
!                         present limits on pressures (from the lbl
!                         computations require that the top level
!                         be 0 mb, and the bottom level pressure
!                         not exceed 1165 mb.
!             pd8, plm8 : same as pd, plm; highest pressure is 0.8*
!                         (highest pressure in plm).
!
!     outputs:
!        trns_interp_lyr_ps_nf : array of interpolated layer transmission
!                             functions for the pressure array (pd).
!        trns_interp_lyr_ps8_nf: array of interpolated layer transmission
!                             functions for the pressure array (pd8).
!        trns_interp_lvl_ps_nf : array of interpolated level transmission
!                             functions for the pressure array (plm).
!        trns_interp_lvl_ps8_nf: array of interpolated level transmission
!                             functions for the pressure array (plm8).
!
!       gasint interpolates carbon dioxide transmission functions
!  from the standard level grid,for which the transmission functions
!  have been pre-calculated, to the grid structure specified by the
!  user.
!
!        method:
!
!      gasint is employable for two purposes: 1) to obtain transmis-
!  sivities between any 2 of an array of user-defined pressures; and
!  2) to obtain layer-mean transmissivities between any 2 of an array
!  of user-defined pressure layers.to clarify these two purposes,see
!  the diagram and discussion below.
!
!     let p be an array of user-defined pressures
!     and plm the array of user-defined level pressures
!     and pd the extent of the interpolation layer.
!       for many purposes,plm will be chosen to be the average
!     pressure in the interpolation layer -i.e.,
!     plm(i) = 0.5*(pd(i-1)+pd(i)).
!
!       - - - - - - - - -   pd(i-1)  -----------!
!                                               !
!       -----------------   plm(i), p(k)-----!  !  interpolation layer
!                                            !  !
!       - - - - - - - - -   pd(i)    -----------!       model layer (i)
!                                            !
!       -----------------   plm(i+1)---------!
!            ...
!            ...                          (the notation used is
!            ...                          consistent with the code)
!            ...
!       - - - - - - - - -   pd(j-1)  -----------!
!                                               !
!       -----------------   plm(j), p(k')----!  !  interpolation layer
!                                            !  !
!       - - - - - - - - -   pd(j)    -----------!       model layer (j)
!                                            !
!       -----------------   plm(j+1)---------!
!
!      purpose 1:   the transmissivity between specific pressures
!      p(k) and p(k') ,tau(p(k),p(k'))  is computed by this program.
!      in this mode,there is no reference to layer pressures pd
!
!      purpose 2:   the transmissivity between a pressure p(k) and
!      the interpolation layer centered at p(k') (taulm(p(k),p(k'))
!      is obtained. it is computed by the integral
!
!                           pd(j)
!                           ----
!             1             !
!        -------------  *   !   tau ( p',plm(i) )  dp'
!        pd(j)-pd(j-1)      !
!                        ----
!                        pd(j-1)
!
!       the level pressures (plm) and layer-mean pressures (pd) are
!       both inputted in for this purpose.
!
!       in general, this integral is done using simpson's rule with
!       7 points. however , when plm(i) = pjm(j) (the nearby layer
!       case) a 51-point quadrature is used for greater accuracy.
!       note that taulm(k,k') is not a symmetrical matrix. also, no
!       quadrature is performed over the layer between the smallest
!       nonzero pressure and zero pressure;
!       taulm is taulm(0,plm(j)) in this case,and taulm(0,0)=1.

!
!            the following paragraphs depict the utilization of this
!       code when used to compute transmissivities between specific
!       pressures. later paragraphs describe additional features needed
!       for layer-mean transmissivities.
!
!          for a given co2 mixing ratio and standard temperature
!      profile,a table of transmission functions for a fixed grid
!     of atmospheric pressures has been pre-calculated.
!      the standard temperature profile is from the us
!     standard atmosphere (1977) table.additionally, the
!     same transmission functions have been pre-calculated for a
!     temperature profile increased and decreased (at all levels)
!     by 25 degrees.
!         this program reads in the prespecified transmission functions
!     and a user-supplied pressure grids (p(k)) and calculates trans-
!     mission functions ,tau(p(k),p(k')), for all (k,k') using one
!     of the above methods, for one of the three temperature profiles.
!     outputs are tables of transmission functions.
!
!     this code may be considered to be version 2 of the
!     interpolater. differences between this code and version 1,
!     written in ~1983, are as follows:
!
!     1) the code is written using arrays (previous code was entirely
!     scalar)
!     2) double precision quantities are removed. it is assumed that
!     this code is being run on a 64-bit machine. if not, the
!     user will have to reinstate double precisioning, or
!     the appropriate KIND statement in Fortran 90.
!     3) many redundant calculations were eliminated
!     4) the error function is defined exactly as in Ref. (2).
!
!     as a result, this version of the code runs 100-200 times as fast
!     as version 1, and is suitable for on-line calculation of the
!     co2 transmission function.
!
!            differences in the answers:
!
!     1) as part of the elimination of redundant calculation, the
!    following algorithmic change was performed:
!    calculation of the error function (error_guess1) at standard
!    pressures is done WITHOUT reference to model (user) pressures.
!    answers (in fractional absorptivity change) differ by 10-3 (max)
!    to 10-5. the new method should be "better", as there is no reason
!    why the error function at standard pressures should care
!    about the pressures where it will be interpolated to.
!
!     2) in the "closely spaced" case (model pressures p,p' both
!    between standard pressures (pa(k),pa(k+1))) the coefficients
!    (c,x,eta,sexp) are interpolated, not the approx function. this
!    is consistent with all other cases. fractional absorptivity changes
!    are ~3x10-5 (max) and only for a(p,p') with p~ or = p'.
!
!    3) double precision to single precision change causes fractional
!    absorptivity changes of < 1x10-6.
!
!             references:
!
!     (1): s.b.fels and m.d.schwarzkopf,"an efficient,accurate
!     algorithm for calculating co2 15 um band cooling rates",journal
!     of geophysical research,vol.86,no. c2, pp.1205-1232,1981.
!     (2): m.d. schwarzkopf and s.b. fels, "Improvements to the
!     algorithm for computing co2 transmissivities and cooling rates",
!     JGR, vol. 90, no. C10, pp10541-10550, 1985.
!
!            author:    m.daniel schwarzkopf
!
!            date:      14 july 1996
!
!            address:
!
!                      GFDL
!                      p.o.box 308
!                      princeton,n.j.08542
!                      u.s.a.
!            telephone:  (609) 452-6521
!
!            e-mail:   ds@gfdl.gov
!
!
!    NOTE: the comment line below is part of the original version
!    of this code, written in the late '70s by Stephen B. Fels.
!    although the code has been extensively rewritten, and might
!    be unrecognizable to Steve, this line is kept as a tribute
!    to him.
!
!      ************   function interpolater routine  *****

      implicit none

!  ---  inputs:
      real (kind=kind_phys), dimension(:), intent(in) :: pd,plm,pd8,plm8

      character(len=*), intent(in) :: gas_type

      real (kind=kind_phys), intent(in) :: gas_vmr, std_lo, std_hi

      logical, intent(in) :: do_lvlcalc,do_lyrcalc,do_ctscalc,callrctrns

      integer, intent(in) :: nf, nt, NLAY, NLP1

!  ---  outputs: (to parent variables)

!  ---  locals:
      logical :: do_triangle
      real (kind=kind_phys) :: wgt_lyr(7), wgt_nearby_lyr(51)

      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS) ::      &
     &       approx_guess1, caintv, sexpintv, xaintv, uexpintv,         &
     &       press_hiv, press_lov, error_guess1, trns_vmr

      real (kind=kind_phys), dimension(NLP1,NLP1) :: pressint_hiv,      &
     &       pressint_lov, caintv2, sexpintv2, xaintv2, uexpintv2,      &
     &       errorint_guess1, approxint_guess1
      real (kind=kind_phys), dimension(51,NLP1)   :: sexpnblv,          &
     &       uexpnblv, canblv, xanblv, pressnbl_lov, pressnbl_hiv,      &
     &       approxnbl_guess1, errornbl_guess1

      integer, dimension(NLP1,NLP1)   :: indx_pressint_hiv,             &
     &                                   indx_pressint_lov
      integer, dimension(51,  NLP1)   :: indx_pressnbl_lov,             &
     &                                   indx_pressnbl_hiv

!      xa, ca, dop_core, uexp, sexp are coefficients for
!      the approximation function (Eq. (4) in Ref. (2)) used in
!      the gas interpolation algorithm. the nomenclature is:
!
!      this code           Ref. (2)
!      ---------           --------
!       xa                  X (see Eq. A1)
!       ca                  C (see Eq. A1)
!       uexp                delta (see Eq. A6b)
!       sexp                gamma (see Eq. A6c)
!       dop_core            core (see Eq. A6a)
!
      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: xa,ca, uexp,sexp
      real (kind=kind_phys) :: dop_core

      integer :: n, k, kp, nklo, nkhi, nkplo, nkphi, nq, nprofile

!
!===> ...  begin here
!

!  ---  compute the layer weights for transmissivities. used only if
!       layer transmissivities are needed (do_lyrcalc = true)

      if ( do_lyrcalc ) then
        wgt_lyr(1) = 1./18.
        wgt_lyr(7) = 1./18.

        do n = 1, 3
          wgt_lyr(2*n) = 4./18.
        enddo

        do n = 1, 2
          wgt_lyr(2*n+1) = 2./18.
        enddo

        wgt_nearby_lyr(1) = 1./150.
        wgt_nearby_lyr(51) = 1./150.

        do n = 1, 25
          wgt_nearby_lyr(2*n) = 4./150.
        enddo

        do n = 1, 24
          wgt_nearby_lyr(2*n+1) = 2./150.
        enddo
      endif

!  ---  define transmission function array for (gas_vmr) over standard
!       pressures (pa), using a call to rctrns if necessary.

      if ( callrctrns ) then

        call rctrns                                                     &
!  ---  inputs:
     &     ( gas_type, std_lo, std_hi, gas_vmr, nf, nt,                 &
!  ---  outputs:
     &       trns_vmr                                                   &
     &     )

      else
        trns_vmr(:,:) = trns_std_hi_nf(:,:,nt)
      endif

      do k = 1, NSTDCO2LVLS
        trns_vmr(k,k) = 1.0
      enddo

!  ---  compute co2 transmission functions for actual co2 concentration
!       using method of section 5, Ref. (2).

      call coeint                                                       &
!  ---  inputs:
     &     ( gas_type, nf, trns_vmr,                                    &
!  ---  outputs:
     &       ca, sexp, xa, uexp, dop_core                               &
     &     )

!  --- ...  compute the interpolation.

!  ---  1) compute approx function at standard (pa) pressures

      do_triangle = .true.

      do k = 1, NSTDCO2LVLS
        do kp = k, NSTDCO2LVLS
          press_hiv(kp,k) = pa  (kp)
          press_lov(kp,k) = pa  (k)
          caintv   (kp,k) = ca  (kp)
          sexpintv (kp,k) = sexp(kp)
          xaintv   (kp,k) = xa  (kp)
          uexpintv (kp,k) = uexp(kp)
        enddo
      enddo

!  ---  the call (and calculations) to pathv2_std has been subsumed
!       into the subroutine approx_fn_std

      call approx_fn_std                                                &
!  --- inputs:
     &    ( press_hiv, press_lov, do_triangle, caintv, sexpintv,        &
     &      xaintv, uexpintv, dop_core,                                 &
!  --- outputs:
     &      approx_guess1                                               &
     &    )

!  ---  2) compute error function at standard (pa) pressures

      do k = 1, NSTDCO2LVLS
        do kp = k+1, NSTDCO2LVLS
          error_guess1(kp,k) = 1.0-trns_vmr(kp,k)-approx_guess1(kp,k)
        enddo

        error_guess1(k,k) = f_zero
      enddo

!  ---  define the actual extents of the level interpolation calculation.
!       this depends on the type of calculation desired (whether
!       do_lvlcalc, do_ctscalc is true).

      if ( do_lvlcalc ) then
        nklo = 1
        nkhi = NLP1
        nkplo = 1
        nkphi = NLP1
      elseif ( do_ctscalc ) then
        nklo = 1
        nkhi = 1
        nkplo = 1
        nkphi = NLP1
      endif

      if ( do_ctscalc .or. do_lvlcalc ) then

        do k = 1, NLP1
          trns_interp_lvl_ps_nf (k,k,nt) = 1.0
          trns_interp_lvl_ps8_nf(k,k,nt) = 1.0
        enddo
        do_triangle = .true.

        do nprofile = 1, 2

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

          if (nprofile == 1) then
            do k = nklo, nkhi
              do kp = k+nkplo, nkphi
                pressint_hiv(kp,k) = plm(kp)
                pressint_lov(kp,k) = plm(k)
              enddo
            enddo
          else
            do k = nklo, nkhi
              do kp = k+nkplo, nkphi
                pressint_hiv(kp,k) = plm8(kp)
                pressint_lov(kp,k) = plm8(k)
              enddo
            enddo
          endif

          call intcoef_2d                                               &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle, ca, xa,           &
     &       uexp, sexp, nklo, nkhi, nkplo, nkphi,                      &
!  ---  outputs:
     &       indx_pressint_hiv, indx_pressint_lov, caintv2,             &
     &       sexpintv2, xaintv2, uexpintv2                              &
     &     )

!  ---  4) interpolate error function to (pressint_hiv, pressint_lov)
!          for relevant (k',k)

          call interp_error                                             &
!  ---  inputs:
     &     ( error_guess1, pressint_hiv, pressint_lov,                  &
     &       indx_pressint_hiv, indx_pressint_lov, do_triangle,         &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       errorint_guess1                                            &
     &     )

!  ---  5) compute approx function for (pressint_hiv, pressint_lov)
!          the call (and calculations) to pathv2 has been subsumed
!          into the subroutine approx_fn

          call approx_fn                                                &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle,                   &
     &       caintv2, sexpintv2, xaintv2, uexpintv2, dop_core,          &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       approxint_guess1                                           &
     &     )

!  ---  6) compute interp transmission function using Eq.(3),
!          Ref.(2).

          if ( nprofile == 1 ) then
            do k = nklo, nkhi
              do kp = k+nkplo, nkphi
                trns_interp_lvl_ps_nf(kp,k,nt) = 1.0 -                  &
     &                 (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                trns_interp_lvl_ps_nf(k,kp,nt) =                        &
     &                 trns_interp_lvl_ps_nf(kp,k,nt)
              enddo
            enddo
          else
            do k = nklo, nkhi
              do kp = k+nkplo, nkphi
                trns_interp_lvl_ps8_nf(kp,k,nt) = 1.0 -                 &
     &                 (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                trns_interp_lvl_ps8_nf(k,kp,nt) =                       &
     &                 trns_interp_lvl_ps8_nf(kp,k,nt)
              enddo
            enddo
          endif

        enddo   ! end do_nprofile_loop
      endif   ! end if_do_ctscalc_block

      if ( do_lyrcalc ) then

!  ---  A): calculate, for (kp,k) pairs with kp > k, a set of 7 transmis-
!       sivities, with the values of p'(kp) encompassing the layer bounded
!       by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!       averaged transmissivity (trns_interp_lyr_ps(8)(kp,k,nt)).
!       B): calculate, for (kp,k) pairs with kp < k, a set of 7 transmis-
!       sivities, with the values of p'(kp) encompassing the layer bounded
!       by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!       averaged transmissivity (trns_interp_lyr_ps(8)(kp,k,nt)).
!       C): calculate, for pairs (kp,kp) with kp > 1, a set of 51 transmis-
!       sivities, with the values of p'(kp) encompassing the layer bounded
!       by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!       averaged transmissivity (trns_interp_lyr_ps(8)(kp,k,nt)).
!
!       note: one of the 7 (or 51) transmissivities equals the level
!       transmissivity (trns_interp_lvl_ps(8))
!
!       initialize the layer transmissivities to zero (summing later)
!       except the (1,1), which are set to 1

        trns_interp_lyr_ps_nf (:,:,nt) = f_zero
        trns_interp_lyr_ps8_nf(:,:,nt) = f_zero
        trns_interp_lyr_ps_nf (1,1,nt) = 1.0
        trns_interp_lyr_ps8_nf(1,1,nt) = 1.0

!  ---  case A): (kp) levels are at higher pressure, hence are used for
!       pressint_hiv. the (fixed) (k) levels are used for pressint_lov

        do_triangle = .true.

        nklo = 1
        nkhi = NLAY
        nkplo = 1
        nkphi = NLP1

        do nprofile = 1, 2

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

          do nq = 1, 7

            if ( nprofile == 1 ) then
              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  pressint_hiv(kp,k) = pd(kp-1)                         &
     &                               + (nq-1)* (pd(kp) - pd(kp-1)) / 6
                  pressint_lov(kp,k) = plm(k)
                enddo
              enddo
            else
              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  pressint_hiv(kp,k) = pd8(kp-1)                        &
     &                               + (nq-1)* (pd8(kp) - pd8(kp-1))/6
                  pressint_lov(kp,k) = plm8(k)
                enddo
              enddo
            endif

            call intcoef_2d                                             &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle, ca, xa,           &
     &       uexp, sexp, nklo, nkhi, nkplo, nkphi,                      &
!  ---  outputs:
     &       indx_pressint_hiv, indx_pressint_lov, caintv2,             &
     &       sexpintv2, xaintv2, uexpintv2                              &
     &     )

!  ---  4) interpolate error function to (pressint_hiv, pressint_lov)
!          for relevant (k',k)

            call interp_error                                           &
!  ---  inputs:
     &     ( error_guess1, pressint_hiv, pressint_lov,                  &
     &       indx_pressint_hiv, indx_pressint_lov, do_triangle,         &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       errorint_guess1                                            &
     &     )

!  ---  5) compute approx function for (pressint_hiv, pressint_lov)
!          the call (and calculations) to pathv2 has been subsumed
!          into the subroutine approx_fn

            call approx_fn                                              &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle,                   &
     &       caintv2, sexpintv2, xaintv2, uexpintv2, dop_core,          &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       approxint_guess1                                           &
     &     )

!  ---  6) compute interp transmission function using Eq.(3), Ref.(2).

            if ( nprofile == 1 ) then

              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  trns_interp_lyr_ps_nf(kp,k,nt) =                      &
     &            trns_interp_lyr_ps_nf(kp,k,nt) + wgt_lyr(nq) * (1.0   &
     &            - (errorint_guess1(kp,k) + approxint_guess1(kp,k)))

!  ---  for the case (nq=4), where  (pressint_hiv(kp,k) = plm(kp)) use
!       the (kp,1) unweighted values (errorint + approxint) for the
!       (1,kp) transmissivity, otherwise uncalculated. (exception:
!       when kp = nkphi, the (nq=7) value must be used)

                  if ( nq == 4 .and. k == nklo ) then
                    trns_interp_lyr_ps_nf(nklo,kp,nt) = 1.0 -           &
     &                (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                  endif
                enddo
              enddo

              if ( nq == 7 ) then
                trns_interp_lyr_ps_nf(nklo,nkphi,nt) = 1.0              &
     &                       - (errorint_guess1(nkphi,nklo)             &
     &                       + approxint_guess1(nkphi,nklo))
              endif

            else

              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  trns_interp_lyr_ps8_nf(kp,k,nt) =                     &
     &            trns_interp_lyr_ps8_nf(kp,k,nt) + wgt_lyr(nq) * (1.0  &
     &            - (errorint_guess1(kp,k) + approxint_guess1(kp,k)))

!  ---  for the case (nq=4), where  (pressint_hiv(kp,k) = plm(kp)) use
!       the (kp,1) unweighted values (errorint + approxint) for the
!       (1,kp) transmissivity, otherwise uncalculated. (exception:
!       when kp = nkphi, the (nq=7) value must be used)

                  if ( nq == 4 .and. k == nklo ) then
                    trns_interp_lyr_ps8_nf(nklo,kp,nt) = 1.0            &
     &                           - (errorint_guess1(kp,k)               &
     &                           + approxint_guess1(kp,k))
                  endif
                enddo
              enddo

              if ( nq == 7 ) then
                trns_interp_lyr_ps8_nf(nklo,nkphi,nt) = 1.0             &
     &                       - (errorint_guess1(nkphi,nklo)             &
     &                       + approxint_guess1(nkphi,nklo))
              endif

            endif   ! end if_nprofile_block
          enddo   ! end do_nq_loop
        enddo   ! end do_nprofile_loop

!  --- case B): (k) levels are at higher pressure, hence are used for
!      pressint_hiv. the (variable) (kp) levels are used for pressint_lov.
!      (kp,k) calculations are loaded into (k,kp) array locations to
!      keep calculations into the "upper sandwich". results are then put
!      into their proper array locations (before weighting function is
!      applied). also, no calculations are made for (1,k). these values
!      are obtained from level calcs for (k,1), nq=4.

        do_triangle = .true.

        nklo = 2
        nkhi = NLAY
        nkplo = 1
        nkphi = NLP1

        do nprofile = 1, 2

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

          do nq = 1, 7

            if ( nprofile == 1 ) then
              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  pressint_hiv(kp,k) = plm(kp)
                  pressint_lov(kp,k) = pd(k-1) + (nq-1)*                &
     &                                 (pd(k) - pd(k-1))/6
                enddo
              enddo

            else

              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  pressint_hiv(kp,k) = plm8(kp)
                  pressint_lov(kp,k) = pd8(k-1) + (nq-1)*               &
     &                                 (pd8(k) - pd8(k-1))/6
                enddo
              enddo
            endif

            call intcoef_2d                                             &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle, ca, xa,           &
     &       uexp, sexp, nklo, nkhi, nkplo, nkphi,                      &
!  ---  outputs:
     &       indx_pressint_hiv, indx_pressint_lov, caintv2,             &
     &       sexpintv2, xaintv2,uexpintv2                               &
     &     )

!  ---  4) interpolate error function to (pressint_hiv, pressint_lov)
!          for relevant (k',k)

            call interp_error                                           &
!  ---  inputs:
     &     ( error_guess1, pressint_hiv, pressint_lov,                  &
     &       indx_pressint_hiv, indx_pressint_lov, do_triangle,         &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       errorint_guess1                                            &
     &     )

!  ---  5) compute approx function for (pressint_hiv, pressint_lov)
!          the call (and calculations) to pathv2 has been subsumed
!          into the subroutine approx_fn

            call approx_fn                                              &
!  ---  inputs:
     &     ( pressint_hiv, pressint_lov, do_triangle,                   &
     &       caintv2, sexpintv2, xaintv2, uexpintv2, dop_core,          &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       approxint_guess1                                           &
     &     )

!  ---  6) compute interp transmission function using Eq.(3), Ref.(2).

            if ( nprofile == 1 ) then
              do k = nklo, nkhi
                do kp = k+nkplo, nkphi
                  trns_interp_lyr_ps_nf(k,kp,nt) =                      &
     &            trns_interp_lyr_ps_nf(k,kp,nt) + wgt_lyr(nq)*(1.0     &
     &            - (errorint_guess1(kp,k) + approxint_guess1(kp,k)))
                enddo
              enddo
            else
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  trns_interp_lyr_ps8_nf(k,kp,nt) =                     &
     &               trns_interp_lyr_ps8_nf(k,kp,nt) + wgt_lyr(nq)*(1.0 &
     &             - (errorint_guess1(kp,k) + approxint_guess1(kp,k)))
                enddo
              enddo
            endif

          enddo   ! end do_nq_loop
        enddo   ! end do_nprofile_loop

!  ---  C): calculate, for pairs (kp,kp) with kp > 1, a set of 51 transmis-
!       sivities, with the values of p'(kp) encompassing the layer bounded
!       by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!       averaged transmissivity (trns_interp_lyr_ps(8)(kp,k,nt)).
!       case C): (kp) levels are at higher pressure, hence are used for
!       pressint_hiv. the (fixed) (k) levels are used for pressint_lov

        do_triangle = .false.

        nklo = 2
        nkhi = NLP1
        nkplo = 1
        nkphi = 51

        do nprofile = 1, 2

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

          if ( nprofile == 1 ) then

            do k = nklo, nkhi-1
              do kp = 1, 25
                pressnbl_lov(kp,k) = pd(k-1) + (kp-1)*                  &
     &                               (pd(k) - pd(k-1))/50
                pressnbl_hiv(kp,k) = plm(k)
              enddo

              pressnbl_lov(26,k) = plm(k)
              pressnbl_hiv(26,k) = plm(k)

              do kp = 27, 51
                pressnbl_hiv(kp,k) = pd(k-1) + (kp-1)*                  &
     &                               (pd(k) - pd(k-1))/50
                pressnbl_lov(kp,k) = plm(k)
              enddo
            enddo

            do kp = 1, 50
              pressnbl_lov(kp,nkhi) = pd(nkhi-1) + (kp-1)*              &
     &                                (pd(nkhi) - pd(nkhi-1))/50
              pressnbl_hiv(kp,nkhi) = plm(nkhi)
            enddo

            pressnbl_lov(51,nkhi) = plm(nkhi)
            pressnbl_hiv(51,nkhi) = plm(nkhi)

          else

            do k = nklo, nkhi-1
              do kp = 1, 25
                pressnbl_lov(kp,k) = pd8(k-1) + (kp-1)*                 &
     &                               (pd8(k) - pd8(k-1))/50
                pressnbl_hiv(kp,k) = plm8(k)
              enddo

              pressnbl_lov(26,k) = plm8(k)
              pressnbl_hiv(26,k) = plm8(k)

              do kp = 27, 51
                pressnbl_hiv(kp,k) = pd8(k-1) + (kp-1)*                 &
     &                               (pd8(k) - pd8(k-1))/50
                pressnbl_lov(kp,k) = plm8(k)
              enddo
            enddo

            do kp = 1, 50
              pressnbl_lov(kp,nkhi) = pd8(nkhi-1) + (kp-1)*             &
     &                                (pd8(nkhi) - pd8(nkhi-1))/50
              pressnbl_hiv(kp,nkhi) = plm8(nkhi)
            enddo

            pressnbl_lov(51,nkhi) = plm8(nkhi)
            pressnbl_hiv(51,nkhi) = plm8(nkhi)

          endif   ! end if_nprofile_block

          call intcoef_2d                                               &
!  ---  inputs:
     &     ( pressnbl_hiv, pressnbl_lov, do_triangle, ca, xa,           &
     &       uexp, sexp, nklo, nkhi, nkplo, nkphi,                      &
!  ---  outputs:
     &       indx_pressnbl_hiv, indx_pressnbl_lov, canblv,              &
     &       sexpnblv, xanblv, uexpnblv                                 &
     &     )

!  ---  4) interpolate error function to (pressnbl_hiv, pressnbl_lov)
!          for relevant (k',k)

          call interp_error                                             &
!  ---  inputs:
     &     ( error_guess1, pressnbl_hiv, pressnbl_lov,                  &
     &       indx_pressnbl_hiv, indx_pressnbl_lov, do_triangle,         &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       errornbl_guess1                                            &
     &     )

!  ---  5) compute approx function for (pressnbl_hiv, pressnbl_lov)
!          the call (and calculations) to pathv2 has been subsumed
!          into the subroutine approx_fn

          call approx_fn                                                &
!  ---  inputs:
     &     ( pressnbl_hiv, pressnbl_lov, do_triangle,                   &
     &       canblv, sexpnblv, xanblv, uexpnblv, dop_core,              &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       approxnbl_guess1                                           &
     &     )

!  ---  6) compute interp transmission function using Eq.(3), Ref.(2).

          if ( nprofile == 1 ) then
            do k = nklo, nkhi
              do kp = 1, 51
                trns_interp_lyr_ps_nf(k,k,nt) =                         &
     &          trns_interp_lyr_ps_nf(k,k,nt) + wgt_nearby_lyr(kp)*     &
     &          (1.0 - (errornbl_guess1(kp,k) + approxnbl_guess1(kp,k)))
              enddo
            enddo
          else
            do k = nklo, nkhi
              do kp = 1, 51
                trns_interp_lyr_ps8_nf(k,k,nt) =                        &
     &          trns_interp_lyr_ps8_nf(k,k,nt) + wgt_nearby_lyr(kp)*    &
     &          (1.0 - (errornbl_guess1(kp,k) + approxnbl_guess1(kp,k)))
              enddo
            enddo
          endif

        enddo   ! end do_nprofile_loop
      endif   ! end if_do_lyrcalc_block
!
!...................................
      end subroutine gasint
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine gasins                                                 &
!...................................
!  ---  inputs:
     &     ( gas_type, do_lvlcalc, do_lvlctscalc, do_lyrcalc,           &
     &       nf, ntbnd, ndimkp, ndimk,                                  &
!  ---  outputs:
     &       dcdt10_lvl, dcdt10_lvlcts, dcdt10_lyr,                     &
     &       co2p10_lvl, co2p10_lvlcts, co2p10_lyr,                     &
     &       d2ct10_lvl, d2ct10_lvlcts, d2ct10_lyr,                     &
     &       dcdt8_lvl,  dcdt8_lvlcts,  dcdt8_lyr ,                     &
     &       co2p8_lvl,  co2p8_lvlcts,  co2p8_lyr ,                     &
     &       d2ct8_lvl,  d2ct8_lvlcts,  d2ct8_lyr                       &
     &     )

!---------------------------------------------------------------------
!        gasins processes transmission functions to produce
!     "consolidated" functions over the specific frequency band
!     ranges needed by the SEA code, and the derivatives needed
!     by the SEA algorithm. writing to a file, formerly done in
!     this module, is now done (if needed) in write_seaco2fcns.F
!---------------------------------------------------------------------

      implicit none
!
!  ---  inputs:
      character(len=*), intent(in) :: gas_type

      logical, intent(in) :: do_lvlcalc, do_lyrcalc, do_lvlctscalc

      integer, intent(in) :: nf, ntbnd, ndimkp, ndimk

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) ::             &
     &       dcdt10_lvl, dcdt10_lyr, co2p10_lvl, co2p10_lyr,            &
     &       d2ct10_lvl, d2ct10_lyr, dcdt8_lvl,  dcdt8_lyr,             &
     &       co2p8_lvl,  co2p8_lyr,  d2ct8_lvl,  d2ct8_lyr

      real (kind=kind_phys), dimension(:),   intent(out) ::             &
     &       dcdt10_lvlcts, co2p10_lvlcts, d2ct10_lvlcts,               &
     &       dcdt8_lvlcts,  co2p8_lvlcts,  d2ct8_lvlcts

!  ---  locals:
      real (kind=kind_phys) :: c1,c2
      integer :: k1, k2, k, kp

!
!===> ...  begin here
!
!  ---  obtain array extents for internal arrays

      k1 = size(trns_interp_lvl_ps_nf,1) ! this corresponds to ndimkp
      k2 = size(trns_interp_lvl_ps_nf,2) ! this corresponds to ndimk

      if ( gas_type == 'co2' ) then

!  ---  the following code is rewritten so that the radiative bands are:
!        nf=1    560-800     (consol.=490-850)
!        nf=2    560-630      consol=490-630
!        nf=3    630-700      consol=630-700
!        nf=4    700-800      consol=700-850
!        nf=5   2270-2380     consol=2270-2380
!       the following loop obtains transmission functions for bands
!       used in radiative model calculations,with the equivalent
!       widths kept from the original consolidated co2 tf's.

        if ( nf == 1 ) then
          c1=1.5
          c2=0.5
        endif

        if ( nf == 2 ) then
          c1=2.0
          c2=1.0
        endif

        if ( nf == 3 ) then
          c1=1.0
          c2=f_zero
        endif

        if ( nf == 4 ) then
          c1=1.5
          c2=0.5
        endif

        if ( nf == 5 ) then
          c1=1.0
          c2=f_zero
        endif
      endif

      if ( gas_type == 'ch4' ) then

!  ---  the following code is rewritten so that the radiative bands are:
!        nf=1    1200-1400    consol=1200-1400
!       the following loop obtains transmission functions for bands
!       used in radiative model calculations,with the equivalent
!       widths kept from the original consolidated co2 tf's.

        if ( nf == 1 ) then
          c1=1.0
          c2=f_zero
        endif
      endif

      if ( gas_type == 'n2o' ) then

!  ---  the following code is rewritten so that the radiative bands are:
!        nf=1    1200-1400    consol=1200-1400
!        nf=2    1070-1200    consol=1070-1200
!        nf=3    560-630    consol=560-630
!       the following loop obtains transmission functions for bands
!       used in radiative model calculations,with the equivalent
!       widths kept from the original consolidated co2 tf's.

        if ( nf == 1 ) then
          c1=1.0
          c2=f_zero
        endif

        if ( nf == 2 ) then
          c1=1.0
          c2=f_zero
        endif

        if ( nf == 3 ) then
          c1=1.0
          c2=f_zero
        endif
      endif

      if ( do_lvlcalc ) then
        do k = 1, k2
          do kp = 1, k1
            co2p10_lvl(kp,k) = c1*trns_interp_lvl_ps_nf (kp,k,1) - c2
            co2p8_lvl (kp,k) = c1*trns_interp_lvl_ps8_nf(kp,k,1) - c2
          enddo
        enddo

        if ( ntbnd == 3 ) then
          do k = 1, k2
            do kp = 1, k1
              dcdt10_lvl(kp,k) = .02 * (trns_interp_lvl_ps_nf(kp,k,2)   &
     &                         - trns_interp_lvl_ps_nf(kp,k,3))* 100.
              dcdt8_lvl(kp,k)  = .02 * (trns_interp_lvl_ps8_nf(kp,k,2)  &
     &                         - trns_interp_lvl_ps8_nf(kp,k,3))* 100.
              d2ct10_lvl(kp,k) =.0016* (trns_interp_lvl_ps_nf(kp,k,2)   &
     &                       + trns_interp_lvl_ps_nf(kp,k,3)            &
     &                       - 2.0*trns_interp_lvl_ps_nf(kp,k,1))* 1000.
              d2ct8_lvl(kp,k)  =.0016* (trns_interp_lvl_ps8_nf(kp,k,2)  &
     &                       + trns_interp_lvl_ps8_nf(kp,k,3)           &
     &                       - 2.0*trns_interp_lvl_ps_nf(kp,k,1))* 1000.
            enddo
          enddo
        endif   ! end if_ntbnd_block
      endif   ! end if_do_lvlcalc_block

      if ( do_lvlctscalc ) then
        do kp = 1, k1
          co2p10_lvlcts(kp) = c1*trns_interp_lvl_ps_nf(kp,1,1) - c2
          co2p8_lvlcts (kp) = c1*trns_interp_lvl_ps8_nf(kp,1,1) - c2
        enddo

        if ( ntbnd == 3 ) then
          do kp = 1, k1
            dcdt10_lvlcts(kp) = .02* (trns_interp_lvl_ps_nf(kp,1,2)     &
     &                        - trns_interp_lvl_ps_nf(kp,1,3))* 100.
            dcdt8_lvlcts(kp)  = .02* (trns_interp_lvl_ps8_nf(kp,1,2)    &
     &                        - trns_interp_lvl_ps8_nf(kp,1,3))* 100.
            d2ct10_lvlcts(kp) =.0016* (trns_interp_lvl_ps_nf(kp,1,2)    &
     &                      + trns_interp_lvl_ps_nf(kp,1,3)             &
     &                      - 2.0*trns_interp_lvl_ps_nf(kp,1,1))* 1000.
            d2ct8_lvlcts(kp)  =.0016* (trns_interp_lvl_ps8_nf(kp,1,2)   &
     &                      + trns_interp_lvl_ps8_nf(kp,1,3)            &
     &                      - 2.0*trns_interp_lvl_ps_nf(kp,1,1))* 1000.
          enddo
        endif   ! end if_ntbnd_block
      endif   ! end if_do_lvlctscalc_block

      if ( do_lyrcalc ) then
        do k = 1, k2
          do kp = 1, k1
            co2p10_lyr(kp,k) = c1*trns_interp_lyr_ps_nf(kp,k,1) - c2
            co2p8_lyr (kp,k) = c1*trns_interp_lyr_ps8_nf(kp,k,1) - c2
          enddo
        enddo

        if ( ntbnd == 3 ) then
          do k = 1, k2
            do kp = 1, k1
              dcdt10_lyr(kp,k) = .02* (trns_interp_lyr_ps_nf(kp,k,2)    &
     &                         - trns_interp_lyr_ps_nf(kp,k,3))* 100.
              dcdt8_lyr(kp,k)  = .02* (trns_interp_lyr_ps8_nf(kp,k,2)   &
     &                         - trns_interp_lyr_ps8_nf(kp,k,3))* 100.
              d2ct10_lyr(kp,k) =.0016* (trns_interp_lyr_ps_nf(kp,k,2)   &
     &                      + trns_interp_lyr_ps_nf(kp,k,3)             &
     &                      - 2.0*trns_interp_lyr_ps_nf(kp,k,1))* 1000.
              d2ct8_lyr(kp,k)  =.0016* (trns_interp_lyr_ps8_nf(kp,k,2)  &
     &                      + trns_interp_lyr_ps8_nf(kp,k,3)            &
     &                      - 2.0*trns_interp_lyr_ps8_nf(kp,k,1))* 1000.
            enddo
          enddo
        endif   ! end if_ntbnd_block
      endif   ! end if_do_lyrcalc_block
!
!...................................
      end subroutine gasins
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine rctrns                                                 &
!...................................
!  ---  inputs:
     &     ( gas_type, co2_std_lo, co2_std_hi, rco2_vmr, nf, nt,        &
!  ---  outputs:
     &       trns_vmr                                                   & 
     &     )

!---------------------------------------------------------------------
!      co2_std_hi = value of higher std co2 concentration in ppmv
!      co2_std_lo = value of lower std co2 concentration in ppmv
!      rco2_vmr   = value of actual co2 concentration in ppmv
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      character(*), intent(in) :: gas_type

      integer,  intent(in) :: nf,nt
      real (kind=kind_phys), intent(in) :: rco2_vmr, co2_std_lo,        &
     &       co2_std_hi

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: trns_vmr

!  --- locals:
      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS) ::      &
     &       approx_guess1, caintv, uexpintv, sexpintv, xaintv,         &
     &       press_hiv, press_lov, error_guess1, approxint_guess1,      &
     &       errorint_guess1, trans_guess1, approxint_guess2,           &
     &       errorint_guess2, trans_guess2

      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: xa,ca, uexp,sexp

      real (kind=kind_phys) :: dop_core
      logical :: do_triangle
      integer :: k, kp

!
!===> ...  begin here
!
!  ---  compute co2 transmission functions for actual co2 concentration
!       using method of section 5, Ref. (2).
!
!       the first part of the method is to obtain a first guess co2
!       transmission function for the desired concentration using only
!       the co2 tf's for the higher standard concentration.

      call coeint                                                       &
!  ---  inputs:
     &     ( gas_type, nf, trns_std_hi_nf(:,:,nt),                      &
!  ---  outputs:
     &       ca, sexp, xa, uexp, dop_core                               &
     &     )

!  --- ...  compute the interpolation.

      do_triangle = .true.

!  ---  1) compute approx function at standard (pa) pressures

      do k = 1, NSTDCO2LVLS
        do kp = k, NSTDCO2LVLS
          press_hiv(kp,k) = pa  (kp)
          press_lov(kp,k) = pa  (k)
          caintv   (kp,k) = ca  (kp)
          sexpintv (kp,k) = sexp(kp)
          xaintv   (kp,k) = xa  (kp)
          uexpintv (kp,k) = uexp(kp)
        enddo
      enddo

!  ---  the call (and calculations) to pathv2_std has been subsumed
!       into the subroutine approx_fn_std

      call approx_fn_std                                                &
!  ---  inputs:
     &     ( press_hiv, press_lov, do_triangle, caintv, sexpintv,       &
     &       xaintv, uexpintv, dop_core,                                &
!  ---  outputs:
     &       approx_guess1                                              &
     &     )

!  ---  2) compute error function at standard (pa) pressures

      do k = 1, NSTDCO2LVLS
        do kp = k+1, NSTDCO2LVLS
          error_guess1(kp,k) = 1.0 - trns_std_hi_nf(kp,k,nt) -          &
     &                         approx_guess1(kp,k)
        enddo

        error_guess1(k,k) = f_zero
      enddo

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

      if ( nf == 1 .and. nt == 1 ) then
        do k = 1, NSTDCO2LVLS
          do kp = k+1, NSTDCO2LVLS
            pressint_hiv_std_pt1(kp,k) = ((rco2_vmr+co2_std_hi)*pa(kp)  &
     &               + (co2_std_hi-rco2_vmr)*pa(k)) / (2.*co2_std_hi)
            pressint_lov_std_pt1(kp,k) = ((co2_std_hi-rco2_vmr)*pa(kp)  &
     &               + (rco2_vmr+co2_std_hi)*pa(k)) / (2.*co2_std_hi)
          enddo
        enddo
      endif

      call intcoef_2d_std                                               &
!  ---  inputs:
     &     ( pressint_hiv_std_pt1, pressint_lov_std_pt1, nf, nt,        &
     &       do_triangle, ca, xa, uexp, sexp,                           &
!  ---  outputs:
     &       indx_pressint_hiv_std_pt1, indx_pressint_lov_std_pt1,      &
     &       caintv, sexpintv, xaintv, uexpintv                         &
     &     )

!  ---  4) interpolate error function to (pressint_hiv, pressint_lov)
!          for all (k,k')

      call interp_error_r                                               &
!  ---  inputs:
     &     ( error_guess1, pressint_hiv_std_pt1,                        &
     &       pressint_lov_std_pt1, indx_pressint_hiv_std_pt1,           &
     &       indx_pressint_lov_std_pt1, do_triangle,                    &
!  ---  outputs:
     &       errorint_guess1                                            &
     &     )

!  ---  5) compute approx function for (pressint_hiv, pressint_lov)
!          the call (and calculations) to pathv2_std has been subsumed
!          into the subroutine approx_fn_std

      call approx_fn_std                                                &
!  ---  inputs:
     &     ( pressint_hiv_std_pt1, pressint_lov_std_pt1, do_triangle,   &
     &       caintv, sexpintv, xaintv, uexpintv, dop_core,              &
!  ---  outputs:
     &       approxint_guess1                                           &
     &     )

!  ---  6) compute first guess transmission function using Eq.(3), Ref.(2).

      do k = 1, NSTDCO2LVLS
        do kp = k+1, NSTDCO2LVLS
          trans_guess1(kp,k) = 1.0 -                                    &
     &         (errorint_guess1(kp,k) + approxint_guess1(kp,k))
        enddo
      enddo

!  ---  the second part of the method is to obtain a second guess co2
!       transmission function for the lower standard  concentration using
!       only the co2 tf's for the higher standard concentration. the
!       coeint call and steps (1-2) of part (1) need not be repeated.

!  ---  3) derive the pressures for interpolation using Eqs. (8a-b)
!          in Ref.(2).

      if ( nf == 1 .and. nt == 1 ) then
        do k = 1, NSTDCO2LVLS
          do kp = k+1, NSTDCO2LVLS
            pressint_hiv_std_pt2(kp,k) =                                &
     &                 ((co2_std_lo+co2_std_hi)*pa(kp) +                &
     &                 (co2_std_hi-co2_std_lo)*pa(k)) / (2.*co2_std_hi)
            pressint_lov_std_pt2(kp,k) =                                &
     &                 ((co2_std_hi-co2_std_lo)*pa(kp) +                &
     &                 (co2_std_lo+co2_std_hi)*pa(k)) / (2.*co2_std_hi)
          enddo
        enddo
      endif

      call intcoef_2d_std                                               &
!  ---  inputs:
     &     ( pressint_hiv_std_pt2, pressint_lov_std_pt2, nf, nt,        &
     &       do_triangle, ca, xa, uexp, sexp,                           &
!  ---  outputs:
     &       indx_pressint_hiv_std_pt2, indx_pressint_lov_std_pt2,      &
     &       caintv, sexpintv, xaintv, uexpintv                         &
     &     )

!  ---  4) interpolate error function to (pressint_hiv, pressint_lov)
!          for all (k,k')

      call interp_error_r                                               &
!  ---  inputs:
     &     ( error_guess1, pressint_hiv_std_pt2,                        &
     &       pressint_lov_std_pt2, indx_pressint_hiv_std_pt2,           &
     &       indx_pressint_lov_std_pt2, do_triangle,                    &
!  ---  outputs:
     &       errorint_guess2                                            &
     &     )

!  ---  5) compute approx function for (pressint_hiv, pressint_lov)
!          the call (and calculations) to pathv2_std has been subsumed
!          into the subroutine approx_fn_std

      call approx_fn_std                                                &
!  ---  inputs:
     &     ( pressint_hiv_std_pt2, pressint_lov_std_pt2, do_triangle,   &
     &       caintv, sexpintv, xaintv, uexpintv, dop_core,              &
!  ---  outputs:
     &       approxint_guess2                                           &
     &     )

!  ---  6) compute second guess transmission function using Eq.(3), Ref.(2).

      do k = 1, NSTDCO2LVLS
        do kp = k+1, NSTDCO2LVLS
          trans_guess2(kp,k) = 1.0 -                                    &
     &          (errorint_guess2(kp,k) + approxint_guess2(kp,k))
        enddo
      enddo

!  ---  finally, obtain transmission function for (co2_vmr) using
!       Eq.(9), Ref. (2).

      do k = 1, NSTDCO2LVLS
        do kp = k+1, NSTDCO2LVLS
          trns_vmr(kp,k) = trans_guess1(kp,k) +                         &
     &          (co2_std_hi - rco2_vmr)/(co2_std_hi - co2_std_lo)*      &
     &          (trns_std_lo_nf(kp,k,nt) - trans_guess2(kp,k))
          trns_vmr(k,kp) = trns_vmr(kp,k)
        enddo
        trns_vmr(k,k) = 1.0
      enddo
!
!...................................
      end subroutine rctrns
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine coeint                                                 &
!...................................
!  ---  inputs:
     &     ( gas_type, nf, trns_val,                                    &
!  ---  outputs:
     &       ca, sexp, xa, uexp, dop_core                               &
     &     )

!---------------------------------------------------------------------
!
!            the transmission function between p1 and p2 is assumed to
!        have the  functional form
!            tau(p1,p2)= 1.0-(c*log(1.0+x*path**delta))**(gamma/delta),
!            where
!               path(p1,p2)=(p1-p2)**2)*(p1+p2+dop_core)
!               and p2 is the larger of the two pressures (p1,p2).
!
!        the coefficients c and x are functions of p2, while dop_core,
!        gamma and delta are predetermined coefficients.
!        (delta,gamma are uexp,sexp in this code).
!            subroutine coeint determines c(i) and x(i) by using actual
!        values of tau(p(i-2),p(i)) and tau(p(i-1),p(i)), obtained
!        from line-by-line calculations.
!             define:
!                patha=(path(p(i),p(i-2),dop_core)**delta
!                pathb=(path(p(i),p(i-1),dop_core)**delta;
!        then
!         r=(1-tau(p(i),p(i-2)))/(1-tau(p(i),p(i-1)))
!          = (log(1+x(p(i))*patha)/log(1+x(p(i))*pathb))**(gamma/delta),
!        since   c(p(i)) cancels out
!        so that
!           r**(delta/gamma)= log(1+x(p(i))*patha)/log(1+x(p(i))*pathb).
!        this equation is solved by newton's method for x and then the
!        result used to find c. this is repeated for each value of i
!        greater than 2 to give the arrays x(i), c(i).
!        there are several possible pitfalls:
!        1) in the course of iteration, x may reach a value which makes
!           1+x*patha negative; in this case the iteration is stopped,
!           and an error message is printed out.
!        2) even if (1) does not occur, it is still possible that x may
!           be negative and large enough to make
!           1+x*path(p(i),0,dop_core) negative. this is checked in
!           a final loop, and if true,a warning is printed out.
!
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      character(len=*), intent(in) :: gas_type

      real (kind=kind_phys), dimension(:,:), intent(in) :: trns_val
      integer,              intent(in) :: nf

!  ---  outputs:
      real (kind=kind_phys), dimension(:),   intent(out) :: ca, xa,     &
     &       sexp, uexp
      real (kind=kind_phys), intent(out) :: dop_core

!  ---  locals:
      real (kind=kind_phys), parameter :: dop_core0 = 25.0

      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: upath0, upatha,  &
     &       upathb, pam1, pam2, pa0, pr_hi, r, rexp, f, f1, f2,        &
     &       fprime, ftest1, ftest2, xx, xxlog, pa2

      integer :: k, ll
      real (kind=kind_phys) :: check

!
!===> ...  begin here
!

!  ---  the following specifications for dop_core, sexp and uexp follow
!       "try9", which has (as of 5/27/97) been found to produce the
!       most accurate co2 40 level 490-850 cm-1 transmissivities, when
!       compared to LBL calculations over the same frequencies and
!       vertical structure.

      if ( gas_type == 'co2' ) then
        if ( nf == 1 ) dop_core = dop_core0
        if ( nf == 2 ) dop_core = dop_core0 * 560.0 / 670.0
        if ( nf == 3 ) dop_core = dop_core0 * 665.0 / 670.0
        if ( nf == 4 ) dop_core = dop_core0 * 775.0 / 670.0
        if ( nf == 5 ) dop_core = dop_core0 * 2325.0/ 670.0
      endif

      if ( gas_type == 'ch4' ) then
        if ( nf == 1 ) dop_core = dop_core0 * 1300.0 / 670.0
      endif

      if ( gas_type == 'n2o' ) then
        if ( nf == 1 ) dop_core = dop_core0 * 1300.0/ 670.0
        if ( nf == 2 ) dop_core = dop_core0 * 1135.0/ 670.0
        if ( nf == 3 ) dop_core = dop_core0 * 595.0 / 670.0
      endif

      do k = 1, NSTDCO2LVLS
        pa2 (k) = pa(k) * pa(k)
        sexp(k) = 0.505 + 2.0e-5*pa(k)                                  &
     &          + 0.035*(pa2(k) - 0.25) / (pa2(k) + 0.25)
        uexp(k) = sexp(k)*(1.0 + 0.33*pa2(k)/(pa2(k) + 40000.0))
      enddo

      do k = 1, NSTDCO2LVLS
        pr_hi(k) = pa(k)
      enddo

      do k = 3, NSTDCO2LVLS
        pam1(k) = pa(k-1)
        pam2(k) = pa(k-2)
        pa0 (k) = f_zero
      enddo

      call pathv1                                                       &
!  ---  inputs:
     &     ( pr_hi, pam1, sexp, dop_core, 3, NSTDCO2LVLS,               &
!  ---  outputs:
     &       upathb                                                     &
     &     )

      call pathv1                                                       &
!  ---  inputs:
     &     ( pr_hi, pam2, sexp, dop_core, 3, NSTDCO2LVLS,               &
!  ---  outputs:
     &       upatha                                                     &
     &     )

      do k = 3, NSTDCO2LVLS
        r(k) = (1.0 -trns_val(k,k-2))/(1.0 -trns_val(k,k-1))
        rexp(k) = r(k)**(uexp(k)/sexp(k))
        upatha(k) = upatha(k)**uexp(k)
        upathb(k) = upathb(k)**uexp(k)
        xx(k) = 2.0*(upathb(k)*rexp(k) - upatha(k))/                    &
     &            (upathb(k)*upathb(k)*rexp(k) - upatha(k)*upatha(k))
      enddo

      do ll = 1, 20
        do k = 3, NSTDCO2LVLS
          ftest1(k) =xx(k)*upatha(k)
          ftest2(k) =xx(k)*upathb(k)

          if ( ftest1(k) <= 1.0e-10 ) then
!  ---  end iteration and solve if ftest1 is small or ftest2 is large
            xa(k)=1.0
            ca(k)=(1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/upatha(k)
          elseif ( ftest2(k) >= 1.0e+8 ) then
            xxlog(k) = (log(upatha(k)) - rexp(k)*LOG(upathb(k)))/       &
     &                 (rexp(k)-1.0 )
            xa(k) = exp(xxlog(k))
            ca(k) = (1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/         &
     &                 (xxlog(k) + log(upatha(k)))
          else
            f1(k) = LOG(1.0 + xx(k)*upatha(k))
            f2(k) = LOG(1.0 + xx(k)*upathb(k))
            f(k) = f1(k)/f2(k) - rexp(k)
            fprime(k) = (f2(k)*upatha(k)/(1.0 + xx(k)*upatha(k)) -      &
     &                   f1(k)*upathb(k)/(1.0 + xx(k)*upathb(k)))/      &
     &                   (f2(k)*f2(k))
            xx(k) = xx(k) - f(k)/fprime(k)
          endif
        enddo
      enddo

      do k = 3, NSTDCO2LVLS
        if ( ftest1(k) > 1.0e-10 .and. ftest2(k) < 1.0e+8 ) then
          ca(k) = (1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/           &
     &            (log(1.0 + xx(k)*upatha(k)) + 1.0e-20)
          xa(k) = xx(k)
        endif
      enddo

!  ---  by assumption, ca, xa for the first two levels  are
!       equal to the values for the third level.

      xa(2)=xa(3)
      xa(1)=xa(3)
      ca(2)=ca(3)
      ca(1)=ca(3)
!
!...................................
      end subroutine coeint
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine approx_fn                                              &
!...................................
!  ---  inputs:
     &     ( press_hi_app, press_lo_app, do_triangle,                   &
     &       ca_app, sexp_app, xa_app, uexp_app, dop_core,              &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       approx                                                     &
     &     )

!---------------------------------------------------------------------
!   approx_fn computes the co2 approximation function
!          A(press_hi_app(i), press_lo_app(j))  (Eq.(4), Ref. (2))
!   for a particular co2 amount and a standard pressure grid (pa).
!   the calculation is performed for all (press_hi(k),press_lo(k')
!   pairs which are possible according to the state of (do_triangle).
!   the path function (upathv) is evaluated using the expression
!   given in Eqs. (5) and (A5) in Ref. (2) for the co2 interpolation
!   program between a lower model pressure (press_lo_app) and
!   a higher model pressure (press_hi_app) using the interpolation
!   coefficients (ca_app, sexp_app, xa_app, uexp_app) computed in
!   subroutine coeint.
!         the output is in (approx).
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: nklo, nkhi, nkplo, nkphi
      logical, intent(in) :: do_triangle

      real (kind=kind_phys), dimension(:,:), intent(in) ::              &
     &       press_hi_app, press_lo_app, ca_app, sexp_app, xa_app,      &
     &       uexp_app

      real (kind=kind_phys), intent(in) :: dop_core

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: approx

!  ---  locals:
      integer :: k, kp, kp0, k1, k2

      real (kind=kind_phys), dimension(:,:), allocatable :: upathv

!
!===> ...  begin here
!

!  ---   obtain array extents for internal arrays

      k1 = size(press_hi_app, 1)       ! this corresponds to ndimkp
      k2 = size(press_hi_app, 2)       ! this corresponds to ndimk

!  ---   and allocate these arrays

      if ( .not. allocated(upathv) ) allocate ( upathv(k1,k2) )

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          upathv(kp,k) = (press_hi_app(kp,k)                            &
     &                 - press_lo_app(kp,k))**(1./sexp_app(kp,k))       &
     &                 * (press_hi_app(kp,k) + press_lo_app(kp,k)       &
     &                 + dop_core)

          upathv(kp,k) = upathv(kp,k)**uexp_app(kp,k)
          approx(kp,k) = (ca_app(kp,k)                                  &
     &                 * LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))**        &
     &                   (sexp_app(kp,k)/uexp_app(kp,k))
        enddo
      enddo
!
!...................................
      end subroutine approx_fn
!----------------------------------- containd in gases_stdtf


!-----------------------------------
      subroutine approx_fn_std                                          &
!...................................
!  ---  inputs:
     &     ( press_hi_pa, press_lo_pa, do_triangle, ca_app,             &
     &       sexp_app, xa_app, uexp_app, dop_core,                      &
!  ---  outputs:
     &       approx                                                     &
     &     )

!---------------------------------------------------------------------
!   approx_fn_std computes the co2 approximation function
!               A(press_hi(i), press_lo(j))  (Eq.(4), Ref. (2))
!   for a particular co2 amount and a standard pressure grid (pa).
!   the calculation is performed for all (press_hi(k),press_lo(k')
!   pairs which are possible according to the state of (do_triangle).
!   the path function (upathv) is evaluated using the expression
!   given in Eqs. (5) and (A5) in Ref. (2) for the co2 interpolation
!   program between a lower standard pressure (press_lo_pa) and
!   a higher standard pressure (press_hi_pa) using the interpolation
!   coefficients (ca_app, sexp_app, xa_app, uexp_app) computed in
!   subroutine coeint.
!         the output is in (approx).
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      real (kind=kind_phys), dimension(:,:), intent(in) :: ca_app,      &
     &       sexp_app, xa_app, uexp_app, press_hi_pa, press_lo_pa

      real (kind=kind_phys), intent(in) :: dop_core
      logical, intent(in) :: do_triangle

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: approx

!  ---  locals:
      integer :: k, kp, kp0
      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS) ::upathv

!
!===> ...  begin here
!
      do k = 1, NSTDCO2LVLS
        if ( do_triangle ) then
          kp0 = k + 1
        else
          kp0 = 1
        endif

        do kp = kp0, NSTDCO2LVLS
          upathv(kp,k) = (press_hi_pa(kp,k)                             &
     &                 - press_lo_pa(kp,k))**(1./sexp_app(kp,k))        &
     &                 * (press_hi_pa(kp,k) + press_lo_pa(kp,k)         &
     &                 + dop_core)

          upathv(kp,k) = upathv(kp,k)**uexp_app(kp,k)
          approx(kp,k) = (ca_app(kp,k)                                  &
     &                 * LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))**        &
     &                   (sexp_app(kp,k)/uexp_app(kp,k))
        enddo
      enddo
!
!...................................
      end subroutine approx_fn_std
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine intcoef_2d                                             &
!...................................
!  ---  inputs:
     &     ( press_hiv, press_lov, do_triangle, ca, xa,                 &
     &       uexp, sexp, nklo, nkhi, nkplo, nkphi,                      &
!  ---  outputs:
     &       indx_hiv, indx_lov, caintv,                                &
     &       sexpintv, xaintv, uexpintv                                 &
     &     )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: nklo, nkhi, nkplo, nkphi
      logical, intent(in) :: do_triangle

      real (kind=kind_phys), dimension(:,:), intent(in) :: press_hiv,   &
     &       press_lov

      real (kind=kind_phys), dimension(:),   intent(in) :: ca, xa,      &
     &       uexp, sexp

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: sexpintv,   &
     &       uexpintv, caintv, xaintv
      integer, dimension(:,:), intent(out) :: indx_hiv, indx_lov

!  ---  locals:
      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: caxa
      real (kind=kind_phys), dimension(:,:), allocatable :: sexp_hiv,   &
     &       uexp_hiv, ca_hiv, prod_hiv, xa_hiv, d1kp, d2kp, bkp, akp,  &
     &       delp_hi

      integer :: k, kp, kp0, kpp, k1, k2

!
!===> ...  begin here
!
!  ---  obtain array extents for internal arrays

      k1 = size(press_hiv, 1)       ! this corresponds to ndimkp
      k2 = size(press_hiv, 2)       ! this corresponds to ndimk

!  ---  and allocate these arrays

      if ( .not. allocated(sexp_hiv) ) then
        allocate (sexp_hiv(k1,k2), uexp_hiv(k1,k2), ca_hiv(k1,k2),      &
     &            prod_hiv(k1,k2), xa_hiv  (k1,k2), d1kp  (k1,k2),      &
     &            d2kp    (k1,k2), bkp     (k1,k2), akp   (k1,k2),      &
     &            delp_hi (k1,k2) )
      endif

!  ---  compute the index of the inputted pressures (press_hiv,
!       press_lov) corresponding to the standard (pa) pressures.

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          if ( press_hiv(kp,k) < pa(1) ) then
            indx_hiv(kp,k) = 1
          endif

          if ( press_hiv(kp,k) >= pa(NSTDCO2LVLS) ) then
            indx_hiv(kp,k) = NSTDCO2LVLS - 1
          endif

          if ( press_lov(kp,k) < pa(1) ) then
            indx_lov(kp,k) = 1
          endif

          if ( press_lov(kp,k) >= pa(NSTDCO2LVLS) ) then
            indx_lov(kp,k) = NSTDCO2LVLS - 1
          endif
        enddo
      enddo

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          do kpp = 1, NSTDCO2LVLS-1
            if ( press_hiv(kp,k) >= pa(kpp) .and.                       &
     &           press_hiv(kp,k) <  pa(kpp+1) ) then
              indx_hiv(kp,k) = kpp
              exit
            endif
          enddo

          do kpp = 1, NSTDCO2LVLS-1
            if ( press_lov(kp,k) >= pa(kpp) .and.                       &
     &           press_lov(kp,k) <  pa(kpp+1) ) then
              indx_lov(kp,k) = kpp
              exit
            endif
          enddo
        enddo
      enddo

!  ---  interpolate values of cint, xint, sexp, uexp for the pressures
!       (press_hiv)

      do k = 1, NSTDCO2LVLS
        caxa(k) = ca(k) * xa(k)
      enddo

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          sexp_hiv(kp,k) = sexp(indx_hiv(kp,k)) +                       &
     &        (sexp(indx_hiv(kp,k)+1) - sexp(indx_hiv(kp,k))) /         &
     &        (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *         &
     &        (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          uexp_hiv(kp,k) = uexp(indx_hiv(kp,k)) +                       &
     &        (uexp(indx_hiv(kp,k)+1) - uexp(indx_hiv(kp,k))) /         &
     &        (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *         &
     &        (press_hiv(kp,k) - pa(indx_hiv(kp,k)))

          if (indx_hiv(kp,k)>2 .and. indx_hiv(kp,k)<NSTDCO2LVLS-1) then

!  ---  use 3-point interpolation: (indx_hiv of 1 or 2 are excluded
!       since ca and xa were arbitrarily set to ca(3),xa(3))

            delp_hi(kp,k) = press_hiv(kp,k) - pa(indx_hiv(kp,k)+1)

!  ---  interpolate xa

            d1kp(kp,k) = (xa(indx_hiv(kp,k)+2) - xa(indx_hiv(kp,k)+1))  &
     &                 / (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp,k) = (xa(indx_hiv(kp,k)+1) -  xa(indx_hiv(kp,k) ))  &
     &                 / (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))                       &
     &                / (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)                          &
     &                * (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            xa_hiv(kp,k) = xa(indx_hiv(kp,k)+1) +                       &
     &            delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))

!  ---  if xa_hiv is negative or zero, the interpolation fails and
!       the model may bomb. to avoid this, use 2-point interpolation
!       in this case. the 3-point interpolation for prod_hiv is
!       stable, so there is no need to change this calculation.

            if ( xa_hiv(kp,k) <= f_zero ) then
              xa_hiv(kp,k) = xa(indx_hiv(kp,k)) +                       &
     &           (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /          &
     &           (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k))) *          &
     &           (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            endif

!  ---  interpolate caxa

            d1kp(kp,k) = (caxa(indx_hiv(kp,k)+2)-caxa(indx_hiv(kp,k)+1))&
     &                 / (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp,k) = (caxa(indx_hiv(kp,k)+1)-caxa(indx_hiv(kp,k) )) &
     &                 / (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))                       &
     &                 / (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)                          &
     &                 * (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)+1) +                   &
     &          delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))

          else
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)) +                     &
     &          (caxa(indx_hiv(kp,k)+1) - caxa(indx_hiv(kp,k))) /       &
     &          (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *       &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            xa_hiv(kp,k) = xa(indx_hiv(kp,k)) +                         &
     &          (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /           &
     &          (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k))) *           &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          endif

          ca_hiv(kp,k) = prod_hiv(kp,k)/xa_hiv(kp,k)
        enddo
      enddo

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          sexpintv(kp,k) = sexp_hiv(kp,k)
          uexpintv(kp,k) = uexp_hiv(kp,k)
          caintv  (kp,k) = ca_hiv(kp,k)
          xaintv  (kp,k) = xa_hiv(kp,k)
        enddo
      enddo
!
!...................................
      end subroutine intcoef_2d
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine intcoef_2d_std                                         &
!...................................
!  ---  inputs:
     &     ( press_hiv, press_lov, nf, nt,                              &
     &       do_triangle, ca, xa, uexp, sexp,                           &
!  ---  outputs:
     &       indx_hiv, indx_lov, caintv,                                &
     &       sexpintv, xaintv, uexpintv                                 &
     &     )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: nf, nt
      logical, intent(in) :: do_triangle

      real (kind=kind_phys), dimension(:,:), intent(in) :: press_hiv,   &
     &       press_lov

      real (kind=kind_phys), dimension(:),   intent(in) :: ca, xa,      &
     &       uexp, sexp

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: sexpintv,   &
     &       uexpintv, caintv, xaintv
      integer, dimension(:,:), intent(out) :: indx_hiv, indx_lov

!  ---  locals:
      real (kind=kind_phys), dimension(NSTDCO2LVLS,NSTDCO2LVLS) ::      &
     &       prod_hiv
      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: d1kp, d2kp,      &
     &       bkp, akp, delp_hi, caxa

      integer :: k, kp, kp0, kpp

!
!===> ...  begin here
!
!  ---  compute the index of the inputted pressures (press_hiv,
!       press_lov) corresponding to the standard (pa) pressures.
!       (only calculate if nf = 1, nt = 1)

      if ( nf == 1 .and. nt == 1 ) then
        do k = 1, NSTDCO2LVLS
          if ( do_triangle ) then
            kp0 = k + 1
          else
            kp0 = 1
          endif

          do kp = kp0, NSTDCO2LVLS
            if ( press_hiv(kp,k) < pa(1) ) then
              indx_hiv(kp,k) = 1
            endif

            if ( press_hiv(kp,k) >= pa(NSTDCO2LVLS) ) then
              indx_hiv(kp,k) = NSTDCO2LVLS - 1
            endif

            if ( press_lov(kp,k) < pa(1) ) then
              indx_lov(kp,k) = 1
            endif

            if ( press_lov(kp,k) >= pa(NSTDCO2LVLS) ) then
              indx_lov(kp,k) = NSTDCO2LVLS - 1
            endif
          enddo
        enddo

        do k = 1, NSTDCO2LVLS
          if ( do_triangle ) then
            kp0 = k + 1
          else
            kp0 = 1
          endif

          do kp = kp0, NSTDCO2LVLS
            do kpp = 1, NSTDCO2LVLS-1
              if ( press_hiv(kp,k) >= pa(kpp) .and.                     &
     &             press_hiv(kp,k) <  pa(kpp+1) ) then
                indx_hiv(kp,k) = kpp
                exit
              endif
            enddo

            do kpp = 1, NSTDCO2LVLS-1
              if ( press_lov(kp,k) >= pa(kpp) .and.                     &
     &             press_lov(kp,k) <  pa(kpp+1) ) then
                indx_lov(kp,k) = kpp
                exit
              endif
            enddo
          enddo
        enddo
      endif

!  ---  interpolate values of cint, xint, sexp, uexp for the pressures
!       (press_hiv) (for all values of nf, nt)

      do k = 1, NSTDCO2LVLS
        caxa(k) = ca(k)*xa(k)
      enddo

      do k = 1, NSTDCO2LVLS
        if ( do_triangle ) then
          kp0 = k + 1
        else
          kp0 = 1
        endif

        do kp = kp0, NSTDCO2LVLS
          sexpintv(kp,k) = sexp(indx_hiv(kp,k)) +                       &
     &          (sexp(indx_hiv(kp,k)+1) - sexp(indx_hiv(kp,k))) /       &
     &          (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *       &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          uexpintv(kp,k) = uexp(indx_hiv(kp,k)) +                       &
     &          (uexp(indx_hiv(kp,k)+1) - uexp(indx_hiv(kp,k))) /       &
     &          (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *       &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))

          if (indx_hiv(kp,k)>2 .and. indx_hiv(kp,k)<NSTDCO2LVLS-1) then

!  ---  use 3-point interpolation: (indx_hiv of 1 or 2 are excluded
!       since ca and xa were arbitrarily set to ca(3),xa(3))

            delp_hi(kp) = press_hiv(kp,k) - pa(indx_hiv(kp,k)+1)

!  ---  interpolate xa

            d1kp(kp) = (xa(indx_hiv(kp,k)+2) - xa(indx_hiv(kp,k)+1)) /  &
     &                 (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp) = (xa(indx_hiv(kp,k)+1) -  xa(indx_hiv(kp,k) )) /  &
     &                 (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/                            &
     &                (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp) = d1kp(kp) - bkp(kp)*                               &
     &                (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            xaintv(kp,k) = xa(indx_hiv(kp,k)+1) +                       &
     &                delp_hi(kp)*(akp(kp) + delp_hi(kp)*bkp(kp))

!  ---  if xaintv is negative or zero, the interpolation fails and
!       the model may bomb. to avoid this, use 2-point interpolation
!       in this case. the 3-point interpolation for prod_hiv is
!       stable, so there is no need to change this calculation.

            if ( xaintv(kp,k) <= f_zero ) then
              xaintv(kp,k) = xa(indx_hiv(kp,k)) +                       &
     &                (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /     &
     &                (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k))) *     &
     &                (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            endif

!  ---  interpolate caxa

            d1kp(kp) = (caxa(indx_hiv(kp,k)+2)-caxa(indx_hiv(kp,k)+1))  &
     &               / (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp) = (caxa(indx_hiv(kp,k)+1)-caxa(indx_hiv(kp,k) ))   &
     &               / (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/                            &
     &                (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp) = d1kp(kp) - bkp(kp)*                               &
     &                (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)+1) +                   &
     &                delp_hi(kp)*(akp(kp) + delp_hi(kp)*bkp(kp))
          else
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)) +                     &
     &          (caxa(indx_hiv(kp,k)+1) - caxa(indx_hiv(kp,k))) /       &
     &          (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *       &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            xaintv(kp,k) = xa(indx_hiv(kp,k)) +                         &
     &          (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /           &
     &          (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k))) *           &
     &          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          endif

          caintv(kp,k) = prod_hiv(kp,k)/xaintv(kp,k)

        enddo
      enddo
!
!...................................
      end subroutine intcoef_2d_std
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine interp_error                                           &
!...................................
!  ---  inputs:
     &     ( error, pressint_hiv, pressint_lov,                         &
     &       indx_press_hiv, indx_press_lov, do_triangle,               &
     &       nklo, nkhi, nkplo, nkphi,                                  &
!  ---  outputs:
     &       errorint                                                   &
     &     )

!---------------------------------------------------------------------
!     press_hiv = pressure on std pa grid of high (kp) pressure
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     error = error ot standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!     errorint = error at interpolated grid
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: nklo, nkhi, nkplo, nkphi
      logical, intent(in) :: do_triangle

      real (kind=kind_phys), dimension(:,:), intent(in) :: error,       &
     &       pressint_hiv, pressint_lov

      integer, dimension(:,:), intent(in) :: indx_press_hiv,            &
     &                                       indx_press_lov

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: errorint

!  ---  locals:
      integer :: k, kp, kp0, k1, k2

      real (kind=kind_phys), dimension(:,:), allocatable :: delp_lo,    &
     &       delp_hi, d1kp, d2kp, bkp, akp, fkp, fkp1, fkp2

!
!===> ... begin here
!
!  ---  obtain array extents for internal arrays

      k1 = size(pressint_hiv, 1)       ! this corresponds to ndimkp
      k2 = size(pressint_hiv, 2)       ! this corresponds to ndimk

!  ---  and allocate these arrays

      if ( .not. allocated(delp_lo) ) then
        allocate( delp_lo(k1,k2), delp_hi(k1,k2), d1kp(k1,k2),          &
     &            d2kp   (k1,k2), bkp    (k1,k2), akp (k1,k2),          &
     &            fkp    (k1,k2), fkp1   (k1,k2), fkp2(k1,k2) )
      endif

      do k = nklo, nkhi
        if ( do_triangle ) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif

        do kp = kp0, nkphi
          if ( indx_press_hiv(kp,k)-indx_press_lov(kp,k) >= 3 .and.     &
     &         indx_press_hiv(kp,k) < NSTDCO2LVLS-1 ) then

!  ---  use quadratic interpolation:

            delp_lo(kp,k) = pressint_lov(kp,k)                          &
     &                    - pa(indx_press_lov(kp,k)+1)

!  ---  1) for fixed (kp), varying (k)

            d1kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+2) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1)  ) /  &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  )  ) /  &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/                      &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*                         &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp(kp,k) =                                                 &
     &        error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) +      &
     &        delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!  ---  2) for fixed (kp+1), varying (k)

            d1kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+2) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) ) / &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/                      &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k) ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*                         &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp1(kp,k) =                                                &
     &        error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) +    &
     &        delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!  ---  3) for fixed (kp+2), varying (k)

            d1kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+2) -   &
     &         error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) ) / &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/                      &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k) ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*                         &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp2(kp,k) =                                                &
     &        error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) +    &
     &        delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!  ---  4) finally, varying (kp) using (fkp,fkp1,fkp2)

            delp_hi(kp,k) =                                             &
     &        pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k)+1)
            d1kp(kp,k) = (fkp2(kp,k) - fkp1(kp,k)) /                    &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            d2kp(kp,k) = (fkp1(kp,k) - fkp (kp,k)) /                    &
     &        (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)+0))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/                      &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*                         &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            errorint(kp,k) = fkp1(kp,k) +                               &
     &        delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))

          elseif (indx_press_hiv(kp,k) .GT. indx_press_lov(kp,k)) then

!  ---  use linear interpolation:

            delp_lo(kp,k) =                                             &
     &        pressint_lov(kp,k) - pa(indx_press_lov(kp,k))

!  ---  1) for fixed (kp), varying (k)

            d2kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  ) ) /   &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            fkp(kp,k) =                                                 &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)) +       &
     &         delp_lo(kp,k)*d2kp(kp,k)

!  ---  2) for fixed (kp+1), varying (k)

            d2kp(kp,k) =                                                &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            fkp1(kp,k) =                                                &
     &        error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)) +      &
     &        delp_lo(kp,k)*d2kp(kp,k)

!  ---  3) linear interpolate (fkp,fkp1):

            errorint(kp,k) = (fkp(kp,k)*                                &
     &        (pa(indx_press_hiv(kp,k)+1) - pressint_hiv(kp,k)) +       &
     &                        fkp1(kp,k)*                               &
     &        (pressint_hiv(kp,k)-pa(indx_press_hiv(kp,k))) ) /         &
     &        (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)))
          else

!  ---  the error function for closely-spaced pressures equals zero
!       (section 3.2, Ref. (2))

            errorint(kp,k) = f_zero
          endif
        enddo
      enddo
!
!...................................
      end subroutine interp_error
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine interp_error_r                                         &
!...................................
!  ---  inputs:
     &     ( error, pressint_hiv, pressint_lov,                         &
     &       indx_press_hiv, indx_press_lov, do_triangle,               &
!  ---  outputs:
     &       errorint                                                   &
     &     )

!---------------------------------------------------------------------
!     press_hiv = pressure on std pa grid of high (kp) pressure
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     error = error at standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!     errorint = error at interpolated grid
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      logical, intent(in) :: do_triangle

      real (kind=kind_phys), dimension(:,:), intent(in) :: error,       &
     &       pressint_hiv, pressint_lov

      integer, dimension(:,:), intent(in) :: indx_press_hiv,            &
     &                                       indx_press_lov

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: errorint

!  ---  locals:
      real (kind=kind_phys), dimension(NSTDCO2LVLS) :: delp_lo, delp_hi,&
     &       d1kp, d2kp, bkp, akp, fkp, d1kp1, d2kp1, bkp1, akp1, fkp1, &
     &       d1kp2, d2kp2, bkp2, akp2, fkp2, d1kpf, d2kpf, bkpf, akpf

      integer :: k, kp, kp0

!
!===> ...  begin here
!
      do k = 1, NSTDCO2LVLS
        if ( do_triangle ) then
          kp0 = k + 1
        else
          kp0 = 1
        endif

        do kp = kp0, NSTDCO2LVLS
          if ( indx_press_hiv(kp,k)-indx_press_lov(kp,k) >= 3 .and.     &
     &         indx_press_hiv(kp,k) < NSTDCO2LVLS-1 ) then

!  ---  use quadratic interpolation:

            delp_lo(kp) = pressint_lov(kp,k)-pa(indx_press_lov(kp,k)+1)

!  ---  1) for fixed (kp), varying (k)

            d1kp(kp) =                                                  &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+2) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) ) /   &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp) =                                                  &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  ) ) /   &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/                            &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k) ))
            akp(kp) = d1kp(kp) - bkp(kp)*                               &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp(kp) =                                                   &
     &        error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) +      &
     &        delp_lo(kp)*(akp(kp) + delp_lo(kp)*bkp(kp))

!  ---  2) for fixed (kp+1), varying (k)

            d1kp1(kp) =                                                 &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+2) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) ) / &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp1(kp) =                                                 &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            bkp1(kp) = (d1kp1(kp) - d2kp1(kp))/                         &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k) ))
            akp1(kp) = d1kp1(kp) - bkp1(kp)*                            &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp1(kp) =                                                  &
     &        error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) +    &
     &        delp_lo(kp)*(akp1(kp) + delp_lo(kp)*bkp1(kp))

!  ---  3) for fixed (kp+2), varying (k)

            d1kp2(kp) =                                                 &
     &        (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+2) -   &
     &         error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) ) / &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp2(kp) =                                                 &
     &        (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            bkp2(kp) = (d1kp2(kp) - d2kp2(kp))/                         &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k) ))
            akp2(kp) = d1kp2(kp) - bkp2(kp)*                            &
     &        (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp2(kp) =                                                  &
     &        error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) +    &
     &        delp_lo(kp)*(akp2(kp) + delp_lo(kp)*bkp2(kp))

!  ---  4) finally, varying (kp) using (fkp,fkp1,fkp2)

            delp_hi(kp) =                                               &
     &        pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k)+1)
            d1kpf(kp) = (fkp2(kp) - fkp1(kp)) /                         &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            d2kpf(kp) = (fkp1(kp) - fkp (kp)) /                         &
     &        (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)+0))
            bkpf(kp) = (d1kpf(kp) - d2kpf(kp))/                         &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)  ))
            akpf(kp) = d1kpf(kp) - bkpf(kp)*                            &
     &        (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            errorint(kp,k) = fkp1(kp) +                                 &
     &        delp_hi(kp)*(akpf(kp) + delp_hi(kp)*bkpf(kp))

          elseif (indx_press_hiv(kp,k) .GT. indx_press_lov(kp,k)) then

!  ---   use linear interpolation:

            delp_lo(kp) = pressint_lov(kp,k)-pa(indx_press_lov(kp,k))

!  ---  1) for fixed (kp), varying (k)

            d2kp(kp) =                                                  &
     &        (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -     &
     &         error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  ) ) /   &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            fkp(kp) =                                                   &
     &        error(indx_press_hiv(kp,k),indx_press_lov(kp,k)) +        &
     &        delp_lo(kp)*d2kp(kp)

!  ---  2) for fixed (kp+1), varying (k)

            d2kp1(kp) =                                                 &
     &        (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -   &
     &         error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  ) ) / &
     &        (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            fkp1(kp) =                                                  &
     &        error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)) +      &
     &        delp_lo(kp)*d2kp1(kp)

!  ---  3) linear interpolate (fkp,fkp1):

            errorint(kp,k) = ( fkp(kp)*                                 &
     &        (pa(indx_press_hiv(kp,k)+1) - pressint_hiv(kp,k))         &
     &                     +   fkp1(kp) *                               &
     &        (pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k))) ) /       &
     &        (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)))
          else

!  ---  the error function for closely-spaced pressures equals zero
!       (section 3.2, Ref. (2))

            errorint(kp,k) = f_zero
          endif
        enddo
      enddo
!
!...................................
      end subroutine interp_error_r
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine pathv1                                                 &
!...................................
!  ---  inputs:
     &     ( press_hi, press_lo, sexp, dop_core,  ndimlo, ndimhi,       &
!  ---  outputs:
     &       upath                                                      &
     &     )

!---------------------------------------------------------------------
!     pathv1 computes the path function given in Eqs. (5) and (A5) in
!     Ref. (2) for the co2 interpolation pgm. between a
!     pressure (press_lo) and a variable pressure (press_hi). This
!     has been modified on 5/27/97.
!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      integer, intent(in) :: ndimlo, ndimhi

      real (kind=kind_phys), dimension(:), intent(in) :: press_hi,      &
     &       press_lo, sexp
      real (kind=kind_phys), intent(in) :: dop_core

!  ---  outputs:
      real (kind=kind_phys), dimension (:), intent(out) :: upath

!  ---  locals:
      integer :: k

!
!===> ...  begin here
!
      do k = ndimlo, ndimhi
        upath(k) = (press_hi(k) - press_lo(k))**(1./sexp(k))*           &
     &             (press_hi(k) + press_lo(k) + dop_core)
      enddo
!
!...................................
      end subroutine pathv1
!----------------------------------- contained in gases_stdtf


!-----------------------------------
      subroutine read_lbltfs                                            &
!...................................
!  ---  inputs:
     &     ( gas_type, callrctrns, nstd_lo, nstd_hi, nf, ntbnd,         &
!  ---  outputs:
     &       trns_std_hi_nf, trns_std_lo_nf                             &
     &     )

!---------------------------------------------------------------------

      implicit none

!  ---  inputs:
      character(len=*), intent(in) :: gas_type

      logical, intent(in) :: callrctrns
      integer, intent(in) :: nstd_lo, nstd_hi, nf
      integer, dimension(:), intent(in) :: ntbnd

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       trns_std_hi_nf, trns_std_lo_nf

!  ---  locals:
      character(len=17) :: input_lblco2name(nfreq_bands_sea_co2,9)
      character(len=17) :: input_lblch4name(nfreq_bands_sea_ch4,5)
      character(len=17) :: input_lbln2oname(nfreq_bands_sea_n2o,4)

      data  input_lblco2name(:,1)  /                                    &
     &  'cns_165_490850   ', 'cns_165_490630   ', 'cns_165_630700   ',  &
     &  'cns_165_700850   ', 'cns_165_43um     ' /
      data  input_lblco2name(:,2)  /                                    &
     &  'cns_300_490850   ', 'cns_300_490630   ', 'cns_300_630700   ',  &
     &  'cns_300_700850   ', 'cns_300_43um     ' /
      data  input_lblco2name(:,3)  /                                    &
     &  'cns_330_490850   ', 'cns_330_490630   ', 'cns_330_630700   ',  &
     &  'cns_330_700850   ', 'cns_330_43um     ' /
      data  input_lblco2name(:,4)  /                                    &
     &  'cns_348_490850   ', 'cns_348_490630   ', 'cns_348_630700   ',  &
     &  'cns_348_700850   ', 'cns_348_43um     ' /
      data  input_lblco2name(:,5)  /                                    &
     &  'cns_356_490850   ', 'cns_356_490630   ', 'cns_356_630700   ',  &
     &  'cns_356_700850   ', 'cns_356_43um     ' /
      data  input_lblco2name(:,6)  /                                    &
     &  'cns_360_490850   ', 'cns_360_490630   ', 'cns_360_630700   ',  &
     &  'cns_360_700850   ', 'cns_360_43um     ' /
      data  input_lblco2name(:,7)  /                                    &
     &  'cns_600_490850   ', 'cns_600_490630   ', 'cns_600_630700   ',  &
     &  'cns_600_700850   ', 'cns_600_43um     ' /
      data  input_lblco2name(:,8)  /                                    &
     &  'cns_660_490850   ', 'cns_660_490630   ', 'cns_660_630700   ',  &
     &  'cns_660_700850   ', 'cns_660_43um     ' /
      data  input_lblco2name(:,9)  /                                    &
     &  'cns_1320_490850  ', 'cns_1320_490630  ', 'cns_1320_630700  ',  &
     &  'cns_1320_700850  ', 'cns_1320_43um    ' /
!
      data  input_lblch4name(:,:)  /                                    &
     &  'cns_700_12001400 ', 'cns_1250_12001400', 'cns_1750_12001400',  &
     &  'cns_2250_12001400', 'cns_2800_12001400' /
!
      data  input_lbln2oname(:,1)  /                                    &
     &  'cns_275_12001400 ', 'cns_275_10701200 ', 'cns_275_560630   ' /
      data  input_lbln2oname(:,2)  /                                    &
     &  'cns_310_12001400 ', 'cns_310_10701200 ', 'cns_310_560630   ' /
      data  input_lbln2oname(:,3)  /                                    &
     &  'cns_340_12001400 ', 'cns_340_10701200 ', 'cns_340_560630   ' /
      data  input_lbln2oname(:,4)  /                                    &
     &  'cns_375_12001400 ', 'cns_375_10701200 ', 'cns_375_560630   ' /

      character(len=17) :: name_lo, name_hi

      integer :: n, nt, inrad, nrec_inhi, nrec_inlo
      real (kind=kind_phys),dimension(NSTDCO2LVLS,NSTDCO2LVLS) ::trns_in

!
!===> ...  begin here
!

!  ---  input format for lbl transmission fctns (4f20.14)

      if ( gas_type == 'co2' ) then
        name_lo = input_lblco2name(nf,nstd_lo)
        name_hi = input_lblco2name(nf,nstd_hi)
      endif

      if ( gas_type == 'ch4' ) then
        name_lo = input_lblch4name(nf,nstd_lo)
        name_hi = input_lblch4name(nf,nstd_hi)
      endif

      if ( gas_type == 'n2o' ) then
        name_lo = input_lbln2oname(nf,nstd_lo)
        name_hi = input_lbln2oname(nf,nstd_hi)
      endif

!  ---  read in tfs of higher std gas concentration

      inrad = 50
      open (inrad, file = name_hi, access = 'direct',                   &
     &      recl = NSTDCO2LVLS*NSTDCO2LVLS*8)
      nrec_inhi = 0

      do nt = 1, ntbnd(nf)
        nrec_inhi = nrec_inhi + 1
        read (inrad, rec = nrec_inhi) trns_in
        trns_std_hi_nf(:,:,nt) = trns_in(:,:)
      enddo

      close (inrad)

!  ---  if necessary, read in tfs of lower standard gas concentration

      if ( callrctrns ) then

        open (inrad, file = name_lo, access = 'direct',                 &
     &        recl = NSTDCO2LVLS*NSTDCO2LVLS*8)
        nrec_inlo = 0

        do nt = 1, ntbnd(nf)
          nrec_inlo = nrec_inlo + 1
          read (inrad, rec = nrec_inlo) trns_in
          trns_std_lo_nf(:,:,nt) = trns_in(:,:)
        enddo

        close (inrad)

      endif
!
!...................................
      end subroutine read_lbltfs
!----------------------------------- contained in gases_stdtf
!
!...................................
      end subroutine gases_stdtf
!-----------------------------------



!  =========================================
!  *****  radiation_utilities section  *****
!  =========================================


!-----------------------------------
      subroutine locate_in_table                                        &
!...................................
!  ---  inputs:
     &     ( table_axis, x, k_min, k_max,                               &
!  ---  outputs:
     &       dx, ix                                                     &
     &     )
! -------------------------------------------------------------------- !
!                                                                      !
!     given array x and an arithmetic sequence of table column headings!
!     tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow,  !
!     ixlow+1, ..., ixupp, Locate returns the array ix is column       !
!     indices and the array dx of residuals.                           !
!                                                                      !
!     author: c. h. goldberg                                           !
!     revised: 1/1/93                                                  !
!     certified:  radiation version 1.0                                !
!                                                                      !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      type (axis_type), intent(in) :: table_axis
      real (kind=kind_phys), dimension(:,:), intent(in) :: x
      integer, intent(in) :: k_min, k_max

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: dx
      integer,dimension(:,:), intent(out) :: ix

!  ---  locals:
      real (kind=kind_phys), dimension(size(x,1), size(x,2)) :: fx
      real (kind=kind_phys) :: table_min, table_inc
      integer :: k, table_col
!
!===> ...  begin here
!
      do k = k_min, k_max
        fx(:,k) = aint( (x(:,k) - table_axis%min_val)                   &
     &          / table_axis%tab_inc )
        dx(:,k) = x(:,k) - fx(:,k)*table_axis%tab_inc                   &
     &          - table_axis%min_val
        ix(:,k) = int( fx(:,k) ) + table_axis%first_col
      enddo
!
      return
!...................................
      end subroutine locate_in_table
!-----------------------------------



!-----------------------------------
      subroutine looktab_type1                                          &
!...................................
!  ---  inputs:
     &     ( tab, ix, iy, dx, dy, k_min, k_max,                         &
!  ---  outputs:
     &       answer                                                     &
     &     )

! -------------------------------------------------------------------- !
!                                                                      !
!     given arrays ix(:,:) and iy(:,:) of integral subscripts and      !
!     arrays dx(:,:) and dy(:,:) of differences from x(:,:,:) and      !
!     y(:,:), calculate answer(:,:) = f(x(:,:), y(:,:))                !
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.       !
!                                                                      !
!     author: c. h. goldberg                                           !
!     revised: 1/1/93                                                  !
!     certified:  radiation version 1.0                                !
!                                                                      !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      type(tab1_type), intent(in) :: tab

      real (kind=kind_phys), dimension(:,:), intent(in) :: dx, dy

      integer, dimension(:,:), intent(in) :: ix, iy
      integer,                 intent(in) :: k_min, k_max

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: answer

!  ---  locals:
      integer :: i_min, i_max, i, k
!
!===> ...  begin here
!
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)

      do k = k_min, k_max
        do i = i_min, i_max
          answer(i,k) = tab%vae(ix(i,k),iy(i,k))                        &
     &                + dx(i,k)*tab%td(ix(i,k),iy(i,k))                 &
     &                + dy(i,k)*tab%md(ix(i,k),iy(i,k))                 &
     &                + dx(i,k)*dy(i,k)*tab%cd(ix(i,k),iy(i,k))
        enddo
      enddo

!
      return
!...................................
      end subroutine looktab_type1
!-----------------------------------



!-----------------------------------
      subroutine looktab_type3                                          &
!...................................
!  ---  inputs:
     &     ( tab, ix, dx,  k_min, k_max, n,                             &
!  ---  outputs:
     &       answer                                                     &
     &     )

! -------------------------------------------------------------------- !
!                                                                      !
!     given arrays ix(:,:) and dx(:,:) of integer subscripts and       !
!     differences from x(:,:) and constant column subscript iyconst,   !
!     calculate answer(:,:) = f(x(:,:), y(:,:)) from four tables       !
!     of values f, df/dx, df/dy, and d2f/dxdy.                         !
!                                                                      !
!     author: c. h. goldberg                                           !
!     revised: 1/1/93                                                  !
!     certified:  radiation version 1.0                                !
!                                                                      !
! -------------------------------------------------------------------- !

      implicit none

!  ---  inputs:
      type(tab3_type), intent(in) :: tab

      real (kind=kind_phys), dimension(:,:), intent(in) :: dx

      integer, dimension(:,:), intent(in) :: ix
      integer,                 intent(in) :: k_min, k_max, n

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: answer

!  ---  locals:
      integer :: i, k, i_min, i_max
!
!===> ...  begin here
!
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)

      do k = k_min, k_max
        do i = i_min, i_max
          answer(i,k) = tab%vae(ix(i,k),n) + dx(i,k)*tab%td(ix(i,k),n)
        enddo
      enddo

!
      return
!...................................
      end subroutine looktab_type3
!-----------------------------------



!-----------------------------------
      subroutine table1_alloc                                           &
!...................................
!  ---  inputs:
     &     ( dim1, dim2,                                                &
!  ---  in/outputs:
     &       tab                                                        &
     &     )
!
      implicit none
!
!  ---  inputs:
      integer, intent(in) :: dim1, dim2

!  ---  in/outputs:
      type(tab1_type), intent(inout) :: tab
!
!===> ...  begin here
!
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))
      allocate (tab%md (dim1, dim2))
      allocate (tab%cd (dim1, dim2))
!
      return
!...................................
      end subroutine table1_alloc
!-----------------------------------



!-----------------------------------
      subroutine table3_alloc                                           &
!...................................
!  ---  inputs:
     &     ( dim1, dim2,                                                &
!  ---  in/outputs:
     &       tab                                                        &
     &     )
!
      implicit none
!
!  ---  inputs:
      integer, intent(in) :: dim1, dim2

!  ---  in/outputs:
      type(tab3_type), intent(inout) :: tab
!
!===> ...  begin here
!
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))

!...................................
      end subroutine table3_alloc
!-----------------------------------


!
!........................................!
      end module module_radlw_main       !
!========================================!

