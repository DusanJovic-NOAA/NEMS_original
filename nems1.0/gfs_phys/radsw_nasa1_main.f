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
!      'module_radsw_cldscal'     -- cloud scaling coefficients table  !
!                                                                      !
!   the 'radsw_nasa1_main.f' contains:                                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main nasa1 sw radiation routine                 !
!         inputs:                                                      !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                     !
!            clouds,iovr,aerosols,sfcalb,                              !
!            cosz,solcon,NPTS,idxday,                                  !
!            IMAX, NLAY, NLP1, iflip, lprnt,                           !
!         outputs:                                                     !
!            hswc,topflx,sfcflx,                                       !
!!        optional outputs:                                            !
!            HSW0,HSWB,FLXPRF,FDNCMP                                   !
!           )                                                          !
!                                                                      !
!      'rswinit'    -- initialization routine                          !
!         inputs:                                                      !
!           ( icwp, me, NLAY )                                         !
!         outputs:                                                     !
!           (none)                                                     !
!                                                                      !
!   all the sw radiation subprograms become contained subprograms      !
!   in module 'module_radsw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!                                                                      !
!                                                                      !
!    derived data type constructs used:                                !
!                                                                      !
!     1. radiation flux at toa: (from module 'module_radsw_parameters')!
!          topfsw_type   -  derived data type for toa rad fluxes       !
!            upfxc              total sky upward flux at toa           !
!            dnfxc              total sky downward flux at toa         !
!            upfx0              clear sky upward flux at toa           !
!                                                                      !
!     2. radiation flux at sfc: (from module 'module_radsw_parameters')!
!          sfcfsw_type   -  derived data type for sfc rad fluxes       !
!            upfxc              total sky upward flux at sfc           !
!            dnfxc              total sky downward flux at sfc         !
!            upfx0              clear sky upward flux at sfc           !
!            dnfx0              clear sky downward flux at sfc         !
!                                                                      !
!     3. radiation flux profiles(from module 'module_radsw_parameters')!
!          profsw_type    -  derived data type for rad vertical prof   !
!            upfxc              total sky level upward flux            !
!            dnfxc              total sky level downward flux          !
!            upfx0              clear sky level upward flux            !
!            dnfx0              clear sky level downward flux          !
!                                                                      !
!     4. surface component fluxes(from module 'module_radsw_parameters'!
!          cmpfsw_type    -  derived data type for component sfc flux  !
!            uvbfc              total sky downward uv-b flux at sfc    !
!            uvbf0              clear sky downward uv-b flux at sfc    !
!            nirbm              surface downward nir direct beam flux  !
!            nirdf              surface downward nir diffused flux     !
!            visbm              surface downward uv+vis direct beam flx!
!            visdf              surface downward uv+vis diffused flux  !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!                                                                      !
!       'module machine'                                               !
!       'module physcons'                                              !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_nasa1_param.f'                                           !
!      'radsw_nasa1_datatb.f'                                          !
!      'radsw_nasa1_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use sw modules     !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!   the original program descriptions:                                 !
!                                                                      !
!   ********************       CLIRAD-SW      *********************    !
!                                                                      !
!  - february 18, 2002                                                 !
!    (1) use new parameterization for ice cloud single-scattering      !
!        albedo by chou, lee, and yang (jgr, 2002).                    !
!    (2) parameterize reff following mcfarquhar (2000) for ice clouds  !
!        and szczodrak et al. (2001) for water clouds.                 !
!                                                                      !
!  - following the nasa technical memorandum (nasa/tm-1999-104606,     !
!    Vol. 15) of chou and suarez (1999), this routine computes solar   !
!    fluxes due to absorption by water vapor, ozone, co2, o2, clouds,  !
!    and aerosols and due to scattering by clouds, aerosols, and gases.!
!                                                                      !
!  - the computer code and documentation are accessible at             !
!        http://climate.gsfc.nasa.gov/~chou/clirad_sw                  !
!                                                                      !
!  - cloud ice, liquid, and rain particles are allowed to co-exist     !
!    in a layer.                                                       !
!                                                                      !
!  - the maximum-random assumption is applied for treating cloud       !
!    overlapping. clouds are grouped into high, middle, and low clouds !
!    separated by the level indices ict and icb.  for detail, see      !
!    subroutine "cldscale". note: ict must be less than icb, and icb   !
!    must be less than np+1.                                           !
!                                                                      !
!  - in a high spatial-resolution atmospheric model, fractional cloud  !
!    cover might be computed to be either 0 or 1.  in such a case,     !
!    scaling of the cloud optical thickness is not necessary, and the  !
!    computation can be made faster by setting overcast=.true.         !
!    otherwise, set the option overcast=.false. (note: for the case    !
!    that fractional cloud cover in a layer is either 0 or 1, the      !
!    results of using either the .true. option or the .false. option   !
!    are identical).                                                   !
!                                                                      !
!  - aerosol optical thickness, single-scattering albedo, and asymmetry!
!    factor can be specified as functions of height and spectral band, !
!    and for various aerosol types. Set aerosol=.true. if aerosols are !
!    included.                                                         !
!                                                                      !
!                                                                      !
!                                                                      !
!   references:                                                        !
!                                                                      !
!       chou, m.-d. and m.j. suarez (1999): a solar radiation          !
!       parameterization for atmospheric studies. nasa/tm-1999-104606, !
!       vol. 15                                                        !
!                                                                      !
!                                                                      !
!                                                                      !
!   ncep modifications history log:                                    !
!                                                                      !
!       jun 2002,  yu-tai hou  -- obtained original code from nasa     !
!       jun 2006,  yu-tai hou                                          !
!                  modified code into standard modular f90 code for    !
!                  ncep model radiation package.                       !
!       apr 2007,  yu-tai hou                                          !
!                  add spectral band heating as optional output        !
!                                                                      !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radsw_main           !
!........................................!
!
      use machine,            only : kind_phys
!     use physcons,           only : con_rd, con_ttp, con_g,            &
!    &                               con_cp, con_amo3, con_amd,         &
!    &                               con_gasv, con_p0
      use module_radsw_parameters
      use module_radsw_cntr_para
!
      implicit   none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGSW='NASA-SW vfeb02  jun 2006'
      character(24), parameter :: VTAGSW='NASA-SW vfeb02  Apr 2007'

!  ---  constant values

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: _zero = 0.0
      real (kind=kind_phys), parameter :: _one  = 1.0
      real (kind=kind_phys), parameter :: ftiny = 1.0e-11
      real (kind=kind_phys), parameter :: fpmin = 1.0e-8
      real (kind=kind_phys), parameter :: fpmax = 0.999999

!! ...  logical flags for optional output fields

      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.

!  ---  those data will be set up only once by "rswinit"

!  ...  heatfac is the factor for heating rates
!       (in k/day, or k/sec set by subroutine 'rlwinit')

      real (kind=kind_phys) :: heatfac

      public  swrad, rswinit


! =================
      contains
! =================


!-----------------------------------
      subroutine swrad                                                  &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,iovr,aerosols,sfcalb,                               &
     &       cosz,solcon,NPTS,idxday,                                   &
     &       IMAX, NLAY, NLP1, iflip, lprnt,                            &
!  ---  outputs:
     &       hswc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HSW0,HSWB,FLXPRF,FDNCMP                                    &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    swrad      computes short-wave radiative heating       !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   plyr (IMAX,NLAY) : model layer mean pressure in mb (not in use)     !
!   plvl (IMAX,NLP1) : model level pressure in mb                       !
!   tlyr (IMAX,NLAY) : model layer mean temperature in k                !
!   tlvl (IMAX,NLP1) : model level temperature in k    (not in use)     !
!   qlyr (IMAX,NLAY) : layer h2o mass mixing ratio (gm/gm) *see inside  !
!   olyr (IMAX,NLAY) : layer ozone mass concentration (gm/gm)           !
!   gasvmr(IMAX,NLAY,:): atmospheric constent gases:                    !
!                      (check module_radiation_gases for definition)    !
!     gasvmr(:,:,1)    -   co2 volume mixing ratio                      !
!     gasvmr(:,:,2)    -   n2o volume mixing ratio     (not used)       !
!     gasvmr(:,:,3)    -   ch4 volume mixing ratio     (not used)       !
!     gasvmr(:,:,4)    -   o2  volume mixing ratio                      !
!     gasvmr(:,:,5)    -   co  volume mixing ratio     (not used)       !
!     gasvmr(:,:,6)    -   cf11 volume mixing ratio    (not used)       !
!     gasvmr(:,:,7)    -   cf12 volume mixing ratio    (not used)       !
!     gasvmr(:,:,8)    -   cf22 volume mixing ratio    (not used)       !
!     gasvmr(:,:,9)    -   ccl4 volume mixing ratio    (not used)       !
!     gasvmr(:,:,10)   -   cf113 volume mixing ratio   (not used)       !
!   clouds(IMAX,NLAY,:): cloud profile                                  !
!                      (check module_raddiation_clouds for definition)  !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud liq water path      (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer cloud ice water path      (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!   iovr             : control flag for cloud overlapping (approxi only)!
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   aerosols(IMAX,NLAY,NBDSW,:) : aerosol optical properties            !
!                      (check module_radiation_aerosols for definition) !
!         (:,:,:,1)   - optical depth                                   !
!         (:,:,:,2)   - single scattering albedo                        !
!         (:,:,:,3)   - asymmetry parameter                             !
!   sfcalb(IMAX, : ) : surface albedo in fraction                       !
!                      (check module_radiation_surface for definition)  !
!         ( :, 1 )    - near ir direct beam albedo                      !
!         ( :, 2 )    - near ir diffused albedo                         !
!         ( :, 3 )    - uv+vis direct beam albedo                       !
!         ( :, 4 )    - uv+vis diffused albedo                          !
!   cosz  (IMAX)     : cosine of solar zenith angle                     !
!   solcon           : solar constant                      (w/m**2)     !
!   NPTS             : num of daytime points                            !
!   idxday(IMAX)     : index array for daytime points                   !
!   IMAX             : number of horizontal points                      !
!   NLAY,NLP1        : vertical layer/lavel numbers                     !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lprnt            : logical check print flag                         !
!                                                                       !
! control variables set in module 'module_radsw_cntr_para':             !
!   iswrate: heating rate unit selections                               !
!            =1: output in k/day                                        !
!            =2: output in k/second                                     !
!   iaersw : flags for aerosols effect                                  !
!            =0: without aerosol effect                                 !
!            >0: include aerosol effect                                 !
!   ioxysw : flags for o2 absorption                                    !
!            =0: without o2 absorption                                  !
!            =1: include o2 absorption                                  !
!   ico2sw : flags for co2 absorption                                   !
!            =0: without co2 absorption                                 !
!            =1: include co2 absorption                                 !
!   iflagc : control flag for cloud scaling method                      !
!            =1; input cloud condensates (pathes and eff radius)        !
!                no scaling, cloud fraction either 0 or 1               !
!            =2; input cloud condensates (pathes and eff radius)        !
!                use fractional cloud amount with scaling               !
!                                                                       !
! output parameter:                                                     !
!   hswc  (IMAX,NLAY): total sky heating rates (k/sec or k/day)         !
!   topflx(IMAX)     : radiation fluxes at toa (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at toa                   !
!     dnfxc            - total sky downward flux at toa                 !
!     upfx0            - clear sky upward flux at toa                   !
!   sfcflx(IMAX)     : radiation fluxes at sfc (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at sfc                   !
!     dnfxc            - total sky downward flux at sfc                 !
!     upfx0            - clear sky upward flux at sfc                   !
!     dnfx0            - clear sky downward flux at sfc                 !
!                                                                       !
!!optional outputs:                                                     !
!   hswb(IMAX,NLAY,NBDSW): spectral band total sky heating rates        !
!   hsw0  (IMAX,NLAY): clear sky heating rates (k/sec or k/day)         !
!   flxprf(IMAX,NLP1): level radiation fluxes (w/m**2), components:     !
!                      (check module_radsw_parameters for definition)   !
!     dnfxc            - total sky downward flux at interface           !
!     upfxc            - total sky upward flux at interface             !
!     dnfx0            - clear sky downward flux at interface           !
!     upfx0            - clear sky upward flux at interface             !
!   fdncmp(IMAX)     : component surface downward fluxes (w/m**2):      !
!                      (check module_radsw_parameters for definition)   !
!     uvbfc            - total sky downward uv-b flux at sfc            !
!     uvbf0            - clear sky downward uv-b flux at sfc            !
!     nirbm            - downward surface nir direct beam flux          !
!     nirdf            - downward surface nir diffused flux             !
!     visbm            - downward surface vis direct beam flux          !
!     visdf            - downward surface vis diffused flux             !
!                                                                       !
! note:                                                                 !
!   for all internal quantities, k=1 is the top level/layer             !
!                                                                       !
! internal variables:                                                   !
!                                                   units      size     !
!  level pressure (pl)                              mb      NPTS*NLP1   !
!  layer temperature (ta)                           k       NPTS*NLAY   !
!  layer specific humidity (wa)                     gm/gm   NPTS*NLAY   !
!  layer ozone concentration (oa)                   gm/gm   NPTS*NLAY   !
!  co2 mixing ratio by volume (co2)                 pppv    NPTS*NLAY   !
!  cloud water mixing ratio (cwc)                  gm/gm    NPTS*NLAY*3 !
!        index 1 for ice particles                                      !
!        index 2 for liquid drops                                       !
!        index 3 for rain drops                                         !
!  cloud amount (fcld)                            fraction  NPTS*NLAY   !
!  level index separating high and middle           n/d         1       !
!        clouds (ict)                                                   !
!  level index separating middle and low            n/d         1       !
!        clouds (icb)                                                   !
!  aerosol optical thickness (taual)                n/d     NPTS*NLAY*11!
!  aerosol single-scattering albedo (ssaal)         n/d     NPTS*NLAY*11!
!  aerosol asymmetry factor (asyal)                 n/d     NPTS*NLAY*11!
!        in the uv region :                                             !
!           index  1 for the 0.175-0.225 micron band                    !
!           index  2 for the 0.225-0.245;0.26-0.28 micron band          !
!           index  3 for the 0.245-0.260 micron band                    !
!           index  4 for the 0.280-0.295 micron band                    !
!           index  5 for the 0.295-0.310 micron band                    !
!           index  6 for the 0.310-0.320 micron band                    !
!           index  7 for the 0.325-0.400 micron band                    !
!        in the par region :                                            !
!           index  8 for the 0.400-0.700 micron band                    !
!        in the infrared region :                                       !
!           index  9 for the 0.700-1.220 micron band                    !
!           index 10 for the 1.220-2.270 micron band                    !
!           index 11 for the 2.270-10.00 micron band                    !
!  reflection of ground/ocean                     fraction    NPTS      !
!   in uv+parc (bands 1-8)                                              !
!        (rgbuv) for beam insolation                                    !
!        (rgfuv) for diffuse insolation                                 !
!  reflection of ground/ocean                     fraction    NPTS      !
!   in ir (bands 9-11)                                                  !
!        (rgbir) for beam insolation                                    !
!        (rgfir) for diffuse insolation                                 !
!                                                                       !
!   all-sky flux (downward minus upward) (flx)     fraction NPTS*NLP1   !
!   clear-sky flux (downward minus upward) (flc)   fraction NPTS*NLP1   !
!   all-sky direct downward uv (0.175-0.4 micron)                       !
!                flux at the surface (fdiruv)      fraction   NPTS      !
!   all-sky diffuse downward uv flux at                                 !
!                the surface (fdifuv)              fraction   NPTS      !
!   all-sky direct downward par (0.4-0.7 micron)                        !
!                flux at the surface (fdirpar)     fraction   NPTS      !
!   all-sky diffuse downward par flux at                                !
!                the surface (fdifpar)             fraction   NPTS      !
!   all-sky direct downward ir (0.7-10 micron)                          !
!                flux at the surface (fdirir)      fraction   NPTS      !
!   all-sky diffuse downward ir flux at                                 !
!                the surface (fdifir)              fraction   NPTS      !
!                                                                       !
!  Notes:                                                               !
!                                                                       !
!  (1) The unit of output fluxes (flx,flc,etc.) is fraction of the total!
!      insolation at the top of the atmosphere.  Therefore, fluxes      !
!      are the output fluxes multiplied by the extra-terrestrial solar  !
!      flux and the cosine of the solar zenith angle.                   !
!  (2) pl(*,1) is the pressure at the top of the model, and             !
!      pl(*,NLP1) is the surface pressure.                              !
!  (3) the pressure levels ict and icb correspond approximately         !
!      to 400 and 700 mb.                                               !
!                                                                       !
!  -- Please notify Ming-Dah Chou at  chou@climate.gsfc.nasa.gov        !
!     for coding errors.                                                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use module_radsw_co2tab
 
      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1, iovr, iflip, NPTS
      integer, intent(in) :: idxday(:)
      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, tlvl,  &
     &       plyr, tlyr, qlyr, olyr, sfcalb

      real (kind=kind_phys), dimension(:,:,:),   intent(in) :: gasvmr,  &
     &       clouds

      real (kind=kind_phys), dimension(:,:,:,:), intent(in) :: aerosols

      real (kind=kind_phys), intent(in) :: cosz(:), solcon

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: hswc
      type (topfsw_type),    dimension(:),   intent(out) :: topflx
      type (sfcfsw_type),    dimension(:),   intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys),dimension(:,:,:),optional,intent(out):: hswb
      real (kind=kind_phys),dimension(:,:),  optional,intent(out):: hsw0
      type (profsw_type), dimension(:,:),optional, intent(out) :: flxprf
      type (cmpfsw_type), dimension(:),  optional, intent(out) :: fdncmp

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLP1,NBDSW) :: fnetb
      real (kind=kind_phys), dimension(NPTS,NLP1) :: pl, flx, flc,      &
     &       flxu, flxd, flcu, flcd, swu, swh, so2, df, dfo2, dfco2a,   &
     &       dfco2b
      real (kind=kind_phys), dimension(NPTS,NLAY) :: pa, ta, wa, oa,    &
     &       co2, o2, fcld, dp, wh, oh, scal
      real (kind=kind_phys), dimension(NPTS)      :: rgbuv, rgfuv,      &
     &       rgbir, rgfir, snt, cnt, cosza, s0cosz
      real (kind=kind_phys), dimension(NLP1)      :: flxk
      real (kind=kind_phys) :: xx, fntop, xo2, xco2a, xco2b

      real (kind=kind_phys), dimension(NPTS,NLAY,3)     :: cwp, reff
      real (kind=kind_phys), dimension(NPTS,NLAY,NBDSW) :: taual,       &
     &       ssaal, asyal

!! ---  for optional surface fluxes components:
      real (kind=kind_phys), dimension(NPTS)      :: fdirpar, fdifpar,  &
     &       fdiruv, fdifuv, fdirir, fdifir, flxuvb, flcuvb

      integer :: i, ib, j, j1, k, k1, in, ntop, mb
      integer, dimension(NPTS) :: nctop, ict, icb

!  ---  parameters for co2 transmission tables
      real (kind=kind_phys) :: uu1,du1,ww1,dw1, uu2,du2,ww2,dw2

      parameter (uu1=-3.0,   ww1=-4.0,  du1=0.15,   dw1=0.15)
      parameter (uu2=.25e-3, ww2=-2.0,  du2=.05e-3, dw2=0.05)

!
!===> ...  begin here
!

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )

!  --- ...  initial output arrays

      hswc(:,:) = _zero
      topflx = topfsw_type ( _zero, _zero, _zero )
      sfcflx = sfcfsw_type ( _zero, _zero, _zero, _zero )

!! ---  optional outputs
      if ( lflxprf ) then
        flxprf = profsw_type ( _zero, _zero, _zero, _zero )
      endif

      if ( lfdncmp ) then
        fdncmp = cmpfsw_type ( _zero, _zero, _zero, _zero, _zero, _zero )
      endif

      if ( lhsw0 ) then
        hsw0(:,:) = _zero
      endif

      if ( lhswb ) then
        hswb(:,:,:) = _zero
      endif

!  --- ...  work on daytime points only
!           note: the internal array is always from top to surface

      if ( iflip == 0 ) then      ! input from toa to sfc

        do i = 1, NPTS
          j1 = idxday(i)
          pl(i,NLP1) = plvl(j1,NLP1)

          do k = 1, NLAY
            pl(i,k)  = plvl(j1,k)
            ta(i,k)  = tlyr(j1,k)
!test use
!           wa(i,k)  = max(_zero, qlyr(j1,k))                   ! input mass mixing ratio
!ncep model use
            wa(i,k)  = max(_zero, qlyr(j1,k)/(_one-qlyr(j1,k))) ! input specific humidity
            oa(i,k)  = max(_zero, olyr(j1,k))                   ! input mass mixing ratio
            co2(i,k) = max(_zero, gasvmr(j1,k,1))               ! co2
!o2vmr      o2 (i,k) = max(_zero, gasvmr(j1,k,4))               ! o2

            fcld(i,k) = clouds(j1,k,1)
            cwp(i,k,1)= clouds(j1,k,4)
            cwp(i,k,2)= clouds(j1,k,2)
            cwp(i,k,3)= clouds(j1,k,6)
            reff(i,k,1) = max( _one, min( 150., clouds(j1,k,5) ) )
            reff(i,k,2) = max( 4.0,  min( 20.0, clouds(j1,k,3) ) )
            reff(i,k,3) = clouds(j1,k,7)
          enddo
        enddo

        if (iaersw > 0) then
          do i = 1, NPTS
            j1 = idxday(i)

            do ib = 1, NBDSW
              do k = 1, NLAY
                taual(i,k,ib) = aerosols(j1,k,ib,1)
                ssaal(i,k,ib) = aerosols(j1,k,ib,2)
                asyal(i,k,ib) = aerosols(j1,k,ib,3)
              enddo
            enddo
          enddo
        else
          taual(:,:,:) = _zero
          ssaal(:,:,:) = _zero
          asyal(:,:,:) = _zero
        endif

      else                        ! input data from sfc to toa

        do i = 1, NPTS
          j1 = idxday(i)
          pl(i,1)  = plvl(j1,NLP1)

          do k = 1, NLAY
            k1 = NLP1 - k
            pl(i,k+1)= plvl(j1,k1)
            ta(i,k)  = tlyr(j1,k1)
!test use
!           wa(i,k)  = max(_zero, qlyr(j1,k1))                   ! input mass mixing ratio
!ncep model use
            wa(i,k)  = max(_zero, qlyr(j1,k1)/(_one-qlyr(j1,k1)))! input specific humidity
            oa(i,k)  = max(_zero, olyr(j1,k1))                   ! input mass mixing ratio
            co2(i,k) = max(_zero, gasvmr(j1,k1,1))               ! co2
!o2vmr      o2 (i,k) = max(_zero, gasvmr(j1,k1,4))               ! o2

            fcld(i,k) = clouds(j1,k1,1)
            cwp(i,k,1)= clouds(j1,k1,4)
            cwp(i,k,2)= clouds(j1,k1,2)
            cwp(i,k,3)= clouds(j1,k1,6)
            reff(i,k,1) = max( _one, min( 150., clouds(j1,k1,5) ) )
            reff(i,k,2) = max( 4.0,  min( 20.0, clouds(j1,k1,3) ) )
            reff(i,k,3) = clouds(j1,k1,7)
          enddo
        enddo

        if ( iaersw > 0 ) then
          do i = 1, NPTS
            j1 = idxday(i)

            do ib = 1, NBDSW
              do k = 1, NLAY
                k1 = NLP1 - k

                taual(i,k,ib) = aerosols(j1,k1,ib,1)
                ssaal(i,k,ib) = aerosols(j1,k1,ib,2)
                asyal(i,k,ib) = aerosols(j1,k1,ib,3)
              enddo
            enddo
          enddo
        else
          taual(:,:,:) = _zero
          ssaal(:,:,:) = _zero
          asyal(:,:,:) = _zero
        endif

      endif                       ! if_iflip

      do i = 1, NPTS
        j1 = idxday(i)
        rgbir(i) = sfcalb(j1,1)
        rgfir(i) = sfcalb(j1,2)
        rgbuv(i) = sfcalb(j1,3)
        rgfuv(i) = sfcalb(j1,4)
        cosza(i) = cosz(j1)
      enddo

!  --- ...  determine boundaries of cloud domains

      ict = 2
      icb = NLAY
      do k = NLAY, 2, -1
        do i = 1, NPTS
          if ( pl(i,k) >= 400.0 ) then
            ict(i) = k - 1
          endif
          if ( pl(i,k) >= 700.0 ) then
            icb(i) = k
          endif
        enddo
      enddo

      do i = 1, NPTS
        if ( ict(i) >= icb(i) ) then
          ict(i) = icb(i) - 1
        endif
      enddo

!  --- ...  initialization and set the secant of the solar zenith angle
      do i = 1, NPTS
        swh(i,1) = _zero
        so2(i,1) = _zero
        snt  (i) = _one / cosza(i)
        s0cosz(i)= solcon * cosza(i)
      enddo

!  ---  ozone unit conversion factor
!o3   xx = 1.0e4*con_rd*(con_amd/con_amo3)*con_ttp/con_p0/con_g

      do k = 1, NLAY
        do i = 1, NPTS

!  ---  compute layer thickness. indices for the surface level and
!       surface layer are NLP1 and np, respectively.
          dp(i,k) = pl(i,k+1) - pl(i,k)

!  ---  compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!       unit is g/cm**2
          pa(i,k) = 0.5 * (pl(i,k) + pl(i,k+1))

          scal(i,k) = dp(i,k) * (pa(i,k)/300.)**0.8
          wh(i,k) = 1.02 * wa(i,k) * scal(i,k)                          &
     &            * (_one + 0.00135*(ta(i,k)-240.)) + ftiny
          swh(i,k+1) = swh(i,k) + wh(i,k)

!  ---  compute ozone amount, unit is (cm-atm)stp the number 466.7 is
!       the unit conversion factor from g/cm**2 to (cm-atm)stp
          oh(i,k) = 1.02*oa(i,k)*dp(i,k)*466.7 + ftiny
!o3       oh(i,k) = xx  *oa(i,k)*dp(i,k)       + ftiny
        enddo
      enddo

!  ---  initialize fluxes for all-sky (flx), clear-sky (flc), and
!       flux reduction (df)

      do k = 1, NLP1
        do i = 1, NPTS
          flxu(i,k) = _zero
          flxd(i,k) = _zero
          flcu(i,k) = _zero
          flcd(i,k) = _zero
        enddo
      enddo

!  ---  compute solar uv and par fluxes

      call soluv                                                        &
!  ---  inputs:
     &     ( cosza,wh,oh,dp,fcld,cwp,reff,ict,icb,                      &
     &       taual,ssaal,asyal,rgbuv,rgfuv,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       flxu,flxd,flcu,flcd                                        &
!! ---  optional outputs:
     &,      fdiruv,fdifuv,fdirpar,fdifpar,flxuvb,flcuvb, fnetb         &
     &     )

!  ---  compute and update solar ir fluxes

      call solir                                                        &
!  ---  inputs:
     &     ( cosza,wh,dp,fcld,cwp,reff,ict,icb,                         &
     &       taual,ssaal,asyal,rgbir,rgfir,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       flxu,flxd,flcu,flcd                                        &
!! ---  optional outputs:
     &,      fdirir,fdifir, fnetb                                       &
     &     )

      dfo2  (:,:) = _zero
      dfco2a(:,:) = _zero
      dfco2b(:,:) = _zero

      if ( ioxysw == 1 ) then

!  ---  compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!       unit is (cm-atm)stp.
!       the constant 165.22 equals (1000/980)*23.14%*(22400/32)
!       for o2 in vmr the constant becomes (1000/980)*o2vmr*(22400/28.97)

!o2vmr  xx = (10.0 / con_g) * (con_gasv*1.0e6 / con_amd)
        do i = 1, NPTS
!o2vmr    cnt(i) = xx    *snt(i)
          cnt(i) = 165.22*snt(i)
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
!o2vmr      so2(i,k+1) = so2(i,k) + scal(i,k)*cnt(i)*o2(i,k)
            so2(i,k+1) = so2(i,k) + scal(i,k)*cnt(i)
          enddo
        enddo

!  ---  compute flux reduction due to oxygen following Eq. (3.18)
!       the constant 0.0633 is the fraction of insolation contained 
!       in the oxygen bands

        do k = 2, NLP1
          do i = 1, NPTS
            dfo2(i,k) = 0.0633*(_one - exp(-0.000145*sqrt(so2(i,k))))
          enddo
        enddo          

      endif    ! end if_ioxysw_block

      if ( ico2sw == 1 ) then

!  ---  for solar heating due to co2
!       the constant 789 equals (1000/980)*(44/28.97)*(22400/44)

!co2    xx = (10.0 / con_g) * (con_gasv*1.0e6 / con_amd)
        do i = 1, NPTS
!co2      cnt(i) = xx    * snt(i)
          cnt(i) = 789.0 * snt(i)
        enddo

!  ---  scale co2 amounts following Eq. (3.5) with f=1. unit is (cm-atm)stp.

        do k = 1, NLAY
          do i = 1, NPTS
            so2(i,k+1) = so2(i,k) + cnt(i)*co2(i,k)*scal(i,k) + ftiny
          enddo
        enddo

!  ---  for co2 absorption in Band 10 where absorption due to water vapor
!       and co2 are both moderate
!       so2 and swh are the co2 and water vapor amounts integrated
!       from the top of the atmosphere

        do k = 2, NLP1
          do i = 1, NPTS
            swu(i,k) = log10( so2(i,k) )
            swh(i,k) = log10( swh(i,k)*snt(i) )
          enddo
        enddo

!  ---  df is the updated flux reduction given by the second term on the
!       right-hand-side of Eq. (3.24) divided by So

        call rflx                                                       &
!  ---  inputs:
     &     ( swu,uu1,du1,NUCAH,swh,ww1,dw1,NWCAH,cah,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       dfco2a                                                     &
     &     )

!  ---  for co2 absorption in Band 11 where the co2 absorption has
!       a large impact on the heating of middle atmosphere.
!       co2 mixing ratio is independent of space

        do k = 2, NLP1
          do i = 1, NPTS
            swu(i,k) = co2(i,k-1) * snt(i)
          enddo
        enddo

!  ---  swh is the logarithm of pressure

        do k = 2, NLP1
          do i = 1, NPTS
            swh(i,k) = log10( pl(i,k) )
          enddo
        enddo

!  ---  df is the updated flux reduction derived from the table given by
!       Eq. (3.19)

        call rflx                                                       &
!  ---  inputs:
     &     ( swu,uu2,du2,NXCOA,swh,ww2,dw2,NYCOA,coa,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       dfco2b                                                     &
     &     )

      endif   ! end if_ico2sw_block

      if ( ioxysw ==1 .or. ico2sw == 1 ) then

        do k = 1, NLP1
          do i = 1, NPTS
            df(i,k) = dfo2(i,k) + dfco2a(i,k) + dfco2b(i,k)
            dfo2(i,k) = dfo2(i,k) / 3.0
          enddo
        enddo

!  ---  identify top cloud-layer

        do i = 1, NPTS
          nctop(i) = NLP1
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
            if ( fcld(i,k) > 0.02 .and. nctop(i) == NLP1 ) then
              nctop(i) = k
            endif
          enddo
        enddo

        do i = 1, NPTS
          ntop  = nctop(i)
          fntop = flxd(i,ntop) - flxu(i,ntop)

!  ---  adjust fluxes below cloud top following Eq. (6.18)
          if ( ntop < NLP1 ) then
            do k = ntop+1, NLP1
              flxk(k) = (flxd(i,k) - flxu(i,k)) / fntop
              xx = df(i,ntop) + (df(i,k)-df(i,ntop))*flxk(k)
              flxd(i,k) = flxd(i,k) - xx
            enddo

!  ---  adjustment for the direct downward flux
            if ( lfdncmp ) then
              fdirir(i) = max(_zero, fdirir(i)-xx)
            endif

!  ---  for spectral band heating
            if ( lhswb ) then
              do k = ntop+1, NLP1
                xo2  = dfo2(i,ntop)                                     &
     &               + (dfo2(i,k) - dfo2(i,ntop)) * flxk(k)
                xco2a= dfco2a(i,ntop)                                   &
     &               + (dfco2a(i,k) - dfco2a(i,ntop)) * flxk(k)
                xco2b= dfco2b(i,ntop)                                   &
     &               + (dfco2b(i,k) - dfco2b(i,ntop)) * flxk(k)

                fnetb(i,k,8) = fnetb(i,k,8)  - xo2
                fnetb(i,k,9) = fnetb(i,k,9)  - xo2
                fnetb(i,k,10)= fnetb(i,k,10) - xo2 - xco2a
                fnetb(i,k,11)= fnetb(i,k,11) - xco2b
              enddo
            endif
          endif    ! end if_ntop

!  ---  adjust fluxes above clouds following Eq. (6.17)
          do k = 1, ntop
            flxd(i,k) = flxd(i,k) - df(i,k)
          enddo

!  ---  for spectral band heating
          if ( lhswb ) then
            do k = 1, ntop
              fnetb(i,k,8) = fnetb(i,k,8)  - dfo2(i,k)
              fnetb(i,k,9) = fnetb(i,k,9)  - dfo2(i,k)
              fnetb(i,k,10)= fnetb(i,k,10) - dfo2(i,k) - dfco2a(i,k)
              fnetb(i,k,11)= fnetb(i,k,11) - dfco2b(i,k)
            enddo
          endif
        enddo


!  ---  adjustment for the effect of o2 cnd co2 on clear-sky fluxes.
!       both flc and df are positive quantities

        do k = 1, NLP1
          do i = 1, NPTS
            flcd(i,k) = flcd(i,k) - df(i,k)
          enddo
        enddo
      endif    ! end if_ioxysw_or_ico2sw_block

!  ---  convert the units of flx and flc from fraction to w/m^2

      do k = 1, NLP1
        do i = 1, NPTS
          flxu(i,k) = flxu(i,k)*s0cosz(i)
          flxd(i,k) = flxd(i,k)*s0cosz(i)
          flcu(i,k) = flcu(i,k)*s0cosz(i)
          flcd(i,k) = flcd(i,k)*s0cosz(i)

          flx (i,k) = flxd(i,k) - flxu(i,k)
        enddo
      enddo

!! ---  optional clear sky heating
      if ( lhsw0 ) then
        flc(:,:) = flcd(:,:) - flcu(:,:)
      endif

!! ---  optional spectral heating
      if ( lhswb ) then
        do mb = 1, NBDSW
          do k = 1, NLP1
            do i = 1, NPTS
              fnetb(i,k,mb) = fnetb(i,k,mb)*s0cosz(i)
            enddo
          enddo
        enddo
      endif

!  ---  compute heating rates

      if ( iflip == 0 ) then      ! output data from toa to sfc

        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            hswc(j1,k) = heatfac * (flx(i,k)-flx(i,k+1))                &
     &                           / (pl(i,k+1)-pl(i,k))
          enddo
        enddo

!! ---  optional clear sky heating
        if ( lhsw0 ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              hsw0(j1,k) = heatfac * (flc(i,k)-flc(i,k+1))              &
     &                             / (pl(i,k+1)-pl(i,k))
            enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhswb ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do mb = 1, NBDSW
            do k = 1, NLAY
              hswb(j1,k,mb) = heatfac * (fnetb(i,k,mb)-fnetb(i,k+1,mb)) &
     &                                / (pl(i,k+1)-pl(i,k))
            enddo
            enddo
          enddo
        endif

!! ---  optional flux profiles
        if ( lflxprf ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              flxprf(j1,k)%upfx0 = flcu(i,k)
              flxprf(j1,k)%dnfx0 = flcd(i,k)
              flxprf(j1,k)%upfxc = flxu(i,k)
              flxprf(j1,k)%dnfxc = flxd(i,k)
            enddo
          enddo
        endif

      else                        ! output data from sfc to toa

        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            k1 = NLP1 - k
            hswc(j1,k1) = heatfac * (flx(i,k)-flx(i,k+1))               &
     &                            / (pl(i,k+1)-pl(i,k))
          enddo
        enddo

!! ---  optional clear sky heating
        if ( lhsw0 ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              k1 = NLP1 - k
              hsw0(j1,k1) = heatfac * (flc(i,k)-flc(i,k+1))             &
     &                              / (pl(i,k+1)-pl(i,k))
            enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhswb ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do mb = 1, NBDSW
            do k = 1, NLAY
              k1 = NLP1 - k
              hswb(j1,k1,mb) = heatfac*(fnetb(i,k,mb)-fnetb(i,k+1,mb))  &
     &                                / (pl(i,k+1)-pl(i,k))
            enddo
            enddo
          enddo
        endif

!! ---  optional flux profiles
        if ( lflxprf ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              k1 = NLP1 - k + 1
              flxprf(j1,k1)%upfx0 = flcu(i,k)
              flxprf(j1,k1)%dnfx0 = flcd(i,k)
              flxprf(j1,k1)%upfxc = flxu(i,k)
              flxprf(j1,k1)%dnfxc = flxd(i,k)
            enddo
          enddo
        endif

      endif                       ! if_iflip

!  ---  toa and sfc fluxes

      do i = 1, NPTS
        j1 = idxday(i)
        topflx(j1)%upfxc = flxu(i,1)
        topflx(j1)%dnfxc = flxd(i,1)
        topflx(j1)%upfx0 = flcu(i,1)
        sfcflx(j1)%upfxc = flxu(i,NLP1)
        sfcflx(j1)%dnfxc = flxd(i,NLP1)
        sfcflx(j1)%upfx0 = flcu(i,NLP1)
        sfcflx(j1)%dnfx0 = flcd(i,NLP1)
      enddo

      if ( lfdncmp ) then
        do i = 1, NPTS
          j1 = idxday(i)
!! ---  optional uv-b surface downward flux
          fdncmp(j1)%uvbfc = flxuvb(i) * s0cosz(i)
          fdncmp(j1)%uvbf0 = flcuvb(i) * s0cosz(i)

!! ---  optional beam and diffused sfc fluxes
          fdncmp(j1)%nirbm = fdirir(i) * s0cosz(i)
          fdncmp(j1)%nirdf = fdifir(i) * s0cosz(i)
          fdncmp(j1)%visbm = (fdirpar(i) + fdiruv(i)) * s0cosz(i)
          fdncmp(j1)%visdf = (fdifpar(i) + fdifuv(i)) * s0cosz(i)
!!        fdncmp(j1)%uvbm  = fdiruv(i) * s0cosz(i)
!!        fdncmp(j1)%uvdf  = fdifuv(i) * s0cosz(i)
!!        fdncmp(j1)%visbm = fdirpar(i) * s0cosz(i)
!!        fdncmp(j1)%visdf = fdifpar(i) * s0cosz(i)
        enddo
      endif
!
      return
!...................................
      end subroutine swrad
!-----------------------------------


!-----------------------------------
      subroutine rswinit                                                &
!...................................

!  ---  inputs:
     &     ( icwp, me, NLAY )
!  ---  outputs: (none)

!  ==================================================================  !
!                                                                      !
!  rswinit is an initialization program for shortwave radiation        !
!  it sets up band lower and upper ranges according to the control     !
!  flags defined in the module "module_radsw_cntr_para".  it also      !
!  check for cloud schemes and the unit of output sw heating rates     !
!                                                                      !
!                                                                      !
!  inputs:                                                             !
!    icwp     -  flag of cloud schemes used by model                   !
!                =0: diagnostic scheme gives cloud tau, omiga, and g   !
!                =1: prognostic scheme gives cloud liq/ice path, etc.  !
!    me       - print control for parallel process                     !
!    NLAY     - number of vertical layers                              !
!                                                                      !
!  outputs: (none)                                                     !
!                                                                      !
!  control flags in module "module_radsw_cntr_para":                   !
!     iswrate - heating rate unit selections                           !
!               =1: output in k/day                                    !
!               =2: output in k/second                                 !
!     iaersw  - flags for aerosols effect                              !
!               =0: without aerosol effect                             !
!               >0: include aerosol effect                             !
!     ioxysw  - flags for o2 absorption                                !
!               =0: without o2 absorption                              !
!               =1: include o2 absorption                              !
!     ico2sw  - flags for co2 absorption                               !
!               =0: without co2 absorption                             !
!               =1: include co2 absorption                             !
!     iflagc  - flag for cloud scaling method                          !
!               =1: input cloud condensates (pathes and eff radius)    !
!                   no scaling, clouds either 0 or 1                   !
!               =2: input cloud condensates (pathes and eff radius)    !
!                   use fractional clouds with scaling                 !
!                                                                      !
!  ==================================================================  !
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: (none)

!  ---  locals:
      integer :: ibnd
!
!===>  ... begin here
!

      if ( me == 0 ) then
        print *,' - Using NASA Shortwave Radiation, Version: ',VTAGSW

        if ( iaersw == 0 ) then
          print *,'   --- Aerosol effect is NOT included in SW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        else
          print *,'   --- Using input aerosol parameters for SW'
        endif

        if ( ioxysw == 0 ) then
          print *,'   --- O2 absorption is NOT included in SW'
        else
          print *,'   --- Include O2 absorption in SW'
        endif

        if ( ico2sw == 0 ) then
          print *,'   --- CO2 absorption is NOT included in SW'
        else
          print *,'   --- Include CO2 absorption in SW'
        endif

        if     ( iflagc == 1 ) then
          print *,'   --- no cloud scaling, use overcast assumption'
        elseif ( iflagc == 2 ) then
          print *,'   --- include cloud scaling for fractional clouds'
          print *,'       clouds are max-ovlp in-domain and ran-ovlp',  &
     &            ' between 3 domains set at 400 and 700 mb'
        else
          print *,'   *** ERROR! check iflagc setting ***'
          stop
        endif

        if ( icwp == 0 ) then
          print *,'   *** ERROR! This version of sw code does not',     &
     &            ' support diagnostic cloud scheme ***'
          stop
        endif
      endif

!  --- ... fheat is the factor for heating rates
!          the 1.0e-2 is to convert pressure from mb to N/m**2

      if ( iswrate == 1 ) then
!       heatfac = con_g * 86400. * 1.e-2 / con_cp    !    (in k/day)
!       heatfac = con_g * 864.0 / con_cp             !    (in k/day)
        heatfac = 8.4410                             !    (in k/day)
      else
!       heatfac = con_g * 1.0e-2 / con_cp            !    (in k/second)
        heatfac = 8.4410 / 86400.0                   !    (in k/second)
      endif

!
      return
!...................................
      end subroutine rswinit
!-----------------------------------


!-----------------------------------
      subroutine soluv                                                  &
!...................................

!  ---  inputs:
     &     ( cosz,wh,oh,dp,fcld,cwp,reff,ict,icb,                       &
     &       taual,ssaal,asyal,rgbuv,rgfuv,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       flxu,flxd,flcu,flcd                                        &
!! ---  optional outputs:
     &,      fdiruv,fdifuv,fdirpar,fdifpar,flxuvb,flcuvb, fnetb         &
     &     )

!******************************************************************
!  compute solar fluxes in the uv+par region. the spectrum is
!  grouped into 8 bands:
!  
!              Band     Micrometer
!
!       UV-C    1.     .175 - .225
!               2.     .225 - .245
!                      .260 - .280
!               3.     .245 - .260
!
!       UV-B    4.     .280 - .295
!               5.     .295 - .310
!               6.     .310 - .320
!      
!       UV-A    7.     .320 - .400
!      
!       PAR     8.     .400 - .700
!
!----- Input parameters:                            units      size
!
!  number of soundings (m)                          n/d         1
!  number of atmospheric layers (np)                n/d         1
!  cosine of solar zenith angle (cosz)              n/d         m
!  layer scaled-water vapor content (wh)          gm/cm^2      m*np
!  layer ozone content (oh)                      (cm-atm)stp   m*np
!  layer pressure thickness (dp)                    mb         m*np
!  option for scaling cloud optical thickness       n/d         1
!        overcast="true" if scaling is NOT required
!        overcast="false" if scaling is required
!  cloud water amount (cwp)                        gm/m**2     m*np*3
!        index 1 for ice particles
!        index 2 for liquid drops
!        index 3 for rain drops
!  effective cloud-particle size (reff)          micrometer    m*np*3
!       index 1 for ice particles
!       index 2 for liquid drops
!       index 3 for rain drops
!  level index separating high and                  n/d         m
!       middle clouds (ict)
!  level index separating middle and                n/d         m
!       low clouds (icb)
!  cloud amount (fcld)                            fraction     m*np
!  aerosol optical thickness (taual)                n/d        m*np*11
!  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
!  aerosol asymmetry factor (asyal)                 n/d        m*np*11
!  reflection of ground/ocean                     fraction      m
!   in uv+parc (bands 1-8)      
!        (rgbuv) for beam insolation
!        (rgfuv) for diffuse insolation
!
!---- temporary array
!
!  scaled cloud optical thickness                   n/d        m*np
!       for beam radiation (tauclb)
!  scaled cloud optical thickness                   n/d        m*np
!       for diffuse radiation  (tauclf)     
!
!----- output (updated) parameters:
!
!  all-sky net downward flux (flx)               fraction      m*(np+1)
!  clear-sky net downward flux (flc)             fraction      m*(np+1)
!  all-sky direct downward uv flux at
!       the surface (fdiruv)                     fraction       m
!  all-sky diffuse downward uv flux at
!       the surface (fdifuv)                     fraction       m
!  all-sky direct downward par flux at
!       the surface (fdirpar)                    fraction       m
!  all-sky diffuse downward par flux at
!       the surface (fdifpar)                    fraction       m
!
!***********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1
      integer, dimension(:), intent(in) :: ict, icb

      real (kind=kind_phys), dimension(:,:),   intent(in) :: wh, oh,    &
     &       dp, fcld
      real (kind=kind_phys), dimension(:),     intent(in) :: cosz,      &
     &       rgbuv, rgfuv
      real (kind=kind_phys), dimension(:,:,:), intent(in) :: cwp, reff, &
     &       taual, ssaal, asyal

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: flxu,     &
     &       flxd, flcu, flcd

!! ---  optional outputs:
      real (kind=kind_phys), dimension(:), intent(out) :: fdiruv,       &
     &       fdifuv, fdirpar, fdifpar, flxuvb, flcuvb
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: fnetb

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) ::  tausto, ssatau,   &
     &       asysto, tautob, ssatob, asytob, tautof, ssatof, asytof,    &
     &       tauclb, tauclf, asycl, rrt, ttt, tdt, rst, tst, dum
      real (kind=kind_phys), dimension(NPTS,NLP1) ::  fxup, fxdn, fcup, &
     &       fcdn, dum1
      real (kind=kind_phys), dimension(NPTS)      ::  fsdir, fsdif,     &
     &       asyclt, dsm, dum2
      real (kind=kind_phys), dimension(NPTS,NLP1,2):: rr, tt, td, rs, ts
      real (kind=kind_phys) :: taucld(NPTS,NLAY,3), cc(NPTS,3), tauc,   &
     &       reff1, reff2, taurs, tauoz, tauwv, g1, g2, g3, xx

      integer :: i, j, k, in, ib, ih1, ih2, im1, im2, is1, is2

!  ---  hk is the fractional extra-terrestrial solar flux in each
!       of the 8 bands. the sum of hk is 0.47074. (table 3)
      real (kind=kind_phys) :: hk(NUVVIS)
      data hk  /  .00057, .00367, .00083, .00417,                       &
     &            .00600, .00556, .05913, .39081  /

!  ---  zk is the ozone absorption coefficient. unit: /(cm-atm)stp
!       (table 3)
      real (kind=kind_phys) :: zk(NUVVIS)
      data zk  /  30.47,187.2,301.9,42.83,7.09,1.25,0.0345,0.0572  /

!  ---  wk is the water vapor absorption coefficient. unit: cm**2/g
!       (table 3)
      real (kind=kind_phys) :: wk(NUVVIS)
      data wk  /  7*0.0, 0.00075  /

!  ---  ry is the extinction coefficient for Rayleigh scattering.
!       unit: /mb. (table 3)
      real (kind=kind_phys) :: ry(NUVVIS)
      data ry  /  .00604, .00170, .00222, .00132,                       &
     &            .00107, .00091, .00055, .00012  /

!  ---  coefficients for computing the extinction coefficients of ice, 
!       water, and rain particles, independent of spectral band. (table 4)
      real (kind=kind_phys) :: aib, awb(2), arb(2)
      data aib  /  1.64  /
      data awb  /  -6.59e-3, 1.65  /
      data arb  /  3.07e-3, 0.00  /

!  ---  coefficients for computing the asymmetry factor of ice, water,
!       and rain particles, independent of spectral band. (table 6)
      real (kind=kind_phys) :: aig(3), awg(3), arg(3)
      data aig  /  .746,   .00282, -.0000230  /
      data awg  /  .82562, .00529, -.0001487  /
      data arg  /  .883,    0.0,     0.0      /
!
!===> ...  begin here
!

!  ---  initialize fdiruv, fdifuv, surface reflectances and transmittances.
!       the reflectance and transmittance of the clear and cloudy portions
!       of a layer are denoted by 1 and 2, respectively.
!       cc is the maximum cloud cover in each of the high, middle, and low
!       cloud groups.
!       1/dsm=1/cos(53) = 1.66
            
      do i = 1, NPTS                
         dsm(i) = 0.602

         rr(i,NLP1,1) = rgbuv(i)
         rr(i,NLP1,2) = rgbuv(i)
         rs(i,NLP1,1) = rgfuv(i)
         rs(i,NLP1,2) = rgfuv(i)
         td(i,NLP1,1) = _zero
         td(i,NLP1,2) = _zero
         tt(i,NLP1,1) = _zero
         tt(i,NLP1,2) = _zero
         ts(i,NLP1,1) = _zero
         ts(i,NLP1,2) = _zero

         cc(i,1) = _zero
         cc(i,2) = _zero
         cc(i,3) = _zero
       enddo

!! ---  for optional fluxes:

      if ( lfdncmp ) then
        do i = 1, NPTS                
           fdiruv(i) = _zero
           fdifuv(i) = _zero
           flxuvb(i) = _zero
           flcuvb(i) = _zero
         enddo
       endif

!  ---  compute cloud optical thickness.  eqs. (4.6) and (4.11)
!       note: the cloud optical properties are assumed to be independent
!       of spectral bands in the UV and PAR regions. the indices 1, 2, 3
!       are for ice, water, and rain particles, respectively.


       do k = 1, NLAY
       do i = 1, NPTS
         taucld(i,k,1) = cwp(i,k,1)*aib/reff(i,k,1)
         taucld(i,k,2) = cwp(i,k,2)*(awb(1)+awb(2)/reff(i,k,2))
         taucld(i,k,3) = cwp(i,k,3)* arb(1)
       enddo
       enddo

!  ---  options for scaling cloud optical thickness

      if ( iflagc == 1 ) then

        do k = 1, NLAY
        do i = 1, NPTS
          tauclb(i,k) = taucld(i,k,1) + taucld(i,k,2) + taucld(i,k,3)
          tauclf(i,k) = tauclb(i,k)
        enddo
        enddo

        do k = 1, 3
        do i = 1, NPTS
          cc(i,k) = _one
        enddo
        enddo

      elseif ( iflagc == 2 ) then

!  ---  scale cloud optical thickness in each layer from taucld (with
!       cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
!       tauclb is the scaled optical thickness for beam radiation and
!       tauclf is for diffuse radiation (see section 7).

        call cldscale                                                   &
!  ---  inputs:
     &     ( cosz, fcld, taucld, ict, icb,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       cc,tauclb,tauclf                                           &
     &     )

      endif

!  ---  cloud asymmetry factor for a mixture of liquid and ice particles.
!       unit of reff is micrometers. eqs. (4.8) and (6.4)

      do k = 1, NLAY
        do i = 1, NPTS
          asyclt(i) = _one
          tauc = taucld(i,k,1) + taucld(i,k,2) + taucld(i,k,3)

          if ( tauc > 0.02 .and. fcld(i,k) > 0.01 ) then
            reff1 = reff(i,k,1)
            reff2 = reff(i,k,2)

            g1 = (aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucld(i,k,1)
            g2 = (awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucld(i,k,2)
            g3 =  arg(1)*taucld(i,k,3)
            asyclt(i) = (g1 + g2 + g3) / tauc
          endif
        enddo

        do i = 1, NPTS
          asycl(i,k) = asyclt(i)
        enddo
      enddo

!  ---  integration over spectral bands

      do 100 ib = 1, NUVVIS

        do k = 1, NLAY
        do i = 1, NPTS

!  ---  compute clear-sky optical thickness, single scattering albedo,
!       and asymmetry factor (eqs. 6.2-6.4)
          taurs = ry(ib)*dp(i,k)
          tauoz = zk(ib)*oh(i,k)
          tauwv = wk(ib)*wh(i,k)
          xx = ssaal(i,k,ib) * taual(i,k,ib)

          tausto(i,k) = taurs + tauoz + tauwv + taual(i,k,ib) + fpmin
          ssatau(i,k) = xx + taurs
          asysto(i,k) = xx * asyal(i,k,ib)

          tautob(i,k)=tausto(i,k)
          asytob(i,k)=asysto(i,k)/ssatau(i,k)
          ssatob(i,k)=ssatau(i,k)/tautob(i,k) + fpmin
          ssatob(i,k)=min(ssatob(i,k), fpmax)
        enddo
        enddo

!  ---  for direct incident radiation

         call deledd                                                    &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, cosz,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rrt, ttt, tdt                                              &
     &     )

!  ---  diffuse incident radiation is approximated by beam radiation with
!       an incident angle of 53 degrees, eqs. (6.5) and (6.6)

         call deledd                                                    &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, dsm,                               &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rst, tst, dum                                              &
     &     )

         do k = 1, NLAY
         do i = 1, NPTS
           rr(i,k,1) = rrt(i,k)
           tt(i,k,1) = ttt(i,k)
           td(i,k,1) = tdt(i,k)
           rs(i,k,1) = rst(i,k)
           ts(i,k,1) = tst(i,k)
         enddo
         enddo

!  ---  compute reflectance and transmittance of the cloudy portion of a layer

        do k = 1, NLAY
        do i = 1, NPTS

!  ---  for direct incident radiation
!       the effective layer optical properties. eqs. (6.2)-(6.4)
          tautob(i,k) = tausto(i,k) + tauclb(i,k)
          ssatob(i,k) = (ssatau(i,k)+ tauclb(i,k))/tautob(i,k) + fpmin
          ssatob(i,k) = min( ssatob(i,k), fpmax )
          asytob(i,k) = (asysto(i,k)+ asycl(i,k)*tauclb(i,k))           &
     &                / (ssatob(i,k)*tautob(i,k))

!  ---  for diffuse incident radiation
          tautof(i,k) = tausto(i,k) + tauclf(i,k)
          ssatof(i,k) = (ssatau(i,k)+ tauclf(i,k))/tautof(i,k) + fpmin
          ssatof(i,k) = min( ssatof(i,k), fpmax )
          asytof(i,k) = (asysto(i,k)+ asycl(i,k)*tauclf(i,k))           &
     &                / (ssatof(i,k)*tautof(i,k))

        enddo
        enddo

!  ---  for direct incident radiation
!       note that the cloud optical thickness is scaled differently for
!       direct and diffuse insolation, eqs. (7.3) and (7.4).

        call deledd                                                     &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, cosz,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rrt, ttt, tdt                                              &
     &     )

!  ---  diffuse incident radiation is approximated by beam radiation with
!       an incident angle of 53 degrees, eqs. (6.5) and (6.6)

        call deledd                                                     &
!  ---  inputs:
     &     ( tautof, ssatof, asytof, dsm,                               &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rst, tst, dum                                              &
     &     )

        do k = 1, NLAY
        do i = 1, NPTS
          rr(i,k,2) = rrt(i,k)
          tt(i,k,2) = ttt(i,k)
          td(i,k,2) = tdt(i,k)
          rs(i,k,2) = rst(i,k)
          ts(i,k,2) = tst(i,k)
        enddo
        enddo

!  ---  flux calculations
!       initialize clear-sky flux (fc), all-sky flux (fx), 
!       and surface downward fluxes (fsdir and fsdif)

        do k = 1, NLP1
        do i = 1, NPTS
          fxup(i,k) = _zero
          fxdn(i,k) = _zero
          fcup(i,k) = _zero
          fcdn(i,k) = _zero
        enddo
        enddo

        do i = 1, NPTS
          fsdir(i) = _zero
          fsdif(i) = _zero
        enddo

!  ---  for clear- and all-sky flux calculations when fractional cloud cover
!       is either 0 or 1.

        if ( iflagc == 1 ) then

          call cldflxy                                                  &
!  ---  inputs:
     &     ( rr, tt, td, rs, ts,                                        &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

        elseif ( iflagc == 2 ) then

!  ---  for clear- and all-sky flux calculations when fractional cloud cover
!       is allowed to be between 0 and 1. the all-sky flux, fall is the
!       summation inside the brackets of eq. (7.11)

          ih1 = 1
          ih2 = 2
          im1 = 1
          im2 = 2
          is1 = 1
          is2 = 2

          call cldflx                                                   &
!  ---  inputs:
     &     ( ict,icb,ih1,ih2,im1,im2,is1,is2,                           &
     &       cc, rr, tt, td, rs, ts,                                    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

        endif

!  ---  flux integration, Eq. (6.1)

        do k = 1, NLP1
        do i = 1, NPTS
          flxu(i,k) = flxu(i,k) + fxup(i,k)*hk(ib)
          flxd(i,k) = flxd(i,k) + fxdn(i,k)*hk(ib)
          flcu(i,k) = flcu(i,k) + fcup(i,k)*hk(ib)
          flcd(i,k) = flcd(i,k) + fcdn(i,k)*hk(ib)
        enddo
        enddo

!! ---  optional spectral band net flux
        if ( lhswb ) then
          do k = 1, NLP1
          do i = 1, NPTS
            fnetb(i,k,ib) = (fxdn(i,k) - fxup(i,k)) * hk(ib)
          enddo
          enddo
        endif

!! ---  optional for computing direct and diffuse downward surface fluxes
!!      in the uv and par regions, and the uv-b spectrum clear and total
!!      sky fluxes

        if ( lfdncmp ) then
          if ( ib < NUVVIS ) then
            do i = 1, NPTS
              fdiruv(i) = fdiruv(i) + fsdir(i)*hk(ib)
              fdifuv(i) = fdifuv(i) + fsdif(i)*hk(ib)
            enddo

            if ( ib >= NUVBS .and. ib <= NUVBE ) then
              do i = 1, NPTS
                flxuvb(i) = flxuvb(i) + fxdn(i,NLP1)*hk(ib)
                flcuvb(i) = flcuvb(i) + fcdn(i,NLP1)*hk(ib)
              enddo
            endif
          else
            do i = 1, NPTS
              fdirpar(i) = fsdir(i)*hk(ib)
              fdifpar(i) = fsdif(i)*hk(ib)
            enddo
          endif
        endif

 100  continue

!
      return
!...................................
      end subroutine soluv
!-----------------------------------


!-----------------------------------
      subroutine solir                                                  &
!...................................

!  ---  inputs:
     &     ( cosz,wh,dp,fcld,cwp,reff,ict,icb,                          &
     &       taual,ssaal,asyal,rgbir,rgfir,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       flxu,flxd,flcu,flcd                                        &
!! ---  optional outputs:
     &,      fdirir,fdifir, fnetb                                       &
     &     )

!************************************************************************
!  compute solar flux in the infrared region. The spectrum is divided
!   into three bands:
!
!          band   wavenumber(/cm)  wavelength (micron)
!          1( 9)    14280-8200         0.70-1.22
!          2(10)     8200-4400         1.22-2.27
!          3(11)     4400-1000         2.27-10.0
!
!----- Input parameters:                            units      size
!
!  number of soundings (m)                          n/d         1
!  number of atmospheric layers (np)                n/d         1
!  cosine of solar zenith angle (cosz)              n/d         m
!  layer scaled-water vapor content (wh)          gm/cm^2      m*np
!  thickness of a layer (dp)                        mb         m*np
!  option for scaling cloud optical thickness       n/d         1
!        overcast="true" if scaling is NOT required
!        overcast="false" if scaling is required
!  cloud water concentration (cwp)                gm/m**2      m*np*3
!        index 1 for ice particles
!        index 2 for liquid drops
!        index 3 for rain drops
!  effective cloud-particle size (reff)           micrometer   m*np*3
!        index 1 for ice particles
!        index 2 for liquid drops
!        index 3 for rain drops
!  level index separating high and                  n/d        m
!        middle clouds (ict)
!  level index separating middle and                n/d        m
!        low clouds (icb)
!  cloud amount (fcld)                            fraction     m*np
!  aerosol optical thickness (taual)                n/d        m*np*11
!  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
!  aerosol asymmetry factor (asyal)                 n/d        m*np*11
!  reflection of ground/ocean in ir               fraction      m
!        (rgbir) for beam insolation
!        (rgfir) for diffuse insolation
!
!---- temporary array
!
!  scaled cloud optical thickness                   n/d        m*np
!          for beam radiation (tauclb)
!  scaled cloud optical thickness                   n/d        m*np
!          for diffuse radiation  (tauclf)     
!
!----- output (updated) parameters:
!
!  all-sky flux (downward-upward) (flx)           fraction     m*(np+1)
!  clear-sky flux (downward-upward) (flc)         fraction     m*(np+1)
!  all-sky direct downward ir flux at
!          the surface (fdirir)                   fraction     m
!  all-sky diffuse downward ir flux at
!          the surface (fdifir)                   fraction     m
!
!**********************************************************************

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1
      integer, dimension(:), intent(in) :: ict, icb

      real (kind=kind_phys), dimension(:,:),  intent(in) :: wh,dp, fcld
      real (kind=kind_phys), dimension(:),    intent(in) :: cosz, rgbir,&
     &       rgfir
      real (kind=kind_phys), dimension(:,:,:),intent(in) :: cwp, reff,  &
     &       taual, ssaal, asyal

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: flxu,     &
     &       flxd, flcu, flcd

!! ---  optional outputs:
      real (kind=kind_phys), dimension(:), intent(out) :: fdirir, fdifir
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: fnetb

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: tauclb, tauclf,    &
     &       ssacl, asycl, tausto, ssatau, asysto, tautob, ssatob,      &
     &       asytob, tautof, ssatof, asytof, rrt, ttt, tdt, rst,        &
     &       tst, dum
      real (kind=kind_phys), dimension(NPTS,NLP1) :: fxup, fxdn, fcup,  &
     &       fcdn, dum1
      real (kind=kind_phys), dimension(NPTS)      :: fsdir, fsdif,      &
     &       ssaclt, asyclt, dsm, dum2
      real (kind=kind_phys), dimension(NPTS,NLP1,2):: rr, tt, td, rs, ts

      real (kind=kind_phys) :: taucld(NPTS,NLAY,3), cc(NPTS,3), taurs,  &
     &       tauwv, tauc, reff1, reff2, w1, w2, w3, g1, g2, g3, xx

      integer :: i, j, k, in, ib, iv, ik, ih1, ih2, im1, im2, is1, is2

!  ---  water vapor absorption coefficient for 10 k-intervals.
!       unit: cm^2/gm (table 2)
      real (kind=kind_phys) :: xk(NK0)
      data xk  /  0.0010, 0.0133, 0.0422, 0.1334, 0.4217,               &
     &            1.334,  5.623,  31.62,  177.8,  1000.0  /

!  ---  water vapor k-distribution function,
!       the sum of hk is 0.52926. unit: fraction (table 2)
      real (kind=kind_phys) :: hk(NIRBND,NK0)
      data hk  /  .20673,.08236,.01074,  .03497,.01157,.00360,          &
     &            .03011,.01133,.00411,  .02260,.01143,.00421,          &
     &            .01336,.01240,.00389,  .00696,.01258,.00326,          &
     &            .00441,.01381,.00499,  .00115,.00650,.00465,          &
     &            .00026,.00244,.00245,  .00000,.00094,.00145  /

!  ---  ry is the extinction coefficient for rayleigh scattering.
!       unit: /mb (table 3)
      real (kind=kind_phys) :: ry(NIRBND)
      data ry  /  .0000156, .0000018, .000000  /

!  ---  coefficients for computing the extinction coefficients of
!       ice, water, and rain particles (table 4)
      real (kind=kind_phys) :: aib, awb(NIRBND,2), arb(NIRBND,2)
      data aib  /  1.64  /
      data awb  /  -0.0101,-0.0166,-0.0339,  1.72, 1.85, 2.16  /
      data arb  /  0.00307,0.00307,0.00307,  0.0 , 0.0 , 0.0   /

!  ---  coefficients for computing the single-scattering co-albedo of
!       ice, water, and rain particles (table 5)
      real (kind=kind_phys), dimension(NIRBND,3) :: aia, awa, ara
      data aia  /   .00000141,   .00112,     .04828,                    &
     &              .00001144,   .001129,    .00547,                    &
     &             -.000000005, -.00000358, -.0000361  /
      data awa  /   .00000007, -.00019934,  .01209318,                  &
     &              .00000845,  .00088757,  .01784739,                  &
     &             -.00000004, -.00000650, -.00036910  /
      data ara  /  .029,.342,.466,  .000,.000,.000,  .000,.000,.000  /

!  ---  coefficients for computing the asymmetry factor of 
!       ice, water, and rain particles (table 6)
      real (kind=kind_phys), dimension(NIRBND,3) :: aig, awg, arg
      data aig  /  .725, .717, .771,    .0037, .00456, .00490,          &
     &            -.0000309,-.00003544,-.0000401  /
      data awg  /  .79375035,  .74513197, .83530748,                    &
     &             .00832441,  .01370071, .00257181,                    &
     &            -.00023263, -.00038203, .00005519  /
      data arg  /  .891,.948,.971,  .000,.000,.000,  .000,.000,.000  /
!
!===> ...  start here
!

!  ---  initialize surface fluxes, reflectances, and transmittances.
!       the reflectance and transmittance of the clear and cloudy portions
!       of a layer are denoted by 1 and 2, respectively. cc is the maximum
!       cloud cover in each of the high, middle, and low cloud groups.
!       1/dsm=1/cos(53)=1.66

      do i = 1, NPTS
        dsm(i) = 0.602

        rr(i,NLP1,1) = rgbir(i)
        rr(i,NLP1,2) = rgbir(i)
        rs(i,NLP1,1) = rgfir(i)
        rs(i,NLP1,2) = rgfir(i)
        td(i,NLP1,1) = _zero
        td(i,NLP1,2) = _zero
        tt(i,NLP1,1) = _zero
        tt(i,NLP1,2) = _zero
        ts(i,NLP1,1) = _zero
        ts(i,NLP1,2) = _zero

        cc(i,1) = _zero
        cc(i,2) = _zero
        cc(i,3) = _zero

!! ---  for optional fluxes:
        fdirir(i) = _zero
        fdifir(i) = _zero
      enddo

!  ---  integration over spectral bands

      do 100 ib = 1, NIRBND

        iv = ib + NUVVIS

!! ---  initial optional spectral net flux output array

        if ( lhswb ) then
          do k = 1, NLP1
          do i = 1, NPTS
            fnetb(i,k,iv) = _zero
          enddo
          enddo
        endif

!  ---  compute cloud optical thickness. eqs. (4.6) and (4.10)
!       the indices 1, 2, 3 are for ice, water, rain particles, respectively.

        do k = 1, NLAY
        do i = 1, NPTS
          taucld(i,k,1) = cwp(i,k,1)*aib/reff(i,k,1)
          taucld(i,k,2) = cwp(i,k,2)*(awb(ib,1)+awb(ib,2)/reff(i,k,2))
          taucld(i,k,3) = cwp(i,k,3)*arb(ib,1)
        enddo
        enddo

!  ---  options for scaling cloud optical thickness

        if ( iflagc == 1 ) then

          do k = 1, NLAY
          do i = 1, NPTS
            tauclb(i,k) = taucld(i,k,1) + taucld(i,k,2) + taucld(i,k,3)
            tauclf(i,k) = tauclb(i,k)
          enddo
          enddo

          do k = 1, 3
          do i = 1, NPTS
            cc(i,k) = _one
          enddo
          enddo

        elseif ( iflagc == 2 ) then

!  ---  scale cloud optical thickness in each layer from taucld (with
!       cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
!       tauclb is the scaled optical thickness for beam radiation and
!       tauclf is for diffuse radiation.

          call cldscale                                                 &
!  ---  inputs:
     &     ( cosz, fcld, taucld, ict, icb,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       cc,tauclb,tauclf                                           &
     &     )

        endif  ! end if_iflagc_block

!  ---  compute cloud single scattering albedo and asymmetry factor
!       for a mixture of ice and liquid particles.
!       eqs.(4.6)-(4.8), (6.2)-(6.4)

        do k = 1, NLAY

          do i = 1, NPTS
            ssaclt(i) = 0.99999
            asyclt(i) = _one
            tauc = taucld(i,k,1) + taucld(i,k,2) + taucld(i,k,3)

            if ( tauc > 0.02 .and. fcld(i,k) > 0.01 ) then
              reff1 = reff(i,k,1)
              reff2 = reff(i,k,2)

              w1 = (_one - ( aia(ib,1) + ( aia(ib,2)                    &
     &           + aia(ib,3)*reff1 )* reff1 )) * taucld(i,k,1)
              w2 = (_one - ( awa(ib,1) + ( awa(ib,2)                    &
     &           + awa(ib,3)*reff2 )* reff2 )) * taucld(i,k,2)
              w3 = (_one - ara(ib,1) ) * taucld(i,k,3)
              ssaclt(i) = (w1 + w2 + w3) / tauc

              g1 = (aig(ib,1) + (aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
              g2 = (awg(ib,1) + (awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
              g3 =  arg(ib,1) * w3
              asyclt(i) = (g1 + g2 + g3) / (w1 + w2 + w3)
            endif
          enddo

          do i = 1, NPTS
            ssacl(i,k) = ssaclt(i)
            asycl(i,k) = asyclt(i)
          enddo
        enddo

!  ---  integration over the k-distribution function

        do 200 ik = 1, NK0

          do k = 1, NLAY
          do i = 1, NPTS
            taurs = ry(ib) * dp(i,k)
            tauwv = xk(ik) * wh(i,k)
            xx = ssaal(i,k,iv) * taual(i,k,iv)
 
!  ---  compute clear-sky optical thickness, single scattering albedo,
!       and asymmetry factor. eqs.(6.2)-(6.4)

            tausto(i,k) = taurs + tauwv + taual(i,k,iv) + fpmin
            ssatau(i,k) = xx + taurs + fpmin
            asysto(i,k) = xx * asyal(i,k,iv)

!  ---  compute reflectance and transmittance of the clear portion of a layer

            tautob(i,k) = tausto(i,k)
            asytob(i,k) = asysto(i,k) / ssatau(i,k)
            ssatob(i,k) = ssatau(i,k) / tautob(i,k) + fpmin
            ssatob(i,k) = min( ssatob(i,k), fpmax )
          enddo
          enddo

!  ---  for direct incident radiation

          call deledd                                                   &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, cosz,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rrt, ttt, tdt                                              &
     &     )

!  ---  diffuse incident radiation is approximated by beam radiation with
!       an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

          call deledd                                                   &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, dsm,                               &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rst, tst, dum                                              &
     &     )

          do k = 1, NLAY
          do i = 1, NPTS
            rr(i,k,1) = rrt(i,k)
            tt(i,k,1) = ttt(i,k)
            td(i,k,1) = tdt(i,k)
            rs(i,k,1) = rst(i,k)
            ts(i,k,1) = tst(i,k)
          enddo
          enddo

!  ---  compute reflectance and transmittance of the cloudy portion of a layer

          do k = 1, NLAY
          do i = 1, NPTS

!  ---  for direct incident radiation. Eqs.(6.2)-(6.4)
            tautob(i,k) = tausto(i,k) + tauclb(i,k)
            ssatob(i,k) = (ssatau(i,k) + ssacl(i,k)*tauclb(i,k))        &
     &                  / tautob(i,k) + fpmin
            ssatob(i,k) = min( ssatob(i,k), fpmax )
            asytob(i,k) = (asysto(i,k)                                  &
     &                  + asycl(i,k)*ssacl(i,k)*tauclb(i,k))            &
     &                  / (ssatob(i,k)*tautob(i,k))

!  ---  for diffuse incident radiation
            tautof(i,k) = tausto(i,k) + tauclf(i,k)
            ssatof(i,k) = (ssatau(i,k) + ssacl(i,k)*tauclf(i,k))        &
     &                  / tautof(i,k) + fpmin
            ssatof(i,k) = min( ssatof(i,k), fpmax )
            asytof(i,k) = (asysto(i,k)                                  &
     &                  + asycl(i,k)*ssacl(i,k)*tauclf(i,k))            &
     &                  / (ssatof(i,k)*tautof(i,k))
          enddo
          enddo

!  ---  for direct incident radiation

          call deledd                                                   &
!  ---  inputs:
     &     ( tautob, ssatob, asytob, cosz,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rrt, ttt, tdt                                              &
     &     )

!  ---  diffuse incident radiation is approximated by beam radiation with
!       an incident angle of 53 degrees, eqs.(6.5) and (6.6)

          call deledd                                                   &
!  ---  inputs:
     &     ( tautof, ssatof, asytof, dsm,                               &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rst, tst, dum                                              &
     &     )

          do k = 1, NLAY
          do i = 1, NPTS
            rr(i,k,2) = rrt(i,k)
            tt(i,k,2) = ttt(i,k)
            td(i,k,2) = tdt(i,k)
            rs(i,k,2) = rst(i,k)
            ts(i,k,2) = tst(i,k)
          enddo
          enddo

!  ---  flux calculations

!  ---  initialize clear-sky flux (fc), all-sky flux (fx), 
!       and surface downward fluxes (fsdir and fsdif)

          do k = 1, NLP1
          do i = 1, NPTS
            fxup(i,k) = _zero
            fxdn(i,k) = _zero
            fcup(i,k) = _zero
            fcdn(i,k) = _zero
          enddo
          enddo

          do i = 1, NPTS
            fsdir(i) = _zero
            fsdif(i) = _zero
          enddo

!  ---  for clear- and all-sky flux calculations when fractional cloud cover
!       is either 0 or 1.

          if ( iflagc == 1 ) then

            call cldflxy                                                &
!  ---  inputs:
     &     ( rr, tt, td, rs, ts,                                        &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

          elseif ( iflagc == 2 ) then

!  ---  for clear- and all-sky flux calculations when fractional cloud cover
!       is allowed to be between 0 and 1. the all-sky flux, fall is the
!       summation inside the brackets of eq. (7.11)

            ih1 = 1
            ih2 = 2
            im1 = 1
            im2 = 2
            is1 = 1
            is2 = 2

            call cldflx                                                 &
!  ---  inputs:
     &     ( ict,icb,ih1,ih2,im1,im2,is1,is2,                           &
     &       cc, rr, tt, td, rs, ts,                                    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

          endif

!  ---  flux integration following eq. (6.1)

          do k = 1, NLP1
          do i = 1, NPTS
            flxu(i,k) = flxu(i,k) + fxup(i,k)*hk(ib,ik)
            flxd(i,k) = flxd(i,k) + fxdn(i,k)*hk(ib,ik)
            flcu(i,k) = flcu(i,k) + fcup(i,k)*hk(ib,ik)
            flcd(i,k) = flcd(i,k) + fcdn(i,k)*hk(ib,ik)
          enddo
          enddo

!! ---  optional spectral band net flux
          if ( lhswb ) then
            do k = 1, NLP1
            do i = 1, NPTS
              fnetb(i,k,iv) = fnetb(i,k,iv)                             &
     &                      + (fxdn(i,k) - fxup(i,k))*hk(ib,ik)
            enddo
            enddo
          endif

!! ---  optional downward surface fluxes in the ir region
          if ( lfdncmp ) then
            do i = 1, NPTS
              fdirir(i) = fdirir(i) + fsdir(i)*hk(ib,ik)
              fdifir(i) = fdifir(i) + fsdif(i)*hk(ib,ik)
            enddo
          endif

  200   continue    ! end do_ik_loop
  100 continue    ! end do_iv_loop

!
      return
!...................................
      end subroutine solir
!-----------------------------------


!-----------------------------------
      subroutine cldscale                                               &
!...................................

!  ---  inputs:
     &     ( cosz, fcld, taucld, ict, icb,                              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       cc,tauclb,tauclf                                           &
     &     )

!********************************************************************
!
!   This subroutine computes the high, middle, and low cloud
!    amounts and scales the cloud optical thickness (Section 7)
!
!   To simplify calculations in a cloudy atmosphere, clouds are
!    grouped into high, middle and low clouds separated by the levels
!    ict and icb (level 1 is the top of the model atmosphere).
!
!   Within each of the three groups, clouds are assumed maximally
!    overlapped, and the cloud cover (cc) of a group is the maximum
!    cloud cover of all the layers in the group.  The optical thickness
!    (taucld) of a given layer is then scaled to new values (tauclb and
!    tauclf) so that the layer reflectance corresponding to the cloud
!    cover cc is the same as the original reflectance with optical
!    thickness taucld and cloud cover fcld.
!
!---input parameters
!
!    number of atmospheric soundings (m)
!    number of atmospheric layers (np)
!    cosine of the solar zenith angle (cosz)
!    fractional cloud cover (fcld)
!    cloud optical thickness (taucld)
!    index separating high and middle clouds (ict)
!    index separating middle and low clouds (icb)
!
!---output parameters
!
!    fractional cover of high, middle, and low cloud groups (cc)
!    scaled cloud optical thickness for direct  radiation (tauclb)
!    scaled cloud optical thickness for diffuse radiation (tauclf)
!
!********************************************************************
!
      use module_radsw_cldsca

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY
      integer, dimension(NPTS), intent(in) :: ict, icb

      real (kind=kind_phys), intent(in) :: cosz(NPTS),                  &
     &       fcld(NPTS,NLAY), taucld(NPTS,NLAY,3)

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: cc(NPTS,3),                 &
     &       tauclb(NPTS,NLAY), tauclf(NPTS,NLAY)

!  ---  locals:
      integer :: i, j, k, in, im, it, ia, kk, ib

      real (kind=kind_phys) :: fm, fn, fa, xai, tauc

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: dm = 0.1      ! size of cosz-interval
      real (kind=kind_phys), parameter :: dn = 0.30103  ! size of taucld-interval
      real (kind=kind_phys), parameter :: da = 0.1      ! size of cld amt-interval
      real (kind=kind_phys), parameter :: t1 = -0.9031

!  ---  include the pre-computed table of mcai for scaling the cloud optical
!       thickness under the assumption that clouds are maximally overlapped
!
!       caib is for scaling the cloud optical thickness for direct radiation
!       caif is for scaling the cloud optical thickness for diffuse radiation

!  ---  clouds within each of the high, middle, and low clouds are assumed
!       to be maximally overlapped, and the cloud cover (cc) for a group
!       (high, middle, or low) is the maximum cloud cover of all the layers
!       within a group

      do i = 1, NPTS
        it = ict(i)
        ib = icb(i)
        cc(i,1) = _zero
        cc(i,2) = _zero
        cc(i,3) = _zero

        do k = 1, NLAY
          if     ( k < it ) then
            cc(i,1) = max( cc(i,1), fcld(i,k) )
          elseif ( k < ib ) then
            cc(i,2) = max( cc(i,2), fcld(i,k) )
          else
            cc(i,3) = max( cc(i,3), fcld(i,k) )
          endif
        enddo
      enddo

!  ---  scale the cloud optical thickness.
!       taucld(i,k,1) is the optical thickness for ice particles
!       taucld(i,k,2) is the optical thickness for liquid particles
!       taucld(i,k,3) is the optical thickness for rain drops
      
      do i = 1, NPTS
        it = ict(i)
        ib = icb(i)

        do k = 1, NLAY

          if     ( k < it ) then
            kk = 1
          elseif ( k >= it .and. k < ib ) then
            kk = 2
          else
            kk = 3
          endif

          tauclb(i,k) = _zero
          tauclf(i,k) = _zero
          tauc = taucld(i,k,1) + taucld(i,k,2) + taucld(i,k,3)

          if ( tauc > 0.02 .and. fcld(i,k) > 0.01 ) then

!  ---  normalize cloud cover following Eq. (7.8)
            fa = fcld(i,k) / cc(i,kk)

!  ---  table look-up
            tauc = min( tauc, 32.0 )

            fm = cosz(i) / dm
            fn = (log10(tauc) - t1) / dn
            fa = fa / da
 
            im = max(2, min(NMCAI-1, int(fm+1.5) ) )
            in = max(2, min(NMCAI-1, int(fn+1.5) ) )
            ia = max(2, min(NMCAI-1, int(fa+1.5) ) )
  
            fm = fm - float(im-1)
            fn = fn - float(in-1)
            fa = fa - float(ia-1)

!  ---  scale cloud optical thickness for beam radiation following
!       eq. (7.3) the scaling factor, xai, is a function of the solar
!       zenith angle, optical thickness, and cloud cover.
 
            xai =       (-caib(ia,in,im-1)*(_one - fm)                  &
     &          +         caib(ia,in,im+1)*(_one + fm))*fm*0.5          &
     &          +         caib(ia,in,im)*(_one - fm*fm)
         
            xai = xai + (-caib(ia,in-1,im)*(_one - fn)                  &
     &          +         caib(ia,in+1,im)*(_one + fn))*fn*0.5          &
     &          +         caib(ia,in,im)*(_one - fn*fn)

            xai = xai + (-caib(ia-1,in,im)*(_one - fa)                  &
     &          +         caib(ia+1,in,im)*(_one + fa))*fa*0.5          &
     &          +         caib(ia,in,im)*(_one - fa*fa)

            xai = xai - 2.0*caib(ia,in,im)
            xai = max(_zero, min(_one, xai) )
     
            tauclb(i,k) = tauc * xai

!  ---  scale cloud optical thickness for diffuse radiation following
!       eq. (7.4) the scaling factor, xai, is a function of the cloud 
!       optical thickness and cover but not the solar zenith angle.

            xai =       (-caif(ia,in-1)*(_one - fn)                     &
     &          +         caif(ia,in+1)*(_one + fn))*fn*0.5             &
     &          +         caif(ia,in)*(_one - fn*fn)

            xai = xai + (-caif(ia-1,in)*(_one - fa)                     &
     &          +         caif(ia+1,in)*(_one + fa))*fa*0.5             &
     &          +         caif(ia,in)*(_one - fa*fa)

            xai = xai - caif(ia,in)
            xai = max(_zero, min(_one, xai) )

            tauclf(i,k) = tauc * xai
          endif

        enddo
      enddo

!
      return
!...................................
      end subroutine cldscale
!-----------------------------------


!-----------------------------------
      subroutine deledd                                                 &
!...................................

!  ---  inputs:
     &     ( tau, ssc, g0, cza,                                         &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       rr, tt, td                                                 &
     &     )

!*********************************************************************
!
!-----uses the delta-eddington approximation to compute the
!     bulk scattering properties of a single layer
!     coded following King and Harshvardhan (JAS, 1986)
!
!  inputs:
!       m:  number of soundings
!      np:  number of atmospheric layers
!     tau:  optical thickness
!     ssc:  single scattering albedo
!     g0:   asymmetry factor
!     cza:  cosine of the zenith angle
!
!  outputs:
!
!     rr:  reflection of the direct beam
!     tt:  total (direct+diffuse) transmission of the direct beam
!     td:  direct transmission of the direct beam
!
!*********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(NPTS,NLAY), intent(in) ::        &
     &       tau, ssc, g0
      real (kind=kind_phys), dimension(NPTS),      intent(in) :: cza

!  ---  outputs:
      real (kind=kind_phys), dimension(NPTS,NLAY), intent(out) ::       &
     &       rr, tt, td

!  ---  locals:
      integer :: i, k
      real (kind=kind_phys) :: zth, ff, xx, taup, sscp, gp, gm1, gm2,   &
     &       gm3, akk, alf1, alf2, all, bll, st7, st8, cll, dll, fll,   &
     &       ell, st1, st2, st3, st4
 
!---------------------------------------------------------------------

      do k = 1, NLAY
      do i = 1, NPTS
        zth = cza(i)
 
!  ---  delta-eddington scaling of single scattering albedo, optical
!       thickness, and asymmetry factor, k & h eqs(27-29)

        ff   = g0(i,k) * g0(i,k)
        xx   = _one - ff*ssc(i,k)
        taup = tau(i,k) * xx
        sscp = ssc(i,k) * (_one - ff) / xx
        gp   = g0(i,k) / (_one + g0(i,k))
 
!  ---  gamma1, gamma2, and gamma3. see table 2 and eq(26) k & h
!       ssc and gp are the d-s single scattering albedo and asymmetry
!       factor.

        xx  = 3.0 * gp 
        gm1 = (7.0  - sscp*(4.0 + xx)) * 0.25
        gm2 =-(_one - sscp*(4.0 - xx)) * 0.25

!  ---  akk is k as defined in eq(25) of k & h
 
        akk = sqrt( (gm1+gm2) * (gm1-gm2) )
 
        xx  = akk * zth
        st7 = _one - xx
        st8 = _one + xx
        st3 = st7 * st8

        if ( abs(st3) < fpmin ) then
          zth = zth + 0.001
          xx  = akk * zth
          st7 = _one - xx
          st8 = _one + xx
          st3 = st7 * st8
        endif

!  ---  extinction of the direct beam transmission
 
        td(i,k) = exp( -taup/zth )

!  ---  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of k & h
 
        gm3 = (2.0 - zth*3.0*gp) * 0.25
        xx   = gm1 - gm2
        alf1 = gm1 - gm3*xx
        alf2 = gm2 + gm3*xx
 
!  ---  all is last term in eq(21) of k & h
!  ---  bll is last term in eq(22) of k & h
 
        xx  = akk * 2.0
        all = (gm3 - alf2*zth) * xx * td(i,k)
        bll = (_one - gm3 + alf1*zth) * xx
 
        xx  = akk * gm3
        cll = (alf2 + xx) * st7
        dll = (alf2 - xx) * st8
 
        xx  = akk * (_one - gm3)
        fll = (alf1 + xx) * st8
        ell = (alf1 - xx) * st7
  
        st2 = exp( -akk*taup )
        st4 = st2*st2
        st1 = sscp / ((akk + gm1 + (akk-gm1)*st4) * st3)
 
!  ---  rr is r-hat of eq(21) of k & h
!  ---  tt is diffuse part of t-hat of eq(22) of k & h
 
        rr(i,k) = ( cll - dll*st4          - all*st2) * st1
        tt(i,k) =-((fll - ell*st4)*td(i,k) - bll*st2) * st1
 
        rr(i,k) = max( rr(i,k), _zero )
        tt(i,k) = max( tt(i,k), _zero )

        tt(i,k) = tt(i,k) + td(i,k)

       enddo
       enddo

!
      return
!...................................
      end subroutine deledd
!-----------------------------------


!-----------------------------------
      subroutine cldflxy                                                &
!...................................

!  ---  inputs:
     &     ( rr, tt, td, rs, ts,                                        &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

!*******************************************************************
!  This subroutine computes fluxes for the case that an atmospehric 
!    layer is either clear or totally cloudy.
!
!  upward and downward fluxes are computed using a two-stream adding 
!    method following equations (6.9)-(6.16).
!
!  input parameters:
!
!   m:   number of soundings
!   np:  number of atmospheric layers
!   rr:  reflection of a layer illuminated by beam radiation
!   tt:  total (direct+diffuse) transmission of a layer illuminated 
!        by beam radiation
!   td:  direct beam transmission
!   rs:  reflection of a layer illuminated by diffuse radiation
!   ts:  transmission of a layer illuminated by diffuse radiation
!
!  output parameters:
!
!     fcdn:  clear-sky downward flux
!     fcup:  clear-sky upward flux
!     fxdn:  all-sky downward flux
!     fxup:  all-sky upward flux
!     fsdir: surface direct downward flux
!     fsdif: surface diffuse downward flux
!
!*********************************************************************c
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: rr, tt,    &
     &       td, rs, ts

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),  intent(out) :: fcup,      &
     &       fcdn, fxup, fxdn
      real (kind=kind_phys), dimension(:),    intent(out) :: fsdir,fsdif

!  ---  locals:
      integer :: i, k, ih

      real (kind=kind_phys), dimension(NPTS,NLP1,2) :: rra, tta, tda,   &
     &       rsa, rxa
      real (kind=kind_phys), dimension(NPTS,NLP1)   :: flxdn, flxup
      real (kind=kind_phys), dimension(NPTS)        :: fdndir, fdndif
      real (kind=kind_phys) :: fupdif, denm, xx, yy
!
!===> ...  begin here
!

!  ---  compute transmittances and reflectances for a composite of
!       layers. layers are added one at a time, going down from the top.
!       tda is the composite direct transmittance illuminated by beam 
!           radiation
!       tta is the composite total transmittance illuminated by beam 
!           radiation
!       rsa is the composite reflectance illuminated from below by 
!           diffuse radiation
!       tta and rsa are computed from eqs. (6.10) and (6.12)

      do ih = 1, 2
        do i = 1, NPTS
          tda(i,1,ih) = td(i,1,ih)
          tta(i,1,ih) = tt(i,1,ih)
          rsa(i,1,ih) = rs(i,1,ih)
        enddo

        do k = 2, NLAY
        do i = 1, NPTS
          denm = ts(i,k,ih) / (_one - rsa(i,k-1,ih)*rs(i,k,ih))

          tda(i,k,ih) = tda(i,k-1,ih)*td(i,k,ih)
          tta(i,k,ih) = tda(i,k-1,ih)*tt(i,k,ih)                        &
     &                + (tda(i,k-1,ih)*rsa(i,k-1,ih)*rr(i,k,ih)         &
     &                + tta(i,k-1,ih) - tda(i,k-1,ih))*denm
          rsa(i,k,ih) = rs(i,k,ih) + ts(i,k,ih)*rsa(i,k-1,ih)*denm
        enddo
        enddo
      enddo

!  ---  layers are added one at a time, going up from the surface.
!       rra is the composite reflectance illuminated by beam radiation
!       rxa is the composite reflectance illuminated from above
!           by diffuse radiation
!       rra and rxa are computed from eqs. (6.9) and (6.11)

      do ih = 1, 2
        do i = 1, NPTS
          rra(i,NLP1,ih) = rr(i,NLP1,ih)
          rxa(i,NLP1,ih) = rs(i,NLP1,ih)
        enddo

        do k = NLAY, 1,-1
        do i = 1, NPTS
          denm = ts(i,k,ih) / (_one - rs(i,k,ih)*rxa(i,k+1,ih))

          rra(i,k,ih) = rr(i,k,ih) + (td(i,k,ih)*rra(i,k+1,ih)          &
     &                + (tt(i,k,ih) - td(i,k,ih))*rxa(i,k+1,ih))*denm
          rxa(i,k,ih) = rs(i,k,ih) + ts(i,k,ih)*rxa(i,k+1,ih)*denm
        enddo
        enddo
      enddo

!  ---  compute fluxes following Eq. (6.15) for fupdif and
!       eq. (6.16) for (fdndir+fdndif)
!       fdndir is the direct  downward flux
!       fdndif is the diffuse downward flux
!       fupdif is the diffuse upward flux

      do ih = 1, 2
        do k = 2, NLP1
        do i = 1, NPTS
          denm = _one / (_one - rsa(i,k-1,ih)*rxa(i,k,ih))

          xx = tda(i,k-1,ih) * rra(i,k,ih)
          yy = tta(i,k-1,ih) - tda(i,k-1,ih)

          fdndir(i) = tda(i,k-1,ih)
          fdndif(i) = (xx*rsa(i,k-1,ih) + yy) * denm
          fupdif    = (xx + yy*rxa(i,k,ih)) * denm

          flxdn(i,k) = fdndir(i) + fdndif(i)
          flxup(i,k) = fupdif
        enddo
        enddo

        do i = 1, NPTS
          flxdn(i,1) = _one
          flxup(i,1) = rra(i,1,ih)
        enddo

!  ---  ih=1 for clear-sky; ih=2 for overcast sky

        do k = 1, NLP1
        do i = 1, NPTS
          if ( ih == 1 ) then
            fcdn(i,k) = flxdn(i,k)
            fcup(i,k) = flxup(i,k)
          else
            fxdn(i,k) = flxdn(i,k)
            fxup(i,k) = flxup(i,k)
          endif
        enddo
        enddo

        if ( ih == 2 ) then
          do i = 1, NPTS
            fsdir(i) = fdndir(i)
            fsdif(i) = fdndif(i)
          enddo
        endif
      enddo 
!
      return
!...................................
      end subroutine cldflxy
!-----------------------------------


!-----------------------------------
      subroutine cldflx                                                 &
!...................................

!  ---  inputs:
     &     ( ict,icb,ih1,ih2,im1,im2,is1,is2,                           &
     &       cc, rr, tt, td, rs, ts,                                    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       fsdir, fsdif, fxup, fxdn, fcup, fcdn                       &
     &     )

!*******************************************************************
!  This subroutine computes fluxes for the case that cloud fraction of
!  an atmospehric can be any values ranging from 0 to 1.

!  compute upward and downward fluxes using a two-stream adding method
!  following equations (6.9)-(6.16).
!
!  clouds are grouped into high, middle, and low clouds which are assumed
!  randomly overlapped. It involves a maximum of 8 sets of calculations.
!  In each set of calculations, each atmospheric layer is homogeneous,
!  either totally filled with clouds or without clouds.

!  input parameters:
!
!   m:   number of soundings
!   np:  number of atmospheric layers
!   ict: the level separating high and middle clouds
!   icb: the level separating middle and low clouds
!   ih1,ih2,im1,im2,is1,is2: indices for three group of clouds
!   cc:  effective cloud covers for high, middle and low clouds
!   rr:  reflection of a layer illuminated by beam radiation
!   tt:  total (direct+diffuse) transmission of a layer illuminated 
!        by beam radiation
!   td:  direct beam transmission
!   rs:  reflection of a layer illuminated by diffuse radiation
!   ts:  transmission of a layer illuminated by diffuse radiation
!
!  output parameters:
!
!     fcdn:  clear-sky downward flux
!     fcup:  clear-sky upward flux
!     fxdn:  all-sky downward flux
!     fxup:  all-sky upward flux
!     fsdir: surface direct downward flux
!     fsdif: surface diffuse downward flux
!
!*********************************************************************c
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1
      integer, intent(in) :: ict(NPTS), icb(NPTS)
      integer, intent(in) :: ih1, ih2, im1, im2, is1, is2

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: rr, tt,    &
     &       td, rs, ts
      real (kind=kind_phys), dimension(:,:),   intent(in) :: cc

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: fxup, fxdn, &
     &       fcup, fcdn
      real (kind=kind_phys), dimension(:),   intent(out) :: fsdir,fsdif

!  ---  locals:
      integer :: i, k, ih, im, is

      real (kind=kind_phys), dimension(NPTS,NLP1,2,2) :: rra, tta,      &
     &       tda, rsa, rxa
      real (kind=kind_phys), dimension(NPTS,NLP1)     :: flxdn, flxup
      real (kind=kind_phys), dimension(NPTS)          :: ch, cm, ct,    &
     &       fdndir, fdndif
      real (kind=kind_phys) :: fupdif, denm, xx, yy
!
!===> ...  begin here
!

!  ---  compute transmittances and reflectances for a composite of
!       layers. layers are added one at a time, going down from the top.
!       tda is the composite direct transmittance illuminated by beam
!           radiation
!       tta is the composite total transmittance illuminated by beam
!           radiation
!       rsa is the composite reflectance illuminated from below by
!           diffuse radiation
!       tta and rsa are computed from eqs. (6.10) and (6.12)

!  ---  to save memory space, tda, tta, and rsa are pre-computed for k<icb.
!       the dimension of these parameters is (m,np,2,2). It would have
!       been (m,np,2,2,2) if these parameters were computed for all k's.

!  ---  for high clouds
!       ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

      do ih = ih1, ih2
        do i = 1, NPTS
          tda(i,1,ih,1) = td(i,1,ih)
          tta(i,1,ih,1) = tt(i,1,ih)
          rsa(i,1,ih,1) = rs(i,1,ih)
          tda(i,1,ih,2) = td(i,1,ih)
          tta(i,1,ih,2) = tt(i,1,ih)
          rsa(i,1,ih,2) = rs(i,1,ih)
        enddo

        do i = 1, NPTS
          do k = 2, ict(i)-1
            denm = ts(i,k,ih) / (_one - rsa(i,k-1,ih,1)*rs(i,k,ih))

            tda(i,k,ih,1) = tda(i,k-1,ih,1)*td(i,k,ih)
            tta(i,k,ih,1) = tda(i,k-1,ih,1)*tt(i,k,ih)                  &
     &                    + (tda(i,k-1,ih,1)*rsa(i,k-1,ih,1)*rr(i,k,ih) &
     &                    + tta(i,k-1,ih,1) - tda(i,k-1,ih,1)) * denm
            rsa(i,k,ih,1) = rs(i,k,ih)+ts(i,k,ih)*rsa(i,k-1,ih,1)*denm
            tda(i,k,ih,2) = tda(i,k,ih,1)
            tta(i,k,ih,2) = tta(i,k,ih,1)
            rsa(i,k,ih,2) = rsa(i,k,ih,1)
          enddo
        enddo

!  ---  for middle clouds
!       im=1 for clear-sky condition, im=2 for cloudy-sky condition

        do im = im1, im2
          do i = 1, NPTS
            do k = ict(i), icb(i)-1
              denm = ts(i,k,im) / (_one - rsa(i,k-1,ih,im)*rs(i,k,im))

              tda(i,k,ih,im) = tda(i,k-1,ih,im)*td(i,k,im)
              tta(i,k,ih,im) = tda(i,k-1,ih,im)*tt(i,k,im)              &
     &                + (tda(i,k-1,ih,im)*rsa(i,k-1,ih,im)*rr(i,k,im)   &
     &                + tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
              rsa(i,k,ih,im) = rs(i,k,im)                               &
     &                + ts(i,k,im) * rsa(i,k-1,ih,im) * denm
            enddo
          enddo
        enddo                 ! end im loop
      enddo                 ! end ih loop

!  ---  layers are added one at a time, going up from the surface.
!       rra is the composite reflectance illuminated by beam radiation
!       rxa is the composite reflectance illuminated from above by
!           diffuse radiation
!       rra and rxa are computed from eqs. (6.9) and (6.11)

!  ---  to save memory space, rra and rxa are pre-computed for k>=icb.
!       the dimension of these parameters is (m,np,2,2). It would have
!       been (m,np,2,2,2) if these parameters were computed for all k's.

!  ---  for the low clouds
!       is=1 for clear-sky condition, is=2 for cloudy-sky condition

      do is = is1, is2
        do i = 1, NPTS
          rra(i,NLP1,1,is) = rr(i,NLP1,is)
          rxa(i,NLP1,1,is) = rs(i,NLP1,is)
          rra(i,NLP1,2,is) = rr(i,NLP1,is)
          rxa(i,NLP1,2,is) = rs(i,NLP1,is)
        enddo

        do i = 1, NPTS
          do k = NLAY, icb(i), -1
            denm = ts(i,k,is) / (_one - rs(i,k,is)*rxa(i,k+1,1,is))

            rra(i,k,1,is) = rr(i,k,is) + (td(i,k,is)*rra(i,k+1,1,is)    &
     &                + (tt(i,k,is)-td(i,k,is))*rxa(i,k+1,1,is))*denm
            rxa(i,k,1,is) = rs(i,k,is) + ts(i,k,is)*rxa(i,k+1,1,is)*denm
            rra(i,k,2,is) = rra(i,k,1,is)
            rxa(i,k,2,is) = rxa(i,k,1,is)
          enddo
        enddo

!  ---  for middle clouds

        do im = im1, im2
          do i = 1, NPTS
          do k = icb(i)-1, ict(i), -1
            denm = ts(i,k,im) / (_one - rs(i,k,im)*rxa(i,k+1,im,is))

            rra(i,k,im,is) = rr(i,k,im) + (td(i,k,im)*rra(i,k+1,im,is)  &
     &                + (tt(i,k,im)-td(i,k,im))*rxa(i,k+1,im,is))*denm
            rxa(i,k,im,is) = rs(i,k,im)+ts(i,k,im)*rxa(i,k+1,im,is)*denm
          enddo
          enddo
        enddo                 ! end im loop
      enddo                 ! end is loop

!  ---  integration over eight sky situations.
!       ih, im, is denote high, middle and low cloud groups.

      do ih = ih1, ih2

        if ( ih == 1 ) then
!  ---  clear portion 
          do i = 1, NPTS
            ch(i) = _one - cc(i,1)
          enddo
        else
!  ---  cloudy portion
          do i = 1, NPTS
            ch(i) = cc(i,1)
          enddo
        endif

        do im = im1, im2

          if ( im == 1 ) then
!  ---  clear portion
            do i = 1, NPTS
              cm(i) = ch(i)*(_one - cc(i,2))
            enddo
          else
!  ---  cloudy portion
            do i = 1, NPTS
              cm(i) = ch(i)*cc(i,2) 
            enddo
          endif

          do is = is1, is2

            if ( is == 1 ) then
!  ---  clear portion
              do i = 1, NPTS
                ct(i) = cm(i)*(_one - cc(i,3)) 
              enddo
            else
!  ---  cloudy portion
              do i = 1, NPTS
                ct(i) = cm(i)*cc(i,3)
              enddo
            endif

!  ---  add one layer at a time, going down.

            do i = 1, NPTS
            do k = icb(i), NLAY
              denm = ts(i,k,is) / (_one - rsa(i,k-1,ih,im)*rs(i,k,is))

              tda(i,k,ih,im) = tda(i,k-1,ih,im)*td(i,k,is)
              tta(i,k,ih,im) = tda(i,k-1,ih,im)*tt(i,k,is)              &
     &                 + (tda(i,k-1,ih,im)*rr(i,k,is)*rsa(i,k-1,ih,im)  &
     &                 + tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
              rsa(i,k,ih,im) = rs(i,k,is)                               &
     &                 + ts(i,k,is)*rsa(i,k-1,ih,im)*denm
            enddo
            enddo

!  ---  add one layer at a time, going up.

            do i = 1, NPTS
            do k = ict(i)-1, 1, -1
              denm = ts(i,k,ih) / (_one - rs(i,k,ih)*rxa(i,k+1,im,is))

              rra(i,k,im,is) = rr(i,k,ih)+(td(i,k,ih)*rra(i,k+1,im,is)  &
     &                 + (tt(i,k,ih)-td(i,k,ih))*rxa(i,k+1,im,is))*denm
              rxa(i,k,im,is) = rs(i,k,ih)                               &
     &                 + ts(i,k,ih)*rxa(i,k+1,im,is)*denm
            enddo
            enddo

!  ---  compute fluxes following Eq. (6.15) for fupdif and
!       eq. (6.16) for (fdndir+fdndif)
!       fdndir is the direct  downward flux
!       fdndif is the diffuse downward flux
!       fupdif is the diffuse upward flux

            do k = 2, NLP1
            do i = 1, NPTS
              denm = _one / (_one - rsa(i,k-1,ih,im)*rxa(i,k,im,is))

              xx = tda(i,k-1,ih,im) * rra(i,k,im,is)
              yy = tta(i,k-1,ih,im) - tda(i,k-1,ih,im)

              fdndir(i) = tda(i,k-1,ih,im)
              fdndif(i) = (xx*rsa(i,k-1,ih,im) + yy) * denm
              fupdif    = (xx + yy*rxa(i,k,im,is)) * denm
              flxdn(i,k)= fdndir(i) + fdndif(i)
              flxup(i,k)= fupdif
            enddo
            enddo

            do i = 1, NPTS
              flxdn(i,1) = _one
              flxup(i,1) = rra(i,1,im,is)
            enddo

!  ---  summation of fluxes over all sky situations;
!       the term in the brackets of Eq. (7.11)

            do k = 1, NLP1
            do i = 1, NPTS
              if ( ih == 1 .and. im == 1 .and. is == 1 ) then
                fcdn(i,k) = flxdn(i,k)
                fcup(i,k) = flxup(i,k)
              endif
              fxdn(i,k) = fxdn(i,k) + flxdn(i,k)*ct(i)
              fxup(i,k) = fxup(i,k) + flxup(i,k)*ct(i)
            enddo
            enddo

            do i = 1, NPTS
              fsdir(i) = fsdir(i) + fdndir(i)*ct(i)
              fsdif(i) = fsdif(i) + fdndif(i)*ct(i)
            enddo
          enddo                 ! end is loop
        enddo                 ! end im loop
      enddo                 ! end ih loop

!
      return
!...................................
      end subroutine cldflx
!-----------------------------------


!-----------------------------------
      subroutine rflx                                                   &
!...................................

!  ---  inputs:
     &     ( swc, u1, du, nu, swh, w1, dw, nw, tbl,                     &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       df                                                         &
     &     )

!*****************************************************************
!     compute the reduction of clear-sky downward solar flux
!     due to co2 absorption.
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, nu, nw

      real (kind=kind_phys), intent(in) :: u1,du,w1,dw
      real (kind=kind_phys), dimension(:,:), intent(in) :: swc,swh,tbl

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: df

!  ---  locals:
      integer :: i, k, ic, iw 
      real (kind=kind_phys) :: clog, wlog, dc, dd, x0,x1,x2, y0,y1,y2
!
!===> ...  begin here
!
!  ---  table look-up for the reduction of clear-sky solar

      x0 = u1 + float(nu)*du
      y0 = w1 + float(nw)*dw

      x1 = u1 - 0.5*du
      y1 = w1 - 0.5*dw

      do k = 2, NLP1
      do i = 1, NPTS
        clog = min( x0, max( (x1-du),     swc(i,k) ))
        wlog = min( y0, max( (y1-3.0*dw), swh(i,k) ))

        ic = max( 2, min( nu, int((clog-x1)/du+_one) ))
        iw = max( 2, min( nw, int((wlog-y1)/dw+_one) ))

        dc = clog - float(ic-2)*du - u1
        dd = wlog - float(iw-2)*dw - w1   

        x2 = tbl(ic-1,iw-1) + (tbl(ic-1,iw)-tbl(ic-1,iw-1))/dw * dd
        y2 = x2 + (tbl(ic,iw-1)-tbl(ic-1,iw-1)) / du * dc

        df(i,k) = df(i,k) + y2
      enddo      
      enddo  
!
      return
!...................................
      end subroutine rflx
!-----------------------------------

!
!........................................!
      end module module_radsw_main       !
!========================================!

