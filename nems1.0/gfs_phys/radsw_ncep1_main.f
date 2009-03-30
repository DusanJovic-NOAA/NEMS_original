!!!!!  ==========================================================  !!!!!
!!!!!            sw-ncep1 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the sw-ncep1 package includes these parts:                         !
!                                                                      !
!      'radsw_ncep1_param.f'                                           !
!      'radsw_ncep1_datatb.f'                                          !
!      'radsw_ncep1_main.f'                                            !
!                                                                      !
!   the 'radsw_ncep1_param.f' contains:                                !
!                                                                      !
!      'module_radsw_parameters'  -- band parameters set up            !
!      'module_radsw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radsw_ncep1_datatb.f' contains:                               !
!                                                                      !
!      'module_radsw_co2tab'      -- co2 lookup tables                 !
!      'module_radsw_cldprtb'     -- cloud property coefficients table !
!                                                                      !
!   the 'radsw_ncep1_main.f' contains:                                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main ncep1 sw radiation routine                 !
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
!      'radsw_ncep1_param.f'                                           !
!      'radsw_ncep1_datatb.f'                                          !
!      'radsw_ncep1_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use sw modules     !
!                                                                      !
!                                                                      !
!                                                                      !
!   program description abstract:                                      !
!                                                                      !
!      this code is a modified version of nasa m.d. chou's sw          !
!      radiation code to fit ncep gfs and climate models.  it computes !
!      sw atmospheric absorption and scattering effects due to o3,     !
!      h2o, co2, o2, clouds, and aerosols, etc.                        !
!      it has 8 uv+vis bands and 3 (or 1) nir bands (10 k-values each) !
!                                                                      !
!                                                                      !
!   references:                                                        !
!                                                                      !
!      chou (1986, j. clim. appl.meteor.)                              !
!      chou (1990, j. clim.), and chou (1992, j. atms. sci.)           !
!      chou and suarez (1999, nasa/tm-1999-104606,vol.15)              !
!                                                                      !
!                                                                      !
!                                                                      !
!   ncep modifications history log:                                    !
!                                                                      !
!      jun 1994,  yu-tai hou  -- obtained original code from  nasa     !
!      feb 1995,  yu-tai hou                                           !
!                 recode for nmc models, change flux adding scheme to  !
!                 coakley et al.(1983) formula, recalculate rayleigh   !
!                 scatt coef for sigma structure.                      !
!      aug 1998,  yu-tai hou                                           !
!                 updated cloud radiative properties calculation. use  !
!                 slingo's method (jas 1989) on water cloud, ebert and !
!                 curry's method (jgr 1992) on ice cloud.              !
!      mar 1999,  yu-tai hou                                           !
!                 updated cloud radiative property calculations, use   !
!                 chou et al. new method (j. clim 1998)                !
!      apr 1999,  yu-tai hou                                           !
!                 updated cloud radiative property calculations, use   !
!                 linear t-adjusted method.                            !
!      sep 1999,  yu-tai hou                                           !
!                 updated to chou's june,99 version                    !
!      jul 2001,  yu-tai hou                                           !
!                 modified code to fortran 90 standard with unified    !
!                 in/out parameters                                    !
!      apr 2005,  yu-tai hou                                           !
!                 modification on module structures                    !
!      apr 2007,  yu-tai hou                                           !
!                 add spectral band heating as optional output         !
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
      use physcons,           only : con_rd, con_ttp, con_g,            &
     &                               con_cp, con_t0c, con_amo3,         &
     &                               con_amd, con_amco2, con_amo2,      &
     &                               con_gasv, con_p0
      use module_radsw_parameters
      use module_radsw_cntr_para
      use module_radsw_co2tab
!
      implicit   none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGSW='NCEP-SW v03.1   jun 2004'
      character(24), parameter :: VTAGSW='NCEP-SW v03.2   Apr 2007'

!  ---  constant values
      real (kind=kind_phys) :: fpmin, fpmax, ftiny, fsmall, fffrs, zero

      parameter (fpmin=1.0e-8,  fpmax=0.999999,  ftiny=1.0e-11)
      parameter (fsmall=1.0e-3, fffrs=0.1, zero=0.0)

!  ---  atomic weights for dry air and ozone (g/mol), and NTP pressure (Pa)
!       faco3 is unit convertion factor from kg/m**2 to m_atm at NTP
!     real (kind=kind_phys) :: amdair, amo3, p0ntp, faco3

!     parameter (amdair=28.9644,  amo3=48.0,  p0ntp=101325.0)
      real (kind=kind_phys), parameter ::                               &
     &       faco3 = con_rd*(con_amd/con_amo3)*con_ttp/con_p0

!  ---  band indices (ih2osw,ioznsw,ibndsw)
      integer, dimension(0:1,0:1,0:1) :: mstr, mend, nstr, nend
      integer       :: mbstr, mbend, nbstr, nbend

      data mstr / 2, 1, 2, 1, 3, 3, 2, 2 /
      data mend / 1, 1, 2, 2, 2, NBD0, 2, NBD0 /
      data nstr / NK0P,  1, NK0P, 1,   NBK0P, NBK0P,NK0P, NK0P /  ! band loop starting index
      data nend / NK0, NK0, NBK0, NBK0, NBK0, NBK1, NBK0, NBK1 /  ! band loop ending index

      integer, parameter :: nuvb1 = NK0 + 4       !starting uv-b band index
!     integer, parameter :: nuvb1 = NK0 + 5       !starting uv-b band index
      integer, parameter :: nuvb2 = NK0 + 6       !ending   uv-b band index

!  ---  spectral bands are arranged as: 1 ir band (10 ks), 7 uv bands, 1 vis band,
!       then follows 3 ir bands (29 ks)

!  ---  indicis look up arrays
      integer, dimension(NBK1) :: jbk, jcc, jmm, jab

!  ---  jbk : index array for k-values
      data jbk /  1, 2, 3, 4, 5, 6, 7, 8, 9,10                          &
     &,           11,12,13,14,15,16,17,18                               &
     &,           1, 2, 3, 4, 5, 6, 7, 8, 9                             &
     &,           1, 2, 3, 4, 5, 6, 7, 8, 9,10                          &
     &,           1, 2, 3, 4, 5, 6, 7, 8, 9,10 /
!  ---  jcc : index array for cloud optical coefficients
      data jcc /  1,1,1,1,1,1,1,1,1,1,  2,2,2,2,2,2,2,2                 &
     &,           3,3,3,3,3,3,3,3,3,    4,4,4,4,4,4,4,4,4,4             &
     &,           5,5,5,5,5,5,5,5,5,5 /
!  ---  jmm : index array for aerosol and Rayleigh scat coefficients
      data jmm /  1,1,1,1,1,1,1,1,1,1,  2,3,4,5,6,7,8,9                 &
     &,           10,10,10,10,10,10,10,10,10                            &
     &,           11,11,11,11,11,11,11,11,11,11                         &
     &,           12,12,12,12,12,12,12,12,12,12 /
!  ---  jab : index array for surface albedo: 1-nir; 2-uv,vis
      data jab /  1,1,1,1,1,1,1,1,1,1,  2,2,2,2,2,2,2,2                 &
     &,           1,1,1,1,1,1,1,1,1,    1,1,1,1,1,1,1,1,1,1             &
     &,           1,1,1,1,1,1,1,1,1,1 /

!  ---  k-value for water vapor
      real (kind=kind_phys), dimension(NBK0) :: wk

      data wk / 0.0010, 0.0133, 0.0422, 0.1334, 0.4217,                 &
     &          1.3340, 5.6230, 31.620, 177.80, 1000.0,                 &
     &          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .75e-3 /

!  ---  k-value for ozone, and rayleigh scattering tau
      real (kind=kind_phys), dimension(NBD1) :: ak, ry

      data ak / 0.0, 30.47, 187.24, 301.92, 42.83, 7.090, 1.250,        &
     &          .0345, .0572,   0.0,  0.0,  0.0 /
!chou99
      data ry / .180e-5, .604e-2, .170e-2, .222e-2, .132e-2, .107e-2,   &
     &          .910e-3, .550e-3, .120e-3, .156e-4, .180e-5, .0      /
!chou99
!mrf95
!     data ry / .128e-1,   7.006,   2.117,   2.453,   1.398,   1.133,   &
!    &           0.9532,  0.6104,  0.1096, .354e-1, .375e-2, .289e-3 /
!mrf95

!  ---  solar energy spectral weights
      real (kind=kind_phys), dimension(NBK1) :: ss

      data ss / .29983, .05014, .04555, .03824, .02965,                 &
     &          .02280, .02321, .01230, .00515, .00239,                 &
     &          .00057, .00367, .00083, .00417,                         &
     &          .00600, .00556, .05913, .39081,                         &
     &          .20673, .03497, .03011, .02260, .01336,                 &
     &          .00696, .00441, .00115, .00026,                         &
     &          .08236, .01157, .01133, .01143, .01240,                 &
     &          .01258, .01381, .00650, .00244, .00094,                 &
     &          .01074, .00360, .00411, .00421, .00389,                 &
     &          .00326, .00499, .00465, .00245, .00145 /

!! ...  logical flags for optional output fields

      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.

!  ---  those data will be set up only once by "rswinit"

!  ...  fheat is the factor for heating rates (in k/day, or k/sec set
!       by subroutine 'rswinit')

      real (kind=kind_phys) :: fheat

      public swrad, rswinit


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
! usage:       call swrad                                               !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   cray c-90, ibm sp, sgi                                   !
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
!   clouds(IMAX,NLAY,:): cloud profile                                  !
!                      (check module_raddiation_clouds for definition)  !
!                ---  for  iflagc > 0  ---                              !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud liq water path      (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer cloud ice water path      (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  iflagc = 0  ---                              !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
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
!   ih2osw : flags for h2o absorption                                   !
!            =0: without h2o absorption                                 !
!            =1: include h2o absorption                                 !
!   ioznsw : flags for o3 absorption                                    !
!            =0: without o3 absorption                                  !
!            =1: include o3 absorption                                  !
!   iswrate: flags for heating rate unit                                !
!            =1: output in k/day                                        !
!            =2: output in k/second                                     !
!   ibndsw : band selection in nir spectrum                             !
!            =0: use 1 nir band                                         !
!            =1: use 3 nir bands                                        !
!   iflagc : control flag for cloud optical properties                  !
!            =0; input cloud optical depth, fixed ssa, asy              !
!            =1: input cwp/cip, t-adj coeff for all bands               !
!            =2: input cwp/cip, chou (1999) coeff for uv                !
!                and 3 nir bands, or gfs's tuned 1 nir                  !
!                band coeff according to ibndsw's setting               !
!            =3: input cwp/cip, chou (2002) coeff for uv                !
!                and 3 nir bands, or gfs's tuned 1 nir                  !
!                band coeff according to ibndsw's setting               !
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
!  ====================    end of description    =====================  !
!
      implicit none


!  ---  constants for co2 look-up table
      real (kind=kind_phys) :: uu1, uu2, ww1, ww2, du1, du2, dw1, dw2

      parameter (uu1=-3.0,   ww1=-4.0,  du1=0.15,   dw1=0.15)
      parameter (uu2=.25e-3, ww2=-2.0,  du2=.05e-3, dw2=0.05)

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

!  ---  local arrays:
      real (kind=kind_phys), dimension(NPTS,NLP1) :: pl1, swu, swh, so2

      real (kind=kind_phys), dimension(NPTS,NLAY) :: ta1, wa1, oa1, dp1,&
     &       rh1, wh1, oh1, scal, dz1, cld1, cwp1, cip1, rew1, rei1,    &
     &       cda1, cda2, cda3, cda4, tauc1, ssat1, asyt1, ffft1, taut1, &
     &       o2vmr, co2vmr

      real (kind=kind_phys) :: s0, tem0d1, tem0d2, tem0d3

      real (kind=kind_phys), dimension(NPTS,NLAY,NBD0) :: taucw, ssacw, &
     &       asycw, fffcw
      real (kind=kind_phys), dimension(NPTS,NLAY,NSWSTR:NSWEND) ::      &
     &       tauae, ssaae, asyae, taurs
      real (kind=kind_phys), dimension(NPTS,NLAY,2)    :: tauto, ssato, &
     &       asyto, fffto

      real (kind=kind_phys), dimension(NPTS)      :: cosz1, scnt1,      &
     &       albb1, albd1, topfd, topfu0, topfuc, sfcfd0, sfcfdc,       &
     &       sfcfu0, sfcfuc, cf0, cf1, cnt, tem1d

      real (kind=kind_phys), dimension(NPTS,2)    :: albbm, albdf,      &
     &       dwsfxb, dwsfxd

      real (kind=kind_phys), dimension(NPTS,NLP1) :: flxdn0, flxdnc,    &
     &       dflx0, dflxc, tem2d

      real (kind=kind_phys), dimension(NPTS,NLP1,2)    :: upflux,dwflux

      real (kind=kind_phys), dimension(NPTS,NLP1,NBDSW):: fnet0, fnetc
      real (kind=kind_phys), dimension(NPTS,NLP1) :: fnet0a, fnetca

!! ---  for optional outputs:
      real (kind=kind_phys), dimension(NPTS,2)    :: sdnbm0, sdndf0,    &
     &       sdnbmc, sdndfc
      real (kind=kind_phys), dimension(NPTS,NLP1) :: flxup0, flxupc,    &
     &       dfo2, dfco2a, dfco2b
      real (kind=kind_phys), dimension(NPTS)      :: suvbf0, suvbfc

      integer, dimension(NPTS) :: kctop

      integer :: i, j1, k, k1
      integer :: ib, icc, ibk, iab, imm, mb

      logical :: lcloudy
!
!===> ... begin here
!

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )

!  ---  s0, the solar constant is scaled to a more current value
!       i.e. if solc=2.0 ly/min then ssolar=1.96 ly/min. then
!       convert unit to w/m**2

!         s0 = solfac * con_solr * 0.980e0 * 6.976670e2    ! old gfdl form
          s0 = solcon

!  ---  initial output arrays

      hswc(:,:) = zero
      topflx = topfsw_type ( zero, zero, zero )
      sfcflx = sfcfsw_type ( zero, zero, zero, zero )

!! ---  initial optional outputs
      if ( lflxprf ) then
        flxprf = profsw_type ( zero, zero, zero, zero )
      endif

      if ( lfdncmp ) then
        fdncmp = cmpfsw_type ( zero, zero, zero, zero, zero, zero )
      endif

      if ( lhsw0 ) then
        hsw0(:,:) = zero
      endif

      if ( lhswb ) then
        hswb(:,:,:) = zero
      endif

!  ---  setup internal arrays

      do i = 1, NPTS
        cf0   (i) = 1.0
        cf1   (i) = 1.0
        kctop (i) = NLP1
      enddo

      do i = 1, NPTS
        j1 = idxday(i)
        topfd(i) = s0 * cosz(j1)
        cosz1(i) = cosz(j1)
        scnt1(i) = 1.0 / cosz(j1)   ! scnt1 = secant of solar zenith angle

        albbm(i,1) = sfcalb(j1,1)   ! nir beam albedo
        albdf(i,1) = sfcalb(j1,2)   ! nir diff albedo
        albbm(i,2) = sfcalb(j1,3)   ! uvv beam albedo
        albdf(i,2) = sfcalb(j1,4)   ! uvv diff albedo
      enddo

!  ---  the internal array is always from top to surface

      if (iflip == 0) then        ! input from toa to sfc

        do i = 1, NPTS
          j1 = idxday(i)
          pl1(i,NLP1)  = plvl (j1,NLP1)

          do k = 1, NLAY
            pl1(i,k)  = plvl (j1,k)
            ta1(i,k)  = tlyr (j1,k)
!test use
!           wa1(i,k)  = max(zero, qlyr(j1,k))                   ! input mass mixing ratio
!ncep model use
            wa1(i,k)  = max(zero, qlyr(j1,k)/(1.0-qlyr(j1,k)))  ! input specific humidity
            oa1(i,k)  = max(zero, olyr(j1,k))                   ! input mass mixing ratio
            co2vmr(i,k)= max(zero, gasvmr(j1,k,1))              ! co2
            o2vmr(i,k) = max(zero, gasvmr(j1,k,4))              ! o2
          enddo
        enddo

        if (iaersw > 0) then
          do i = 1, NPTS
            j1 = idxday(i)

            do ib = 1, NBDSW
              mb = NSWSTR + ib - 1

              do k = 1, NLAY
                tauae(i,k,mb) = aerosols(j1,k,ib,1)
                ssaae(i,k,mb) = aerosols(j1,k,ib,2)
                asyae(i,k,mb) = aerosols(j1,k,ib,3)
              enddo
            enddo
          enddo
        else
          tauae(:,:,:) = zero
          ssaae(:,:,:) = zero
          asyae(:,:,:) = zero
        endif

        if (iflagc > 0) then      ! use prognostic cloud method
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              cld1(i,k) = clouds(j1,k,1)
              cwp1(i,k) = clouds(j1,k,2)
              rew1(i,k) = clouds(j1,k,3)
              cip1(i,k) = clouds(j1,k,4)
              rei1(i,k) = clouds(j1,k,5)
              cda1(i,k) = clouds(j1,k,6)
              cda2(i,k) = clouds(j1,k,7)
              cda3(i,k) = clouds(j1,k,8)
              cda4(i,k) = clouds(j1,k,9)
            enddo
          enddo
        else                      ! use diagnostic cloud method
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              cld1(i,k) = clouds(j1,k,1)
              cda1(i,k) = clouds(j1,k,2)
              cda2(i,k) = clouds(j1,k,3)
              cda3(i,k) = clouds(j1,k,4)
            enddo
          enddo
        endif                     ! end_if_iflagc

      else                        ! input data from sfc to toa

        do i = 1, NPTS
          j1 = idxday(i)
          pl1(i,1)  = plvl (j1,NLP1)

          do k = 1, NLAY
            k1 = NLP1 - k
            pl1(i,k+1)= plvl (j1,k1)
            ta1(i,k)  = tlyr (j1,k1)
!test use
!           wa1(i,k)  = max(zero, qlyr(j1,k1))                  ! input mass mixing ratio
!ncep model use
            wa1(i,k)  = max(zero, qlyr(j1,k1)/(1.0-qlyr(j1,k1)))! input specific humidity
            oa1(i,k)  = max(zero, olyr(j1,k1))                  ! input mass mixing ratio
            co2vmr(i,k)= max(zero, gasvmr(j1,k1,1))             ! co2
            o2vmr(i,k) = max(zero, gasvmr(j1,k1,4))             ! o2
          enddo
        enddo

        if (iaersw > 0) then
          do i = 1, NPTS
            j1 = idxday(i)

            do ib = 1, NBDSW
              mb = NSWSTR + ib - 1

              do k = 1, NLAY
                k1 = NLP1 - k

                tauae(i,k,mb) = aerosols(j1,k1,ib,1)
                ssaae(i,k,mb) = aerosols(j1,k1,ib,2)
                asyae(i,k,mb) = aerosols(j1,k1,ib,3)
              enddo
            enddo
          enddo
        else
          tauae(:,:,:) = zero
          ssaae(:,:,:) = zero
          asyae(:,:,:) = zero
        endif

        if (iflagc > 0) then      ! use prognostic cloud method
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              k1 = NLP1 - k
              cld1(i,k) = clouds(j1,k1,1)
              cwp1(i,k) = clouds(j1,k1,2)
              rew1(i,k) = clouds(j1,k1,3)
              cip1(i,k) = clouds(j1,k1,4)
              rei1(i,k) = clouds(j1,k1,5)
              cda1(i,k) = clouds(j1,k1,6)
              cda2(i,k) = clouds(j1,k1,7)
              cda3(i,k) = clouds(j1,k1,8)
              cda4(i,k) = clouds(j1,k1,9)
            enddo
          enddo
        else                      ! use diagnostic cloud method
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              j1 = idxday(i)
              cld1(i,k) = clouds(j1,k1,1)
              cda1(i,k) = clouds(j1,k1,2)
              cda2(i,k) = clouds(j1,k1,3)
              cda3(i,k) = clouds(j1,k1,4)
            enddo
          enddo
        endif                     ! end_if_iflagc

      endif                       ! if_iflip

      do k = 1, NLAY
        do i = 1, NPTS
          dp1(i,k) = pl1(i,k+1) - pl1(i,k)
        enddo
      enddo

!  ---  find top layer of cloud

      do i = 1, NPTS
        lab_do_k_cld1 : do k = 1, NLAY
          if (cld1(i,k) > zero) then         ! cloudy layer
            kctop(i) = k
            exit lab_do_k_cld1
          endif
        enddo  lab_do_k_cld1
      enddo

!  ---  compute fractions of clear sky view
!       (cf0 and cf1 are initialized to 1.0)

      if (iovr == 0) then                 ! random overlap
        do k = 1, NLAY
          do i = 1, NPTS
            cf0(i) = cf0(i) * (1.0 - cld1(i,k))
          enddo
        enddo
      else                                ! max/ran overlap
        do k = 1, NLAY
          do i = 1, NPTS
            if (cld1(i,k) > zero) then       ! cloudy layer
              cf1(i) = min ( cf1(i), 1.0-cld1(i,k) )
            elseif (cf1(i) < 1.0) then
              cf0(i) = cf0(i) * cf1(i)
              cf1(i) = 1.0
            endif
          enddo
        enddo

        do i = 1, NPTS
          cf0(i) = cf0(i) * cf1(i)
        enddo
      endif

      tem0d1 = 1.0 - fsmall
      lcloudy = .false.
      do i = 1, NPTS
        if (cf0(i) < fsmall) cf0(i) = zero
        if (cf0(i) > tem0d1) cf0(i) = 1.0
        cf1(i) = 1.0 - cf0(i)
        lcloudy = lcloudy .or. cf0(i)<1.0
      enddo
!
!===> ...  compute cloud optical properties
!
      if ( lcloudy ) then

        do i = 1, NPTS
          if (cf1(i) > zero) then
            tem0d1 = 1.0 / cf1(i)
            do k = 1, NLAY
              cld1(i,k) = cld1(i,k) * tem0d1
            enddo
          endif
        enddo

        call cldprsw                                                    &
!  ---  inputs:
     &     ( ta1,cld1,cf1,cwp1,cip1,rew1,rei1,cda1,cda2,cda3,cda4,      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       taucw,ssacw,asycw,fffcw                                    &
     &     )

      endif       ! end_if_lcloudy

!
!===> ... compute rayleigh optical thickness
!

!chou99
      do j1 = 1, NBDSW
        mb = NSWSTR + j1 - 1
        do k = 1, NLAY
          do i = 1, NPTS
            taurs(i,k,mb) = ry(mb) * dp1(i,k)
          enddo
        enddo
      enddo
!chou99
!mrf95
!     do j1 = 1, NBDSW
!       mb = NSWSTR + j1 - 1
!       do k = 1, NLAY
!         do i = 1, NPTS
!           taurs(i,k,mb) = ry(mb) * dp1(i,k) / pl1(i,NLP1)
!         enddo
!       enddo
!     enddo
!mrf95
!
!===> ...  pressure scaling function for water vapor absorption, and
!          scaled layer water vapor and ozone amount, factor 466.7 is
!          the unit conversion factor from g/cm**2 to (cm-atm)stp, here
!          faco3 and con_g are in mks units, and factor 1.e4 change to
!          cgs units.
      tem0d1 = 0.5 / 300.0
!     tem0d2 = 1.02e0 * 466.7e0
      tem0d2 = 1.0e4 * faco3 / con_g

      do k = 1, NLAY
        do i = 1, NPTS
          scal(i,k) = dp1(i,k) * (tem0d1 * (pl1(i,k)+pl1(i,k+1)))**0.8

!  ---  scaled absorber amounts for h2o(wh,swh), unit : g/cm**2
          tem0d3 = 0.00135 * (ta1(i,k) - 240.0)
          wh1(i,k) = 1.02 * wa1(i,k) * scal(i,k)                        &
!    &             * exp( 0.00135*(ta1(i,k) - 240.0))
     &             * (1.0 + tem0d3 + 0.5*tem0d3*tem0d3) + ftiny

!  ---  scaled absorber amounts for ozoen,       unit : cm-atm
          oh1(i,k) = tem0d2 * oa1(i,k) * dp1(i,k) + ftiny
        enddo
      enddo

!
!===> ...  initialize local arrays for radiation flux calculations
!

      do i = 1, NPTS
        topfu0(i)  = zero
        sfcfu0(i)  = zero
        sfcfd0(i)  = zero
        dflx0(i,1) = zero
        topfuc(i)  = zero
        sfcfuc(i)  = zero
        sfcfdc(i)  = zero
        dflxc(i,1) = zero
      enddo

      do mb = 1, NBDSW
        do k = 1, NLP1
        do i = 1, NPTS
          fnet0(i,k,mb) = zero
          fnetc(i,k,mb) = zero
        enddo
        enddo
      enddo

      do k = 1, NLP1
        do i = 1, NPTS
          fnet0a(i,k) = zero
          fnetca(i,k) = zero
          flxdn0(i,k) = zero
          flxdnc(i,k) = zero
        enddo
      enddo

!! ---  optional fields

      if ( lfdncmp ) then
        do i = 1, NPTS
!! ---  optional uv-b surface downward flux
          suvbfc(i)  = zero
          suvbf0(i)  = zero

!! ---  optional beam and difussed fluxes
          sdnbm0(i,1) = zero
          sdnbm0(i,2) = zero
          sdndf0(i,1) = zero
          sdndf0(i,2) = zero
          sdnbmc(i,1) = zero
          sdnbmc(i,2) = zero
          sdndfc(i,1) = zero
          sdndfc(i,2) = zero
        enddo
      endif

      if ( lflxprf ) then
        do k = 1, NLP1
          do i = 1, NPTS
            flxup0(i,k) = zero
            flxupc(i,k) = zero
          enddo
        enddo
      endif

!
!===> ... loop over all k's to compute fluxes
!

      lab_do_ib : do ib = nbstr, nbend

        ibk = jbk(ib)                ! index for k-values
        icc = jcc(ib)                ! index for cloud properties
        imm = jmm(ib)                ! index for aerosol and Rayleigh scat
        iab = jab(ib)                ! index for surface albedo

!  ---  compute clear sky optical thickness, single scattering albedo,
!       asymmetry factor, and rayleigh optical thickness

        do k = 1, NLAY
          do i = 1, NPTS
            tem0d1     = ssaae(i,k,imm) * tauae(i,k,imm)

            taut1(i,k) = wk(ibk)*wh1(i,k) + ak(imm)*oh1(i,k)            &
     &                 + tauae(i,k,imm) + taurs(i,k,imm)
            ssat1(i,k) = tem0d1 + taurs(i,k,imm)
            asyt1(i,k) = tem0d1 * asyae(i,k,imm)
            ffft1(i,k) = asyt1(i,k)*asyae(i,k,imm)+fffrs*taurs(i,k,imm)

            tem2d(i,k) = 1.0 / max(fpmin, ssat1(i,k))
          enddo
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
            tauto(i,k,1) = max(fpmin, taut1(i,k))
            ssato(i,k,1) = min(fpmax, ssat1(i,k) / tauto(i,k,1))
            asyto(i,k,1) = asyt1(i,k) * tem2d(i,k)
            fffto(i,k,1) = ffft1(i,k) * tem2d(i,k)
          enddo
        enddo

!  ---  compute cloudy sky optical thickness, single scattering albedo,
!       asymmetry factor, and rayleigh optical thickness

        if ( lcloudy ) then

          do k = 1, NLAY
            do i = 1, NPTS
              if (cld1(i,k) > zero) then
                tem0d1 = ssat1(i,k) + ssacw(i,k,icc)
                tem0d2 = 1.0 / max(fpmin, tem0d1)
                tem0d3 = max(fpmin, taut1(i,k)+taucw(i,k,icc))

                tauto(i,k,2) = tem0d3
                ssato(i,k,2) = min(fpmax, tem0d1 / tem0d3)
                asyto(i,k,2) = (asyt1(i,k)+asycw(i,k,icc)) * tem0d2
                fffto(i,k,2) = (ffft1(i,k)+fffcw(i,k,icc)) * tem0d2
              else
                tauto(i,k,2) = tauto(i,k,1)
                ssato(i,k,2) = ssato(i,k,1)
                asyto(i,k,2) = asyto(i,k,1)
                fffto(i,k,2) = fffto(i,k,1)
              endif
            enddo
          enddo

!  ---  delta scaling layer optical parameters

          do k = 1, NLAY
            do i = 1, NPTS
              if (cld1(i,k) > zero) then
                tem0d1 = 1.0 - fffto(i,k,2)*ssato(i,k,2)
                tem0d2 = 1.0 - fffto(i,k,2)

                tauto(i,k,2) = tauto(i,k,2) * tem0d1
                ssato(i,k,2) = ssato(i,k,2) * tem0d2 / tem0d1
                asyto(i,k,2) = (asyto(i,k,2)-fffto(i,k,2)) / tem0d2
              endif
            enddo
          enddo

        else

          do k = 1, NLAY
            do i = 1, NPTS
              tauto(i,k,2) = tauto(i,k,1)
              ssato(i,k,2) = ssato(i,k,1)
              asyto(i,k,2) = asyto(i,k,1)
              fffto(i,k,2) = fffto(i,k,1)
            enddo
          enddo

        endif       ! end_if_lcloudy

!  --  set up surface albedo

        do i = 1, NPTS
          albb1(i) = albbm(i,iab)
          albd1(i) = albdf(i,iab)
        enddo
!
!===> ... compute radiation fluxes for clear sky and total sky
!
        call swflux                                                     &
!  ---  inputs:
     &     ( tauto,ssato,asyto,scnt1,cosz1,albb1,albd1,cld1,lcloudy,    &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       upflux,dwflux,dwsfxb,dwsfxd                                &
     &     )

!  ---  clear sky net fluxes

        mb = imm - NSWSTR + 1
        do k = 1, NLP1
          do i = 1, NPTS
            fnet0(i,k,mb) = fnet0(i,k,mb)                               &
     &                    + (dwflux(i,k,1) - upflux(i,k,1)) * ss(ib)
            flxdn0(i,k) = flxdn0(i,k) + dwflux(i,k,1)*ss(ib)
          enddo
        enddo

!  ---  toa and sfc fluxes
        do i = 1, NPTS
          topfu0(i) = topfu0(i) + upflux(i,1,1)   *ss(ib)
          sfcfu0(i) = sfcfu0(i) + upflux(i,NLP1,1)*ss(ib)
          sfcfd0(i) = sfcfd0(i) + dwflux(i,NLP1,1)*ss(ib)
        enddo

!! ---  optional flux profile
        if ( lflxprf ) then
          do k = 1, NLP1
          do i = 1, NPTS
            flxup0(i,k) = flxup0(i,k) + upflux(i,k,1)*ss(ib)
          enddo
          enddo
        endif    ! end if_lflxprf

        if ( lfdncmp ) then
!! ---  optional uv-b surface downward flux
          if (ib >= nuvb1 .and. ib <= nuvb2) then
            do i = 1, NPTS
              suvbf0(i) = suvbf0(i) + dwflux(i,NLP1,1)*ss(ib)
            enddo
          endif

!! ---  optional beam and diffused sfc fluxes
          do i = 1, NPTS
            sdnbm0(i,iab) = sdnbm0(i,iab) + dwsfxb(i,1)*ss(ib)
            sdndf0(i,iab) = sdndf0(i,iab) + dwsfxd(i,1)*ss(ib)
          enddo
        endif    ! end if_lfdncmp

!  ---  total sky net fluxes

        if (lcloudy) then

          do k = 1, NLP1
            do i = 1, NPTS
              fnetc(i,k,mb) = fnetc(i,k,mb)                             &
     &                      + (dwflux(i,k,2) - upflux(i,k,2)) * ss(ib)
              flxdnc(i,k) = flxdnc(i,k) + dwflux(i,k,2)*ss(ib)
            enddo
          enddo

!  ---  toa and sfc fluxes
          do i = 1, NPTS
            topfuc(i) = topfuc(i) + upflux(i,1,2)   *ss(ib)
            sfcfuc(i) = sfcfuc(i) + upflux(i,NLP1,2)*ss(ib)
            sfcfdc(i) = sfcfdc(i) + dwflux(i,NLP1,2)*ss(ib)
          enddo

!! ---  optional flux profile
          if ( lflxprf ) then
            do k = 1, NLP1
            do i = 1, NPTS
              flxupc(i,k) = flxupc(i,k) + upflux(i,k,2)*ss(ib)
            enddo
            enddo
          endif     ! end if_lflxprf

          if ( lfdncmp ) then
!! ---  optional uv-b surface downward flux
            if (ib >= nuvb1 .and. ib <= nuvb2) then
              do i = 1, NPTS
                suvbfc(i) = suvbfc(i) + dwflux(i,NLP1,2)*ss(ib)
              enddo
            endif

!! ---  optional beam and diffused sfc fluxes
            do i = 1, NPTS
              sdnbmc(i,iab) = sdnbmc(i,iab) + dwsfxb(i,2)*ss(ib)
              sdndfc(i,iab) = sdndfc(i,iab) + dwsfxd(i,2)*ss(ib)
            enddo
          endif     ! end if_lfdncmp

        endif       ! end_if_lcloudy
 
      enddo lab_do_ib

!  ---  if no cloudy point, copy clear sky fluxes to total sky fluxes

      if (.not. lcloudy ) then
        fnetc(:,:,:) = fnet0(:,:,:)
        flxdnc(:,:) = flxdn0(:,:)

!! ---  optional
        if ( lflxprf ) then
          flxupc(:,:) = flxup0(:,:)
        endif

!  ---  toa and sfc fluxes
        do i = 1, NPTS
          topfuc(i) = topfu0(i)
          sfcfuc(i) = sfcfu0(i)
          sfcfdc(i) = sfcfd0(i)
        enddo

        if ( lfdncmp ) then
!! ---  optional uv-b surface downward flux
          do i = 1, NPTS
            suvbfc(i) = suvbf0(i)

!! ---  optional beam and diffused sfc fluxes
            sdnbmc(i,iab) = sdnbm0(i,iab)
            sdndfc(i,iab) = sdndf0(i,iab)
          enddo
        endif    ! end if_lfdncmp
      endif    ! end if_not_lcloudy

!
!===> ... compute the absorption due to oxygen,chou(1990,j.climate,209-217)
!         pressure scaled amounts for o2(o2,so2), unit is (cm-atm)stp for o2.
!         the constant 165.22=(1000/980)*23.14%*(22400/32). note: for o2 in vmr
!         the constant will becomes (1000/980)*o2vmr*(22400/28.97)
 
      dfo2  (:,:) = zero
      dfco2a(:,:) = zero
      dfco2b(:,:) = zero

      lab_if_o2 : if (ioxysw == 1) then

        tem0d1 = (10.0 / con_g) * (con_gasv*1.0e6 / con_amd)
        do i = 1, NPTS
          so2(i,1) = zero
          cnt(i) = tem0d1 * scnt1(i)
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
            so2(i,k+1) = so2(i,k) + cnt(i)*o2vmr(i,k)*scal(i,k)
          enddo
        enddo
!
!  ---  compute flux reduction due to oxygen, the constant 0.0633 is
!       the fraction of insolation contained in the oxygen bands.
!       here tem0d1 is the broadband transmission function for oxygen.

        do k = 2, NLP1
          do i = 1, NPTS
            tem0d1    = exp(-0.145e-3 * sqrt(so2(i,k)) )
            dfo2(i,k) = 0.0633 * (1.0 - tem0d1)
          enddo
        enddo

      endif lab_if_o2

!
!===> ... table look-up for the absorption due to co2
!
      lab_if_co2 : if (ico2sw == 1) then

!
!  ---  compute scaled amounts for co2(wc,so2).
!       the constant 789=(1000/980)*(44/28.97)*(22400/44)
!       the factor (44/28.97 is to convert co2 vmr to mmr)

        tem0d1 = (10.0 / con_g) * (con_gasv*1.0e6 / con_amd)
        do i = 1, NPTS
!old      cnt(i) = 789.0  * scnt1(i)
          cnt(i) = tem0d1 * scnt1(i)
          so2(i,1) = ftiny
          swh(i,1) = zero
          swu(i,1) = zero
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
            so2(i,k+1) = so2(i,k) + cnt(i)*co2vmr(i,k)*scal(i,k) ! in cm-atm-stp
            swh(i,k+1) = swh(i,k) + wh1(i,k)
          enddo
        enddo
!
!  ---  for co2 absorption in spectrum 1.220-2.270 micron
!       both water vapor and co2 absorptions are moderate
!        (so2 and swh are the co2 and water vapor amounts integrated
!         from the top of the atmosphere)

        do k = 2, NLP1
          do i = 1, NPTS
            swu(i,k) = log10( so2(i,k) )
            swh(i,k) = log10( swh(i,k)*scnt1(i) )
          enddo
        enddo

        call flxco2                                                     &
!  ---  inputs:
     &     ( swu,uu1,du1,NUCO2,swh,ww1,dw1,NWCO2,cah,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outs:
     &       dfco2a                                                     &
     &     )
!
!  ---  for co2 absorption in spectrum 2.270-10.00 micron where co2
!       absorption has a large impact on the heating of middle atmos
!        (co2 mixing ratio is independent of space and swh is the 
!         logarithm of pressure)

        do k = 2, NLP1
          do i = 1, NPTS
            swu(i,k) = co2vmr(i,k-1) * scnt1(i)
            swh(i,k) = log10( pl1(i,k) )
          enddo
        enddo

        call flxco2                                                     &
!  ---  inputs:
     &     ( swu,uu2,du2,NXCO2,swh,ww2,dw2,NYCO2,coa,                   &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outs:
     &       dfco2b                                                     &
     &     )

      endif lab_if_co2

!
!===> ... adjust for the effect of o2 and co2 on net fluxes
!         and convert them to w/m**2

      lab_if_o2co2 : if (ioxysw == 1 .or. ico2sw == 1) then

        do k = 1, NLP1
          do i = 1, NPTS
            dflx0(i,k) = dfo2(i,k) + dfco2a(i,k) + dfco2b(i,k)
            dflxc(i,k) = dflx0(i,k)
            tem2d(i,k) = 1.0
          enddo
        enddo

        do i = 1, NPTS
          j1 =  kctop(i) + 1
          do k = j1, NLP1
            tem2d(i,k) = flxdnc(i,k) / flxdn0(i,k)
            dflxc(i,k) = dflx0(i,k) * tem2d(i,k)
          enddo
        enddo

        do i = 1, NPTS
          tem1d (i) = 0.5 * (albbm(i,1) + albbm(i,2))
          sfcfd0(i) = sfcfd0(i) - dflx0(i,NLP1)
          sfcfdc(i) = sfcfdc(i) - dflxc(i,NLP1)
        enddo

!! ---  optional fluxe profile
        if ( lflxprf ) then
          do k = 1, NLP1
            do i = 1, NPTS
              flxdn0(i,k) = flxdn0(i,k) - dflx0(i,k)
              flxdnc(i,k) = flxdnc(i,k) - dflxc(i,k)
              flxup0(i,k) = flxup0(i,k) - dflx0(i,k)*tem1d(i)
              flxupc(i,k) = flxupc(i,k) - dflxc(i,k)*tem1d(i)
            enddo
          enddo
        endif    ! end if_lflxprf

!! ---  optional surface component fluxes
        if ( lfdncmp ) then
          do i = 1, NPTS
            sdnbm0(i,1) = sdnbm0(i,1) - dflx0(i,NLP1)
            sdnbmc(i,1) = sdnbmc(i,1) - dflxc(i,NLP1)
          enddo
        endif

        if ( NSWSTR == 1 ) then       ! use 1 nir band

          tem1d(:) = 1.0 - tem1d(:)
          do k = 1, NLP1
            do i = 1, NPTS
              fnet0(i,k,1) = fnet0(i,k,1) - dflx0(i,k)*tem1d(i)
              fnetc(i,k,1) = fnetc(i,k,1) - dflxc(i,k)*tem1d(i)
            enddo
          enddo

        else                          ! use 3 nir bands

!  ---  o2 absorption

          dflx0(:,:) = dfo2(:,:) / 3.0
          dflxc(:,:) = dflx0(:,:) * tem2d(:,:)
          tem1d(:) = 1.0 - tem1d(:)

          do mb = 8, 10
            do k = 1, NLP1
              do i = 1, NPTS
                fnet0(i,k,mb) = fnet0(i,k,mb) - dflx0(i,k)*tem1d(i)
                fnetc(i,k,mb) = fnetc(i,k,mb) - dflxc(i,k)*tem1d(i)
              enddo
            enddo
          enddo

!  ---  co2 absorption for 1.22-2.27 micron band

          dflxc(:,:) = dfco2a(:,:) * tem2d(:,:)
          do k = 1, NLP1
            do i = 1, NPTS
              fnet0(i,k,10) = fnet0(i,k,10) - dfco2a(i,k)*tem1d(i)
              fnetc(i,k,10) = fnetc(i,k,10) - dflxc (i,k)*tem1d(i)
            enddo
          enddo

!  ---  co2 absorption for 2.27-10.0 micron band

          dflxc(:,:) = dfco2b(:,:) * tem2d(:,:)
          do i = 1, NPTS
            do k = 1, NLP1
              fnet0(i,k,11) = fnet0(i,k,11) - dfco2b(i,k)*tem1d(i)
              fnetc(i,k,11) = fnetc(i,k,11) - dflxc (i,k)*tem1d(i)
            enddo
          enddo

        endif    ! end if_NSWSTR

      endif lab_if_o2co2

!
!===> ... compute final fluxes and convert to w/m**2
!

      do k = 1, NLP1
        do i = 1, NPTS
          do mb = 1, NBDSW
            fnet0a(i,k) = fnet0a(i,k) + fnet0(i,k,mb)
            fnetca(i,k) = fnetca(i,k) + fnetc(i,k,mb)
          enddo
        enddo
      enddo

      do k = 1, NLP1
        do i = 1, NPTS
          fnetca(i,k) = (cf0(i)*fnet0a(i,k)+cf1(i)*fnetca(i,k))*topfd(i)
        enddo
      enddo

!  ---  toa and sfc fluxes
      do i = 1, NPTS
        topfuc(i) = (cf0(i)*topfu0(i) + cf1(i)*topfuc(i)) * topfd(i)
        sfcfuc(i) = (cf0(i)*sfcfu0(i) + cf1(i)*sfcfuc(i)) * topfd(i)
        sfcfdc(i) = (cf0(i)*sfcfd0(i) + cf1(i)*sfcfdc(i)) * topfd(i)
      enddo

      do i = 1, NPTS
        topfu0(i) = topfu0(i) * topfd(i)
        sfcfu0(i) = sfcfu0(i) * topfd(i)
        sfcfd0(i) = sfcfd0(i) * topfd(i)
      enddo

!
!===> ... fheat is the factor for heating rates (in k/day, or k/sec)
!

      if (iflip == 0) then        ! input data from toa to sfc
        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            hswc(j1,k) = (fnetca(i,k)-fnetca(i,k+1))*fheat / dp1(i,k)
          enddo
        enddo
      else                        ! input data from sfc to toa
        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            k1 = NLP1 - k
            hswc(j1,k1) = (fnetca(i,k)-fnetca(i,k+1))*fheat / dp1(i,k)
          enddo
        enddo
      endif                       ! if_flip

!  ---  toa and sfc fluxes

      do i = 1, NPTS
        j1 = idxday(i)
        topflx(j1)%upfxc = topfuc(i)
        topflx(j1)%dnfxc = topfd (i)
        topflx(j1)%upfx0 = topfu0(i)
        sfcflx(j1)%upfxc = sfcfuc(i)
        sfcflx(j1)%dnfxc = sfcfdc(i)
        sfcflx(j1)%upfx0 = sfcfu0(i)
        sfcflx(j1)%dnfx0 = sfcfd0(i)
      enddo

!! ---  optional heating rates

      if ( lhswb ) then
        do mb = 1, NBDSW
          do k = 1, NLP1
          do i = 1, NPTS
            fnetc(i,k,mb) = (cf0(i)*fnet0(i,k,mb)                       &
     &                    +  cf1(i)*fnetc(i,k,mb)) * topfd(i)
          enddo
          enddo
        enddo

        if (iflip == 0) then        ! input data from toa to sfc
          do i = 1, NPTS
            j1 = idxday(i)
            do mb = 1, NBDSW
            do k = 1, NLAY
              hswb(j1,k,mb) = (fnetc(i,k,mb)-fnetc(i,k+1,mb))*fheat     &
     &                      / dp1(i,k)
            enddo
            enddo
          enddo
        else                        ! input data from sfc to toa
          do i = 1, NPTS
            j1 = idxday(i)
            do mb = 1, NBDSW
            do k = 1, NLAY
              k1 = NLP1 - k
              hswb(j1,k1,mb) = (fnetc(i,k,mb)-fnetc(i,k+1,mb))*fheat    &
     &                       / dp1(i,k)
            enddo
            enddo
          enddo
        endif                       ! if_flip
      endif    ! end if_lhswb

      if ( lhsw0 ) then
        do k = 1, NLP1
        do i = 1, NPTS
          fnet0a(i,k) = fnet0a(i,k) * topfd(i)
        end do
        end do

        if (iflip == 0) then        ! input data from toa to sfc
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              hsw0(j1,k) = (fnet0a(i,k)-fnet0a(i,k+1))*fheat/dp1(i,k)
            enddo
          enddo
        else                        ! input data from sfc to toa
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              k1 = NLP1 - k
              hsw0(j1,k1) = (fnet0a(i,k)-fnet0a(i,k+1))*fheat/dp1(i,k)
            enddo
          enddo
        endif                       ! end if_iflip
      endif    ! end if_lhsw0

!! ---  optional fluxes

      if ( lflxprf ) then
        do k = 1, NLP1
        do i = 1, NPTS
          flxdnc(i,k) = cf0(i)*flxdn0(i,k) + cf1(i)*flxdnc(i,k)
          flxupc(i,k) = cf0(i)*flxup0(i,k) + cf1(i)*flxupc(i,k)
        enddo
        enddo

        if (iflip == 0) then        ! input data from toa to sfc
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              flxprf(j1,k)%upfx0 = flxup0(i,k) * topfd(i)
              flxprf(j1,k)%dnfx0 = flxdn0(i,k) * topfd(i)
              flxprf(j1,k)%upfxc = flxupc(i,k) * topfd(i)
              flxprf(j1,k)%dnfxc = flxdnc(i,k) * topfd(i)
            enddo
          enddo
        else                        ! input data from sfc to toa
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              k1 = NLP1 - k + 1
              flxprf(j1,k1)%upfx0 = flxup0(i,k) * topfd(i)
              flxprf(j1,k1)%dnfx0 = flxdn0(i,k) * topfd(i)
              flxprf(j1,k1)%upfxc = flxupc(i,k) * topfd(i)
              flxprf(j1,k1)%dnfxc = flxdnc(i,k) * topfd(i)
            enddo
          enddo
        endif                       ! if_iflip
      endif    ! end if_lflxprf

!! ---  optional surface downward component fluxes

      if ( lfdncmp ) then
        do i = 1, NPTS
          sdnbmc(i,1) = cf0(i)*sdnbm0(i,1) + cf1(i)*sdnbmc(i,1)
          sdnbmc(i,2) = cf0(i)*sdnbm0(i,2) + cf1(i)*sdnbmc(i,2)
          sdndfc(i,1) = cf0(i)*sdndf0(i,1) + cf1(i)*sdndfc(i,1)
          sdndfc(i,2) = cf0(i)*sdndf0(i,2) + cf1(i)*sdndfc(i,2)
        enddo

        do i = 1, NPTS
          j1 = idxday(i)

!! ---  optional uv-b surface downward flux
          fdncmp(j1)%uvbfc = suvbfc(i) * topfd(i)
          fdncmp(j1)%uvbf0 = suvbf0(i) * topfd(i)

!! ---  optional beam and diffused sfc fluxes
          fdncmp(j1)%nirbm = sdnbmc(i,1) * topfd(i)
          fdncmp(j1)%nirdf = sdndfc(i,1) * topfd(i)
          fdncmp(j1)%visbm = sdnbmc(i,2) * topfd(i)
          fdncmp(j1)%visdf = sdndfc(i,2) * topfd(i)
        enddo
      endif    ! end if_lfdncmp

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
!     ih2osw  - h2o absorption control flag                            !
!               =0: without h2o absorption                             !
!               =1: include h2o absorption                             !
!     ioznsw  - ozone absorption control flag                          !
!               =0: without ozone absorption                           !
!               =1: include ozone absorption                           !
!     ibndsw  - band selections in nir spectrum                        !
!               =0: use 1 nir band                                     !
!               =1: use 3 nir bands                                    !
!     iflagc  - flag for cloud property method                         !
!               =0: input cloud opt depth, fixed ssa, asy              !
!               =1: input cwp/cip, t-adj coeff for all bands           !
!               =2: input cwp/cip, chou (1999) coeff for uv and 3 nir  !
!                   bands, or gfs's tuned 1 nir band coefficents       !
!               =3: input cwp/cip, chou (2002) coeff for uv and 3 nir  !
!                   bands, or averaged 1 nir band coefficents          !
!                                                                      !
!  ==================================================================  !
!

      implicit none

!  ---  inputs: 
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: (none)

!  ---  locals: (none)

!
!===>  ... begin here
!

      if (me == 0) then
        print *,' - Using NCEP Shortwave Radiation, Version: ',VTAGSW

        if (iaersw == 0) then
          print *,'   --- Aerosol effect is NOT included in SW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        else
          print *,'   --- Using input aerosol parameters for SW'
        endif

        if (ioznsw == 0) then
          print *,'   --- O3 absorption is NOT included in SW'
        else
          print *,'   --- Include O3 absorption in SW'
        endif

        if (ih2osw == 0) then
          print *,'   --- H2O absorption is NOT included in SW'
        else
          print *,'   --- Include W.V. absorption in SW'
        endif

        if (ioxysw == 0) then
          print *,'   --- O2 absorption is NOT included in SW'
        else
          print *,'   --- Include O2 absorption in SW'
        endif

        if (ico2sw == 0) then
          print *,'   --- CO2 absorption is NOT included in SW'
        else
          print *,'   --- Include CO2 absorption in SW'
        endif

        if (ibndsw == 0) then
          print *,'   --- use 1 nir band, for faster computation'
        else
          print *,'   --- use 3 nir bands, for better accuracy'
        endif
      endif

!  --- ...  check cloud flags for consistency

      if ((icwp == 0 .and. iflagc /= 0) .or.                            &
     &    (icwp == 1 .and. iflagc == 0)) then
        print *, ' *** Model cloud scheme inconsistent with SW',        &
     &           ' radiation cloud radiative property setup !!'
        stop
      endif

!  --- ... setup starting/ending band parameters
!

      if (ih2osw < 0 .or. ih2osw > 1 .or.                               &
     &    ioznsw < 0 .or. ioznsw > 1 .or.                               &
     &    ibndsw < 0 .or. ibndsw > 1) then
        write(6,11) ih2osw,ioznsw,ibndsw
  11    format(/' *** Error setting: ih2osw,ioznsw,ibndsw =',3i4)
        stop 11
      endif

      mbstr = mstr(ih2osw,ioznsw,ibndsw)
      mbend = mend(ih2osw,ioznsw,ibndsw)
      nbstr = nstr(ih2osw,ioznsw,ibndsw)
      nbend = nend(ih2osw,ioznsw,ibndsw)


!  --- ... fheat is the factor for heating rates
!          the 1.0e-2 is to convert pressure from mb to N/m**2

      if (iswrate == 1) then
!       fheat = con_g * 86400. * 1.e-2 / con_cp    !    (in k/day)
        fheat = con_g * 864.0 / con_cp             !    (in k/day)
!       fheat = 8.4410328                          !    (in k/day)
      else
        fheat = con_g * 1.0e-2 / con_cp            !    (in k/second)
!       fheat = 8.4410328 / 86400.0                !    (in k/second)
      endif

!
      return
!...................................
      end subroutine rswinit
!-----------------------------------



!-----------------------------------
      subroutine cldprsw                                                &
!...................................

!  ---  inputs:
     &     ( tasw,clsw,ctot,cwp,cip,rew,rei,cda1,cda2,cda3,cda4,        &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       taucw,ssacw,asycw,fffcw                                    &
     &     )

!  ==================================================================  !
!                                                                      !
!  compute the cloud optical property functions for each cloudy layer  !
!                                                                      !
!                                                                      !
!  inputs:                                                             !
!     tasw    - layer mean temperature (k)                    NPTS*NLAY!
!     clsw    - layer cloud fraction                          NPTS*NLAY!
!     ctot    - column total cloud in fraction                NPTS     !
!     - - -  for prognostic cloud scheme  (iflagc > 0)  - - -          !
!     cwp,cip - layer cloud water/ice path  (g/m**2)          NPTS*NLAY!
!     rew,rei - effective radius for water/ice cloud (micron) NPTS*NLAY!
!     cda1    - layer rain water path  (g/m**2)               NPTS*NLAY!
!     cda2    - effective radius for rain drop (micron)       NPTS*NLAY!
!     cda3    - layer snow water path  (g/m**2)               NPTS*NLAY!
!               (if use fu's snow formula, snow path needs to be       !
!                normalized by snow density (g/m**3/1.0e6) to get micron
!     cda4    - effective radius for snow flakes     (micron) NPTS*NLAY!
!     - - -  for diagnostic cloud scheme  (iflagc = 0)  - - -          !
!     cda1    - input cloud optical depth                     NPTS*NLAY!
!     cda2    - input cloud single scattering albedo          NPTS*NLAY!
!     cda3    - input cloud asymmetry factor                  NPTS*NLAY!
!     cda4    -  optional use                                 NPTS*NLAY!
!     NPTS    - horizontal points                             1        !
!     NLAY,NLP1-vertical layer/level numbers                  1        !
!                                                                      !
!  control flag set up in module "module_radsw_cntr_para"              !
!     iflagc  - control flag for different cloud property method       !
!               =0: input cloud opt depth, fixed ssa, asy              !
!               =1: input cwp/cip, chou's coeff for uv and 3 nir bands !
!               =2: input cwp/cip, chou (1999) coeff for uv bands      !
!                   and 3 nir bands, or gfs's tuned 1 nir band coeff   !
!               =3: input cwp/cip, chou (2002) coeff for uv bands      !
!                   and 3 nir bands, or gfs's tuned 1 nir band coeff   !
!                                                                      !
!  outputs:                                                            !
!     taucw   - cloud optical depth weighting function      NPTS*L*NBD0!
!     ssacw   - cloud single scattering weighting function  NPTS*L*NBD0!
!     asycw   - cloud asymmetry coeff weighting function    NPTS*L*NBD0!
!     fffcw   - cloud forward scattering weighting function NPTS*L*NBD0!
!                                                                      !
!  ==================================================================  !
!
      use module_radsw_cldprtb

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: tasw, clsw,  &
     &       cwp, cip, rew, rei, cda1, cda2, cda3, cda4

      real (kind=kind_phys), dimension(:),   intent(in) :: ctot

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: taucw,    &
     &       ssacw, asycw, fffcw

!  ---  locals:
      real (kind=kind_phys),dimension(NLAY) :: cwp1, cip1, crp1,        &
     &       csp1, rew1, rei1, res1, rewi, reii, resi,                  &
     &       tem1d1, tem1d2, tem1d3, tem1d4

      real (kind=kind_phys) :: tauc, fice, taur, taus

      real (kind=kind_phys) :: tau1, tau2, ssa1, ssa2, ssa3, ssa4,      &
     &       asy1, asy2, asy3, asy4, tem0d1, tem0d2, tem0d3, tem0d4
      integer :: i, k, m

!
!===> ...  begin here
!
      do m = 1, NBD0
        do k = 1, NLAY
          do i = 1, NPTS
            taucw(i,k,m) = zero
            ssacw(i,k,m) = zero
            asycw(i,k,m) = zero
            fffcw(i,k,m) = zero
          enddo
        enddo
      enddo

!  --- ... compute only cloudy sky points

      lab_do_NPTS : do i = 1, NPTS

        lab_if_ctot : if (ctot(i) > zero) then

          if (iflagc <= 1) then
            tem0d1 = con_t0c
            tem0d2 = con_t0c - 10.0
            tem0d3 = 1.0 / 30.0

            do k = 1, NLAY
              tem1d1(k) = max(zero, min(10.0, tem0d1-tasw(i,k)))*0.1
              tem1d2(k) = 1.0 - tem1d1(k)
              tem1d3(k) = max(zero, min(30.0, tem0d2-tasw(i,k) ))*tem0d3
              tem1d4(k) = 1.0 - tem1d3(k)
            enddo
          endif

          if (iflagc >= 1) then
            do k = 1, NLAY
              cwp1(k) = cwp(i,k)                                ! liq water
              cip1(k) = cip(i,k)                                ! ice water
              rew1(k) = rew(i,k)
              rei1(k) = rei(i,k)     
              crp1(k) = cda1(i,k)                               ! rain
              csp1(k) = cda3(i,k)                               ! snow
              res1(k) = cda4(i,k)
            enddo

            where (rew1 > fsmall)
              rewi = 1.0 / rew1
            elsewhere
              rewi = zero
            end where

            where (rei1 > fsmall)
              reii = 1.0 / rei1
            elsewhere
              reii = zero
            end where

            where (res1 > fsmall)
              resi = 1.0 / res1
            elsewhere
              resi = zero
            end where
          endif

!  --- ... compute cloud optical properties

          lab_if_iflg : if (iflagc == 0) then
!  --- ... use input optical depth, no cwp/cip given

            tem0d1 = 20.0
            tem0d2 = 1.0 / tem0d1
            do k = 1, NLAY
              if (clsw(i,k) > zero) then
                tauc = cda1(i,k)
                fice = max(zero, min(tem0d1, con_ttp-tasw(i,k)))*tem0d2

                tau2 = fice * tauc
                tau1 = tauc - tau2

                do m = mbstr, mbend       ! loop for spectral bands
                  taucw(i,k,m) = tauc

                  ssa1 = tau1                                           &
     &                 * (tem1d1(k)*ssaw0(m,1) + tem1d2(k)*ssaw0(m,2))
                  ssa2 = tau2                                           &
     &                 * (tem1d3(k)*ssai0(m,1) + tem1d4(k)*ssai0(m,2))
                  ssacw(i,k,m) = ssa1 + ssa2

                  asy1 = tem1d1(k)*asyw0(m,1) + tem1d2(k)*asyw0(m,2)
                  asy2 = tem1d3(k)*asyi0(m,1) + tem1d4(k)*asyi0(m,2)
                  tem0d1 = asy1 * ssa1
                  tem0d2 = asy2 * ssa2
                  asycw(i,k,m) = tem0d1 + tem0d2
                  fffcw(i,k,m) = asy1*tem0d1 + asy2*tem0d2
                enddo
              endif

            enddo   ! end_do_k_loop

          elseif (iflagc == 1) then       lab_if_iflg
!  --- ... use cwp/cip to compute tau, ssa, asy, with t-adjust
!          rain drops follows chou, snow follows ncar or fu, their optical
!          depth are independent of frequency bands.  for fu's formula, csp1
!          needs to be normalized by snow density.

            do k = 1, NLAY
              if (clsw(i,k) > zero) then
                taur = crp1(k) * a0r
                taus = csp1(k) * (a0s + a1s*resi(k))

                do m = mbstr, mbend       ! loop for spectral bands
                  tau1 = cwp1(k)                                        &
     &               * (tem1d1(k) * (a0w1(m,1) + a1w1(m,1)*rewi(k))     &
     &               +  tem1d2(k) * (a0w1(m,2) + a1w1(m,2)*rewi(k)))
                  tau2 = cip1(k)                                        &
     &               * (tem1d3(k) * (a0i1(m,1) + a1i1(m,1)*reii(k))     &
     &               +  tem1d4(k) * (a0i1(m,2) + a1i1(m,2)*reii(k)))
                  taucw(i,k,m) = tau1 + tau2 + taur + taus

                  ssa1 = tau1 * (1.0                                    &
     &               - (tem1d1(k) * (b0w1(m,1) + b1w1(m,1)*rew1(k))     &
     &               +  tem1d2(k) * (b0w1(m,2) + b1w1(m,2)*rew1(k))))
                  ssa2 = tau2 * (1.0 - (tem1d3(k) * (b0i1(m,1)          &
     &               + (b1i1(m,1) + b2i1(m,1) * rei1(k)) * rei1(k))     &
     &                            +  tem1d4(k) * (b0i1(m,2)             &
     &               + (b1i1(m,2) + b2i1(m,2) * rei1(k)) * rei1(k))))
                  ssa3 = taur * (1.0 - b0r(m))
                  ssa4 = taus * (1.0 - (b0s(m) + b1s(m)*res1(k)))
                  ssacw(i,k,m) = ssa1 + ssa2 + ssa3 + ssa4

                  asy1 = tem1d1(k) * (c0w1(m,1) + c1w1(m,1)*rew1(k))    &
     &                 + tem1d2(k) * (c0w1(m,2) + c1w1(m,2)*rew1(k))
                  asy2 = tem1d3(k) * (c0i1(m,1) + (c1i1(m,1)            &
     &                             +  c2i1(m,1) * rei1(k)) * rei1(k))   &
     &                 + tem1d4(k) * (c0i1(m,2) + (c1i1(m,2)            &
     &                             +  c2i1(m,2) * rei1(k)) * rei1(k))
                  asy3 = c0r(m)
                  asy4 = c0s(m)
                  tem0d1 = asy1 * ssa1
                  tem0d2 = asy2 * ssa2
                  tem0d3 = asy3 * ssa3
                  tem0d4 = asy4 * ssa4
                  asycw(i,k,m) = tem0d1 + tem0d2 + tem0d3 + tem0d4
                  fffcw(i,k,m) = asy1*tem0d1 + asy2*tem0d2              &
     &                         + asy3*tem0d3 + asy4*tem0d4
                enddo
              endif

            enddo   ! end_do_k_loop

          elseif (iflagc == 2) then       lab_if_iflg
!  --- ... use cwp/cip to compute tau, ssa, asy from chou (1999)
!          rain drops follows chou, snow follows ncar or fu, their optical
!          depth are independent of frequency bands.  for fu's formula, csp1
!          needs to be normalized by snow density.

            do k = 1, NLAY
              if (clsw(i,k) > zero) then
                taur = crp1(k) * a0r
                taus = csp1(k) * (a0s + a1s*resi(k))

                do m = mbstr, mbend       ! loop for spectral bands
                  tau1 =  cwp1(k) * (a0w2(m) + a1w2(m)*rewi(k))
                  tau2 =  cip1(k) * (a0i2(m) + a1i2(m)*reii(k))
                  taucw(i,k,m) = tau1 + tau2 + taur + taus

                  ssa1 =  tau1 * (1.0 - (b0w2(m) + (b1w2(m)             &
     &                 +  b2w2(m) * rew1(k)) * rew1(k)))
                  ssa2 =  tau2 * (1.0 - (b0i2(m) + (b1i2(m)             &
     &                 +  b2i2(m) * rei1(k)) * rei1(k)))
                  ssa3 = taur * (1.0 - b0r(m))
                  ssa4 = taus * (1.0 - (b0s(m) + b1s(m)*res1(k)))
                  ssacw(i,k,m) = ssa1 + ssa2 + ssa3 + ssa4

                  asy1 = c0w2(m) + (c1w2(m) + c2w2(m)*rew1(k)) * rew1(k)
                  asy2 = c0i2(m) + (c1i2(m) + c2i2(m)*rei1(k)) * rei1(k)
                  asy3 = c0r(m)
                  asy4 = c0s(m)
                  tem0d1 = asy1 * ssa1
                  tem0d2 = asy2 * ssa2
                  tem0d3 = asy3 * ssa3
                  tem0d4 = asy4 * ssa4
                  asycw(i,k,m) = tem0d1 + tem0d2 + tem0d3 + tem0d4
                  fffcw(i,k,m) = asy1*tem0d1 + asy2*tem0d2                &
     &                       + asy3*tem0d3 + asy4*tem0d4
                enddo
              endif

            enddo   ! end_do_k_loop

          elseif (iflagc == 3) then       lab_if_iflg
!  --- ... use cwp/cip to compute tau, ssa, asy from chou (2002)
!          rain drops follows chou, snow follows ncar or fu, their optical
!          depth are independent of frequency bands.  for fu's formula, csp1
!          needs to be normalized by snow density.

            do k = 1, NLAY
              if (clsw(i,k) > zero) then
                taur = crp1(k) * a0r
                taus = csp1(k) * (a0s + a1s*resi(k))

                do m = mbstr, mbend       ! loop for spectral bands
                  tau1 =  cwp1(k) * (a0w3(m) + a1w3(m)*rewi(k))
                  tau2 =  cip1(k) * (a0i3(m) + a1i3(m)*reii(k))
                  taucw(i,k,m) = tau1 + tau2 + taur + taus

                  ssa1 =  tau1 * (1.0 - (b0w3(m) + (b1w3(m)             &
     &                 +  b2w3(m) * rew1(k)) * rew1(k)))
                  ssa2 =  tau2 * (1.0 - (b0i3(m) + (b1i3(m)             &
     &                 +  b2i3(m) * rei1(k)) * rei1(k)))
                  ssa3 = taur * (1.0 - b0r(m))
                  ssa4 = taus * (1.0 - (b0s(m) + b1s(m)*res1(k)))
                  ssacw(i,k,m) = ssa1 + ssa2 + ssa3 + ssa4

                  asy1 = c0w3(m) + (c1w3(m) + c2w3(m)*rew1(k)) * rew1(k)
                  asy2 = c0i3(m) + (c1i3(m) + c2i3(m)*rei1(k)) * rei1(k)
                  asy3 = c0r(m)
                  asy4 = c0s(m)
                  tem0d1 = asy1 * ssa1
                  tem0d2 = asy2 * ssa2
                  tem0d3 = asy3 * ssa3
                  tem0d4 = asy4 * ssa4
                  asycw(i,k,m) = tem0d1 + tem0d2 + tem0d3 + tem0d4
                  fffcw(i,k,m) = asy1*tem0d1 + asy2*tem0d2              &
     &                         + asy3*tem0d3 + asy4*tem0d4
                enddo
              endif

            enddo   ! end_do_k_loop

          else   lab_if_iflg

            print *,'  error setting of iflagc! stop in cldprsw'
            stop

          endif  lab_if_iflg

        endif  lab_if_ctot

      enddo  lab_do_NPTS

!
      return
!...................................
      end subroutine cldprsw
!-----------------------------------



!-----------------------------------
      subroutine swflux                                                 &
!...................................

!  ---  inputs:
     &     ( tau,ssc,g0,csm,zth,alb,ald,cld1,lcld,                      &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  outputs:
     &       upflux,dwflux,dwsfcb,dwsfcd                                &
     &     )

!  ==================================================================  !
!                                                                      !
!  uses the delta-eddington approximation to compute the bulk          !
!  scattering properties of a single layer coded following             !
!  coakley et al.  (jas, 1982)                                         !
!                                                                      !
!  inputs:                                                             !
!    tau  : the effective optical thickness                NPTS*NLAY*2 !
!    ssc  : the effective single scattering albedo         NPTS*NLAY*2 !
!    g0   : the effective asymmetry factor                 NPTS*NLAY*2 !
!    csm  : the effective secant of the zenith angle       NPTS        !
!    zth  : the cosin of the zenith angle                  NPTS        !
!    alb  : surface albedo for direct radiation            NPTS        !
!    ald  : surface albedo for diffused radiation          NPTS        !
!    cld1 : layer cloud fraction normalized by total cloud NPTS*NLAY   !
!    lcld : logical flag: t=cld exist; f=all clear points   1          !
!    NPTS : number of horizontal points                     1          !
!    NLAY : vertical number of layers                       1          !
!    NLP1 : vertical number of levels (=NLAY+1)             1          !
!                                                                      !
!  outputs:                                                            !
!    upflux: upward fluxes                                 NPTS*NLP1*2 !
!    dwflux: downward fluxes                               NPTS*NLP1*2 !
!    dwsfcb: downward surface flux direct component        NPTS*2      !
!    dwsfcd: downward surface flux diffused component      NPTS*2      !
!                                                                      !
!  ==================================================================  !
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1

      logical, intent(in) :: lcld

      real (kind=kind_phys),dimension(:,:,:),intent(in) :: tau, ssc, g0

      real (kind=kind_phys),dimension(:,:),  intent(in) :: cld1

      real (kind=kind_phys),dimension(:), intent(in) :: csm,zth,alb,ald

!  ---  outputs:
      real (kind=kind_phys),dimension(:,:,:),intent(out):: upflux,dwflux

      real (kind=kind_phys),dimension(:,:),  intent(out):: dwsfcb,dwsfcd


!  ---  diffuse incident radiation is approximated by beam radiation
!       with an incident angle of 53 degrees. cos(53) = 0.602

      real (kind=kind_phys) :: zthdf, csmdf
      parameter (zthdf=0.602,  csmdf=1.0/zthdf)        ! csmdf=1/zthdf

!  ---  local arrays:
      real (kind=kind_phys), dimension(NPTS,NLP1,2) :: ttb0, ttdf, ttbm,&
     &       rrdf, rrbm
      real (kind=kind_phys), dimension(NPTS,NLP1)   :: tbm, tdn, rup,   &
     &       rfu, rfd

      real (kind=kind_phys) :: taup, sscp, gp, oms1, ogs1, tlam, den1,  &
     &       rf1, tf1, goms, gogs, slam, u1, e1, u1p1, u1m1, u1e, u1eme,&
     &       gama, alfa, amg, apg, za, zb, zc, zd, zthf, zthd, zzf, zzd,&
     &       tbf, tbd, denf, dend

      integer :: i, k, kp1, km1, ip
!
!===> ... begin
!
      lab_do_ip1 : do ip = 1, 2       ! ip=1: clear-sky, =2: cloudy-sky

!  ---  use delta-eddington method

        lab_do_NLAY : do k = 1, NLAY
          lab_do_NPTS : do i = 1, NPTS

            taup = tau(i,k,ip)
            sscp = ssc(i,k,ip)
            gp   = g0(i,k,ip)

            oms1 = 1.0 - sscp
            ogs1 = 1.0 - sscp*gp
            tlam = 3.0 * oms1*ogs1
            slam = sqrt( tlam )

            u1   = 1.5 * ogs1 / slam
            u1p1 = u1 + 1.0
            u1m1 = u1 - 1.0
            e1   = exp(-min(100.0e0, taup*slam))
            u1e  = u1 * e1
            u1eme= u1e - e1

            den1 = 1.0 / ((u1p1 + u1eme) * (u1p1 - u1eme))
            rf1  = (u1p1 + u1e + e1) * (u1m1 - u1eme) * den1
            tf1  = 4.0 * u1e * den1
!
!===> ... compute layer transmissions and reflections
!         rrdf :  layer reflection for diffused radiation
!         rrbm :  layer reflection for beam radiation
!         ttdf :  layer total transmission for diffused radiation
!         ttbm :  layer total transmission for beam radiation
!         ttb0 :  layer direct transmission for beam radiation
 
!  ---  diffuse radiation
!       -----------------

            zthf = zthdf
            zzf  = zthf * zthf
            denf = 1.0 - tlam*zzf
            if (abs(denf) < fpmin) then         !===> ... safety check
              zthf = zthf + fsmall
              zzf  = zthf * zthf
              denf = 1.0 - tlam*zzf
            endif
            denf = sscp / denf

            goms = 3.00 * gp * oms1
            gogs = 0.75 * (gp + ogs1)

            gama = 0.50 * (1.0 + goms*zzf) * denf
            alfa = gogs * zthdf * denf
            amg = alfa - gama
            apg = alfa + gama

            tbf = exp(-min(100.0e0, taup*csmdf))
            zc  = amg * tbf
            rrdf(i,k,ip) = max( zero, zc*tf1 + apg*rf1 - amg )
            ttdf(i,k,ip) = max( zero, zc*rf1 + apg*tf1 + (1.0-apg)*tbf )
!
!  ---  direct radiation
!       ----------------
            zthd = zth(i)
            zzd  = zthd * zthd
            dend = 1.0 - tlam*zzd
            if (abs(dend) < fpmin) then         !===> ... safety check
              zthd = zthd + fsmall
              zzd  = zthd * zthd
              dend = 1.0 - tlam*zzd
            endif
            dend = sscp / dend

            gama = 0.50 * (1.0 + goms*zzd) * dend
            alfa = gogs * zthd * dend
            amg = alfa - gama
            apg = alfa + gama
 
            tbd = exp( -min(100.0e0, taup*csm(i)) )
            zd  = amg * tbd
            ttb0(i,k,ip) = max( zero, tbd )
            rrbm(i,k,ip) = max( zero, zd*tf1 + apg*rf1 - amg )
            ttbm(i,k,ip) = max( zero, zd*rf1 + apg*tf1 + (1.0-apg)*tbd )

          enddo  lab_do_NPTS
        enddo  lab_do_NLAY

        do i = 1, NPTS
          ttb0(i,NLP1,ip) = zero
          rrbm(i,NLP1,ip) = alb(i)
          ttbm(i,NLP1,ip) = zero
          rrdf(i,NLP1,ip) = ald(i)
          ttdf(i,NLP1,ip) = zero
        enddo

        if ( ip==1 .and. .not.lcld ) then
          do k = 1, NLP1
          do i = 1, NPTS
            ttb0(i,k,2) = ttb0(i,k,1)
            rrbm(i,k,2) = rrbm(i,k,1)
            ttbm(i,k,2) = ttbm(i,k,1)
            rrdf(i,k,2) = rrdf(i,k,1)
            ttdf(i,k,2) = ttdf(i,k,1)
          enddo
          enddo

          exit  lab_do_ip1
        endif

      enddo  lab_do_ip1

!  ---  combine clear-sky and cloudy-sky parameters

      if ( lcld ) then
        do k = 1, NLAY
        do i = 1, NPTS
          za = 1.0 - cld1(i,k)
          zb = cld1(i,k)

          ttb0(i,k,2) = za*ttb0(i,k,1) + zb*ttb0(i,k,2)
          rrbm(i,k,2) = za*rrbm(i,k,1) + zb*rrbm(i,k,2)
          ttbm(i,k,2) = za*ttbm(i,k,1) + zb*ttbm(i,k,2)
          rrdf(i,k,2) = za*rrdf(i,k,1) + zb*rrdf(i,k,2)
          ttdf(i,k,2) = za*ttdf(i,k,1) + zb*ttdf(i,k,2)
        enddo
        enddo
      endif
!
!===> ... compute upward and downward fluxes
!

      lab_do_ip2 : do ip = 1, 2         ! ip=1: clear-sky; =2: cloudy-sky

!  ---  set values at the top and surface

        do i = 1, NPTS
          tbm(i,1) = ttb0(i,1,ip)
          tdn(i,1) = ttbm(i,1,ip)
          rfd(i,1) = rrdf(i,1,ip)

          rfu(i,NLP1) = ald(i)
          rup(i,NLP1) = alb(i)
        enddo

!  ---  layers added downward starting from top

        do k = 1, NLAY
          kp1 = k + 1
          do i = 1, NPTS
            tbm(i,kp1) = tbm(i,k) * ttb0(i,kp1,ip)
          enddo
        enddo

        where (tbm < 1.0e-30) tbm = zero

        lab_do_topdown : do k = 2, NLP1
          km1 = k - 1
          do i = 1, NPTS
            den1 = ttdf(i,k,ip) / (1.0 - rfd(i,km1) * rrdf(i,k,ip))
            tdn(i,k) = tbm(i,km1)*ttbm(i,k,ip) + (tdn(i,km1)-tbm(i,km1) &
     &               + tbm(i,km1)*rrbm(i,k,ip)*rfd(i,km1)) * den1
            rfd(i,k) = rrdf(i,k,ip) + ttdf(i,k,ip)*rfd(i,km1)*den1
          enddo
        enddo  lab_do_topdown

!  ---  layers added upward starting from surface

        lab_do_bottomup : do k = NLAY, 1, -1
          kp1 = k + 1
          do i = 1, NPTS
            den1 = ttdf(i,k,ip) / (1.0 - rfu(i,kp1) * rrdf(i,k,ip))
            rup(i,k) = rrbm(i,k,ip) + ((ttbm(i,k,ip)-ttb0(i,k,ip))      &
     &               * rfu(i,kp1) + ttb0(i,k,ip)*rup(i,kp1))*den1
            rfu(i,k) = rrdf(i,k,ip) + ttdf(i,k,ip)*rfu(i,kp1)*den1
          enddo
        enddo  lab_do_bottomup

!  ---  find upward and downward fluxes

        do k = 2, NLP1
          km1 = k - 1
          do i = 1, NPTS
            den1 = 1.0 / (1.0 - rfd(i,km1)*rfu(i,k))
            za = tbm(i,km1) * rup(i,k)
            zb = tdn(i,km1) - tbm(i,km1)
            upflux(i,k,ip) = (za + zb*rfu(i,k)) * den1
            dwflux(i,k,ip) = tbm(i,km1) + (za*rfd(i,km1) + zb) * den1
          enddo
        enddo

!  ---  surface downward fluxes

        do i = 1, NPTS
          upflux(i,1,ip) = rup(i,1)
          dwflux(i,1,ip) = 1.0
          dwsfcb(i,ip) = tbm(i,NLAY)
          dwsfcd(i,ip) = dwflux(i,NLP1,ip) - tbm(i,NLAY)
        enddo

        if (ip==1 .and. .not.lcld) then
          do k = 1, NLP1
          do i = 1, NPTS
            upflux(i,k,2) = upflux(i,k,1)
            dwflux(i,k,2) = dwflux(i,k,1)
          enddo
          enddo

          do i = 1, NPTS
            dwsfcb(i,2) = dwsfcb(i,1)
            dwsfcd(i,2) = dwsfcd(i,1)
          enddo

          exit  lab_do_ip2
        endif

      enddo  lab_do_ip2
 
!
      return
!...................................
      end subroutine swflux
!-----------------------------------



!-----------------------------------
      subroutine flxco2                                                 &
!...................................

!  ---  inputs:
     &     ( swc,u1,du,NU,swh,w1,dw,NW,tbl,                             &
     &       NPTS, NLAY, NLP1,                                          &
!  ---  in/outputs:
     &       dflx                                                       &
     &     )

!  ==================================================================  !
!  compute the absorption due to co2. ref: chou (j. climate, 1990,     !
!     209-217)                                                         !
!     updated sep. 1999 based on nasa/tm-1999-104606, vol 15.          !
!                                                                      !
!  the effect of co2 absorption below the cloud top is neglected.      !
!                                                                      !
!  input variables:                                                    !
!     swc         : column amount of co2                 NPTS*NLP1     !
!     swh         : column amount of water vapor         NPTS*NLP1     !
!     u1,du,w1,dw : coefficients                             1         !
!     tbl         : look up co2 absorption table           NU*NW       !
!     NU,NW       : table dimensions                         1         !
!     NPTS        : horizontal points                        1         !
!     NLAY,NLP1   : vertical layer/level numbers             1         !
!                                                                      !
!  output variables:                                                   !
!     dflx        : additional flux reduction due to co2 for clear sky !
!                                                        NPTS*NLP1     !
!  ==================================================================  !
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: NU, NW, NPTS, NLAY, NLP1

      real (kind=kind_phys),dimension(:,:),intent(in) :: tbl, swc, swh

      real (kind=kind_phys), intent(in) :: u1, du, w1, dw

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: dflx

!  ---  locals:
      real (kind=kind_phys) :: x1, y1, x2, y2, tbl0, dc, dd, clog, wlog
      integer :: i, k, ic0, iw0, ic1, iw1
!
!===>  ... table look-up for the reduction of clear-sky solar
!
      x1 = u1 - 0.5*du
      y1 = w1 - 0.5*dw

      do k = 2, NLP1
        do i = 1, NPTS

          clog = swc(i,k)
          wlog = swh(i,k)

          ic0 = max(2, min(NU, int((clog - x1)/du + 1.0) ))
          iw0 = max(2, min(NW, int((wlog - y1)/dw + 1.0) ))
          ic1 = ic0 - 1
          iw1 = iw0 - 1

          dc = clog - float(ic0-2)*du - u1
          dd = wlog - float(iw0-2)*dw - w1

          tbl0 = tbl(ic1,iw1)

          x2 = tbl0 + (tbl(ic1,iw0) - tbl0)/dw * dd
          y2 =        (tbl(ic0,iw1) - tbl0)/du * dc
          dflx(i,k) = dflx(i,k) + x2 + y2

        enddo
      enddo
!
      return
!...................................
      end subroutine flxco2
!-----------------------------------

!
!........................................!
      end module module_radsw_main       !
!========================================!

