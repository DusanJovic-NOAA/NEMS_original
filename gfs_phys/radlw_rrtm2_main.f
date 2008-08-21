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
!          inputs:                                                     !
!           (pmid,pint,tmid,tint,qnm,o3mr,gasvmr,                      !
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
!            upfxc              level upward flux for total sky        !
!            dnfxc              level downward flux for total sky      !
!            upfx0              level upward flux for clear sky        !
!            dnfx0              level downward flux for clear sky      !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'                                               !
!       'module physcons'                                              !
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
!    ncep modifications history log:                                   !
!                                                                      !
!       aug 2002,  yu-tai hou -- received v3.0 version from aer        !
!       dec 2004,  yu-tai hou                                          !
!                  rewritten code into fortran 90 (w/o optim codes)    !
!       apr 2005,  yu-tai hou                                          !
!                  minor modifications on module structures            !
!       mar 2007,  yu-tai hou                                          !
!                  add aerosol effect for lw radiation                 !
!       apr 2007,  yu-tai hou                                          !
!                  add spectral band heating as optional output        !
!                                                                      !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radlw_main           !
!........................................!
!
      use machine,             only : kind_phys
      use physcons,            only : con_g, con_cp, con_avgd, con_amd, &
     &                                con_amw, con_amo3

      use module_radlw_parameters
      use module_radlw_cntr_para

      use module_radlw_refprof
      use module_radlw_avplank
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGLW='RRTM-LW v3.0    Oct 2002'
!     character(24), parameter :: VTAGLW='RRTM-LW v3.0    Mar 2007'
      character(24), parameter :: VTAGLW='RRTM-LW v3.0    Apr 2007'

!  ---  constant values
      real (kind=kind_phys) :: eps, oneminus, bpade, stpfac, wtdiff     &
     &,      co2fac, tblint, secdiff, rec_6, f_zero

      parameter (eps=1.0e-6,  oneminus=1.0-eps)
      parameter (bpade=1.0/0.278)      ! pade approximation constant
      parameter (stpfac=296./1013.)
      parameter (wtdiff=0.5)
!     parameter (avgdro=6.02214199e+23)! avogadro constant  (1/mol)
      parameter (tblint=10000.0)
      parameter (secdiff=1.66)
      parameter (rec_6=0.166667)
      parameter (f_zero=0.0)

!  ...  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys) :: amdw, amdo3

      parameter (amdw =con_amd/con_amw)
      parameter (amdo3=con_amd/con_amo3)

!  ...  band indices
      integer :: nspa(NBANDS), nspb(NBANDS), ngb(NBANDS)

      data nspa / 1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9 /
      data nspb / 1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0 /
      data ngb  / NG01, NG02, NG03, NG04, NG05, NG06, NG07, NG08,       &
     &            NG09, NG10, NG11, NG12, NG13, NG14, NG15, NG16 /

!  ...  band wavenumber intervals
!     real (kind=kind_phys) :: wavenum1(NBANDS), wavenum2(NBANDS)
!     data wavenum1/                                                    &
!    &         10.,  350.,  500.,  630.,  700.,  820.,  980., 1080.,    &
!    &       1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600. /
!     data wavenum2/                                                    &
!    &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
!    &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /
      real (kind=kind_phys) :: delwave(NBANDS)
      data delwave / 340., 150., 130.,  70., 120., 160., 100., 100.,    &
     &               210.,  90., 320., 280., 170., 130., 220., 650. /

!  ...  when standard first-order gaussian quadrature is chosen as the
!       method to approximate the integral over angles that yields flux
!       from radiances, then secreg(i,j) is the secant of the ith (out of
!       a total of j angles) and wtreg(i,j) is the corresponding weight.

      real (kind=kind_phys), dimension(MAXANG,MAXANG) :: secreg, wtreg

      data  secreg(1:4,1:4) /                                           &
     &         1.50000000,  0.0       ,  0.0       ,  0.0       ,       &
     &         1.18350343,  2.81649655,  0.0       ,  0.0       ,       &
     &         1.09719858,  1.69338507,  4.70941630,  0.0       ,       &
     &         1.06056257,  1.38282560,  2.40148179,  7.15513024 /

      data  wtreg(1:4,1:4)  /                                           &
     &         0.5         , 0.0         , 0.0         , 0.0         ,  &
     &         0.3180413817, 0.1819586183, 0.0000000000, 0.0         ,  &
     &         0.2009319137, 0.2292411064, 0.0698269799, 0.0         ,  &
     &         0.1355069134, 0.2034645680, 0.1298475476, 0.0311809710 /

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

!  ---  those data will be set up only once by "rlwinit"

!  ...  fluxfac, heatfac are factors for fluxes (in w/m**2) and heating
!       rates (in k/day, or k/sec set by subroutine 'rlwinit')
!       semiss0 are default surface emissivity for each bands

      real (kind=kind_phys) :: fluxfac, heatfac, semiss0(NBANDS)

      real (kind=kind_phys), dimension(0:NTBL) :: tautb, tf, trans

      public lwrad, rlwinit


! =================
      contains
! =================

!-----------------------------------
      subroutine lwrad                                                  &
!...................................

!  ---  inputs:
     &     ( pmid,pint,tmid,tint,qnm,o3mr,gasvmr,                       &
     &       clouds,iovr,aerosols,sfemis,                               &
     &       NPTS, NLAY, NLP1, iflip, lprnt,                            &
!  ---  outputs:
     &       hlwc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HLW0,HLWB,FLXPRF                                           &
     &     )

!  ====================  defination of variables  ===================  !
!                                                                      !
!  input variables:                                                    !
!     pmid   (NPTS,NLAY)    - layer pressures (mb)                     !
!     pint   (NPTS,NLP1)    - interface pressures (mb)                 !
!     tmid   (NPTS,NLAY)    - layer temperature (k)                    !
!     tint   (NPTS,NLP1)    - interface temperatures (k)               !
!     qnm    (NPTS,NLAY)    - layer h2o mixing ratio (gm/gm)*see inside!
!     o3mr   (NPTS,NLAY)    - layer o3 mixing ratio (gm/gm) *see inside!
!     gasvmr (NPTS,NLAY,:)  -  atmospheric gases amount:               !
!                       (check module_radiation_gases for definition)  !
!       gasvmr(:,:,1)   -      co2 volume mixing ratio                 !
!       gasvmr(:,:,2)   -      n2o volume mixing ratio                 !
!       gasvmr(:,:,3)   -      ch4 volume mixing ratio                 !
!       gasvmr(:,:,4)   -      o2  volume mixing ratio                 !
!       gasvmr(:,:,5)   -      co  volume mixing ratio                 !
!       gasvmr(:,:,6)   -      cfc11 volume mixing ratio               !
!       gasvmr(:,:,7)   -      cfc12 volume mixing ratio               !
!       gasvmr(:,:,8)   -      cfc22 volume mixing ratio               !
!       gasvmr(:,:,9)   -      ccl4  volume mixing ratio               !
!     clouds (NPTS,NLAY,:)  - layer cloud profiles:                    !
!                       (check module_radiation_clouds for definition) !
!                ---  for  iflagliq > 0  ---                           !
!       clouds(:,:,1)  -   layer total cloud fraction                  !
!       clouds(:,:,2)  -   layer cloud liq water path      (g/m**2)    !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)    !
!       clouds(:,:,4)  -   layer cloud ice water path      (g/m**2)    !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)    !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)    !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)    !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)    !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)    !
!                ---  for  iflagliq = 0  ---                           !
!       clouds(:,:,1)  -   layer total cloud fraction                  !
!       clouds(:,:,2)  -   layer cloud optical depth                   !
!       clouds(:,:,3)  -   layer cloud single scattering albedo        !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                !
!     iovr                  - control flag for cloud overlapping       !
!                             =0: random overlapping clouds            !
!                             =1: max/ran overlapping clouds           !
!     aerosols(NPTS,NLAY,NBANDS,:) - aerosol optical properties:       !
!                       (check module_radiation_aerosols for definition!
!        (:,:,:,1)     - optical depth                                 !
!        (:,:,:,2)     - single scattering albedo                      !
!        (:,:,:,3)     - asymmetry parameter                           !
!     sfemis (NPTS)         - surface emissivity                       !
!     NPTS                  - total number of horizontal points        !
!     NLAY,NLP1             - total number of vertical layers, levels  !
!     iflip                 - control flag for in/out vertical index   !
!                             =0: index from toa to surface            !
!                             =1: index from surface to toa            !
!     lprnt                 - cntl flag for diagnostic print out       !
!                                                                      !
!  control parameters in module "module_radlw_cntr_para":              !
!     ilwrate               - heating rate unit selections             !
!                             =1: output in k/day                      !
!                             =2: output in k/second                   !
!     iaerlw                - control flag for aerosols (not yet)      !
!                             =0: do not include aerosol effect        !
!                             >0: include aerosol effect               !
!     irgaslw               - control flag for rare gases              !
!                             (ch4,n2o,o2, etc.)                       !
!                             =0: do not include rare gases            !
!                             =1: include all rare gases               !
!     icfclw                - control flag for cfc gases               !
!                             =0: do not include cfc gases             !
!                             =1: include all cfc gases                !
!     iflagliq              - liq-cloud optical properties contrl flag !
!                             =0: input cld opt dep, ignor iflagice    !
!                             =1: input cwp,cip, (ccm2) ignor iflagice !
!                             =2: input cwp rew, (ccm3 method)         !
!                             =3: input cwp rew, hu and stamnes (1993) !
!     iflagice              - ice-cloud optical properties contrl flag !
!                       * * * if iflagliq .lt. 2, iflafice is ignored  !
!                             =0: input cip rei, (ccm3 method)         !
!                             =1: input cip rei, ebert and curry (1997)!
!                             =2: input cip rei, streamer (1996)       !
!     numangs               - num of angles for radiance calculation   !
!                             =0: diffusivity approx, cos(ang)=1/1.66  !
!                             =n: n=1-4, cosine of the angles have the !
!                                 standard x-axis values in first-moment
!                                 gaussian quadrature                  !
!                                                                      !
!  output variables:                                                   !
!     hlwc   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)  !
!     topflx (NPTS)         - radiation fluxes at top, component:      !
!                       (check module_radlw_parameters for definition) !
!        upfxc                 total sky upward flux at top (w/m2)     !
!        upfx0                 clear sky upward flux at top (w/m2)     !
!     sfcflx (NPTS)         - radiation fluxes at sfc, component:      !
!                       (check module_radlw_parameters for definition) !
!        upfxc                 total sky upward flux at sfc (w/m2)     !
!        dnfxc                 total sky downward flux at sfc (w/m2)   !
!        dnfx0                 clear sky downward flux at sfc (w/m2)   !
!                                                                      !
!! optional output variables:                                          !
!     hlwb(NPTS,NLAY,NBANDS)- spectral band total sky heating rates    !
!     hlw0   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)  !
!     flxprf (NPTS,NLP1)    - level radiative fluxes (w/m2), components!
!                       (check module_radlw_parameters for definition) !
!        upfxc                 total sky upward flux                   !
!        dnfxc                 total sky dnward flux                   !
!        upfx0                 clear sky upward flux                   !
!        dnfx0                 clear sky dnward flux                   !
!                                                                      !
!  module parameters, control and local variables:                     !
!     NBANDS                - number of longwave spectral bands        !
!     MAXGAS                - maximum number of absorbing gaseous      !
!     MAXXSEC               - maximum number of cross-sections         !
!     NGnn   (nn=1-16)      - number of g-points in band nn            !
!     nspa,nspb(NBANDS)     - number of lower/upper ref atm's per band !
!     delwave(NBANDS)       - longwave band width (wavenumbers)        !
!     bpade                 - pade approximation constant (1/0.278)    !
!     pavel  (NLAY)         - layer pressures (mb)                     !
!     delp   (NLAY)         - layer pressure thickness (mb)            !
!     tavel  (NLAY)         - layer temperatures (k)                   !
!     tz     (0:NLAY)       - level (interface) temperatures (k)       !
!     semiss (NBANDS)       - surface emissivity for each band         !
!     wx     (NLAY,MAXXSEC) - cross-section molecules concentration    !
!     coldry (NLAY)         - dry air column amount                    !
!                                   (1.e-20*molecules/cm**2)           !
!     colbrd (NLAY)         - column density for broadening gases      !
!                                   (1.e-20*molecules/cm**2)           !
!     cldfrac(0:NLP1)       - layer cloud fraction                     !
!     taucloud(NLAY,NBANDS) - layer cloud optical depth for each band  !
!     taug   (NLAY,NGMX)    - gaseous optical depths                   !
!     colamt (NLAY,MAXGAS)  - column amounts of absorbing gases        !
!                             1-MAXGAS are for watervapor, carbon      !
!                             dioxide, ozone, nitrous oxide, methane,  !
!                             oxigen, carbon monoxide, respectively    !
!                             (molecules/cm**2)                        !
!     facij  (NLAY)         - indicator of interpolation factors       !
!                             =0/1: indicate lower/higher temp & height!
!     selffac(NLAY)         - scale factor for self-continuum, equals  !
!                          (w.v. density)/(atm density at 296K,1013 mb)!
!     selffrac(NLAY)        - factor for temp interpolation of ref     !
!                             self-continuum data                      !
!     indself(NLAY)         - index of the lower two appropriate ref   !
!                             temp for the self-continuum interpolation!
!     laytrop               - layer at which switch is made from one   !
!                             combination of key species to another    !
!     totuflux(0:NLAY)      - upward longwave flux (w/m2)              !
!     totdflux(0:NLAY)      - downward longwave flux (w/m2)            !
!     totuclfl(0:NLAY)      - clear-sky upward longwave flux (w/m2)    !
!     totdclfl(0:NLAY)      - clear-sky downward longwave flux (w/m2)  !
!     fnet    (0:NLAY)      - net longwave flux (w/m2)                 !
!     fnetc   (0:NLAY)      - clear-sky net longwave flux (w/m2)       !
!                                                                      !
!                                                                      !
!  =====================    end of definitions    ===================  !
!
      implicit none

!  ---  inputs:
      integer,  intent(in) :: NPTS, NLAY, NLP1, iovr, iflip

      logical,  intent(in) :: lprnt

      real (kind=kind_phys), dimension(:,:), intent(in) :: pint, tint,  &
     &       pmid, tmid, qnm, o3mr

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: gasvmr,    &
     &       clouds

      real (kind=kind_phys), dimension(:,:,:,:), intent(in) :: aerosols

      real (kind=kind_phys), dimension(:),     intent(in) :: sfemis

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: hlwc

      type (topflw_type),    dimension(:),   intent(out) :: topflx
      type (sfcflw_type),    dimension(:),   intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys),dimension(:,:,:),optional,intent(out):: hlwb
      real (kind=kind_phys),dimension(:,:),optional,intent(out):: hlw0
      type (proflw_type),   dimension(:,:),optional,intent(out):: flxprf

!  ---  locals:
      real (kind=kind_phys), dimension(0:NLP1) :: cldfrac

      real (kind=kind_phys), dimension(0:NLAY) :: totuflux, totdflux,   &
     &       totuclfl, totdclfl, htr, htrcl, tz

      real (kind=kind_phys), dimension(NLAY)   :: pavel, tavel, delp,   &
     &       taucl, cwp1, cip1, rew1, rei1, cda1, cda2, cda3, cda4,     &
     &       coldry, colbrd, h2ovmr, o3vmr, fac00, fac01, fac10, fac11, &
     &       forfac, forfrac, selffac, selffrac, minorfrac, scaleminor, &
     &       scaleminorn2, temcol

      real (kind=kind_phys), dimension(NLAY,2) :: rat_h2oco2, rat_h2oo3,&
     &       rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2

      real (kind=kind_phys) :: colamt(NLAY,MAXGAS), wx(NLAY,MAXXSEC),   &
     &       taucloud(NLAY,NBANDS), semiss(NBANDS), plankbnd(NBANDS),   &
     &       planklay(NLAY,NBANDS), planklev(0:NLAY,NBANDS),            &
     &       tauaer(NLAY,NBANDS), htrb(NLAY,NBANDS)

      real (kind=kind_phys) :: fp, ft, ft1, tbound, tem0, tem1, tem2

!aer  real (kind=kind_phys) :: rh1(NLAY),dz1(NLAY),cmix1(NXC),denn1(NXC)
!aer  integer              :: idxc1(NXC),idm1(NLAY)

      integer, dimension(NLAY) :: jp, jt, jt1, indself, indfor, indminor
      integer  :: jp1, j, k, k1, iplon, laytrop, ireflect, icldatm

!
!===> ... begin here
!

!     ireflect = 1       ! --- specular surface reflection
      ireflect = 0       ! --- lambertian surface reflection

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

!  ---  loop over horizontal NPTS profiles

      lab_do_iplon : do iplon = 1, NPTS

        if (sfemis(iplon) > eps .and. sfemis(iplon) <= 1.0) then  ! input surface emissivity
          do j = 1, NBANDS
            semiss(j) = sfemis(iplon)
          enddo
        else                                                      ! use default values
          do j = 1, NBANDS
            semiss(j) = semiss0(j)
          enddo
        endif

!  ---  prepare atmospheric profile for use in rrtm
!       the vertical index of internal array is from surface to top

        if (iflip == 0) then        ! input from toa to sfc

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tint(iplon,NLP1)
          tbound= tint(iplon,NLP1)

          do k = 1, NLAY
            k1 = NLP1 - k
            pavel(k)= pmid(iplon,k1)
            delp(k) = pint(iplon,k1+1) - pint(iplon,k1)
            tavel(k)= tmid(iplon,k1)
            tz(k)   = tint(iplon,k1)

!  ---  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qnm(iplon,k1)*amdw)                   ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qnm(iplon,k1))                        ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,                                      &
     &                     qnm(iplon,k1)*amdw/(1.0-qnm(iplon,k1)))      ! input specific humidity
            o3vmr (k)= max(f_zero,o3mr(iplon,k1)*amdo3)                 ! input mass mixing ratio
!test use   o3vmr (k)= max(f_zero,o3mr(iplon,k1))                       ! input vol mixing ratio

            tem0 = (1.0 - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2 * delp(k) / (tem1*tem0*(1.0 + h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) =                coldry(k)*h2ovmr(k)            ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr(iplon,k1,1))  ! co2
            colamt(k,3) =                coldry(k)*o3vmr(k)             ! o3
          enddo

!  ---  set aerosol optical properties

          if (iaerlw > 0) then
            do j = 1, NBANDS
              do k = 1, NLAY
                k1 = NLP1 - k
                tauaer(k,j) = aerosols(iplon,k1,j,1)                    &
     &                      * (1.0 - aerosols(iplon,k1,j,2))
              enddo
            enddo
          else
            tauaer(:,:) = f_zero
          endif

          if (iflagliq > 0) then   ! use prognostic cloud method
            do k = 1, NLAY
              k1 = NLP1 - k
              cldfrac(k)= clouds(iplon,k1,1)
              cwp1 (k)  = clouds(iplon,k1,2)
              rew1 (k)  = clouds(iplon,k1,3)
              cip1 (k)  = clouds(iplon,k1,4)
              rei1 (k)  = clouds(iplon,k1,5)
              cda1 (k)  = clouds(iplon,k1,6)
              cda2 (k)  = clouds(iplon,k1,7)
              cda3 (k)  = clouds(iplon,k1,8)
              cda4 (k)  = clouds(iplon,k1,9)
            enddo
          else                        ! use diagnostic cloud method
            do k = 1, NLAY
              k1 = NLP1 - k
              cldfrac(k)= clouds(iplon,k1,1)
              cda1(k)   = clouds(iplon,k1,2)
            enddo
          endif                       ! end if_iflagliq

          cldfrac(0)    = 1.0         ! padding value only
          cldfrac(NLP1) = f_zero      ! padding value only

        else                        ! input from sfc to toa

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tint(iplon,1)
          tbound= tint(iplon,1)

          do k = 1, NLAY
            pavel(k)= pmid(iplon,k)
            delp(k) = pint(iplon,k) - pint(iplon,k+1)
            tavel(k)= tmid(iplon,k)
            tz(k)   = tint(iplon,k+1)

!  ---  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qnm(iplon,k)*amdw)                    ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qnm(iplon,k))                         ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qnm(iplon,k)*amdw/(1.0-qnm(iplon,k))) ! input specific humidity
            o3vmr (k)= max(f_zero,o3mr(iplon,k)*amdo3)                  ! input mass mixing ratio
!test use   o3vmr (k)= max(f_zero,o3mr(iplon,k))                        ! input vol mixing ratio

            tem0 = (1.0 - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2 * delp(k) / (tem1*tem0*(1.0 + h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) =                coldry(k)*h2ovmr(k)           ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr(iplon,k,1))  ! co2
            colamt(k,3) =                coldry(k)*o3vmr(k)            ! o3
          enddo

!  ---  set aerosol optical properties

          if (iaerlw > 0) then
            do j = 1, NBANDS
              do k = 1, NLAY
                tauaer(k,j) = aerosols(iplon,k,j,1)                     &
     &                      * (1.0 - aerosols(iplon,k,j,2))
              enddo
            enddo
          else
            tauaer(:,:) = f_zero
          endif

          if (iflagliq > 0) then   ! use prognostic cloud method
            do k = 1, NLAY
              cldfrac(k)= clouds(iplon,k,1)
              cwp1 (k)  = clouds(iplon,k,2)
              rew1 (k)  = clouds(iplon,k,3)
              cip1 (k)  = clouds(iplon,k,4)
              rei1 (k)  = clouds(iplon,k,5)
              cda1 (k)  = clouds(iplon,k,6)
              cda2 (k)  = clouds(iplon,k,7)
              cda3 (k)  = clouds(iplon,k,8)
              cda4 (k)  = clouds(iplon,k,9)
            enddo
          else                        ! use diagnostic cloud method
            do k = 1, NLAY
              cldfrac(k)= clouds(iplon,k,1)
              cda1(k)   = clouds(iplon,k,2)
            enddo
          endif

          cldfrac(0)    = 1.0         ! padding value only
          cldfrac(NLP1) = f_zero      ! padding value only

        endif                       ! if_iflip

!  ---  set up col amount for rare gases, convert from volume mixing ratio to
!       molec/cm2 based on coldry (scaled to 1.0e-20) for use in rrtm

        if (iflip == 0) then        ! input from toa to sfc

          if (irgaslw == 1) then
            do k = 1, NLAY
              k1 = NLP1 - k
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr(iplon,k1,2))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr(iplon,k1,3))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr(iplon,k1,4))  ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr(iplon,k1,5))  ! co
            enddo
          else
            do k = 1, NLAY
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co
            enddo
          endif

          if (icfclw == 1) then
            do k = 1, NLAY
              k1 = NLP1 - k
              wx(k,1) = max( f_zero, coldry(k)*gasvmr(iplon,k1,9) )   ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr(iplon,k1,6) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr(iplon,k1,7) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr(iplon,k1,8) )   ! cf22
            enddo
          else
            wx(:,:) = f_zero
          endif

        else                        ! input from sfc to toa

          if (irgaslw == 1) then
            do k = 1, NLAY
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr(iplon,k,2))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr(iplon,k,3))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr(iplon,k,4))  ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr(iplon,k,5))  ! co
            enddo
          else
            do k = 1, NLAY
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co
            enddo
          endif

          if (icfclw == 1) then
            do k = 1, NLAY
              wx(k,1) = max( f_zero, coldry(k)*gasvmr(iplon,k,9) )   ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr(iplon,k,6) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr(iplon,k,7) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr(iplon,k,8) )   ! cf22
            enddo
          else
            wx(:,:) = f_zero
          endif

        endif                       ! if_iflip

!  ---  get column density for broadening gases (molecules/cm**2)

        do k = 1, NLAY
          colbrd(k) = coldry(k) - colamt(k,2)-colamt(k,3)-colamt(k,4)   &
     &              - colamt(k,5)-colamt(k,6)-colamt(k,7)
        enddo

!     if (lprnt) then
!     print *,'  coldry',coldry
!     print *,' wx(*,1) ',(wx(k,1),k=1,NLAY)
!     print *,' wx(*,2) ',(wx(k,2),k=1,NLAY)
!     print *,' wx(*,3) ',(wx(k,3),k=1,NLAY)
!     print *,' wx(*,4) ',(wx(k,4),k=1,NLAY)
!     print *,' iplon ',iplon
!     print *,'  pavel ',pavel
!     print *,'  delp ',delp
!     print *,'  tavel ',tavel
!     print *,'  tz ',tz
!     print *,' h2ovmr ',h2ovmr
!     print *,' o3vmr ',o3vmr
!     endif

!  ---  calculate cloud optical properties

        call cldprop                                                    &
!  ---  inputs:
     &     ( cldfrac, cwp1, cip1, rew1, rei1, cda1, cda2, cda3, cda4,   &
     &       NLAY,                                                      &
!  ---  output:
     &       taucloud, icldatm                                          &
     &     )

!     if (lprnt) then
!     print *,' after cldprop'
!     print *,' cwp1',cwp1
!     print *,' cip1',cip1
!     print *,' rew1',rew1
!     print *,' rei1',rei1
!     print *,' taucl',cda1
!     print *,' cldfrac',cldfrac
!     print *,' taucloud',taucloud
!     endif

        call setcoef                                                    &
!  ---  inputs:
     &     ( pavel,tavel,tz,tbound,semiss,                              &
     &       h2ovmr,colamt,colbrd,coldry,                               &
     &       NLAY,                                                      &
!  ---  outputs:
     &       laytrop,fac00,fac01,fac10,fac11,jp,jt,jt1,                 &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd            &
     &     )


!     if (lprnt) then
!      print *,'laytrop',laytrop
!      print *,'colh2o',(colamt(k,1),k=1,NLAY)
!      print *,'colco2',(colamt(k,2),k=1,NLAY)
!      print *,'colo3', (colamt(k,3),k=1,NLAY)
!      print *,'coln2o',(colamt(k,4),k=1,NLAY)
!      print *,'colch4',(colamt(k,5),k=1,NLAY)
!      print *,'fac00',fac00
!      print *,'fac01',fac01
!      print *,'fac10',fac10
!      print *,'fac11',fac11
!      print *,'jp',jp
!      print *,'jt',jt
!      print *,'jt1',jt1
!      print *,'selffac',selffac
!      print *,'selffrac',selffrac
!      print *,'indself',indself
!      print *,'forfac',forfac
!     endif

!  ---  call the radiative transfer routine.

        if (numangs == 0) then

          call rtr                                                      &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       tauaer, NLAY, icldatm,                                     &
!  ---  outputs:
     &       totuclfl,totdclfl,htrcl,htrb                               &
     &     )

          if (icldatm == 1) then
            if (iovr == 0) then

              call rtrcld                                               &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr,htrb                                 &
     &     )

            else

              call rtrcldmr                                             &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr,htrb                                 &
     &     )

            endif

          else

            do k = 0, NLAY
              totuflux(k) = totuclfl(k)
              totdflux(k) = totdclfl(k)
              htr     (k) = htrcl   (k)
            enddo

          endif

        else

          call rtreg                                                    &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       tauaer, NLAY, icldatm,                                     &
!  ---  outputs:
     &       totuclfl,totdclfl,htrcl,htrb                               &
     &     )

          if (icldatm == 1) then

            if (iovr == 0) then

              call rtregcld                                             &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr,htrb                                 &
     &     )

            else

              call rtregcldmr                                           &
!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr,htrb                                 &
     &     )

            endif

          else

            do k = 0, NLAY
              totuflux(k) = totuclfl(k)
              totdflux(k) = totdclfl(k)
              htr     (k) = htrcl   (k)
            enddo

          endif

        endif


!  ---  output total-sky and clear-sky fluxes and heating rates.

        topflx(iplon)%upfxc = totuflux(NLAY)
        topflx(iplon)%upfx0 = totuclfl(NLAY)

        sfcflx(iplon)%upfxc = totuflux(0)
        sfcflx(iplon)%dnfxc = totdflux(0)
        sfcflx(iplon)%dnfx0 = totdclfl(0)

        if (iflip == 0) then        ! output from toa to sfc

          do k = 0, NLAY-1
            k1 = NLAY - k
            hlwc(iplon,k1) = htr(k)
          enddo

!! ---  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 0, NLAY-1
              k1 = NLAY - k
              hlw0(iplon,k1) = htrcl(k)
            enddo
          endif

!! ---  optional spectral band heating rate
          if ( lhlwb ) then
            do k = 1, NLAY
              k1 = NLP1 - k
              do j = 1, NBANDS
                hlwb(iplon,k1,j) = htrb(k,j)
              enddo
            enddo
          endif

!! ---  optional fluxes
          if ( lflxprf ) then
            do k = 0, NLAY
              k1 = NLP1 - k
              flxprf(iplon,k1)%upfxc = totuflux(k)
              flxprf(iplon,k1)%dnfxc = totdflux(k)
              flxprf(iplon,k1)%upfx0 = totuclfl(k)
              flxprf(iplon,k1)%dnfx0 = totdclfl(k)
            enddo
          endif

        else                        ! output from sfc to toa

          do k = 0, NLAY-1
            hlwc(iplon,k+1) = htr(k)
          enddo

!! ---  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 0, NLAY-1
              hlw0(iplon,k+1) = htrcl(k)
            enddo
          endif

!! ---  optional spectral band heating rate
          if ( lhlwb ) then
            do k = 1, NLAY
              do j = 1, NBANDS
                hlwb(iplon,k,j) = htrb(k,j)
              enddo
            enddo
          endif

!! ---  optional fluxes
          if ( lflxprf ) then
            do k = 0, NLAY
              flxprf(iplon,k+1)%upfxc = totuflux(k)
              flxprf(iplon,k+1)%dnfxc = totdflux(k)
              flxprf(iplon,k+1)%upfx0 = totuclfl(k)
              flxprf(iplon,k+1)%dnfx0 = totdclfl(k)
            enddo
          endif

        endif                       ! if_iflip

      enddo  lab_do_iplon

      return
!...................................
      end subroutine lwrad
!-----------------------------------



!-----------------------------------
      subroutine rlwinit                                                &
!...................................

!  ---  inputs:
     &     ( icwp, me, NLAY )
!  ---  outputs: (none)

!  *******************************************************************  !
!                                                                       !
!  rrtm longwave radiative transfer model                               !
!  atmospheric and environmental research, inc., cambridge, ma          !
!                                                                       !
!                                                                       !
!  original version:       michael j. iacono; july, 1998                !
!  revision for ncar ccm:  michael j. iacono; september, 1998           !
!                                                                       !
!  this subroutine performs calculations necessary for the initialization
!  of the lw model, rrtm.  lookup tables are computed for use in the lw !
!  radiative transfer, and input absorption coefficient data for each   !
!  spectral band are reduced from 256 g-points to 140 for use in rrtm.  !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!     arrays for 5000-point look-up tables:                             !
!     tautb- clear-sky optical depth (used in cloudy radiative transfer)!
!     tf     tautb transition function; i.e. the transition of the planck
!            function from that for the mean layer temperature to that  !
!            for the layer boundary temperature as a function of optical!
!            depth. the "linear in tau" method is used to make the table!
!     trans- transmittance                                              !
!                                                                       !
!  *******************************************************************  !
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
!  control flags in module "module_radlw_cntr_para":                    !
!     ilwrate - heating rate unit selections                            !
!               =1: output in k/day                                     !
!               =2: output in k/second                                  !
!     iaerlw  - control flag for aerosols                               !
!               =0: do not include aerosol effect                       !
!               >0: include aerosol effect                              !
!     irgaslw - control flag for rare gases (ch4,n2o,o2, etc.)          !
!               =0: do not include rare gases                           !
!               =1: include all rare gases                              !
!     icfclw  - control flag for cfc gases                              !
!               =0: do not include cfc gases                            !
!               =1: include all cfc gases                               !
!     iflagliq- cloud optical properties contrl flag                    !
!               =0: input cloud opt depth from diagnostic scheme        !
!               >0: input cwp,cip, and other cloud content parameters   !
!                                                                       !
!  *******************************************************************  !
!
      implicit none
!
!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys) :: tfn, fp, rtfp, pival, explimit
      integer               :: i
!
!===> ... begin here
!

      if (me == 0) then
        print *,' - Using AER Longwave Radiation, Version: ', VTAGLW

        if (iaerlw > 0) then
          print *,'   --- Using input aerosol parameters for LW'
        else
          print *,'   --- Aerosol effect is NOT included in LW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        endif

        if (irgaslw == 1) then
          print *,'   --- Include rare gases N2O, CH4, O2, absorptions',&
     &            ' in LW'
        else
          print *,'   --- Rare gases effect is NOT included in LW'
        endif

        if (icfclw == 1) then
          print *,'   --- Include CFC gases absorptions in LW'
        else
          print *,'   --- CFC gases effect is NOT included in LW'
        endif

        print *,'   --- number of angles used in integration:', numangs
      endif

!  --- ...  check cloud flags for consistency

      if ((icwp == 0 .and. iflagliq /= 0) .or.                          &
     &    (icwp == 1 .and. iflagliq == 0)) then
        print *, ' *** Model cloud scheme inconsistent with LW',        &
     &           ' radiation cloud radiative property setup !!'
        stop
      endif

!  --- ...  setup default surface emissivity for each band here

      semiss0(:) = 1.0

!  --- ...  setup constant factors for flux and heating rate
!           the 1.0e-2 is to convert pressure from mb to N/m**2

      pival = 2.0*asin(1.0)
      fluxfac = pival * 2.0d4
!     fluxfac = 3.1415927410125732 * 2.0d4
!     fluxfac = 62831.85307179586                   ! = 2 * pi * 1.0e4

      if (ilwrate == 1) then
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
!       heatfac = 9.80665 * 864.0 / 1.0046e3
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!  --- ...  compute lookup tables for transmittance, tau transition
!           function, and clear sky tau (for the cloudy sky radiative
!           transfer).  tau is computed as a function of the tau
!           transition function, transmittance is calculated as a 
!           function of tau, and the tau transition function is 
!           calculated using the linear in tau formulation at values of
!           tau above 0.01.  tf is approximated as tau/6 for tau < 0.01.
!           all tables are computed at intervals of 0.001.  the inverse
!           of the constant used in the pade approximation to the tau
!           transition function is set to b.

      tautb(0) = f_zero
      tf   (0) = f_zero
      trans(0) = 1.0

      tautb(NTBL) = 1.e10
      tf   (NTBL) = 1.0
      trans(NTBL) = f_zero

      explimit = aint( -log(tiny(trans(0))) )

      do i = 1, NTBL-1
         tfn = real(i, kind_phys) / real(NTBL-i, kind_phys)
         tautb(i) = bpade * tfn
         if (tautb(i) >= explimit) then
           trans(i) = f_zero
         else
           trans(i) = exp(-tautb(i))
         endif

         if (tautb(i) < 0.06) then
            tf(i) = tautb(i) / 6.0
         else
            tf(i) = 1. - 2.*( (1./tautb(i)) - (trans(i)/(1.-trans(i))) )
         endif
      enddo

!...................................
      end subroutine rlwinit
!-----------------------------------



!-----------------------------------
      subroutine cldprop                                                &
!...................................

!  ---  inputs:
     &     ( cldfrac,cliqp,cicep,reliq,reice,cdat1,cdat2,cdat3,cdat4,   &
     &       NLAY,                                                      &
!  ---  output:
     &       taucloud,icldatm                                           &
     &     )

!  *******************************************************************  !
!                                                                       !
!    purpose:  compute the cloud optical depth(s) for each cloudy layer.!
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!     cldfrac - layer cloud fraction                               0:L+1!
!        - - -  for iflagliq > 0  (prognostic cloud sckeme)  - - -      !
!     cliqp   - layer cloud liquid water path  (g/m**2)              L  !
!     reliq   - effective radius for water cloud (micron)            L  !
!     cicep   - layer cloud ice water path  (g/m**2)                 L  !
!     reice   - effective radius for ice cloud (micron)              L  !
!     cdat1   - layer rain drop water path  (g/m**2)                 L  !
!     cdat2   - effective radius for rain drop (microm)              L  !
!     cdat3   - layer snow flake water path (g/m**2)                 L  !
!               (if use fu's formula it needs to be normalized by       !
!                snow density (g/m**3/1.0e6) to get unit of micron)     !
!     cdat4   - effective radius for snow flakes (micron)            L  !
!        - - -  for iflagliq = 0  (diagnostic cloud sckeme)  - - -      !
!     cdat1   - input cloud optical depth                            L  !
!     cdat2   - optional use                                         L  !
!     cdat3   - optional use                                         L  !
!     cdat4   - optional use                                         L  !
!     cliqp   - not used                                             L  !
!     reliq   - not used                                             L  !
!     cicep   - not used                                             L  !
!     reice   - not used                                             L  !
!                                                                       !
!     NLAY    - vertical layer/level numbers                         1  !
!                                                                       !
!    explanation of the method for each value of iflagliq, and iflagice.!
!    set up in module "module_radlw_cntr_para"                          !
!                                                                       !
!     iflagliq=0 and =1 do not distingish being liquid and ice clouds.  !
!     iflagliq=2 and =3 does distinguish between liquid and ice clouds, !
!                  and requires further user input (iflagice) to specify!
!                  the method to be used to compute the aborption due to!
!                  liquid and ice parts.                                !
!  ...................................................................  !
!                                                                       !
!     iflagliq=0:  for each cloudy layer, the cloud fraction and (gray) !
!                  optical depth are input.                             !
!     iflagliq=1:  for each cloudy layer, the cloud fraction and cloud  !
!                  water path (g/m2) are input.  using clp only. the    !
!                  (gray) cloud optical depth is computed as in ccm2.   !
!     iflagliq=2:  the optical depths due to water clouds are computed  !
!                  as in ccm3.                                          !
!     iflagliq=3:  the water droplet effective radius (microns) is input!
!                  and the opt depths due to water clouds are computed  !
!                  as in hu and stamnes, j., clim., 6, 728-742, (1993). !
!                  the values for absorption coefficients appropriate for
!                  the spectral bands in rrtm have been obtained for a  !
!                  range of effective radii by an averaging procedure   !
!                  based on the work of j. pinto (private communication).
!                  linear interpolation is used to get the absorption   !
!                  coefficients for the input effective radius.         !
!                                                                       !
!     iflagice=0:  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in ccm3.                      !
!     iflagice=1:  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in ebert and curry, jgr, 97,  !
!                  3831-3836 (1992).  the spectral regions in this work !
!                  have been matched with the spectral bands in rrtm to !
!                  as great an extent as possible:                      !
!                     e&c 1      ib = 5      rrtm bands 9-16            !
!                     e&c 2      ib = 4      rrtm bands 6-8             !
!                     e&c 3      ib = 3      rrtm bands 3-5             !
!                     e&c 4      ib = 2      rrtm band 2                !
!                     e&c 5      ib = 1      rrtm band 1                !
!     iflagice=2:  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in streamer (reference: j. key,
!                  streamer user's guide, technical report 96-01,       !
!                  department of geography, boston university, 85 pp.   !
!                  (1996)).  the values of absorption coefficients      !
!                  appropriate for the spectral bands of rrtm were      !
!                  obtained by an averaging procedure based on the work !
!                  of j. pinto (private communication).                 !
!                                                                       !
!  outputs:                                                             !
!     taucloud - cloud optical depth                        NBANDS*L    !
!                                                                       !
!  *******************************************************************  !
!
      use module_radlw_cldprlw

      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY

      real (kind=kind_phys), dimension(0:), intent(in) :: cldfrac

      real (kind=kind_phys), dimension(:),  intent(in) :: cliqp, cicep, &
     &       reliq, reice, cdat1, cdat2, cdat3, cdat4

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: taucloud

      integer, intent(out) :: icldatm

!  ---  locals:
      real (kind=kind_phys) :: cliq, cice, radliq, radice, factor, fint
      real (kind=kind_phys) :: taurain, tausnow
      integer               :: j, k, index

!
!===> ... begin here
!
      do j = 1, NBANDS
        do k = 1, NLAY
          taucloud(k,j) = f_zero
        enddo
      enddo

      icldatm = 0

      lab_do_k : do k = 1, NLAY

        lab_if_cld : if (cldfrac(k) > eps) then

          icldatm = 1

!  ---  ice clouds and water clouds combined.
          lab_if_liq : if (iflagliq == 0) then

            do j = 1, NBANDS
              taucloud(k,j) = cdat1(k)
            enddo

          elseif (iflagliq == 1) then  lab_if_liq

            taurain = absrain * cdat1(k)                 ! ncar formula
            tausnow = abssnow0 * cdat3(k)                ! ncar formula
!           tausnow = abssnow1 * cdat3(k) / cdat4(k)     ! fu's formula

!           taucloud(k,1) = absliq1 * (cliqp(k) + cicep(k))
!           taucloud(k,1) = absliq1 * cliqp(k)
            taucloud(k,1) = absliq1*cliqp(k) + taurain + tausnow
            do j = 2,NBANDS
              taucloud(k,j) = taucloud(k,1)
            enddo

!  ---  separate treatement of ice clouds and water clouds.
          else  lab_if_liq

            taurain = absrain * cdat1(k)                 ! ncar formula
            tausnow = abssnow0 * cdat3(k)                ! ncar formula
!           tausnow = abssnow1 * cdat3(k) / cdat4(k)     ! fu's formula

            cliq = cliqp(k)
            cice = cicep(k)
            radliq = reliq(k)
            radice = reice(k)

!  ---  calculation of absorption coefficients due to liquid clouds.
            if (cliq == f_zero) then
              do j = 1, NBANDS
                abscoliq(j) = f_zero
              enddo
            elseif (iflagliq == 2) then
              abscoliq(1) = cliq * absliq2
              do j = 2, NBANDS
                abscoliq(j) = abscoliq(1)
              enddo
            elseif (iflagliq == 3) then
              radliq = max(1.5e0, min(60.0e0, real(reliq(k)) ))
              factor = radliq - 1.5
              index = max(1, min(57, int(factor)))
              fint = radliq - 1.5 - index
              do j = 1, NBANDS
                abscoliq(j) = cliq * (absliq3(index,j) + fint *         &
     &                    (absliq3(index+1,j) - (absliq3(index,j))))
              enddo
            endif

!  ---  calculation of absorption coefficients due to ice clouds.
            if (cice == f_zero) then
              do j = 1, NBANDS
                abscoice(j) = f_zero
              enddo
            elseif (iflagice == 0) then
              radice = max(10.e0, real(reice(k)) )
              abscoice(1) = cice * (absice0(1) + absice0(2)/radice)
              do j = 2, NBANDS
                abscoice(j) = abscoice(1)
              enddo
            elseif (iflagice == 1) then
              radice = max(13.e0, min(130.e0, real(reice(k)) ))
              do j = 1, NBANDS
                index = ipat(j)
                abscoice(j) = cice * (absice1(1,index)                  &
     &                      + absice1(2,index)/radice)
              enddo
            elseif (iflagice == 2) then
              radice = max(13.e0, min(130.e0, real(reice(k)) ))
              factor = (radice - 10.0) / 3.0
              index = min(39, int(factor))
              fint = factor - index
              do j = 1, NBANDS
                abscoice(j) = cice * (absice2(index,j) + fint *         &
     &                    (absice2(index+1,j) - (absice2(index,j))))
              enddo
            elseif (iflagice == 3) then
              radice = max(10.0, min(140.0, real(reice(k)) ))
              factor = (radice - 5.0) / 5.0
              index = min(26, int(factor))
              fint = factor - index
              do j = 1, NBANDS
                abscoice(j) = cice * (absice3(index,j) + fint *         &
     &                    (absice3(index+1,j) - (absice3(index,j))))
              enddo
            endif

            do j = 1, NBANDS
!             taucloud(k,j) = abscoice(j) + abscoliq(j)
              taucloud(k,j) = abscoice(j) + abscoliq(j)                 &
     &                      + taurain + tausnow
            enddo

          endif  lab_if_liq

        endif  lab_if_cld

      enddo  lab_do_k

      return
!...................................
      end subroutine cldprop
!-----------------------------------



!-----------------------------------
      subroutine setcoef                                                &
!...................................

!  ---  inputs:
     &     ( pavel,tavel,tz,tbound,semiss,                              &
     &       h2ovmr,colamt,colbrd,coldry,                               &
     &       NLAY,                                                      &
!  ---  outputs:
     &       laytrop,fac00,fac01,fac10,fac11,jp,jt,jt1,                 &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd            &
     &     )

!  *******************************************************************  !
!                                                                       !
!   for a given atmosphere, this program calculates the indices and     !
!   fractions related to the pressure and temperature interpolations.   !
!   also calculate the values of the integrated planck functions for    !
!   each band at the level and layer temperatures.                      !
!                                                                       !
!     revision:  3.1       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in)  ::  NLAY

      real (kind=kind_phys), dimension(0:), intent(in) :: tz

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel,        &
     &       tavel, h2ovmr, colbrd, coldry, semiss

      real (kind=kind_phys), dimension(:,:),intent(in) :: colamt

      real (kind=kind_phys),                intent(in) :: tbound

!  ---  outputs:
      real (kind=kind_phys), dimension(:),  intent(out) :: fac00,       &
     &       fac01, fac10, fac11, selffac, selffrac, forfac, forfrac,   &
     &       minorfrac, scaleminor, scaleminorn2, plankbnd

      real (kind=kind_phys), dimension(:,:),intent(out) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       planklay

      real (kind=kind_phys), dimension(0:,:),intent(out) :: planklev

      integer, dimension(:), intent(out) :: jp, jt, jt1, indself,       &
     &       indfor, indminor

      integer,               intent(out) ::  laytrop

!  ---  locals:
      real (kind=kind_phys) :: tbndfrac, t0frac, tlayfrac, tlevfrac,    &
     &       dbdtlev, dbdtlay, plog, fp, compfp, factor, scalefac,      &
     &       ft, ft1, tem1, tem2

      integer  ::  j, k, iband, indbound, indlev0, indlay, indlev
      integer  ::  jp1
!
!===> ... begin here
!
      rat_h2oco2(:,:) = f_zero
      rat_h2oo3 (:,:) = f_zero
      rat_h2on2o(:,:) = f_zero
      rat_h2och4(:,:) = f_zero
      rat_n2oco2(:,:) = f_zero
      rat_o3co2 (:,:) = f_zero

      laytrop = 0

      indbound = max(1, min(180, int(tbound-159.0) ))
      indlev0  = max(1, min(180, int(tz(0) -159.0) ))

      tbndfrac = tbound - 159.0 - float(indbound)
      t0frac   = tz(0)  - 159.0 - float(indlev0)

      lab_do_k : do k = 1, NLAY

!  ---  calculate the integrated planck functions for each band at the
!       surface, level, and layer temperatures.

        indlay = max(1, min(180, int(tavel(k) - 159.0) ))
        indlev = max(1, min(180, int(tz   (k) - 159.0) ))

        tlayfrac = tavel(k) - 159.0 - float(indlay)
        tlevfrac = tz   (k) - 159.0 - float(indlev)

        do iband = 1, NBANDS-1
          if (k == 1) then
            dbdtlev = totplnk(indbound+1,iband)-totplnk(indbound,iband)
            plankbnd(iband) = semiss(iband)                             &
     &           * (totplnk(indbound,iband) + tbndfrac * dbdtlev)

            dbdtlev = totplnk(indlev0+1,iband)-totplnk(indlev0,iband)
            planklev(0,iband) = totplnk(indlev0,iband) + t0frac*dbdtlev
          endif

          dbdtlev = totplnk(indlev+1,iband) - totplnk(indlev,iband)
          dbdtlay = totplnk(indlay+1,iband) - totplnk(indlay,iband)
          planklay(k,iband) = totplnk(indlay,iband) + tlayfrac*dbdtlay
          planklev(k,iband) = totplnk(indlev,iband) + tlevfrac*dbdtlev
        enddo

!  ---  for band 16, if radiative transfer will be performed on just
!       this band, use integrated planck values up to 3250 cm-1.
!       if radiative transfer will be performed across all 16 bands,
!       then include in the integrated planck values for this band
!       contributions from 2600 cm-1 to infinity.

!       if (istart .eq. 16) then
!         if (k == 1) then
!           dbdtlev = totplk16(indbound+1) - totplk16(indbound)
!           plankbnd(NBANDS) = semiss(NBANDS)                           &
!    &           * (totplk16(indbound) + tbndfrac*dbdtlev)
!
!           dbdtlev = totplnk(indlev0+1,NBANDS)-totplnk(indlev0,NBANDS)
!           planklev(0,NBANDS) = totplk16(indlev0) + t0frac*dbdtlev
!         endif
!
!         dbdtlev = totplk16(indlev+1) - totplk16(indlev)
!         dbdtlay = totplk16(indlay+1) - totplk16(indlay)
!         planklay(k,NBANDS) = totplk16(indlay) + tlayfrac*dbdtlay
!         planklev(k,NBANDS) = totplk16(indlev) + tlevfrac*dbdtlev
!       else
          if (k == 1) then
            dbdtlev =totplnk(indbound+1,NBANDS)-totplnk(indbound,NBANDS)
            plankbnd(NBANDS) = semiss(NBANDS)                           &
     &           * (totplnk(indbound,NBANDS) + tbndfrac*dbdtlev)

            dbdtlev = totplnk(indlev0+1,NBANDS)-totplnk(indlev0,NBANDS)
            planklev(0,NBANDS) = totplnk(indlev0,NBANDS)+t0frac*dbdtlev
          endif

          dbdtlev = totplnk(indlev+1,NBANDS) - totplnk(indlev,NBANDS)
          dbdtlay = totplnk(indlay+1,NBANDS) - totplnk(indlay,NBANDS)
          planklay(k,NBANDS) = totplnk(indlay,NBANDS)+tlayfrac*dbdtlay
          planklev(k,NBANDS) = totplnk(indlev,NBANDS)+tlevfrac*dbdtlev
!       endif

!  ---  find the two reference pressures on either side of the layer
!       pressure.  store them in jp and jp1.  store in fp the fraction
!       of the difference (in ln(pressure)) between these two values
!       that the layer pressure lies.

        plog = log(pavel(k))
        jp(k) = max(1, min(58, int(36.0 - 5.0*(plog + 0.04)) ))
        jp1 = jp(k) + 1
!  ---  limit pressure extrapolation at the top
        fp  = max(f_zero, min(1.0, 5.0*(preflog(jp(k))-plog) ))
!org    fp  = 5.0 * (preflog(jp(k)) - plog)

!  ---  determine, for each reference pressure (jp and jp1), which
!       reference temperature (these are different for each reference
!       pressure) is nearest the layer temperature but does not exceed
!       it.  store these indices in jt and jt1, resp. store in ft (resp.
!       ft1) the fraction of the way between jt (jt1) and the next
!       highest reference temperature that the layer temperature falls.

        tem1 = (tavel(k) - tref(jp(k))) / 15.0
        tem2 = (tavel(k) - tref(jp1  )) / 15.0
        jt (k) = max(1, min(4, int(3.0 + tem1) ))
        jt1(k) = max(1, min(4, int(3.0 + tem2) ))
!  ---  restrict extrapolation ranges by limiting abs(det t) < 37.5 deg
        ft  = max(-0.5, min(1.5, tem1 - float(jt (k) - 3) ))
        ft1 = max(-0.5, min(1.5, tem2 - float(jt1(k) - 3) ))
!org    ft  = tem1 - float(jt (k) - 3)
!org    ft1 = tem2 - float(jt1(k) - 3)

!  ---  set up factors needed to separately include the minor gases
!       in the calculation of absorption coefficient

        scaleminor(k) = pavel(k) / tavel(k)
        scaleminorn2(k) = (pavel(k) / tavel(k))                         &
     &                  * (colbrd(k) / (coldry(k) + colamt(k,1)))
        factor = (tavel(k) - 180.8) / 7.2
        indminor(k) = min(18, max(1, int(factor)))
        minorfrac(k) = factor - float(indminor(k))

!  ---  if the pressure is less than ~100mb, perform a different
!       set of species interpolations.

        scalefac = pavel(k) * stpfac / tavel(k)
        if (plog > 4.56) then
          laytrop =  laytrop + 1

          forfac(k) = scalefac / (1.0 + h2ovmr(k))
          factor = (332.0 - tavel(k)) / 36.0
          indfor(k) = min(2, max(1, int(factor)))
          forfrac(k) = factor - float(indfor(k))

!  ---  set up factors needed to separately include the water vapor
!       self-continuum in the calculation of absorption coefficient.

          selffac(k) = h2ovmr(k) * forfac(k)
          factor = (tavel(k) - 188.0) / 7.2
          indself(k) = min(9, max(1, int(factor)-7))
          selffrac(k) = factor - float(indself(k) + 7)

!  ---  setup reference ratio to be used in calculation of binary
!       species parameter in lower atmosphere.

          rat_h2oco2(k,1) = chi_mls(1,jp(k)  ) / chi_mls(2,jp(k)  )
          rat_h2oco2(k,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rat_h2oo3(k,1) = chi_mls(1,jp(k)  ) / chi_mls(3,jp(k)  )
          rat_h2oo3(k,2) = chi_mls(1,jp(k)+1) / chi_mls(3,jp(k)+1)

          rat_h2on2o(k,1) = chi_mls(1,jp(k)  ) / chi_mls(4,jp(k)  )
          rat_h2on2o(k,2) = chi_mls(1,jp(k)+1) / chi_mls(4,jp(k)+1)

          rat_h2och4(k,1) = chi_mls(1,jp(k)  ) / chi_mls(6,jp(k)  )
          rat_h2och4(k,2) = chi_mls(1,jp(k)+1) / chi_mls(6,jp(k)+1)

          rat_n2oco2(k,1) = chi_mls(4,jp(k)  ) / chi_mls(2,jp(k)  )
          rat_n2oco2(k,2) = chi_mls(4,jp(k)+1) / chi_mls(2,jp(k)+1)

!  ---  above laytrop.
        else

          forfac(k) = scalefac / (1.0 + h2ovmr(k))
          factor = (tavel(k) - 188.0) / 36.0
          indfor(k) = 3
          forfrac(k) = factor - 1.0

          selffac(k)  = f_zero
          selffrac(k) = f_zero
          indself(k)  = 0

!  ---  setup reference ratio to be used in calculation of binary
!       species parameter in upper atmosphere.

          rat_h2oco2(k,1) = chi_mls(1,jp(k)  ) / chi_mls(2,jp(k)  )
          rat_h2oco2(k,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rat_o3co2(k,1) = chi_mls(3,jp(k)  ) / chi_mls(2,jp(k)  )
          rat_o3co2(k,2) = chi_mls(3,jp(k)+1) / chi_mls(2,jp(k)+1)
        endif

!  ---  we have now isolated the layer ln pressure and temperature,
!       between two reference pressures and two reference temperatures
!       (for each reference pressure).  we multiply the pressure
!       fraction fp with the appropriate temperature fractions to get
!       the factors that will be needed for the interpolation that yields
!       the optical depths (performed in routines taugbn for band n).`

        compfp = 1.0 - fp
        fac10(k) = compfp * ft
        fac00(k) = compfp * (1.0 - ft)
        fac11(k) = fp * ft1
        fac01(k) = fp * (1.0 - ft1)

!  ---  rescale selffac and forfac for use in taumol

        selffac(k) = colamt(k,1)*selffac(k)
        forfac(k) = colamt(k,1)*forfac(k)
      enddo  lab_do_k

      return
!...................................
      end subroutine setcoef
!-----------------------------------



!-----------------------------------
      subroutine rtr                                                    &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       tauaer, NLAY, icldatm,                                     &
!  ---  outputs:
     &       totuclfl,totdclfl,htrcl, htrb                              &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary atmosphere.  the input to this       !
!   program is the atmospheric profile and all planck function          !
!   information.  only one angle is used from standard gaussian         !
!   quadrature.                                                         !
!                                                                       !
!     revision:  3.1       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, laytrop, ireflect, icldatm
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel, delp,  &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       plankbnd, semiss

      real (kind=kind_phys), dimension(:,:),intent(in) :: rat_h2oco2,   &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htrcl,       &
     &       totuclfl, totdclfl
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX) :: taug, fracs

      real (kind=kind_phys), dimension(0:NLAY)    :: urad1, drad1,      &
     &       dflux, uflux, fnetc

      real (kind=kind_phys), dimension(NLAY)      :: bbu1, atrans1

      real (kind=kind_phys) :: rad, rad0, radlu1, radld1, blay, bbd1,   &
     &       dplankup, dplankdn, odepth1, tblind, trans1, tausfac1,     &
     &       plfrac, factot, reflect

      integer :: j, k, ig, iband, itr1

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuclfl(0:NLAY)     ! clear sky upward longwave flux (w/m2)       !
!    totdclfl(0:NLAY)     ! clear sky downward longwave flux (w/m2)     !
!    htrcl   (0:NLAY)     ! clear sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnetc   (0:NLAY)     ! clear sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!

      do k = 0, NLAY
         urad1(k) = f_zero
         drad1(k) = f_zero
         totuclfl(k) = f_zero
         totdclfl(k) = f_zero
      enddo

!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)

          radld1 = f_zero

!  ---  downward radiative transfer.  due to the simple form taken by
!       certain equations when the optical depth is small, these
!       conditions are tested for.

          do k = NLAY, 1, -1
            blay = planklay(k,iband)
            plfrac = fracs(k,ig)
            dplankup = planklev(k,iband) - blay
            dplankdn = planklev(k-1,iband) - blay
            odepth1 = secdiff * taug(k,ig)

            if (odepth1 <= 0.06) then
              if (odepth1 < f_zero) odepth1 = f_zero
              atrans1(k) = odepth1 - 0.5*odepth1*odepth1
              odepth1 = rec_6 * odepth1
              bbd1 =  plfrac * (blay + dplankdn*odepth1)
              radld1 = radld1 + (bbd1 - radld1)*atrans1(k)
              bbu1(k) =  plfrac * (blay + dplankup*odepth1)
            else
              tblind = odepth1 / (bpade + odepth1)
              itr1 = tblint*tblind + 0.5
              trans1 = trans(itr1)
              atrans1(k) = 1.0 - trans1
              tausfac1 = tf(itr1)
              bbd1 = plfrac * (blay + tausfac1*dplankdn)
              radld1 = radld1 + (bbd1 - radld1)*atrans1(k)
              bbu1(k) = plfrac * (blay + tausfac1*dplankup)
            endif

            drad1(k-1) = drad1(k-1) + radld1
          enddo


!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)

!  ---  add in reflection of surface downward radiance.
          reflect = 1.0 - semiss(iband)
          if (ireflect == 1) then         !  specular reflection.
            radlu1 = rad0 + reflect * radld1
          else                            !  lambertian reflection.
            rad = 2.0 * radld1 * wtdiff
            radlu1 = rad0 + reflect * rad
          endif

          urad1(0) = urad1(0) + radlu1
          do k = 1, NLAY
            radlu1 = radlu1 + (bbu1(k) - radlu1)*atrans1(k)
            urad1(k) = urad1(k) + radlu1
          enddo

        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

!       do k = NLAY, 0, -1
        do k = 0, NLAY
          uflux(k) = urad1(k) * wtdiff
          dflux(k) = drad1(k) * wtdiff
          urad1(k) = f_zero
          drad1(k) = f_zero
          totuclfl(k) = totuclfl(k) + uflux(k) * delwave(iband)
          totdclfl(k) = totdclfl(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb .and. icldatm==0 ) then
          do k = 0, NLAY
            fnetc(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnetc(k-1) - fnetc(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!  ---  convert radiances to fluxes and heating rates for total sky.
!       calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuclfl(k) = totuclfl(k) * fluxfac
        totdclfl(k) = totdclfl(k) * fluxfac

        fnetc(k) = totuclfl(k) - totdclfl(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htrcl(j) = heatfac * (fnetc(j) - fnetc(k)) / delp(k)
      enddo

      htrcl(NLAY) = f_zero

      return
!...................................
      end subroutine rtr
!-----------------------------------



!-----------------------------------
      subroutine rtrcld                                                 &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr, htrb                                &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary cloudy atmosphere.  the input to     !
!   this program is the atmospheric profile, including cloud properties,!
!   and all planck function information.  only one angle is used from   !
!   standard gaussian quadrature.                                       !
!                                                                       !
!     revision:  3.2       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer,               intent(in) :: NLAY, laytrop, ireflect
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(0:), intent(in) :: cldfrac

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel, delp,  &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       semiss, plankbnd

      real (kind=kind_phys), dimension(:,:), intent(in) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, taucloud, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htr,         &
     &       totuflux, totdflux
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX) :: taug, fracs

      real (kind=kind_phys), dimension(0:NLAY)    :: urad, drad,        &
     &       dflux, uflux, fnet

      real (kind=kind_phys), dimension(NLAY)      :: bbugas, bbutot,    &
     &       atrans, atot

      real (kind=kind_phys), dimension(NLAY,NBANDS):: odcld, abscld,    &
     &       efclfrac

      real (kind=kind_phys) :: rad, rad0, radlu, radld, blay, bbd,      &
     &       bbdtot, dplankup, dplankdn, odepth, odtot, odepth_rec,     &
     &       odtot_rec, gassrc, tblind, transc, transcld, tausfac,      &
     &       tfactot, tfacgas, plfrac, reflect

      integer :: j, k, ig, iband, icldlyr(NLAY), ittot, itgas, itr

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    cldfrac (0:NLP1)     ! layer cloud fraction (padded at 2 ends)     !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuflux(0:NLAY)     ! total sky upward longwave flux (w/m2)       !
!    totdflux(0:NLAY)     ! total sky downward longwave flux (w/m2)     !
!    htr     (0:NLAY)     ! total sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnet    (0:NLAY)     ! total sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!

      do k = 0, NLAY
         urad(k) = f_zero
         drad(k) = f_zero
         totuflux(k) = f_zero
         totdflux(k) = f_zero
      enddo

      do k = 1, NLAY
        do iband = 1, NBANDS
          if (cldfrac(k) >= 1.e-6) then
            odcld(k,iband) = secdiff * taucloud(k,iband)
            transcld = exp(-odcld(k,iband))
            abscld(k,iband) = 1.0 - transcld
            efclfrac(k,iband) = abscld(k,iband) * cldfrac(k)
            icldlyr(k) = 1
          else
            odcld(k,iband) = f_zero
            abscld(k,iband) = f_zero
            efclfrac(k,iband) = f_zero
            icldlyr(k) = 0
          endif
        enddo
      enddo


!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)
          radld = f_zero

!  ---  downward radiative transfer loop.
!           here are some variable definitions:
!             odtot      optical depth of gas and cloud
!             atrans     absorptivity for gas only
!             atot       absorptivity for gas and cloud
!             tfacgas    gas-only pade factor, used for planck fn
!             tfactot    gas and cloud pade factor, used for planck fn
!             bbdgas     gas-only planck function for downward rt
!             bbdtot     gas and cloud planck function for downward rt
!             bbutot     gas and cloud planck function for upward calc.
!             gassrc     source radiance due to gas only

          do k = NLAY, 1, -1
            plfrac = fracs(k,ig)
            blay = planklay(k,iband)
            dplankup = planklev(k,iband) - blay
            dplankdn = planklev(k-1,iband) - blay
            odepth = secdiff * taug(k,ig)
            if (odepth < f_zero) odepth = f_zero

            if (icldlyr(k) == 1) then
              odtot = odepth + odcld(k,iband)

              if (odtot < 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth_rec = rec_6 * odepth
                gassrc = plfrac*(blay + dplankdn*odepth_rec)*atrans(k)

                atot(k) =  odtot - 0.5*odtot*odtot
                odtot_rec = rec_6 * odtot
                bbdtot =  plfrac * (blay + dplankdn*odtot_rec)
                radld = radld - radld * ( atrans(k)                     &
     &                + efclfrac(k,iband)*(1.0 - atrans(k)) ) + gassrc  &
     &                + cldfrac(k) * (bbdtot*atot(k) - gassrc)
                drad(k-1) = drad(k-1) + radld

                bbugas(k) =  plfrac * (blay + dplankup*odepth_rec)
                bbutot(k) =  plfrac * (blay + dplankup*odtot_rec)

              elseif (odepth <= 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth_rec = rec_6 * odepth
                gassrc = plfrac*(blay + dplankdn*odepth_rec)*atrans(k)

                odtot = odepth + odcld(k,iband)
                tblind = odtot / (bpade+odtot)
                ittot = tblint*tblind + 0.5
                tfactot = tf(ittot)
                bbdtot = plfrac * (blay + tfactot*dplankdn)
                atot(k) = 1.0 - trans(ittot)

                radld = radld - radld * ( atrans(k)                     &
     &                + efclfrac(k,iband)*(1.0 - atrans(k)) ) + gassrc  &
     &                + cldfrac(k)*(bbdtot*atot(k) - gassrc)
                drad(k-1) = drad(k-1) + radld

                bbugas(k) = plfrac * (blay + dplankup*odepth_rec)
                bbutot(k) = plfrac * (blay + tfactot*dplankup)
              else
                tblind = odepth / (bpade + odepth)
                itgas = tblint*tblind + 0.5
                odepth = tautb(itgas)
                atrans(k) = 1.0 - trans(itgas)
                tfacgas = tf(itgas)
                gassrc = atrans(k)*plfrac*(blay + tfacgas*dplankdn)

                odtot = odepth + odcld(k,iband)
                tblind = odtot / (bpade + odtot)
                ittot = tblint*tblind + 0.5
                tfactot = tf(ittot)
                bbdtot = plfrac * (blay + tfactot*dplankdn)
                atot(k) = 1.0 - trans(ittot)

                radld = radld - radld * ( atrans(k)                     &
     &                + efclfrac(k,iband)*(1.0 - atrans(k)) ) + gassrc  &
     &                + cldfrac(k)*(bbdtot*atot(k) - gassrc)
                drad(k-1) = drad(k-1) + radld

                bbugas(k) = plfrac * (blay + tfacgas*dplankup)
                bbutot(k) = plfrac * (blay + tfactot*dplankup)
              endif
            else
              if (odepth <= 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth = rec_6 * odepth
                bbd       =  plfrac * (blay + dplankdn*odepth)
                bbugas(k) =  plfrac * (blay + dplankup*odepth)
              else
                tblind = odepth / (bpade + odepth)
                itr = tblint*tblind + 0.5
                transc = trans(itr)
                atrans(k) = 1.0 - transc
                tausfac = tf(itr)
                bbd       = plfrac * (blay + tausfac*dplankdn)
                bbugas(k) = plfrac * (blay + tausfac*dplankup)
              endif

              radld = radld + (bbd - radld)*atrans(k)
              drad(k-1) = drad(k-1) + radld
            endif
          enddo

!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)

!  ---  add in reflection of surface downward radiance.
          reflect = 1.0 - semiss(iband)
          if (ireflect == 1) then         !  specular reflection.
            radlu = rad0 + reflect*radld
          else                            !  lambertian reflection.
            rad = 2.0 * radld * wtdiff
            radlu = rad0 + reflect*rad
          endif

          urad(0) = urad(0) + radlu
          do k = 1, NLAY
            if (icldlyr(k) == 1) then
              gassrc = bbugas(k) * atrans(k)
              radlu = radlu - radlu*( atrans(k)                         &
     &              + efclfrac(k,iband)*(1.0 - atrans(k)) ) + gassrc    &
     &              + cldfrac(k)*(bbutot(k)*atot(k) - gassrc)
              urad(k) = urad(k) + radlu
            else
              radlu = radlu + (bbugas(k) - radlu)*atrans(k)
              urad(k) = urad(k) + radlu
            endif
          enddo

        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

        do k = NLAY, 0, -1
          uflux(k) = urad(k) * wtdiff
          dflux(k) = drad(k) * wtdiff
          urad(k)  = f_zero
          drad(k)  = f_zero
          totuflux(k) = totuflux(k) + uflux(k) * delwave(iband)
          totdflux(k) = totdflux(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb ) then
          do k = 0, NLAY
            fnet(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnet(k-1) - fnet(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!  ---  convert radiances to fluxes and heating rates for total sky.
!       calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuflux(k) = totuflux(k) * fluxfac
        totdflux(k) = totdflux(k) * fluxfac

        fnet(k) = totuflux(k) - totdflux(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htr(j) = heatfac * (fnet(j) - fnet(k)) / delp(k)
      enddo

      htr(NLAY) = f_zero

      return
!...................................
      end subroutine rtrcld
!-----------------------------------



!-----------------------------------
      subroutine rtrcldmr                                               &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr, htrb                                &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary cloudy atmosphere.  the input to     !
!   this program is the atmospheric profile, including cloud properties,!
!   and all planck function information.  only one angle is used from   !
!   standard gaussian quadrature.                                       !
!                                                                       !
!     revision:  3.2       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer,               intent(in) :: NLAY, laytrop, ireflect
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(0:), intent(in) :: cldfrac

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel, delp,  &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       semiss, plankbnd

      real (kind=kind_phys), dimension(:,:),intent(in) :: rat_h2oco2,   &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, taucloud, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htr,         &
     &       totuflux, totdflux
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX) :: taug, fracs

      real (kind=kind_phys), dimension(0:NLAY)    :: urad, drad,        &
     &       dflux, uflux, fnet

      real (kind=kind_phys), dimension(NLAY)      :: bbugas, bbutot,    &
     &       atrans, atot

      real (kind=kind_phys), dimension(NLAY,NBANDS):: odcld, abscld

      real (kind=kind_phys), dimension(NLAY+1)     :: faccld1, faccld2, &
     &       facclr1, facclr2, faccmb1, faccmb2

      real (kind=kind_phys), dimension(0:NLAY)     :: faccld1d,         &
     &       faccld2d, facclr1d, facclr2d, faccmb1d, faccmb2d

      integer :: icldlyr(NLAY), istcld(NLAY+1), istcldd(0:NLAY)

      real (kind=kind_phys) :: rad, rad0, radlu, radld, blay, bbd,      &
     &       bbdtot, dplankup, dplankdn, odepth, odtot, odepth_rec,     &
     &       odtot_rec, gassrc, tblind, transc, transcld, tausfac,      &
     &       tfactot, tfacgas, plfrac, reflect, cldradd, clrradd,       &
     &       cldradu, clrradu, cldsrc, radmod, oldcld, oldclr, fmax,    &
     &       rat1, rat2, fmin, ttot

      integer :: j, k, ig, iband, ittot, itgas, itr, iclddn

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    cldfrac (0:NLP1)     ! layer cloud fraction (padded at 2 ends)     !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuflux(0:NLAY)     ! total sky upward longwave flux (w/m2)       !
!    totdflux(0:NLAY)     ! total sky downward longwave flux (w/m2)     !
!    htr     (0:NLAY)     ! total sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnet    (0:NLAY)     ! total sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!

      do k = 0, NLAY
         urad(k) = f_zero
         drad(k) = f_zero
         totuflux(k) = f_zero
         totdflux(k) = f_zero
      enddo

      do k = 1, NLAY
        do iband = 1, NBANDS
          if (cldfrac(k) >= 1.e-6) then
            odcld(k,iband) = secdiff * taucloud(k,iband)
            transcld = exp(-odcld(k,iband))
            abscld(k,iband) = 1.0 - transcld
            icldlyr(k) = 1
          else
            odcld(k,iband) = f_zero
            abscld(k,iband) = f_zero
            icldlyr(k) = 0
          endif
        enddo
      enddo

!  ---  maximum/random cloud overlap parameter

      istcld(1) = 1
      istcldd(NLAY) = 1

      do k = 1, NLAY

        if (icldlyr(k) == 1) then
          istcld(k+1) = 0

          if (k == NLAY) then
            faccld1(k+1) = f_zero
            faccld2(k+1) = f_zero
            facclr1(k+1) = f_zero
            facclr2(k+1) = f_zero
            faccmb1(k+1) = f_zero
            faccmb2(k+1) = f_zero
          elseif (cldfrac(k+1) >= cldfrac(k)) then
            faccld1(k+1) = f_zero
            faccld2(k+1) = f_zero

            if (istcld(k) == 1) then
              facclr1(k+1) = f_zero
              facclr2(k+1) = f_zero

              if (cldfrac(k) < 1.0) then
                facclr2(k+1) = (cldfrac(k+1) - cldfrac(k))              &
     &                       / (1.0 - cldfrac(k))
              endif

              facclr2(k) = f_zero
              faccld2(k) = f_zero
            else
              fmax = max(cldfrac(k), cldfrac(k-1))

              if (cldfrac(k+1) > fmax) then
                facclr1(k+1) = rat2
                facclr2(k+1) = (cldfrac(k+1) - fmax) / (1.0 - fmax)
              elseif (cldfrac(k+1) < fmax) then
                facclr1(k+1) = (cldfrac(k+1) - cldfrac(k))              &
     &                       / (cldfrac(k-1) - cldfrac(k))
                facclr2(k+1) = f_zero
              else
                facclr1(k+1) = rat2
                facclr2(k+1) = f_zero
              endif
            endif

            if (facclr1(k+1)>f_zero .or. facclr2(k+1)>f_zero) then
              rat1 = 1.0
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            facclr1(k+1) = f_zero
            facclr2(k+1) = f_zero
            if (istcld(k) == 1) then
              faccld1(k+1) = f_zero
              faccld2(k+1) = (cldfrac(k) - cldfrac(k+1)) / cldfrac(k)
              facclr2(k) = f_zero
              faccld2(k) = f_zero
            else
              fmin = min(cldfrac(k), cldfrac(k-1))

              if (cldfrac(k+1) <= fmin) then
                faccld1(k+1) = rat1
                faccld2(k+1) = (fmin - cldfrac(k+1)) / fmin
              else
                faccld1(k+1) = (cldfrac(k) - cldfrac(k+1))              &
     &                       / (cldfrac(k) - fmin)
                faccld2(k+1) = f_zero
              endif
            endif

            if (faccld1(k+1)>f_zero .or. faccld2(k+1)>f_zero) then
              rat1 = f_zero
              rat2 = 1.0
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1(k+1) = facclr1(k+1)*faccld2(k)*cldfrac(k-1)
          faccmb2(k+1) = faccld1(k+1)*facclr2(k)*(1.0-cldfrac(k-1))
        else
          istcld(k+1) = 1
        endif
      enddo

      do k = NLAY, 1, -1
        if (icldlyr(k) == 1) then
          istcldd(k-1) = 0

          if (k == 1) then
            faccld1d(k-1) = f_zero
            faccld2d(k-1) = f_zero
            facclr1d(k-1) = f_zero
            facclr2d(k-1) = f_zero
            faccmb1d(k-1) = f_zero
            faccmb2d(k-1) = f_zero
          elseif (cldfrac(k-1) >= cldfrac(k)) then
            faccld1d(k-1) = f_zero
            faccld2d(k-1) = f_zero

            if (istcldd(k) == 1) then
              facclr1d(k-1) = f_zero
              facclr2d(k-1) = f_zero
              if (cldfrac(k) < 1.0) then
                facclr2d(k-1) = (cldfrac(k-1) - cldfrac(k))             &
     &                        / (1.0 - cldfrac(k))
              endif
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmax = max(cldfrac(k), cldfrac(k+1))

              if (cldfrac(k-1) > fmax) then
                facclr1d(k-1) = rat2
                facclr2d(k-1) = (cldfrac(k-1) - fmax) / (1.0 - fmax)
              elseif (cldfrac(k-1) < fmax) then
                facclr1d(k-1) = (cldfrac(k-1) - cldfrac(k))             &
     &                        / (cldfrac(k+1) - cldfrac(k))
                facclr2d(k-1) = f_zero
              else
                facclr1d(k-1) = rat2
                facclr2d(k-1) = f_zero
              endif
            endif

            if (facclr1d(k-1)>f_zero .or. facclr2d(k-1)>f_zero) then
              rat1 = 1.0
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            facclr1d(k-1) = f_zero
            facclr2d(k-1) = f_zero

            if (istcldd(k) == 1) then
              faccld1d(k-1) = f_zero
              faccld2d(k-1) = (cldfrac(k) - cldfrac(k-1)) / cldfrac(k)
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmin = min(cldfrac(k), cldfrac(k+1))

              if (cldfrac(k-1) <= fmin) then
                faccld1d(k-1) = rat1
                faccld2d(k-1) = (fmin - cldfrac(k-1)) / fmin
              else
                faccld1d(k-1) = (cldfrac(k) - cldfrac(k-1))             &
     &                        / (cldfrac(k) - fmin)
                faccld2d(k-1) = f_zero
              endif
            endif

            if (faccld1d(k-1)>f_zero .or. faccld2d(k-1)>f_zero) then
              rat1 = f_zero
              rat2 = 1.0
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1d(k-1) = facclr1d(k-1)*faccld2d(k)*cldfrac(k+1)
          faccmb2d(k-1) = faccld1d(k-1)*facclr2d(k)*(1.0-cldfrac(k+1))
        else
          istcldd(k-1) = 1
        endif
      enddo

!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)
          radld = f_zero
          iclddn = 0

!  ---  downward radiative transfer loop.
!           here are some variable definitions:
!             odtot      optical depth of gas and cloud
!             atrans     absorptivity for gas only
!             atot       absorptivity for gas and cloud
!             tfacgas    gas-only pade factor, used for planck fn
!             tfactot    gas and cloud pade factor, used for planck fn
!             bbdgas     gas-only planck function for downward rt
!             bbdtot     gas and cloud planck function for downward rt
!             bbutot     gas and cloud planck function for upward calc.
!             gassrc     source radiance due to gas only

          do k = NLAY, 1, -1
            plfrac = fracs(k,ig)
            blay = planklay(k,iband)
            dplankup = planklev(k,iband) - blay
            dplankdn = planklev(k-1,iband) - blay
            odepth = secdiff * taug(k,ig)
            if (odepth < f_zero) odepth = f_zero

            if (icldlyr(k) == 1) then
              odtot = odepth + odcld(k,iband)

              if (odtot < 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth_rec = rec_6 * odepth
                gassrc = plfrac*(blay + dplankdn*odepth_rec)*atrans(k)

                atot(k) =  odtot - 0.5*odtot*odtot
                odtot_rec = rec_6 * odtot
                bbdtot =  plfrac * (blay + dplankdn*odtot_rec)

                bbugas(k) =  plfrac * (blay + dplankup*odepth_rec)
                bbutot(k) =  plfrac * (blay + dplankup*odtot_rec)
              elseif (odepth <= 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth_rec = rec_6 * odepth
                gassrc = plfrac*(blay + dplankdn*odepth_rec)*atrans(k)

                odtot = odepth + odcld(k,iband)
                tblind = odtot / (bpade+odtot)
                ittot = tblint*tblind + 0.5
                tfactot = tf(ittot)
                bbdtot = plfrac * (blay + tfactot*dplankdn)
                atot(k) = 1.0 - trans(ittot)

                bbugas(k) = plfrac * (blay + dplankup*odepth_rec)
                bbutot(k) = plfrac * (blay + tfactot*dplankup)
              else
                tblind = odepth / (bpade + odepth)
                itgas = tblint*tblind + 0.5
                odepth = tautb(itgas)
                atrans(k) = 1.0 - trans(itgas)
                tfacgas = tf(itgas)
                gassrc = atrans(k)*plfrac*(blay + tfacgas*dplankdn)

                odtot = odepth + odcld(k,iband)
                tblind = odtot / (bpade + odtot)
                ittot = tblint*tblind + 0.5
                tfactot = tf(ittot)
                bbdtot = plfrac * (blay + tfactot*dplankdn)
                atot(k) = 1.0 - trans(ittot)

                bbugas(k) = plfrac * (blay + tfacgas*dplankup)
                bbutot(k) = plfrac * (blay + tfactot*dplankup)
              endif

              if (istcldd(k) == 1) then
                cldradd = cldfrac(k) * radld
                clrradd = radld - cldradd
                oldcld = cldradd
                oldclr = clrradd
                rad = f_zero
              endif

              ttot = 1.0 - atot(k)
              cldsrc = bbdtot * atot(k)
              cldradd = cldradd*ttot + cldfrac(k)*cldsrc
              clrradd = clrradd*(1.0 - atrans(k))                       &
     &                + gassrc*(1.0 - cldfrac(k))
              radld = cldradd + clrradd
              drad(k-1) = drad(k-1) + radld

              radmod = rad * ( facclr1d(k-1)*(1.0 - atrans(k))          &
     &               + faccld1d(k-1)* ttot ) - faccmb1d(k-1)*gassrc     &
     &               + faccmb2d(k-1)*cldsrc

              oldcld = cldradd - radmod
              oldclr = clrradd + radmod
              rad = -radmod + facclr2d(k-1)*oldclr-faccld2d(k-1)*oldcld
              cldradd = cldradd + rad
              clrradd = clrradd - rad
            else
              if (odepth <= 0.06) then
                atrans(k) = odepth - 0.5*odepth*odepth
                odepth = rec_6 * odepth
                bbd       =  plfrac * (blay + dplankdn*odepth)
                bbugas(k) =  plfrac * (blay + dplankup*odepth)
              else
                tblind = odepth / (bpade + odepth)
                itr = tblint*tblind + 0.5
                transc = trans(itr)
                atrans(k) = 1.0 - transc
                tausfac = tf(itr)
                bbd       = plfrac * (blay + tausfac*dplankdn)
                bbugas(k) = plfrac * (blay + tausfac*dplankup)
              endif

              radld = radld + (bbd - radld)*atrans(k)
              drad(k-1) = drad(k-1) + radld
            endif
          enddo

!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)

!  ---  add in reflection of surface downward radiance.
          reflect = 1.0 - semiss(iband)
          if (ireflect == 1) then         !  specular reflection.
            radlu = rad0 + reflect*radld
          else                            !  lambertian reflection.
            rad = 2.0 * radld * wtdiff
            radlu = rad0 + reflect*rad
          endif

          urad(0) = urad(0) + radlu
          do k = 1, NLAY
            if (icldlyr(k) == 1) then
              gassrc = bbugas(k) * atrans(k)

              if (istcld(k) == 1) then
                cldradu = cldfrac(k) * radlu
                clrradu = radlu - cldradu
                oldcld = cldradu
                oldclr = clrradu
                rad = f_zero
              endif

              ttot = 1.0 - atot(k)
              cldsrc = bbutot(k) * atot(k)
              cldradu = cldradu*ttot + cldfrac(k)*cldsrc
              clrradu = clrradu*(1.0 - atrans(k))                       &
     &                + gassrc*(1.0 - cldfrac(k))

!  ---  total sky radiance

              radlu = cldradu + clrradu
              urad(k) = urad(k) + radlu
              radmod = rad*( facclr1(k+1)*(1.0 - atrans(k))             &
     &               + faccld1(k+1)*ttot ) - faccmb1(k+1)*gassrc        &
     &               + faccmb2(k+1)*cldsrc
              oldcld = cldradu - radmod
              oldclr = clrradu + radmod
              rad = -radmod + facclr2(k+1)*oldclr - faccld2(k+1)*oldcld
              cldradu = cldradu + rad
              clrradu = clrradu - rad
            else
              radlu = radlu + (bbugas(k) - radlu)*atrans(k)
              urad(k) = urad(k) + radlu
            endif
          enddo

        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

        do k = NLAY, 0, -1
          uflux(k) = urad(k) * wtdiff
          dflux(k) = drad(k) * wtdiff
          urad(k)  = f_zero
          drad(k)  = f_zero
          totuflux(k) = totuflux(k) + uflux(k) * delwave(iband)
          totdflux(k) = totdflux(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb ) then
          do k = 0, NLAY
            fnet(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnet(k-1) - fnet(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!  ---  convert radiances to fluxes and heating rates for total sky.
!       calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuflux(k) = totuflux(k) * fluxfac
        totdflux(k) = totdflux(k) * fluxfac

        fnet(k) = totuflux(k) - totdflux(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htr(j) = heatfac * (fnet(j) - fnet(k)) / delp(k)
      enddo

      htr(NLAY) = f_zero

      return
!...................................
      end subroutine rtrcldmr
!-----------------------------------



!-----------------------------------
      subroutine rtreg                                                  &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       tauaer, NLAY, icldatm,                                     &
!  ---  outputs:
     &       totuclfl,totdclfl,htrcl, htrb                              &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary atmosphere.  the input to this       !
!   program is the atmospheric profile and all planck function          !
!   information.  only one angle is used from standard gaussian         !
!   quadrature.                                                         !
!                                                                       !
!     revision:  3.2       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, laytrop, ireflect, icldatm
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(:), intent(in) :: pavel, delp,   &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       semiss, plankbnd

      real (kind=kind_phys), dimension(:,:), intent(in) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htrcl,       &
     &       totuclfl, totdclfl
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX):: taug, fracs
      real (kind=kind_phys), dimension(0:NLAY)   :: dflux, uflux, fnetc

      real (kind=kind_phys), dimension(0:NLAY,MAXANG) :: urad, drad
      real (kind=kind_phys), dimension(NLAY,MAXANG)   :: bbu, atrans
      real (kind=kind_phys), dimension(MAXANG) :: rad, secang, angweigh

      real (kind=kind_phys) :: rad0, radlu, radld, blay, bbd, radsum,   &
     &       dplankup, dplankdn, odepth, tblind, tausfac, plfrac,       &
     &       reflect

      integer :: j, k, ig, iband, iang, numang, itr


!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuclfl(0:NLAY)     ! clear sky upward longwave flux (w/m2)       !
!    totdclfl(0:NLAY)     ! clear sky downward longwave flux (w/m2)     !
!    htrcl   (0:NLAY)     ! clear sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnetc   (0:NLAY)     ! clear sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!
      radsum = f_zero
      numang = numangs

!  ---  load angle data in arrays depending on angular quadrature scheme.

      do iang = 1, numang
        secang(iang) = secreg(iang,numang)
        angweigh(iang) = wtreg(iang,numang)
      enddo

      do k = 0, NLAY
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

      do iang = 1, numang
        do k = 0, NLAY
          urad(k,iang) = f_zero
          drad(k,iang) = f_zero
        enddo
      enddo

!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)

!  ---  loop over each angle for which the radiance is to be computed.

          lab_iang : do iang = 1, numang
            radld = f_zero

!  ---  downward radiative transfer.

            do k = NLAY, 1, -1
              blay = planklay(k,iband)
              plfrac = fracs(k,ig)
              dplankup = planklev(k,iband) - blay
              dplankdn = planklev(k-1,iband) - blay
              odepth = secang(iang) * taug(k,ig)

              if (odepth <= 0.06) then
                if (odepth < f_zero) odepth = f_zero
                atrans(k,iang) = odepth - 0.5*odepth*odepth
                odepth = rec_6 * odepth
                bbd =  plfrac * (blay + dplankdn*odepth)
                radld = radld + (bbd - radld)*atrans(k,iang)
                drad(k-1,iang) = drad(k-1,iang) + radld
                bbu(k,iang) =  plfrac * (blay + dplankup*odepth)
              else
                tblind = odepth / (bpade + odepth)
                itr = tblint*tblind + 0.5
                atrans(k,iang) = 1.0 - trans(itr)
                tausfac = tf(itr)
                bbd = plfrac * (blay + tausfac*dplankdn)
                radld = radld + (bbd - radld)*atrans(k,iang)
                drad(k-1,iang) = drad(k-1,iang) + radld
                bbu(k,iang) = plfrac * (blay + tausfac*dplankup)
              endif
            enddo

            rad(iang) = radld
            radsum = radsum + angweigh(iang) * radld
          enddo  lab_iang

!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)
          reflect = 1.0 - semiss(iband)

!  ---  add in reflection of surface downward radiance.

          do iang = 1, numang
            if (ireflect == 1) then
              radlu = rad0 + reflect*rad(iang)      ! --- specular reflection
            else
              radlu = rad0 + 2.0*reflect*radsum     ! --- lambertian reflection
            endif

            urad(0,iang) = urad(0,iang) + radlu
            do k = 1, NLAY
              radlu = radlu + (bbu(k,iang) - radlu)*atrans(k,iang)
              urad(k,iang) = urad(k,iang) + radlu
            enddo
          enddo

          radsum = f_zero
        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

        do k = NLAY, 0, -1
          uflux(k) = f_zero
          dflux(k) = f_zero

          do iang = 1, numang
            uflux(k) = uflux(k) + urad(k,iang)*angweigh(iang)
            dflux(k) = dflux(k) + drad(k,iang)*angweigh(iang)
            urad(k,iang) = f_zero
            drad(k,iang) = f_zero
          enddo

          totuclfl(k) = totuclfl(k) + uflux(k) * delwave(iband)
          totdclfl(k) = totdclfl(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb .and. icldatm==0 ) then
          do k = 0, NLAY
            fnetc(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnetc(k-1) - fnetc(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!===> ...  convert radiances to fluxes and heating rates for total sky.
!          calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuclfl(k) = totuclfl(k) * fluxfac
        totdclfl(k) = totdclfl(k) * fluxfac

        fnetc(k) = totuclfl(k) - totdclfl(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htrcl(j) = heatfac * (fnetc(j) - fnetc(k)) / delp(k)
      enddo

      htrcl(NLAY) = f_zero

      return
!...................................
      end subroutine rtreg
!-----------------------------------



!-----------------------------------
      subroutine rtregcld                                               &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr, htrb                                &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary cloudy atmosphere.  the input to     !
!   this program is the atmospheric profile, including cloud properties,!
!   and all planck function information.  only one angle is used from   !
!   standard gaussian quadrature.                                       !
!                                                                       !
!     revision:  3.2       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer,               intent(in) :: NLAY, laytrop, ireflect
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(0:), intent(in) :: cldfrac

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel, delp,  &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       semiss, plankbnd

      real (kind=kind_phys), dimension(:,:), intent(in) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, taucloud, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htr,         &
     &       totuflux, totdflux
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX) :: taug, fracs
      real (kind=kind_phys), dimension(0:NLAY)    :: dflux, uflux, fnet

      real (kind=kind_phys), dimension(0:NLAY,MAXANG) :: urad, drad
      real (kind=kind_phys), dimension(NLAY,MAXANG)   :: bbugas, atrans,&
     &       bbutot, atot
      real (kind=kind_phys), dimension(MAXANG) :: rad, secang, angweigh

      real (kind=kind_phys), dimension(NLAY,NBANDS,MAXANG) :: odcld,    &
     &       abscld, efclfrac

      real (kind=kind_phys) ::      rad0, radlu, radld, blay, bbd,      &
     &       bbdtot, dplankup, dplankdn, odepth, odtot, odepth_rec,     &
     &       odtot_rec, gassrc, tblind, transc, transcld, tausfac,      &
     &       tfactot, tfacgas, plfrac, reflect, radsum

      integer :: j, k, ig, iband, icldlyr(NLAY), ittot, itgas, itr,     &
     &       iang, numang

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    cldfrac (0:NLP1)     ! layer cloud fraction (padded at 2 ends)     !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuflux(0:NLAY)     ! total sky upward longwave flux (w/m2)       !
!    totdflux(0:NLAY)     ! total sky downward longwave flux (w/m2)     !
!    htr     (0:NLAY)     ! total sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnet    (0:NLAY)     ! total sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!
      radsum = f_zero
      numang = numangs

!  ---  load angle data in arrays depending on angular quadrature scheme.

      do iang = 1, numang
        secang(iang) = secreg(iang,numang)
        angweigh(iang) = wtreg(iang,numang)
      enddo

      totuflux(0) = f_zero
      totdflux(0) = f_zero
      do iang = 1, numang
        urad(0,iang) = f_zero
        drad(0,iang) = f_zero
      enddo

      do k = 1, NLAY
        totuflux(k) = f_zero
        totdflux(k) = f_zero

        do iang = 1, numang
          urad(k,iang) = f_zero
          drad(k,iang) = f_zero

          do iband = 1, NBANDS
            if (cldfrac(k) >= 1.e-6) then
              odcld(k,iband,iang) = secang(iang) * taucloud(k,iband)
              transcld = exp( -odcld(k,iband,iang) )
              abscld(k,iband,iang) = 1.0 - transcld
              efclfrac(k,iband,iang) = abscld(k,iband,iang)             &
     &                               * cldfrac(k)
              icldlyr(k) = 1
            else
              odcld(k,iband,iang) = f_zero
              abscld(k,iband,iang) = f_zero
              efclfrac(k,iband,iang) = f_zero
              icldlyr(k) = 0
            endif
          enddo
        enddo
      enddo

!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)

!  ---  downward radiative transfer loop.
!           here are some variable definitions:
!             odtot      optical depth of gas and cloud
!             atrans     absorptivity for gas only
!             atot       absorptivity for gas and cloud
!             tfacgas    gas-only pade factor, used for planck fn
!             tfactot    gas and cloud pade factor, used for planck fn
!             bbdgas     gas-only planck function for downward rt
!             bbdtot     gas and cloud planck function for downward rt
!             bbutot     gas and cloud planck function for upward calc.
!             gassrc     source radiance due to gas only

!  ---  loop over each angle for which the radiance is to be computed.

          do iang = 1, numang
            radld = f_zero

            do k = NLAY, 1, -1
              plfrac = fracs(k,ig)
              blay = planklay(k,iband)
              dplankup = planklev(k,iband) - blay
              dplankdn = planklev(k-1,iband) - blay
              odepth = secang(iang) * taug(k,ig)
              if (odepth < f_zero) odepth = f_zero

              if (icldlyr(k) == 1) then
                odtot = odepth + odcld(k,iband,iang)

                if (odtot < 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth_rec = rec_6 * odepth
                  gassrc = plfrac*(blay + dplankdn*odepth_rec)          &
     &                   * atrans(k,iang)

                  atot(k,iang) = odtot - 0.5*odtot*odtot
                  odtot_rec = rec_6 * odtot
                  bbdtot =  plfrac * (blay + dplankdn*odtot_rec)
                  radld = radld - radld * ( atrans(k,iang)              &
     &                  + efclfrac(k,iband,iang)*(1.0 - atrans(k,iang)))&
     &                  + gassrc + cldfrac(k)                           &
     &                  * (bbdtot*atot(k,iang) - gassrc)
                  drad(k-1,iang) = drad(k-1,iang) + radld

                  bbugas(k,iang) =  plfrac*(blay + dplankup*odepth_rec)
                  bbutot(k,iang) =  plfrac*(blay + dplankup*odtot_rec)
                elseif (odepth <= 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth_rec = rec_6 * odepth
                  gassrc = plfrac*(blay + dplankdn*odepth_rec)          &
     &                   * atrans(k,iang)

                  odtot = odepth + odcld(k,iband,iang)
                  tblind = odtot / (bpade + odtot)
                  ittot = tblint*tblind + 0.5
                  tfactot = tf(ittot)
                  bbdtot = plfrac * (blay + tfactot*dplankdn)
                  atot(k,iang) = 1.0 - trans(ittot)

                  radld = radld - radld * ( atrans(k,iang)              &
     &                  + efclfrac(k,iband,iang)*(1.0 - atrans(k,iang)))&
     &                  + gassrc + cldfrac(k)                           &
     &                  * (bbdtot*atot(k,iang) - gassrc)
                  drad(k-1,iang) = drad(k-1,iang) + radld

                  bbugas(k,iang) = plfrac*(blay + dplankup*odepth_rec)
                  bbutot(k,iang) = plfrac*(blay + tfactot*dplankup)
                else
                  tblind = odepth / (bpade + odepth)
                  itgas = tblint*tblind + 0.5
                  odepth = tautb(itgas)
                  atrans(k,iang) = 1.0 - trans(itgas)
                  tfacgas = tf(itgas)
                  gassrc = atrans(k,iang)*plfrac*(blay+tfacgas*dplankdn)

                  odtot = odepth + odcld(k,iband,iang)
                  tblind = odtot / (bpade + odtot)
                  ittot = tblint*tblind + 0.5
                  tfactot = tf(ittot)
                  bbdtot = plfrac * (blay + tfactot*dplankdn)
                  atot(k,iang) = 1.0 - trans(ittot)

                  radld = radld - radld * ( atrans(k,iang)              &
     &                  + efclfrac(k,iband,iang)*(1.0 - atrans(k,iang)))&
     &                  + gassrc + cldfrac(k)                           &
     &                  * (bbdtot*atot(k,iang) - gassrc)
                  drad(k-1,iang) = drad(k-1,iang) + radld

                  bbugas(k,iang) = plfrac * (blay + tfacgas*dplankup)
                  bbutot(k,iang) = plfrac * (blay + tfactot*dplankup)
                endif
              else
                if (odepth <= 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth = rec_6 * odepth
                  bbd            =  plfrac * (blay + dplankdn*odepth)
                  bbugas(k,iang) =  plfrac * (blay + dplankup*odepth)
                else
                  tblind = odepth / (bpade + odepth)
                  itr = tblint*tblind + 0.5
                  transc = trans(itr)
                  atrans(k,iang) = 1.0 - transc
                  tausfac = tf(itr)
                  bbd            = plfrac * (blay + tausfac*dplankdn)
                  bbugas(k,iang) = plfrac * (blay + tausfac*dplankup)
                endif

                radld = radld + (bbd - radld)*atrans(k,iang)
                drad(k-1,iang) = drad(k-1,iang) + radld
              endif
            enddo

            rad(iang) = radld
            radsum = radsum + angweigh(iang)*radld
          enddo

!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)

!  ---  add in reflection of surface downward radiance.
          reflect = 1.0 - semiss(iband)

          do iang = 1, numang
            if (ireflect == 1) then         !  specular reflection.
              radlu = rad0 + reflect*rad(iang)
            else                            !  lambertian reflection.
              radlu = rad0 + 2.0*reflect*radsum
            endif

            urad(0,iang) = urad(0,iang) + radlu
            do k = 1, NLAY
              if (icldlyr(k) == 1) then
                gassrc = bbugas(k,iang) * atrans(k,iang)
                radlu = radlu - radlu*( atrans(k,iang)                  &
     &                + efclfrac(k,iband,iang)*(1.0 - atrans(k,iang)) ) &
     &                + gassrc + cldfrac(k)                             &
     &                * (bbutot(k,iang)*atot(k,iang) - gassrc)
                urad(k,iang) = urad(k,iang) + radlu
              else
                radlu = radlu + (bbugas(k,iang) - radlu)*atrans(k,iang)
                urad(k,iang) = urad(k,iang) + radlu
              endif
            enddo
          enddo

          radsum = f_zero
        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

        do k = NLAY, 0, -1
          uflux(k) = f_zero
          dflux(k) = f_zero

          do iang = 1, numang
            uflux(k) = uflux(k) + urad(k,iang)*angweigh(iang)
            dflux(k) = dflux(k) + drad(k,iang)*angweigh(iang)
            urad(k,iang)  = f_zero
            drad(k,iang)  = f_zero
          enddo

          totuflux(k) = totuflux(k) + uflux(k) * delwave(iband)
          totdflux(k) = totdflux(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb ) then
          do k = 0, NLAY
            fnet(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnet(k-1) - fnet(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!  ---  convert radiances to fluxes and heating rates for total sky.
!       calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuflux(k) = totuflux(k) * fluxfac
        totdflux(k) = totdflux(k) * fluxfac

        fnet(k) = totuflux(k) - totdflux(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htr(j) = heatfac * (fnet(j) - fnet(k)) / delp(k)
      enddo

      htr(NLAY) = f_zero

      return
!...................................
      end subroutine rtregcld
!-----------------------------------



!-----------------------------------
      subroutine rtregcldmr                                             &
!...................................

!  ---  inputs:
     &     ( pavel,delp,semiss,ireflect,laytrop,colamt,wx,colbrd,       &
     &       coldry,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,planklay,planklev,plankbnd,           &
     &       cldfrac,taucloud,                                          &
     &       tauaer, NLAY,                                              &
!  ---  outputs:
     &       totuflux,totdflux,htr, htrb                                &
     &     )

!  *******************************************************************  !
!                                                                       !
!   this program calculates the upward fluxes, downward fluxes, and     !
!   heating rates for an arbitrary cloudy atmosphere.  the input to     !
!   this program is the atmospheric profile, including cloud properties,!
!   and all planck function information.   first order standard gaussian!
!   quadrature is used for the angle integration.  clouds are treated   !
!   with maximum/random overlap scheme.                                 !
!                                                                       !
!     revision:  3.2       created:  2002/08/15                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer,               intent(in) :: NLAY, laytrop, ireflect
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(0:), intent(in) :: cldfrac

      real (kind=kind_phys), dimension(:),  intent(in) :: pavel, delp,  &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       semiss, plankbnd

      real (kind=kind_phys), dimension(:,:), intent(in) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, planklay, taucloud, tauaer

      real (kind=kind_phys), dimension(0:,:),intent(in) :: planklev

!  ---  outputs:
      real (kind=kind_phys), dimension(0:), intent(out) :: htr,         &
     &       totuflux, totdflux
      real (kind=kind_phys), dimension(:,:), intent(out) :: htrb

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY,NGMX) :: taug, fracs
      real (kind=kind_phys), dimension(0:NLAY)    :: dflux, uflux, fnet

      real (kind=kind_phys), dimension(0:NLAY,MAXANG) :: urad, drad
      real (kind=kind_phys), dimension(NLAY,MAXANG)   :: bbugas, atrans,&
     &       bbutot, atot
      real (kind=kind_phys), dimension(MAXANG) :: rad, secang, angweigh

      real (kind=kind_phys), dimension(NLAY,NBANDS,MAXANG)::odcld,abscld

      real (kind=kind_phys), dimension(NLAY+1)    :: faccld1, faccld2,  &
     &       facclr1, facclr2, faccmb1, faccmb2

      real (kind=kind_phys), dimension(0:NLAY)    :: faccld1d,          &
     &       faccld2d, facclr1d, facclr2d, faccmb1d, faccmb2d

      integer :: icldlyr(NLAY), istcld(NLAY+1), istcldd(0:NLAY)

      real (kind=kind_phys) :: radt, rad0, radlu, radld, blay, bbd,     &
     &       bbdtot, dplankup, dplankdn, odepth, odtot, odepth_rec,     &
     &       odtot_rec, gassrc, tblind, transc, transcld, tausfac,      &
     &       tfactot, tfacgas, plfrac, reflect, radsum, cldradd,        &
     &       clrradd, cldradu, clrradu, cldsrc, radmod, oldcld, oldclr, &
     &       fmax, rat1, rat2, fmin, ttot

      integer :: j, k, ig, iband, ittot, itgas, itr, numang, iang


!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    pavel   (NLAY)       ! layer pressure (mb)                         !
!    delp    (NLAY)       ! layer pressure thickness (mb)               !
!    cldfrac (0:NLP1)     ! layer cloud fraction (padded at 2 ends)     !
!    NLAY                 ! number of model layers/levels               !
!                                                                       !
!  constants or shared variables:                                       !
!    NBANDS               ! number of longwave spectral bands           !
!    wtdiff               ! weight for radiance to flux conversion      !
!    bpade                ! pade constant                               !
!    tf                   ! tautb transition function look-up table     !
!    trans                ! clear sky transmittance look-up table       !
!                                                                       !
!  output variables:                                                    !
!    totuflux(0:NLAY)     ! total sky upward longwave flux (w/m2)       !
!    totdflux(0:NLAY)     ! total sky downward longwave flux (w/m2)     !
!    htr     (0:NLAY)     ! total sky longwave heating rate (k/day)     !
!                                                                       !
!  local variables:                                                     !
!                                                                       !
!    fnet    (0:NLAY)     ! total sky net longwave flux (w/m2)          !
!                                                                       !
!  =====================    end of definitions    ====================  !

!
!===> ... begin here
!
      radsum = f_zero
      numang = numangs

!  ---  load angle data in arrays depending on angular quadrature scheme.

      do iang = 1, numang
        secang(iang) = secreg(iang,numang)
        angweigh(iang) = wtreg(iang,numang)
      enddo

      totuflux(0) = f_zero
      totdflux(0) = f_zero

      do iang = 1, numang
        urad(0,iang) = f_zero
        drad(0,iang) = f_zero
      enddo

      do k = 1, NLAY
        totuflux(k) = f_zero
        totdflux(k) = f_zero

        do iang = 1, numang
          urad(k,iang) = f_zero
          drad(k,iang) = f_zero

          do iband = 1, NBANDS
            if (cldfrac(k) >= 1.e-6) then
              odcld(k,iband,iang) = secang(iang)*taucloud(k,iband)
              transcld = exp(-odcld(k,iband,iang))
              abscld(k,iband,iang) = 1.0 - transcld
              icldlyr(k) = 1
            else
              odcld(k,iband,iang) = f_zero
              abscld(k,iband,iang) = f_zero
              icldlyr(k) = 0
            endif
          enddo
        enddo
      enddo

!  ---  maximum/random cloud overlap parameter

      istcld(1) = 1
      istcldd(NLAY) = 1

      do k = 1, NLAY

        if (icldlyr(k) == 1) then
          istcld(k+1) = 0

          if (k == NLAY) then
            faccld1(k+1) = f_zero
            faccld2(k+1) = f_zero
            facclr1(k+1) = f_zero
            facclr2(k+1) = f_zero
            faccmb1(k+1) = f_zero
            faccmb2(k+1) = f_zero
          elseif (cldfrac(k+1) >= cldfrac(k)) then
            faccld1(k+1) = f_zero
            faccld2(k+1) = f_zero

            if (istcld(k) == 1) then
              facclr1(k+1) = f_zero
              facclr2(k+1) = f_zero

              if (cldfrac(k) < 1.0) then
                facclr2(k+1) = (cldfrac(k+1) - cldfrac(k))              &
     &                       / (1.0 - cldfrac(k))
              endif

              facclr2(k) = f_zero
              faccld2(k) = f_zero
            else
              fmax = max(cldfrac(k), cldfrac(k-1))

              if (cldfrac(k+1) > fmax) then
                facclr1(k+1) = rat2
                facclr2(k+1) = (cldfrac(k+1) - fmax) / (1.0 - fmax)
              elseif (cldfrac(k+1) < fmax) then
                facclr1(k+1) = (cldfrac(k+1) - cldfrac(k))              &
     &                       / (cldfrac(k-1) - cldfrac(k))
                facclr2(k+1) = f_zero
              else
                facclr1(k+1) = rat2
                facclr2(k+1) = f_zero
              endif
            endif

            if (facclr1(k+1)>f_zero .or. facclr2(k+1)>f_zero) then
              rat1 = 1.0
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            facclr1(k+1) = f_zero
            facclr2(k+1) = f_zero
            if (istcld(k) == 1) then
              faccld1(k+1) = f_zero
              faccld2(k+1) = (cldfrac(k) - cldfrac(k+1)) / cldfrac(k)
              facclr2(k) = f_zero
              faccld2(k) = f_zero
            else
              fmin = min(cldfrac(k), cldfrac(k-1))

              if (cldfrac(k+1) <= fmin) then
                faccld1(k+1) = rat1
                faccld2(k+1) = (fmin - cldfrac(k+1)) / fmin
              else
                faccld1(k+1) = (cldfrac(k) - cldfrac(k+1))              &
     &                       / (cldfrac(k) - fmin)
                faccld2(k+1) = f_zero
              endif
            endif

            if (faccld1(k+1)>f_zero .or. faccld2(k+1)>f_zero) then
              rat1 = f_zero
              rat2 = 1.0
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1(k+1) = facclr1(k+1)*faccld2(k)*cldfrac(k-1)
          faccmb2(k+1) = faccld1(k+1)*facclr2(k)*(1.0-cldfrac(k-1))
        else
          istcld(k+1) = 1
        endif
      enddo

      do k = NLAY, 1, -1
        if (icldlyr(k) == 1) then
          istcldd(k-1) = 0

          if (k == 1) then
            faccld1d(k-1) = f_zero
            faccld2d(k-1) = f_zero
            facclr1d(k-1) = f_zero
            facclr2d(k-1) = f_zero
            faccmb1d(k-1) = f_zero
            faccmb2d(k-1) = f_zero
          elseif (cldfrac(k-1) >= cldfrac(k)) then
            faccld1d(k-1) = f_zero
            faccld2d(k-1) = f_zero

            if (istcldd(k) == 1) then
              facclr1d(k-1) = f_zero
              facclr2d(k-1) = f_zero
              if (cldfrac(k) < 1.0) then
                facclr2d(k-1) = (cldfrac(k-1) - cldfrac(k))             &
     &                        / (1.0 - cldfrac(k))
              endif
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmax = max(cldfrac(k), cldfrac(k+1))

              if (cldfrac(k-1) > fmax) then
                facclr1d(k-1) = rat2
                facclr2d(k-1) = (cldfrac(k-1) - fmax) / (1.0 - fmax)
              elseif (cldfrac(k-1) < fmax) then
                facclr1d(k-1) = (cldfrac(k-1) - cldfrac(k))             &
     &                        / (cldfrac(k+1) - cldfrac(k))
                facclr2d(k-1) = f_zero
              else
                facclr1d(k-1) = rat2
                facclr2d(k-1) = f_zero
              endif
            endif

            if (facclr1d(k-1)>f_zero .or. facclr2d(k-1)>f_zero) then
              rat1 = 1.0
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            facclr1d(k-1) = f_zero
            facclr2d(k-1) = f_zero

            if (istcldd(k) == 1) then
              faccld1d(k-1) = f_zero
              faccld2d(k-1) = (cldfrac(k) - cldfrac(k-1)) / cldfrac(k)
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmin = min(cldfrac(k), cldfrac(k+1))

              if (cldfrac(k-1) <= fmin) then
                faccld1d(k-1) = rat1
                faccld2d(k-1) = (fmin - cldfrac(k-1)) / fmin
              else
                faccld1d(k-1) = (cldfrac(k) - cldfrac(k-1))             &
     &                        / (cldfrac(k) - fmin)
                faccld2d(k-1) = f_zero
              endif
            endif

            if (faccld1d(k-1)>f_zero .or. faccld2d(k-1)>f_zero) then
              rat1 = f_zero
              rat2 = 1.0
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1d(k-1) = facclr1d(k-1)*faccld2d(k)*cldfrac(k+1)
          faccmb2d(k-1) = faccld1d(k-1)*facclr2d(k)*(1.0-cldfrac(k+1))
        else
          istcldd(k-1) = 1
        endif
      enddo

!  ---  loop over frequency bands.

      lab_iband : do iband = 1, NBANDS

        call taumol                                                     &
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       iband, NLAY,                                               &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ---  radiative transfer starts here.  loop over g-channels.

        lab_ig : do ig = 1, ngb(iband)

          do iang = 1, numang
            radld = f_zero
            cldradd = f_zero
            clrradd = f_zero
            oldcld = f_zero
            oldclr = f_zero
            radt = f_zero

!  ---  downward radiative transfer loop.
!           here are some variable definitions:
!             odtot      optical depth of gas and cloud
!             atrans     absorptivity for gas only
!             atot       absorptivity for gas and cloud
!             tfacgas    gas-only pade factor, used for planck fn
!             tfactot    gas and cloud pade factor, used for planck fn
!             bbdgas     gas-only planck function for downward rt
!             bbdtot     gas and cloud planck function for downward rt
!             bbutot     gas and cloud planck function for upward calc.
!             gassrc     source radiance due to gas only

            do k = NLAY, 1, -1
              plfrac = fracs(k,ig)
              blay = planklay(k,iband)
              dplankup = planklev(k,  iband) - blay
              dplankdn = planklev(k-1,iband) - blay
              odepth = secang(iang) * taug(k,ig)
              if (odepth < f_zero) odepth = f_zero

              if (icldlyr(k) == 1) then
                odtot = odepth + odcld(k,iband,iang)

                if (odtot < 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth_rec = rec_6 * odepth
                  gassrc = plfrac*(blay + dplankdn*odepth_rec)          &
     &                   * atrans(k,iang)

                  atot(k,iang) =  odtot - 0.5*odtot*odtot
                  odtot_rec = rec_6 * odtot
                  bbdtot =  plfrac*(blay + dplankdn*odtot_rec)

                  bbugas(k,iang) =  plfrac*(blay + dplankup*odepth_rec)
                  bbutot(k,iang) =  plfrac*(blay + dplankup*odtot_rec)
                elseif (odepth <= 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth_rec = rec_6 * odepth
                  gassrc = plfrac*(blay + dplankdn*odepth_rec)          &
     &                   * atrans(k,iang)

                  odtot = odepth + odcld(k,iband,iang)
                  tblind = odtot / (bpade+odtot)
                  ittot = tblint*tblind + 0.5
                  tfactot = tf(ittot)
                  bbdtot = plfrac*(blay + tfactot*dplankdn)
                  atot(k,iang) = 1.0 - trans(ittot)

                  bbugas(k,iang) = plfrac*(blay + dplankup*odepth_rec)
                  bbutot(k,iang) = plfrac*(blay + tfactot*dplankup)
                else
                  tblind = odepth / (bpade + odepth)
                  itgas = tblint*tblind + 0.5
                  odepth = tautb(itgas)
                  atrans(k,iang) = 1.0 - trans(itgas)
                  tfacgas = tf(itgas)
                  gassrc = atrans(k,iang)*plfrac                        &
     &                   * (blay + tfacgas*dplankdn)

                  odtot = odepth + odcld(k,iband,iang)
                  tblind = odtot / (bpade + odtot)
                  ittot = tblint*tblind + 0.5
                  tfactot = tf(ittot)
                  bbdtot = plfrac*(blay + tfactot*dplankdn)
                  atot(k,iang) = 1.0 - trans(ittot)

                  bbugas(k,iang) = plfrac*(blay + tfacgas*dplankup)
                  bbutot(k,iang) = plfrac*(blay + tfactot*dplankup)
                endif

                if (istcldd(k) == 1) then
                  cldradd = cldfrac(k) * radld
                  clrradd = radld - cldradd
                  oldcld = cldradd
                  oldclr = clrradd
                  radt = f_zero
                endif

                ttot = 1.0 - atot(k,iang)
                cldsrc = bbdtot * atot(k,iang)
                cldradd = cldradd*ttot + cldfrac(k)*cldsrc
                clrradd = clrradd*(1.0 - atrans(k,iang))                &
     &                  + gassrc*(1.0 - cldfrac(k))
                radld = cldradd + clrradd
                drad(k-1,iang) = drad(k-1,iang) + radld

                radmod = radt * ( facclr1d(k-1)*(1.0 - atrans(k,iang))  &
     &                 + faccld1d(k-1)* ttot ) - faccmb1d(k-1)*gassrc   &
     &                 + faccmb2d(k-1)*cldsrc

                oldcld = cldradd - radmod
                oldclr = clrradd + radmod
                radt = -radmod + facclr2d(k-1)*oldclr                   &
     &               - faccld2d(k-1)*oldcld
                cldradd = cldradd + radt
                clrradd = clrradd - radt
              else
                if (odepth <= 0.06) then
                  atrans(k,iang) = odepth - 0.5*odepth*odepth
                  odepth = rec_6 * odepth
                  bbd            =  plfrac*(blay + dplankdn*odepth)
                  bbugas(k,iang) =  plfrac*(blay + dplankup*odepth)
                else
                  tblind = odepth / (bpade + odepth)
                  itr = tblint*tblind + 0.5
                  transc = trans(itr)
                  atrans(k,iang) = 1.0 - transc
                  tausfac = tf(itr)
                  bbd            = plfrac*(blay + tausfac*dplankdn)
                  bbugas(k,iang) = plfrac*(blay + tausfac*dplankup)
                endif

                radld = radld + (bbd - radld)*atrans(k,iang)
                drad(k-1,iang) = drad(k-1,iang) + radld
              endif
            enddo

            rad(iang) = radld
            radsum = radsum + angweigh(iang) * radld
          enddo

!  ---  upward radiative transfer

          rad0 = fracs(1,ig) * plankbnd(iband)

!  ---  add in reflection of surface downward radiance.
          reflect = 1.0 - semiss(iband)

          do iang = 1, numang
            radld = f_zero
            cldradd = f_zero
            clrradd = f_zero
            oldcld = f_zero
            oldclr = f_zero
            radt = f_zero

            if (ireflect == 1) then         !  specular reflection.
              radlu = rad0 + reflect*rad(iang)
            else                            !  lambertian reflection.
              radlu = rad0 + 2.0*reflect*radsum
            endif

            urad(0,iang) = urad(0,iang) + radlu
            do k = 1, NLAY
              if (icldlyr(k) == 1) then
                gassrc = bbugas(k,iang) * atrans(k,iang)

                if (istcld(k) == 1) then
                  cldradu = cldfrac(k) * radlu
                  clrradu = radlu - cldradu
                  oldcld = cldradu
                  oldclr = clrradu
                  radt = f_zero
                endif

                ttot = 1.0 - atot(k,iang)
                cldsrc = bbutot(k,iang) * atot(k,iang)
                cldradu = cldradu*ttot + cldfrac(k)*cldsrc
                clrradu = clrradu*(1.0 - atrans(k,iang))                &
     &                  + gassrc*(1.0 - cldfrac(k))

!  ---  total sky radiance

                radlu = cldradu + clrradu
                urad(k,iang) = urad(k,iang) + radlu
                radmod = radt*( facclr1(k+1)*(1.0 - atrans(k,iang))     &
     &                 + faccld1(k+1)*ttot ) - faccmb1(k+1)*gassrc      &
     &                 + faccmb2(k+1)*cldsrc
                oldcld = cldradu - radmod
                oldclr = clrradu + radmod
                radt = -radmod + facclr2(k+1)*oldclr                    &
     &               - faccld2(k+1)*oldcld
                cldradu = cldradu + radt
                clrradu = clrradu - radt
              else
                radlu = radlu + (bbugas(k,iang) - radlu)*atrans(k,iang)
                urad(k,iang) = urad(k,iang) + radlu
              endif
            enddo
          enddo

          radsum = f_zero
        enddo  lab_ig

!  ---  calculate upward, downward, and net flux.

        do k = NLAY, 0, -1
          uflux(k) = f_zero
          dflux(k) = f_zero

          do iang = 1, numang
            uflux(k) = uflux(k) + urad(k,iang)*angweigh(iang)
            dflux(k) = dflux(k) + drad(k,iang)*angweigh(iang)
            urad(k,iang)  = f_zero
            drad(k,iang)  = f_zero
          enddo

          totuflux(k) = totuflux(k) + uflux(k) * delwave(iband)
          totdflux(k) = totdflux(k) + dflux(k) * delwave(iband)
        enddo

!! ---  optional spectral band heating

        if ( lhlwb ) then
          do k = 0, NLAY
            fnet(k) = fluxfac * delwave(iband) * (uflux(k) - dflux(k))
          enddo

          do k = 1, NLAY
            htrb(k,iband) = heatfac * (fnet(k-1) - fnet(k)) / delp(k)
          enddo
        endif

      enddo  lab_iband

!  ---  convert radiances to fluxes and heating rates for total sky.
!       calculates clear sky surface and toa values.

      do k = 0, NLAY
        totuflux(k) = totuflux(k) * fluxfac
        totdflux(k) = totdflux(k) * fluxfac

        fnet(k) = totuflux(k) - totdflux(k)
      enddo

!  --- ...  calculate heating rates.
      do k = 1, NLAY
         j = k - 1
         htr(j) = heatfac * (fnet(j) - fnet(k)) / delp(k)
      enddo

      htr(NLAY) = f_zero

      return
!...................................
      end subroutine rtregcldmr
!-----------------------------------



!-----------------------------------
      subroutine taumol                                                 &
!...................................
!  ---  inputs:
     &     ( pavel,laytrop,colamt,wx,colbrd,coldry,                     &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                         &
     &       indself,selffac,selffrac,indfor,forfac,forfrac,            &
     &       indminor,scaleminor,scaleminorn2,minorfrac,                &
     &       rat_h2oco2,rat_h2oo3,rat_h2on2o,rat_h2och4,                &
     &       rat_n2oco2,rat_o3co2,tauaer,                               &
     &       IB, NLAY,                                                  &
!  ---  outputs:
     &       taug, fracs                                                &
     &     )

!  ************    original subprogram description    ***************  *
!                                                                      *
!                  optical depths developed for the                    *
!                                                                      *
!                rapid radiative transfer model (rrtm)                 *
!                                                                      *
!            atmospheric and environmental research, inc.              *
!                        131 hartwell avenue                           *
!                        lexington, ma 02421                           *
!                                                                      *
!                           eli j. mlawer                              *
!                         jennifer delamere                            *
!                         steven j. taubman                            *
!                         shepard a. clough                            *
!                                                                      *
!                       email:  mlawer@aer.com                         *
!                       email:  jdelamer@aer.com                       *
!                                                                      *
!        the authors wish to acknowledge the contributions of the      *
!        following people:  karen cady-pereira, patrick d. brown,      *
!        michael j. iacono, ronald e. farren, luke chen,               *
!        robert bergstrom.                                             *
!                                                                      *
!                                                                      *
!     taumol                                                           *
!                                                                      *
!     this file contains the subroutines taugbn (where n goes from     *
!     1 to 16).  taugbn calculates the optical depths and planck       *
!     fractions per g-value and layer for band n.                      *
!                                                                      *
!  output:  optical depths (unitless)                                  *
!           fractions needed to compute planck functions at every layer*
!           and g-value                                                *
!                                                                      !
!  description:                                                        !
!     NG##        - number of g-values in band ## (##=01-16)           !
!     nspa(iband) - for the lower atmosphere, the number of reference  !
!                   atmospheres that are stored for band iband per     !
!                   pressure level and temperature.  each of these     !
!                   atmospheres has different relative amounts of the  !
!                   key species for the band (i.e. different binary    !
!                   species parameters).                               !
!     nspb(iband) - same for upper atmosphere                          !
!     oneminus    - since problems are caused in some cases by         !
!                   interpolation parameters equal to or greater than  !
!                   1, for these cases these parameters are set to this!
!                   value, slightly < 1.                               !
!     laytrop     - layer at which switch is made from one combination !
!                   of key species to another                          !
!     colamt(NLAY,MAXGAS)                                              !
!                 - column amounts of water vapor,carbon dioxide,      !
!                   ozone, nitrous oxide, methane, o2, co, respectively!
!                   (molecules/cm**2)                                  !
!     coldry(NLAY)- dry air column amount (1.e-20*molecules/cm**2)     !
!     facij (NLAY)- for layer lay, these are factors that are needed to!
!                   compute the interpolation factors that multiply the!
!                   appropriate reference k-values.  a value of 0 (1)  !
!                   for i,j indicates that the corresponding factor    !
!                   multiplies reference k-value for the lower (higher)!
!                   of the two appropriate temperatures, and altitudes,!
!                   respectively.                                      !
!     jp (NLAY)   - the index of the lower (in altitude) of the two    !
!                   appropriate reference pressure levels needed for   !
!                   interpolation.                                     !
!     jt, jt1(NLAY)-the indices of the lower of the two appropriate    !
!                   reference temperatures needed for interpolation    !
!                   (for pressure levels jp and jp+1, respectively)    !
!     selffac(NLAY)-scale factor needed to water vapor self-continuum, !
!                   equals (water vapor density)/(atmospheric density  !
!                   at 296k and 1013 mb)                               !
!     selffrac(NLAY)factor needed for temperature interpolation of     !
!                   reference water vapor self-continuum data          !
!     indself(NLAY)-index of the lower of the two appropriate reference!
!                   temperatures needed for the self-continuum         !
!                   interpolation                                      !
!     forfac       -scale factor needed for water vapor foreign-       !
!                   continuum.                                         !
!     forfrac      -factor needed for temperature interpolation of     !
!                   reference water vapor foreign-continuum data       !
!     indfor       -index of the lower of the two appropriate reference!
!                   temperatures needed for the foreign-continuum      !
!                   interpolation                                      !
!                                                                      !
!  data input                                                          !
!     absa(nspa(nn),5,13,NGnn), absb(nspb(nn),5,13:59,NGnn),           !
!     selfref(10,NGnn)                                                 !
!                      (note:  nn is the band number)                  !
!                                                                      !
!     absa    - k-values for low reference atmospheres (no water vapor !
!               self-continuum) (units: cm**2/molecule)                !
!     absb    - k-values for high reference atmospheres (all sources)  !
!               (units: cm**2/molecule)                                !
!     ama'gas'- k-values for low reference atmosphere minor species    !
!               (units: cm**2/molecule)                                !
!     amb'gas'- k-values for high reference atmosphere minor species   !
!               (units: cm**2/molecule)                                !
!     selfref - k-values for water vapor self-continuum for reference  !
!               atmospheres (used below laytrop)                       !
!               (units: cm**2/molecule)                                !
!     forref  - k-values for water vapor foreign-continuum for         !
!               reference atmospheres (used below/above laytrop)       !
!               (units: cm**2/molecule)                                !
!                                                                      !
!     dimension absa(65*nspa(n),NG##), absb(235*nspb(n),NG##)          !
!                                                                      !
!  ******************************************************************  !
!
      implicit none
!
!  ---  inputs:
      integer,               intent(in) :: laytrop, IB, NLAY
      integer, dimension(:), intent(in) :: jp, jt, jt1, indself,        &
     &       indfor, indminor

      real (kind=kind_phys), dimension(:), intent(in) :: pavel,         &
     &       colbrd, coldry, fac00, fac01, fac10, fac11, scaleminorn2,  &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor

      real (kind=kind_phys), dimension(:,:), intent(in) :: rat_h2oco2,  &
     &       rat_h2oo3, rat_h2on2o, rat_h2och4, rat_n2oco2, rat_o3co2,  &
     &       colamt, wx, tauaer

!  ---  outputs:
      real (kind=kind_phys),dimension(:,:), intent(out) :: taug, fracs

!  ---  locals:
      real (kind=kind_phys)     :: tem1, tem2
      integer :: k0, id0(NLAY), id1(NLAY)
!
!===> ... begin here
!
      do k0 = 1, laytrop
        id0(k0) = ((jp(k0)-1) *5 + jt (k0) - 1) * nspa(IB)
        id1(k0) = ( jp(k0)    *5 + jt1(k0) - 1) * nspa(IB)
      enddo

      do k0 = laytrop+1, NLAY
        id0(k0) = ((jp(k0)-13)*5 + jt (k0) - 1) * nspb(IB)
        id1(k0) = ((jp(k0)-12)*5 + jt1(k0) - 1) * nspb(IB)
      enddo

      if (IB == 1) then
        call taugb01
      else if (IB ==  2) then
        call taugb02
      else if (IB ==  3) then
        call taugb03
      else if (IB ==  4) then
        call taugb04
      else if (IB ==  5) then
        call taugb05
      else if (IB ==  6) then
        call taugb06
      else if (IB ==  7) then
        call taugb07
      else if (IB ==  8) then
        call taugb08
      else if (IB ==  9) then
        call taugb09
      else if (IB == 10) then
        call taugb10
      else if (IB == 11) then
        call taugb11
      else if (IB == 12) then
        call taugb12
      else if (IB == 13) then
        call taugb13
      else if (IB == 14) then
        call taugb14
      else if (IB == 15) then
        call taugb15
      else if (IB == 16) then
        call taugb16
      endif


! =================
      contains
! =================

!-----------------------------------
      subroutine taugb01
!...................................

!  ------------------------------------------------------------------  !
!     written by eli j. mlawer, atmospheric & environmental research.  !
!                                                                      !
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)             !
!                          (high key - h2o; high minor - n2)           !
!     note: previous versions of rrtm band 1:                          !
!              10-250 cm-1 (low - h2o; high - h2o)                     !
!                                                                      !
!  ------------------------------------------------------------------  !
!
      use module_radlw_kgb01
!
      implicit none

!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf, indm
!
      real (kind=kind_phys) :: corradj, scalen2o, tauself, taufor,      &
     &      taun2o

!   ---  compute the optical depth by interpolating in ln(pressure) and
!        temperature.  below laytrop, the water vapor self-continuum and
!        foreign continuum is interpolated (in temperature) separately.

      do k = 1, laytrop
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA01, id1(k) + 1 )
        ind12 = min(MSA01, ind11 + 1)
        inds  = indself(k)
        indf  = indfor(k)
        indm  = indminor(k)

        if (pavel(k) < 250.0) then
          corradj = 1.0 - 0.15*(250.0 - pavel(k)) / 154.4
        else
          corradj = 1.0
        endif
        scalen2o = colbrd(k) * scaleminorn2(k)

        do j = 1, NG01
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          taun2o  = scalen2o   * ( aman2o (indm,j) + minorfrac(k)       &
     &            * (aman2o (indm+1,j) - aman2o (indm,j)) )

          taug(k,j) = corradj * ( colamt(k,1)                           &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + taun2o ) + tauaer(k,1)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1
        indf  = indfor(k)
        indm  = indminor(k)

        corradj =  1.0 - 0.15*( pavel(k)/95.6 )
        scalen2o = colbrd(k) * scaleminorn2(k)

        do j = 1, NG01
          taufor = forfac(k) * ( forref(indf,j) + forfrac(k)            &
     &           * (forref(indf+1,j) - forref(indf,j)) )
          taun2o = scalen2o  * ( ambn2o(indm,j) + minorfrac(k)          &
     &           * (ambn2o(indm+1,j) - ambn2o(indm,j)) )

          taug(k,j) = corradj * ( colamt(k,1)                           &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + taufor + taun2o ) + tauaer(k,1)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb01
!-----------------------------------


!-----------------------------------
      subroutine taugb02
!...................................

!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)

!     note: previous version of rrtm band 2:
!              250-500 cm-1 (low - h2o; high - h2o)
!
      use module_radlw_kgb02
!
      implicit none

!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf
!
      real (kind=kind_phys) :: corradj, tauself, taufor

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  below laytrop, the water vapor self-continuum is 
!     interpolated (in temperature) separately.

      do k = 1, laytrop
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA02, id1(k) + 1 )
        ind12 = min(MSA02, ind11 + 1 )
        inds  = indself(k)
        indf  = indfor(k)

        corradj = 1.0 - 0.05*(pavel(k) - 100.0)/900.0

        do j = 1, NG02
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = corradj * ( colamt(k,1)                           &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor ) + tauaer(k,2)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1
        indf  = indfor(k)

        do j = 1, NG02
          taufor = forfac(k) * ( forref(indf,j) + forfrac(k)            &
     &           * (forref(indf+1,j) - forref(indf,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + taufor + tauaer(k,2)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb02
!-----------------------------------


!-----------------------------------
      subroutine taugb03
!...................................

!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!
      use module_radlw_kgb03
!
      implicit none

!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jmn2o, jpl, jpp1
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_mn2o, specmult_mn2o, speccomb_planck,  &
     &       specmult_planck, refrat_planck_a, refrat_planck_b,         &
     &       refrat_m_a, refrat_m_b, fs, fs1, fmn2o, fmn2omf, fpl,      &
     &       fk00, fk10, fk20, fk01, fk11, fk21, n2om1, n2om2, absn2o,  &
     &       ratn2o, adjcoln2o, tauself, taufor, tem1, tem2

!  ---  p = 212.725 mb
      refrat_planck_a = chi_mls(1,9)  / chi_mls(2,9)

!  ---  p = 95.58 mb
      refrat_planck_b = chi_mls(1,13) / chi_mls(2,13)

!  ---  p = 706.270mb
      refrat_m_a = chi_mls(1,3)  / chi_mls(2,3)

!  ---  p = 95.58 mb
      refrat_m_b = chi_mls(1,13) / chi_mls(2,13)

!     compute the optical depth by interpolating in ln(pressure) and
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated (in
!     temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,2)
        specmult_mn2o = 8.0 * min(oneminus, colamt(k,1)/speccomb_mn2o)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - int(specmult_mn2o)

        fmn2omf = minorfrac(k) * fmn2o

!     in atmospheres where the amount of n2o is too great to be considered
!     a minor species, adjust the column amount of n2o by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratn2o  = colamt(k,4) / (coldry(k) * chi_mls(4,jpp1))
        if (ratn2o > 1.5) then
          adjcoln2o = (0.5 + (ratn2o - 0.5)**0.65)                      &
     &              * chi_mls(4,jpp1) * coldry(k)
        else
          adjcoln2o = colamt(k,4)
        endif

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA03, id1(k) + js1 )
          ind12 = min(MSA03, ind11 + 1 )
          ind13 = min(MSA03, ind11 + 2 )
          ind14 = min(MSA03, ind11 + 9 )
          ind15 = min(MSA03, ind11 + 10)
          ind16 = min(MSA03, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA03, id1(k) + js1 )
          ind11 = min(MSA03, ind12 + 1 )
          ind13 = min(MSA03, ind12 - 1 )
          ind14 = min(MSA03, ind12 + 10)
          ind15 = min(MSA03, ind12 + 9 )
          ind16 = min(MSA03, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA03, id1(k) + js1 )
          ind12 = min(MSA03, ind11 + 1 )
          ind13 = min(MSA03, ind11     )
          ind14 = min(MSA03, ind11 + 9 )
          ind15 = min(MSA03, ind11 + 10)
          ind16 = min(MSA03, ind11     )
        endif

        do j = 1, NG03
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          n2om1   = aman2o(jmn2o,indm  ,j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm  ,j) - aman2o(jmn2o,indm  ,j))
          n2om2   = aman2o(jmn2o,indm+1,j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm+1,j) - aman2o(jmn2o,indm+1,j))
          absn2o  = n2om1 + minorfrac(k) *(n2om2 - n2om1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + adjcoln2o*absn2o + tauaer(k,3)

          fracs(k,j) = fracrefa(j,jpl)                                  &
     &               + fpl * (fracrefa(j,jpl+1)-fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        speccomb = colamt(k,1) + rat_h2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 4.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 4.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        fac000 = (1. - fs) * fac00(k)
        fac010 = (1. - fs) * fac10(k)
        fac100 =       fs  * fac00(k)
        fac110 =       fs  * fac10(k)

        fac001 = (1. - fs1) * fac01(k)
        fac011 = (1. - fs1) * fac11(k)
        fac101 =       fs1  * fac01(k)
        fac111 =       fs1  * fac11(k)

        speccomb_mn2o = colamt(k,1) + refrat_m_b*colamt(k,2)
        specmult_mn2o = 4.0 * min(oneminus, colamt(k,1)/speccomb_mn2o)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - int(specmult_mn2o)

        fmn2omf = minorfrac(k) * fmn2o

!     in atmospheres where the amount of n2o is too great to be considered
!     a minor species, adjust the column amount of n2o by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratn2o  = colamt(k,4) / (coldry(k) * chi_mls(4,jpp1))
        if (ratn2o > 1.5) then
          adjcoln2o = (0.5 + (ratn2o - 0.5)**0.65)                      &
     &              * chi_mls(4,jpp1) * coldry(k)
        else
          adjcoln2o = colamt(k,4)
        endif

        speccomb_planck = colamt(k,1) + refrat_planck_b*colamt(k,2)
        specmult_planck = 4.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        indf = indfor(k)
        indm = indminor(k)

        ind01 = id0(k) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k) + js1
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        do j = 1, NG03
          taufor = forfac(k) * ( forref(indf,j) + forfrac(k)            &
     &           * (forref(indf+1,j) - forref(indf,j)) )
          n2om1  = ambn2o(jmn2o,indm  ,j) + fmn2o                       &
     &           * (ambn2o(jmn2o+1,indm  ,j) - ambn2o(jmn2o,indm  ,j))
          n2om2  = ambn2o(jmn2o,indm+1,j) + fmn2o                       &
     &           * (ambn2o(jmn2o+1,indm+1,j) - ambn2o(jmn2o,indm+1,j))
          absn2o = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absb(ind01,j) + fac100*absb(ind02,j)           &
     &         +  fac010*absb(ind03,j) + fac110*absb(ind04,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absb(ind11,j) + fac101*absb(ind12,j)           &
     &         +  fac011*absb(ind13,j) + fac111*absb(ind14,j))          &
     &         + taufor + adjcoln2o*absn2o + tauaer(k,3)

          fracs(k,j) = fracrefb(j,jpl)                                  &
     &         + fpl * (fracrefb(j,jpl+1) - fracrefb(j,jpl))
        enddo
      enddo

      return
!...................................
      end subroutine taugb03
!-----------------------------------


!-----------------------------------
      subroutine taugb04
!...................................

!     band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
      use module_radlw_kgb04
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, js, js1, jpl
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck,               &
     &       refrat_planck_a, refrat_planck_b, fs, fs1, fpl,            &
     &       fk00, fk10, fk20, fk01, fk11, fk21, tauself, taufor,       &
     &       tem1, tem2

!  ---  p = 142.5940 mb
      refrat_planck_a = chi_mls(1,11) / chi_mls(2,11)

!  ---  p = 95.58350 mb
      refrat_planck_b = chi_mls(3,13) / chi_mls(2,13)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA04, id1(k) + js1 )
          ind12 = min(MSA04, ind11 + 1 )
          ind13 = min(MSA04, ind11 + 2 )
          ind14 = min(MSA04, ind11 + 9 )
          ind15 = min(MSA04, ind11 + 10)
          ind16 = min(MSA04, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA04, id1(k) + js1 )
          ind11 = min(MSA04, ind12 + 1 )
          ind13 = min(MSA04, ind12 - 1 )
          ind14 = min(MSA04, ind12 + 10)
          ind15 = min(MSA04, ind12 + 9 )
          ind16 = min(MSA04, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA04, id1(k) + js1 )
          ind12 = min(MSA04, ind11 + 1 )
          ind13 = min(MSA04, ind11     )
          ind14 = min(MSA04, ind11 + 9 )
          ind15 = min(MSA04, ind11 + 10)
          ind16 = min(MSA04, ind11     )
        endif


        do j = 1, NG04
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + tauaer(k,4)

          fracs(k,j) = fracrefa(j,jpl)                                  &
     &               + fpl * (fracrefa(j,jpl+1)-fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        speccomb = colamt(k,3) + rat_o3co2(k,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,3) + rat_o3co2(k,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        fac000 = (1. - fs) * fac00(k)
        fac010 = (1. - fs) * fac10(k)
        fac100 =       fs  * fac00(k)
        fac110 =       fs  * fac10(k)

        fac001 = (1. - fs1) * fac01(k)
        fac011 = (1. - fs1) * fac11(k)
        fac101 =       fs1  * fac01(k)
        fac111 =       fs1  * fac11(k)

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specmult_planck = 4.0*min(oneminus,colamt(k,3)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        ind01 = id0(k) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k) + js1
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        do j = 1, NG04
          taug(k,j) = speccomb                                          &
     &         * (fac000*absb(ind01,j) + fac100*absb(ind02,j)           &
     &         +  fac010*absb(ind03,j) + fac110*absb(ind04,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absb(ind11,j) + fac101*absb(ind12,j)           &
     &         +  fac011*absb(ind13,j) + fac111*absb(ind14,j))          &
     &         + tauaer(k,4)

          fracs(k,j) = fracrefb(j,jpl)                                  &
     &         + fpl * (fracrefb(j,jpl+1) - fracrefb(j,jpl))
        enddo

!     empirical modification to code to improve stratospheric cooling rates
!     for co2.

        taug(k, 8) = taug(k, 8) * 0.92
        taug(k, 9) = taug(k, 9) * 0.88
        taug(k,10) = taug(k,10) * 1.07
        taug(k,11) = taug(k,11) * 1.1
        taug(k,12) = taug(k,12) * 0.99
        taug(k,13) = taug(k,13) * 0.88
        taug(k,14) = taug(k,14) * 0.83

      enddo

      return
!...................................
      end subroutine taugb04
!-----------------------------------


!-----------------------------------
      subroutine taugb05
!...................................

!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!
      use module_radlw_kgb05
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jpl, jmo3
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck,               &
     &       refrat_planck_a, refrat_planck_b, fs, fs1, fpl,            &
     &       fk00, fk10, fk20, fk01, fk11, fk21, fmo3, refrat_m_a,      &
     &       speccomb_mo3, specmult_mo3, o3m1, o3m2, abso3, tauself,    &
     &       taufor, tem1, tem2

!  ---  calculate reference ratio to be used in calculation of planck
!       fraction in lower/upper atmosphere.

!  ---  p = 473.420 mb
      refrat_planck_a = chi_mls(1,5) / chi_mls(2,5)

!  ---  p = 0.2369 mb
      refrat_planck_b = chi_mls(3,43) / chi_mls(2,43)

!  ---  p = 317.3480
      refrat_m_a = chi_mls(1,7) / chi_mls(2,7)

!  ---  compute the optical depth by interpolating in ln(pressure), 
!       temperature, and appropriate species.  below laytrop, the water
!       vapor self-continuum and foreign continuum is interpolated
!       (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mo3 = colamt(k,1) + refrat_m_a*colamt(k,2)
        specmult_mo3 = 8.0*min(oneminus, colamt(k,1)/speccomb_mo3)
        jmo3 = 1 + int(specmult_mo3)
        fmo3 = specmult_mo3 - int(specmult_mo3)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA05, id1(k) + js1 )
          ind12 = min(MSA05, ind11 + 1 )
          ind13 = min(MSA05, ind11 + 2 )
          ind14 = min(MSA05, ind11 + 9 )
          ind15 = min(MSA05, ind11 + 10)
          ind16 = min(MSA05, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA05, id1(k) + js1 )
          ind11 = min(MSA05, ind12 + 1 )
          ind13 = min(MSA05, ind12 - 1 )
          ind14 = min(MSA05, ind12 + 10)
          ind15 = min(MSA05, ind12 + 9 )
          ind16 = min(MSA05, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA05, id1(k) + js1 )
          ind12 = min(MSA05, ind11 + 1 )
          ind13 = min(MSA05, ind11     )
          ind14 = min(MSA05, ind11 + 9 )
          ind15 = min(MSA05, ind11 + 10)
          ind16 = min(MSA05, ind11     )
        endif

        do j = 1, NG05
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          o3m1    = amao3(jmo3,indm  ,j) + fmo3                         &
     &            * (amao3(jmo3+1,indm  ,j) - amao3(jmo3,indm  ,j))
          o3m2    = amao3(jmo3,indm+1,j) + fmo3                         &
     &            * (amao3(jmo3+1,indm+1,j) - amao3(jmo3,indm+1,j))
          abso3   = o3m1 + minorfrac(k) *(o3m2 - o3m1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + abso3*colamt(k,3) + wx(k,1)*ccl4(j) &
     &         + tauaer(k,5)

          fracs(k,j) = fracrefa(j,jpl)                                  &
     &               + fpl * (fracrefa(j,jpl+1)-fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        speccomb = colamt(k,3) + rat_o3co2(k,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,3) + rat_o3co2(k,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        fac000 = (1. - fs) * fac00(k)
        fac010 = (1. - fs) * fac10(k)
        fac100 =       fs  * fac00(k)
        fac110 =       fs  * fac10(k)

        fac001 = (1. - fs1) * fac01(k)
        fac011 = (1. - fs1) * fac11(k)
        fac101 =       fs1  * fac01(k)
        fac111 =       fs1  * fac11(k)

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specmult_planck = 4.0*min(oneminus,colamt(k,3)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        ind01 = id0(k) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k) + js1
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        do j = 1, NG05
          taug(k,j) = speccomb                                          &
     &         * (fac000*absb(ind01,j) + fac100*absb(ind02,j)           &
     &         +  fac010*absb(ind03,j) + fac110*absb(ind04,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absb(ind11,j) + fac101*absb(ind12,j)           &
     &         +  fac011*absb(ind13,j) + fac111*absb(ind14,j))          &
     &         + wx(k,1)*ccl4(j) + tauaer(k,5)

          fracs(k,j) = fracrefb(j,jpl)                                  &
     &         + fpl * (fracrefb(j,jpl+1) - fracrefb(j,jpl))
        enddo
      enddo

      return
!...................................
      end subroutine taugb05
!-----------------------------------


!-----------------------------------
      subroutine taugb06
!...................................

!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!               (high key - nothing; high minor - cfc11, cfc12)
!
      use module_radlw_kgb06
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf, indm,    &
     &       jpp1
!
      real (kind=kind_phys) :: chi_co2, ratco2, adjfac, adjcolco2,      &
     &       tauself, taufor, absco2

!     compute the optical depth by interpolating in ln(pressure) and
!     temperature. the water vapor self-continuumand foreign continuum
!     is interpolated (in temperature) separately.  

      do k = 1, laytrop

!     in atmospheres where the amount of co2 is too great to be
!     considered a minor species, adjust the column amount of co2
!     by an empirical factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratco2 = colamt(k,2) / (coldry(k) * chi_mls(2,jpp1))
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2 - 2.0)**0.77
          adjcolco2 = adjfac * chi_mls(2,jpp1) * coldry(k)
        else
          adjcolco2 = colamt(k,2)
        endif

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA06, id1(k) + 1 )
        ind12 = min(MSA06, ind11 + 1 )

        do j = 1, NG06
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          absco2  =              ( amaco2 (indm,j) + minorfrac(k)       &
     &            * (amaco2 (indm+1,j) - amaco2 (indm,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + adjcolco2*absco2                    &
     &         + wx(k,2)*cfc11adj(j) + wx(k,3)*cfc12(j) + tauaer(k,6)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

!     nothing important goes on above laytrop in this band.

      do k = laytrop+1, NLAY
        do j = 1, NG06
          taug(k,j) = wx(k,2)*cfc11adj(j) + wx(k,3)*cfc12(j)            &
     &              + tauaer(k,6)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb06
!-----------------------------------


!-----------------------------------
      subroutine taugb07
!...................................

!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!
      use module_radlw_kgb07
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jpl, jmco2, jpp1
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck, speccomb_mco2,&
     &       specmult_mco2, adjfac, adjcolco2, ratco2, refrat_planck_a, &
     &       refrat_m_a, fs, fs1, fpl, fk00, fk10, fk20, fk01, fk11,    &
     &       fk21, fmco2, co2m1, co2m2, absco2, tauself, taufor,        &
     &       tem1, tem2

!  ---  minor gas mapping level :
!        lower - co2, p = 706.2620 mbar, t= 278.94 k
!        upper - co2, p = 12.9350 mbar, t = 234.01 k

!     calculate reference ratio to be used in calculation of planck
!     fraction in lower atmosphere.

!  ---  p = 706.2620 mb
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

!  ---  p = 706.2720 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2oo3(k,1)*colamt(k,3)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oo3(k,2)*colamt(k,3)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,3)
        specmult_mco2 = 8.0 * min(oneminus, colamt(k,1)/speccomb_mco2)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = specmult_mco2 - int(specmult_mco2)

!     in atmospheres where the amount of co2 is too great to be considered
!     a minor species, adjust the column amount of co2 by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratco2 = colamt(k,2) / (coldry(k) * chi_mls(2,jpp1))
        if (ratco2 > 3.0) then
          adjfac = 3.0 + (ratco2 - 3.0)**0.79
          adjcolco2 = adjfac * chi_mls(2,jpp1) * coldry(k)
        else
          adjcolco2 = colamt(k,2)
        endif

        speccomb_planck = colamt(k,1) + refrat_planck_a * colamt(k,3)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl= 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA07, id1(k) + js1 )
          ind12 = min(MSA07, ind11 + 1 )
          ind13 = min(MSA07, ind11 + 2 )
          ind14 = min(MSA07, ind11 + 9 )
          ind15 = min(MSA07, ind11 + 10)
          ind16 = min(MSA07, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA07, id1(k) + js1 )
          ind11 = min(MSA07, ind12 + 1 )
          ind13 = min(MSA07, ind12 - 1 )
          ind14 = min(MSA07, ind12 + 10)
          ind15 = min(MSA07, ind12 + 9 )
          ind16 = min(MSA07, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA07, id1(k) + js1 )
          ind12 = min(MSA07, ind11 + 1 )
          ind13 = min(MSA07, ind11     )
          ind14 = min(MSA07, ind11 + 9 )
          ind15 = min(MSA07, ind11 + 10)
          ind16 = min(MSA07, ind11     )
        endif

        do j = 1, NG07
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          co2m1   = amaco2(jmco2,indm  ,j) + fmco2                      &
     &            * (amaco2(jmco2+1,indm  ,j) - amaco2(jmco2,indm  ,j))
          co2m2   = amaco2(jmco2,indm+1,j) + fmco2                      &
     &            * (amaco2(jmco2+1,indm+1,j) - amaco2(jmco2,indm+1,j))
          absco2  = co2m1 + minorfrac(k) *(co2m2 - co2m1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + adjcolco2*absco2 + tauaer(k,7)

          fracs(k,j) = fracrefa(j,jpl)                                  &
     &               + fpl * (fracrefa(j,jpl+1)-fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY

!     in atmospheres where the amount of co2 is too great to be considered
!     a minor species, adjust the column amount of co2 by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratco2 = colamt(k,2) / (coldry(k) * chi_mls(2,jpp1))
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2 - 2.0)**0.79
          adjcolco2 = adjfac * chi_mls(2,jpp1) * coldry(k)
        else
          adjcolco2 = colamt(k,2)
        endif

        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG07
          absco2 = ambco2(indm,j) + minorfrac(k)                        &
     &           * (ambco2(indm+1,j) - ambco2(indm,j))

          taug(k,j) = colamt(k,3)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + adjcolco2*absco2 + tauaer(k,7)

          fracs(k,j) = fracrefb(j)
        enddo

!     empirical modification to code to improve stratospheric cooling rates
!     for o3.

        taug(k, 8) = taug(k, 8) * 0.92
        taug(k, 9) = taug(k, 9) * 0.88
        taug(k,10) = taug(k,10) * 1.07
        taug(k,11) = taug(k,11) * 1.1
        taug(k,12) = taug(k,12) * 0.99
        taug(k,13) = taug(k,13) * 0.88
        taug(k,14) = taug(k,14) * 0.83

      enddo

      return
!...................................
      end subroutine taugb07
!-----------------------------------


!-----------------------------------
      subroutine taugb08
!...................................

!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!
      use module_radlw_kgb08
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf, indm,    &
     &       jpp1
!
      real (kind=kind_phys) :: tauself, taufor, absco2, abso3, absn2o,  &
     &       ratco2, adjfac, adjcolco2

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated (in
!     temperature) separately.

      do k = 1, laytrop

!     in atmospheres where the amount of co2 is too great to be considered
!     a minor species, adjust the column amount of co2 by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratco2 = colamt(k,2) / (coldry(k) * chi_mls(2,jpp1))
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2 - 2.0)**0.65
          adjcolco2 = adjfac * chi_mls(2,jpp1) * coldry(k)
        else
          adjcolco2 = colamt(k,2)
        endif

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA08, id1(k) + 1 )
        ind12 = min(MSA08, ind11 + 1 )

        do j = 1, NG08
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          absco2  =              ( amaco2 (indm,j) + minorfrac(k)       &
     &            * (amaco2 (indm+1,j) - amaco2 (indm,j)) )
          abso3   =              ( amao3  (indm,j) + minorfrac(k)       &
     &            * (amao3  (indm+1,j) - amao3  (indm,j)) )
          absn2o  =              ( aman2o (indm,j) + minorfrac(k)       &
     &            * (aman2o (indm+1,j) - aman2o (indm,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + adjcolco2*absco2                    &
     &         + colamt(k,3)*abso3 + colamt(k,4)*absn2o                 &
     &         + wx(k,3)*cfc12(j) + wx(k,4)*cfc22adj(j) + tauaer(k,8)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY

!     in atmospheres where the amount of co2 is too great to be considered
!     a minor species, adjust the column amount of co2 by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratco2 = colamt(k,2) / (coldry(k) * chi_mls(2,jpp1))
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2 - 2.0)**0.65
          adjcolco2 = adjfac * chi_mls(2,jpp1) * coldry(k)
        else
          adjcolco2 = colamt(k,2)
        endif

        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG08
          absco2 = (ambco2(indm,j) + minorfrac(k)                       &
     &           * (ambco2(indm+1,j) - ambco2(indm,j)))
          absn2o = (ambn2o(indm,j) + minorfrac(k)                       &
     &           * (ambn2o(indm+1,j) - ambn2o(indm,j)))

          taug(k,j) = colamt(k,3)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + adjcolco2*absco2 + colamt(k,4)*absn2o                  &
     &         + wx(k,3)*cfc12(j) + wx(k,4)*cfc22adj(j) + tauaer(k,8)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb08
!-----------------------------------


!-----------------------------------
      subroutine taugb09
!...................................

!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
!
      use module_radlw_kgb09
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jpl, jmn2o, jpp1
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck, speccomb_mn2o,&
     &       specmult_mn2o, adjfac, adjcoln2o, ratn2o, refrat_planck_a, &
     &       refrat_m_a, fs, fs1, fpl, fk00, fk10, fk20, fk01, fk11,    &
     &       fk21, fmn2o, n2om1, n2om2, absn2o, tauself, taufor,        &
     &       tem1, tem2

!  ---  minor gas mapping level :
!         lower - n2o, p = 706.272 mbar, t = 278.94 k
!         upper - n2o, p = 95.58 mbar, t = 215.7 k

!     calculate reference ratio to be used in calculation of planck
!     fraction in lower/upper atmosphere.

!  ---  p = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

!  ---  p = 706.272 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)


!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2och4(k,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2och4(k,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,5)
        specmult_mn2o = 8.0*min(oneminus,colamt(k,1)/speccomb_mn2o)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - int(specmult_mn2o)

!     in atmospheres where the amount of n2o is too great to be considered
!     a minor species, adjust the column amount of n2o by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratn2o = colamt(k,4) / (coldry(k) * chi_mls(4,jpp1))
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * chi_mls(4,jpp1) * coldry(k)
        else
          adjcoln2o = colamt(k,4)
        endif

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,5)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA09, id1(k) + js1 )
          ind12 = min(MSA09, ind11 + 1 )
          ind13 = min(MSA09, ind11 + 2 )
          ind14 = min(MSA09, ind11 + 9 )
          ind15 = min(MSA09, ind11 + 10)
          ind16 = min(MSA09, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA09, id1(k) + js1 )
          ind11 = min(MSA09, ind12 + 1 )
          ind13 = min(MSA09, ind12 - 1 )
          ind14 = min(MSA09, ind12 + 10)
          ind15 = min(MSA09, ind12 + 9 )
          ind16 = min(MSA09, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA09, id1(k) + js1 )
          ind12 = min(MSA09, ind11 + 1 )
          ind13 = min(MSA09, ind11     )
          ind14 = min(MSA09, ind11 + 9 )
          ind15 = min(MSA09, ind11 + 10)
          ind16 = min(MSA09, ind11     )
        endif

        do j = 1, NG09
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          n2om1   = aman2o(jmn2o,indm,  j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm,  j) - aman2o(jmn2o,indm,  j))
          n2om2   = aman2o(jmn2o,indm+1,j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm+1,j) - aman2o(jmn2o,indm+1,j))
          absn2o  = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + adjcoln2o*absn2o + tauaer(k,9)

          fracs(k,j) = fracrefa(j,jpl) + fpl                            &
     &               * (fracrefa(j,jpl+1) - fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY

!     in atmospheres where the amount of n2o is too great to be considered
!     a minor species, adjust the column amount of n2o by an empirical
!     factor to obtain the proper contribution.

        jpp1 = jp(k) + 1
        ratn2o = colamt(k,4) / (coldry(k) * chi_mls(4,jpp1))
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * chi_mls(4,jpp1) * coldry(k)
        else
          adjcoln2o = colamt(k,4)
        endif

        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG09
          absn2o = ambn2o(indm,j) + minorfrac(k)                        &
     &           * (ambn2o(indm+1,j) - ambn2o(indm,j))

          taug(k,j) = colamt(k,5)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + adjcoln2o*absn2o + tauaer(k,9)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb09
!-----------------------------------


!-----------------------------------
      subroutine taugb10
!...................................

!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!
      use module_radlw_kgb10
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf
!
      real (kind=kind_phys) :: tauself, taufor

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  below laytrop, the water vapor self-continuum and
!     foreign continuum is interpolated (in temperature) separately. 

      do k = 1, laytrop
        inds = indself(k)
        indf = indfor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA10, id1(k) + 1 )
        ind12 = min(MSA10, ind11 + 1 )

        do j = 1, NG10
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + tauaer(k,10)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY
        indf = indfor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG10
          taufor = forfac(k) * ( forref(indf,j) + forfrac(k)            &
     &           * (forref(indf+1,j) - forref(indf,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + taufor + tauaer(k,10)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb10
!-----------------------------------


!-----------------------------------
      subroutine taugb11
!...................................

!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!
      use module_radlw_kgb11
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf, indm
!
      real (kind=kind_phys) :: tauself, taufor, tauo2, scaleo2
!
!  ---  minor gas mapping level :
!          lower - o2, p = 706.2720 mbar, t = 278.94 k
!          upper - o2, p = 4.758820 mbarm t = 250.85 k

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  below laytrop, the water vapor self-continuum and
!     foreign continuum is interpolated (in temperature) separately.  

      do k = 1, laytrop
        scaleo2 = colamt(k,6)*scaleminor(k)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA11, id1(k) + 1 )
        ind12 = min(MSA11, ind11 + 1 )

        do j = 1, NG11
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          tauo2   = scaleo2    * ( amao2  (indm,j) + minorfrac(k)       &
     &            * (amao2  (indm+1,j) - amao2  (indm,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + tauo2 + tauaer(k,11)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY
        scaleo2 = colamt(k,6)*scaleminor(k)

        indf = indfor(k)
        indm = indminor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG11
          taufor = forfac(k) * ( forref(indf,j) + forfrac(k)            &
     &           * (forref(indf+1,j) - forref(indf,j)) )
          tauo2  = scaleo2   * ( ambo2 (indm,j) +  minorfrac(k)         &
     &           * (ambo2 (indm+1,j) - ambo2 (indm,j)) )

          taug(k,j) = colamt(k,1)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + taufor + tauo2 + tauaer(k,11)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb11
!-----------------------------------


!-----------------------------------
      subroutine taugb12
!...................................

!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
      use module_radlw_kgb12
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, js, js1, jpl
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck,               &
     &       refrat_planck_a, fs, fs1, fpl, fk00, fk10, fk20, fk01,     &
     &       fk11, fk21, tauself, taufor, tem1, tem2

!  ---  calculate reference ratio to be used in calculation of planck
!       fraction in lower/upper atmosphere.

!  ---  p = 174.164 mb
      refrat_planck_a = chi_mls(1,10) / chi_mls(2,10)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA12, id1(k) + js1 )
          ind12 = min(MSA12, ind11 + 1 )
          ind13 = min(MSA12, ind11 + 2 )
          ind14 = min(MSA12, ind11 + 9 )
          ind15 = min(MSA12, ind11 + 10)
          ind16 = min(MSA12, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA12, id1(k) + js1 )
          ind11 = min(MSA12, ind12 + 1 )
          ind13 = min(MSA12, ind12 - 1 )
          ind14 = min(MSA12, ind12 + 10)
          ind15 = min(MSA12, ind12 + 9 )
          ind16 = min(MSA12, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA12, id1(k) + js1 )
          ind12 = min(MSA12, ind11 + 1 )
          ind13 = min(MSA12, ind11     )
          ind14 = min(MSA12, ind11 + 9 )
          ind15 = min(MSA12, ind11 + 10)
          ind16 = min(MSA12, ind11     )
        endif

        do j = 1, NG12
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + tauaer(k,12)

          fracs(k,j) = fracrefa(j,jpl) + fpl                            &
     &               * (fracrefa(j,jpl+1) - fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        do j = 1, NG12
          taug(k,j)  = tauaer(k,12)
          fracs(k,j) = f_zero
        enddo
      enddo

      return
!...................................
      end subroutine taugb12
!-----------------------------------


!-----------------------------------
      subroutine taugb13
!...................................

!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!
      use module_radlw_kgb13
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jpl, jmco2, jmco
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck, speccomb_mco2,&
     &       specmult_mco2, adjfac, adjcolco2, ratco2, refrat_planck_a, &
     &       refrat_m_a, refrat_m_a3, fs, fs1, fpl, fk00, fk10, fk20,   &
     &       fk01, fk11, fk21, fmco2, co2m1, co2m2, absco2, tauself,    &
     &       taufor, speccomb_mco, specmult_mco, fmco, com1, com2,      &
     &       absco, abso3, tem0, tem1, tem2

!  ---  minor gas mapping levels :
!          lower - co2, p = 1053.63 mb, t = 294.2 k
!          lower - co,  p = 706 mb,     t = 278.94 k
!          upper - o3,  p = 95.5835 mb, t = 215.7 k

!  ---  calculate reference ratio to be used in calculation of planck
!       fraction in lower/upper atmosphere.

!  ---  p = 473.420 mb (level 5)
      refrat_planck_a = chi_mls(1,5) / chi_mls(4,5)

!  ---  p = 1053. (level 1)
      refrat_m_a = chi_mls(1,1) / chi_mls(4,1)

!  ---  p = 706. (level 3)
      refrat_m_a3 = chi_mls(1,3) / chi_mls(4,3)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2on2o(k,1)*colamt(k,4)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2on2o(k,2)*colamt(k,4)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,4)
        specmult_mco2 = 8.0 * min(oneminus, colamt(k,1)/speccomb_mco2)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = specmult_mco2 - int(specmult_mco2)

!  ---  in atmospheres where the amount of co2 is too great to be
!       considered a minor species, adjust the column amount of co2 by
!       an empirical factor to obtain the proper contribution.

        tem0 = coldry(k) * 3.55e-4
        ratco2 = colamt(k,2) / tem0
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2 - 2.0)**0.68
          adjcolco2 = adjfac * tem0
        else
          adjcolco2 = colamt(k,2)
        endif

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,4)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)


        speccomb_mco = colamt(k,1) + refrat_m_a3*colamt(k,4)
        specmult_mco = 8.0 * min(oneminus, colamt(k,1)/speccomb_mco)
        jmco = 1 + int(specmult_mco)
        fmco = specmult_mco - int(specmult_mco)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA13, id1(k) + js1 )
          ind12 = min(MSA13, ind11 + 1 )
          ind13 = min(MSA13, ind11 + 2 )
          ind14 = min(MSA13, ind11 + 9 )
          ind15 = min(MSA13, ind11 + 10)
          ind16 = min(MSA13, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA13, id1(k) + js1 )
          ind11 = min(MSA13, ind12 + 1 )
          ind13 = min(MSA13, ind12 - 1 )
          ind14 = min(MSA13, ind12 + 10)
          ind15 = min(MSA13, ind12 + 9 )
          ind16 = min(MSA13, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA13, id1(k) + js1 )
          ind12 = min(MSA13, ind11 + 1 )
          ind13 = min(MSA13, ind11     )
          ind14 = min(MSA13, ind11 + 9 )
          ind15 = min(MSA13, ind11 + 10)
          ind16 = min(MSA13, ind11     )
        endif

        do j = 1, NG13
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )
          co2m1   = amaco2(jmco2,indm,  j) + fmco2                      &
     &            * (amaco2(jmco2+1,indm,  j) - amaco2(jmco2,indm,  j))
          co2m2   = amaco2(jmco2,indm+1,j) + fmco2                      &
     &            * (amaco2(jmco2+1,indm+1,j) - amaco2(jmco2,indm+1,j))
          absco2  = co2m1 + minorfrac(k) * (co2m2 - co2m1)
          com1    = amaco(jmco,indm,  j) + fmco                         &
     &            * (amaco(jmco+1,indm,  j) - amaco(jmco,indm,  j))
          com2    = amaco(jmco,indm+1,j) + fmco                         &
     &            * (amaco(jmco+1,indm+1,j) - amaco(jmco,indm+1,j))
          absco   = com1 + minorfrac(k) * (com2 - com1)

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor                                       &
     &         + adjcolco2*absco2 + colamt(k,7)*absco + tauaer(k,13)

          fracs(k,j) = fracrefa(j,jpl) + fpl                            &
     &               * (fracrefa(j,jpl+1) - fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        do j = 1, NG13
          abso3 = ambo3(indm,j) + minorfrac(k)                          &
     &          * (ambo3(indm+1,j) - ambo3(indm,j))

          taug(k,j) = colamt(k,3) * abso3 + tauaer(k,13)

          fracs(k,j) =  fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb13
!-----------------------------------


!-----------------------------------
      subroutine taugb14
!...................................

!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
      use module_radlw_kgb14
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind11, ind12, inds, indf
!
      real (kind=kind_phys) :: tauself, taufor

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  below laytrop, the water vapor self-continuum 
!     and foreign continuum is interpolated (in temperature) separately.  

      do k = 1, laytrop
        inds = indself(k)
        indf = indfor(k)

        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA14, id1(k) + 1 )
        ind12 = min(MSA14, ind11 + 1 )

        do j = 1, NG14
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = colamt(k,2)                                       &
     &         * (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &         +  fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))      &
     &         + tauself + taufor + tauaer(k,14)

          fracs(k,j) = fracrefa(j)
        enddo
      enddo

      do k = laytrop+1, NLAY
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = id1(k) + 1
        ind12 = ind11 + 1

        do j = 1, NG14
          taug(k,j) = colamt(k,2)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + tauaer(k,14)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb14
!-----------------------------------


!-----------------------------------
      subroutine taugb15
!...................................

!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!
      use module_radlw_kgb15
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, indm,       &
     &       js, js1, jpl, jmn2o
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck, speccomb_mn2o,&
     &       specmult_mn2o, adjfac, adjcoln2o, ratn2o, refrat_planck_a, &
     &       refrat_m_a, scalen2o, fs, fs1, fpl, fk00, fk10, fk20, fk01,&
     &       fk11, fk21, fmn2o, n2om1, n2om2, absn2o, tauself, taufor,  &
     &       taun2o, tem0, tem1, tem2

!  ---  minor gas mapping level :
!          lower - nitrogen continuum, p = 1053., t = 294.

!  ---  calculate reference ratio to be used in calculation of planck
!       fraction in lower atmosphere.

!  ---  p = 1053. mb (level 1)
      refrat_planck_a = chi_mls(4,1) / chi_mls(2,1)

!  ---  p = 1053.
      refrat_m_a = chi_mls(4,1) / chi_mls(2,1)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,4) + rat_n2oco2(k,1)*colamt(k,2)
        specparm = colamt(k,4) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,4) + rat_n2oco2(k,2)*colamt(k,2)
        specparm1 = colamt(k,4) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_mn2o = colamt(k,4) + refrat_m_a*colamt(k,2)
        specmult_mn2o = 8.0 * min(oneminus, colamt(k,4)/speccomb_mn2o)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - int(specmult_mn2o)

        speccomb_planck = colamt(k,4) + refrat_planck_a*colamt(k,2)
        specmult_planck = 8.0*min(oneminus,colamt(k,4)/speccomb_planck)
        jpl= 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        scalen2o = colbrd(k) * scaleminor(k)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA15, id1(k) + js1 )
          ind12 = min(MSA15, ind11 + 1 )
          ind13 = min(MSA15, ind11 + 2 )
          ind14 = min(MSA15, ind11 + 9 )
          ind15 = min(MSA15, ind11 + 10)
          ind16 = min(MSA15, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA15, id1(k) + js1 )
          ind11 = min(MSA15, ind12 + 1 )
          ind13 = min(MSA15, ind12 - 1 )
          ind14 = min(MSA15, ind12 + 10)
          ind15 = min(MSA15, ind12 + 9 )
          ind16 = min(MSA15, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA15, id1(k) + js1 )
          ind12 = min(MSA15, ind11 + 1 )
          ind13 = min(MSA15, ind11     )
          ind14 = min(MSA15, ind11 + 9 )
          ind15 = min(MSA15, ind11 + 10)
          ind16 = min(MSA15, ind11     )
        endif

        do j = 1, NG15
!err      tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
!    &            * (selfref(inds+1,j) - selfref(inds,j)) )
!         taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
!    &            * (forref (indf+1,j) - forref (indf,j)) )
          n2om1   = aman2o(jmn2o,indm,  j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm,  j) - aman2o(jmn2o,indm,  j))
          n2om2   = aman2o(jmn2o,indm+1,j) + fmn2o                      &
     &            * (aman2o(jmn2o+1,indm+1,j) - aman2o(jmn2o,indm+1,j))
          taun2o  = scalen2o * (n2om1 + minorfrac(k) * (n2om2 - n2om1))

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
!err &         + tauself + taufor
     &         + taun2o + tauaer(k,15)

          fracs(k,j) = fracrefa(j,jpl) + fpl                            &
     &               * (fracrefa(j,jpl+1)-fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        do j = 1, NG15
          taug(k,j)  = tauaer(k,15)

          fracs(k,j) = f_zero
        enddo
      enddo

      return
!...................................
      end subroutine taugb15
!-----------------------------------


!-----------------------------------
      subroutine taugb16
!...................................

!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!
      use module_radlw_kgb16
!
      implicit none
!
!  ---  local variables:
      integer :: j, k, ind01, ind02, ind03, ind04, ind05, ind06, ind11, &
     &       ind12, ind13, ind14, ind15, ind16, inds, indf, js, js1, jpl
!
      real (kind=kind_phys) :: fac000, fac010, fac100, fac110, fac001,  &
     &       fac011, fac101, fac111, fac200, fac201, fac210, fac211,    &
     &       speccomb, specmult, specparm, speccomb1, specmult1,        &
     &       specparm1, speccomb_planck, specmult_planck,               &
     &       refrat_planck_a, fs, fs1, fpl, fk00, fk10, fk20, fk01,     &
     &       fk11, fk21, tauself, taufor, tem1, tem2

!  ---  calculate reference ratio to be used in calculation of planck
!       fraction in lower atmosphere.

!  ---  p = 387. mb (level 6)
      refrat_planck_a = chi_mls(1,6) / chi_mls(6,6)

!     compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  below laytrop, the water
!     vapor self-continuum and foreign continuum is interpolated
!     (in temperature) separately.  

      do k = 1, laytrop
        speccomb = colamt(k,1) + rat_h2och4(k,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(oneminus, specparm)
        js = 1 + int(specmult)
        fs = specmult - int(specmult)

        speccomb1 = colamt(k,1) + rat_h2och4(k,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(oneminus, specparm1)
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - int(specmult1)

        speccomb_planck = colamt(k,1)+refrat_planck_a*colamt(k,5)
        specmult_planck = 8.0*min(oneminus,colamt(k,1)/speccomb_planck)
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - int(specmult_planck)

        inds = indself(k)
        indf = indfor(k)

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          tem1 = fs - 1.0
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = fs1 - 1.0
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01 + 2
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01 + 11

          ind11 = min(MSA16, id1(k) + js1 )
          ind12 = min(MSA16, ind11 + 1 )
          ind13 = min(MSA16, ind11 + 2 )
          ind14 = min(MSA16, ind11 + 9 )
          ind15 = min(MSA16, ind11 + 10)
          ind16 = min(MSA16, ind11 + 11)
        else if (specparm>0.875 .and. specparm1>0.875) then
          tem1 = -fs
          fk00 = tem1**4
          fk10 = 1.0 - tem1 - 2.0*fk00
          fk20 = tem1 + fk00

          fac000 = fk00*fac00(k)
          fac100 = fk10*fac00(k)
          fac200 = fk20*fac00(k)
          fac010 = fk00*fac10(k)
          fac110 = fk10*fac10(k)
          fac210 = fk20*fac10(k)

          tem2 = -fs1
          fk01 = tem2**4
          fk11 = 1.0 - tem2 - 2.0*fk01
          fk21 = tem2 + fk01

          fac001 = fk01*fac01(k)
          fac101 = fk11*fac01(k)
          fac201 = fk21*fac01(k)
          fac011 = fk01*fac11(k)
          fac111 = fk11*fac11(k)
          fac211 = fk21*fac11(k)

          ind02 = id0(k) + js
          ind01 = ind02 + 1
          ind03 = ind02 - 1
          ind04 = ind02 + 10
          ind05 = ind02 + 9
          ind06 = ind02 + 8

          ind12 = min(MSA16, id1(k) + js1 )
          ind11 = min(MSA16, ind12 + 1 )
          ind13 = min(MSA16, ind12 - 1 )
          ind14 = min(MSA16, ind12 + 10)
          ind15 = min(MSA16, ind12 + 9 )
          ind16 = min(MSA16, ind12 + 8 )
        else
          fac000 = (1. - fs) * fac00(k)
          fac100 =       fs  * fac00(k)
          fac200 = f_zero
          fac010 = (1. - fs) * fac10(k)
          fac110 =       fs  * fac10(k)
          fac210 = f_zero

          fac001 = (1. - fs1) * fac01(k)
          fac101 =       fs1  * fac01(k)
          fac201 = f_zero
          fac011 = (1. - fs1) * fac11(k)
          fac111 =       fs1  * fac11(k)
          fac211 = f_zero

          ind01 = id0(k) + js
          ind02 = ind01 + 1
          ind03 = ind01
          ind04 = ind01 + 9
          ind05 = ind01 + 10
          ind06 = ind01

          ind11 = min(MSA16, id1(k) + js1 )
          ind12 = min(MSA16, ind11 + 1 )
          ind13 = min(MSA16, ind11     )
          ind14 = min(MSA16, ind11 + 9 )
          ind15 = min(MSA16, ind11 + 10)
          ind16 = min(MSA16, ind11     )
        endif

        do j = 1, NG16
          tauself = selffac(k) * ( selfref(inds,j) + selffrac(k)        &
     &            * (selfref(inds+1,j) - selfref(inds,j)) )
          taufor  = forfac (k) * ( forref (indf,j) + forfrac (k)        &
     &            * (forref (indf+1,j) - forref (indf,j)) )

          taug(k,j) = speccomb                                          &
     &         * (fac000*absa(ind01,j) + fac100*absa(ind02,j)           &
     &         +  fac200*absa(ind03,j) + fac010*absa(ind04,j)           &
     &         +  fac110*absa(ind05,j) + fac210*absa(ind06,j))          &
     &         +      speccomb1                                         &
     &         * (fac001*absa(ind11,j) + fac101*absa(ind12,j)           &
     &         +  fac201*absa(ind13,j) + fac011*absa(ind14,j)           &
     &         +  fac111*absa(ind15,j) + fac211*absa(ind16,j))          &
     &         + tauself + taufor + tauaer(k,16)

          fracs(k,j) = fracrefa(j,jpl) + fpl                            &
     &         * (fracrefa(j,jpl+1) - fracrefa(j,jpl))
        enddo
      enddo

      do k = laytrop+1, NLAY
        ind01 = id0(k) + 1
        ind02 = ind01 + 1
        ind11 = min(MSA16, id1(k) + 1 )
        ind12 = min(MSA16, ind11 + 1 )

        do j = 1, NG16
          taug(k,j) = colamt(k,5)                                       &
     &         * (fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &         +  fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j))      &
     &         + tauaer(k,16)

          fracs(k,j) = fracrefb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taugb16
!-----------------------------------


!...................................
      end subroutine taumol
!-----------------------------------


!
!........................................!
      end module module_radlw_main       !
!========================================!

