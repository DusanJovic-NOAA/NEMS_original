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
!      'swrad'      -- main gfdl1 sw radiation routine                 !
!         inputs:                                                      !
!           ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                    !
!             clouds,iovr,aerosols,sfcalb,                             !
!             cosz,solcon,NPTS,idxday,                                 !
!             IMAX, NLAY, NLP1, iflip, lprnt,                          !
!         outputs:                                                     !
!             hswc,topflx,sfcflx,                                      !
!!        optional outputs:                                            !
!             HSW0,HSWB,FLXPRF,FDNCMP                                  !
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
!!    3. radiation flux profiles(from module 'module_radsw_parameters')!
!!         profsw_type    -  derived data type for rad vertical prof   !
!!           upfxc              total sky level upward flux            !
!!           dnfxc              total sky level downward flux          !
!!           upfx0              clear sky level upward flux            !
!!           dnfx0              clear sky level downward flux          !
!                                                                      !
!!    4. surface component fluxes(from module 'module_radsw_parameters'!
!!         cmpfsw_type    -  derived data type for component sfc flux  !
!!           uvbfc              total sky downward uv-b flux at sfc    !
!!           uvbf0              clear sky downward uv-b flux at sfc    !
!!           nirbm              surface downward nir direct beam flux  !
!!           nirdf              surface downward nir diffused flux     !
!!           visbm              surface downward uv+vis direct beam flx!
!!           visdf              surface downward uv+vis diffused flux  !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!                                                                      !
!       'module machine'                                               !
!       'module physcons'                                              !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_gfdl1_param.f'                                           !
!      'radsw_gfdl1_datatb.f'                                          !
!      'radsw_gfdl1_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use sw modules     !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!   descriptions of the program:                                       !
!                                                                      !
!     esf -- 'exponential-sum fit' sw radiation scheme                 !
!                                                                      !
!   swresf uses the delta-eddington technique in conjunction with a    !
!   multiple-band parameterization for h2o+co2+o2+o3 absorption to     !
!   derive solar fluxes and heating rates.                             !
!                                                                      !
!                                                                      !
!   reference:                                                         !
!                                                                      !
!     freidenreich and ramaswamy, 1999: a new multiple-band solar      !
!     radiative parameterization for general circulation models.       !
!     j. geophys. res., 104, 31389-31409.                              !
!                                                                      !
!                                                                      !
!                                                                      !
!   ncep modifications history log:                                    !
!                                                                      !
!       sep 2000,  ken campana  -- received original code from gfdl    !
!       apr 2005,  yu-tai hou                                          !
!                  modified for ncep model applications                !
!       may 2006,  yu-tai hou                                          !
!                  recoded to fit in ncep unified radiation package    !
!       apr 2007,  yu-tai hou                                          !
!                  add spectral band heating as optional output        !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radsw_main           !
!........................................!
!
      use machine,                 only : kind_phys
      use physcons,                only : con_pi, con_rd, con_g,        &
     &                                    con_cp, con_p0

      use module_radsw_parameters
      use module_radsw_cntr_para

      use module_radsw_bandtbl
      use module_radsw_cldprtb
!
      implicit   none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGSW='GFDL-SW RESF.02 mar 2005'
      character(24), parameter :: VTAGSW='GFDL-SW RESF.02 Apr 2007'

!  --- ...  constant values
      real (kind=kind_phys), parameter :: rhoair = 1.292269
      real (kind=kind_phys), parameter :: _zero  = 0.0
      real (kind=kind_phys), parameter :: _one   = 1.0

!  --- ...  parameters for determining rayleigh optical depth
 
      real (kind=kind_phys), parameter :: temprefray = 288.15
      real (kind=kind_phys), parameter :: densmolref = 2.54743e+19
      real (kind=kind_phys), parameter :: refquanray = densmolref       &
     &                                      * temprefray / con_p0
      real (kind=kind_phys), parameter :: convfac    = 1.0e+18
      real (kind=kind_phys), parameter :: depfac     = 1.39e-02      ! depolorization factor

!     real (kind=kind_phys), parameter :: vers_num = 0.02

!  --- ...  the following data will be set up only once by "rswinit"

!  ---  bdens   is the quantity which multiples the molecular density to
!               yield the rayleigh scattering coefficient
!       coszstr is the gaussian angles for evaluation of the diffuse beam
!       radcon  is the factor for heating rates (in k/day, or k/sec)

      real (kind=kind_phys) :: bdens(NBANDS), coszstr(NSTREAMS),        &
     &       gausswt(NSOLWG), radcon

!  ---  c4co2.. and c4o2.. are coefficients for co2, and o2 optical depth
!       solflxtotal is the total solar flux

      real (kind=kind_phys), dimension(NH2OBANDS)  :: c4co2, c4co2str,  &
     &       c4o2, c4o2str
      real (kind=kind_phys), dimension(TOT_WVNUMS) :: solarfluxtoa

      real (kind=kind_phys) :: c4o2strschrun, solflxtotal

!  ---  define the solar weights and interval counters that are used to 
!       determine the single-scattering properties for the parameterization 
!       band spectral intervals, from the specified spectral intervals for
!       cloud water drops, ice particles, rain, and snow.

      integer, dimension(NBANDS)  ::  nv1liq, nv2liq,  nv1ice, nv2ice,  &
     &       nv1rain, nv2rain, nv1snow, nv2snow

      real (kind=kind_phys), dimension(NBANDS,NLIQCLDV)  :: solvliq
      real (kind=kind_phys), dimension(NBANDS,NICECLDV)  :: solvice
      real (kind=kind_phys), dimension(NBANDS,NRAINCLDV) :: solvrain
      real (kind=kind_phys), dimension(NBANDS,NSNOWCLDV) :: solvsnow

!! ...  logical flags for optional output fields

      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.

!  --- ...  interfaces

      public  swrad, rswinit
 

! =================
      contains
! =================


!-----------------------------------
      subroutine swrad                                                  &
!...................................

!  ---  inputs:
     &    ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                       &
     &      clouds,IOVR,aerosols,sfcalb,                                &
     &      cosz,solcon,NPTS,idxday,                                    &
     &      IMAX, NLAY, NLP1, iflip, lprnt,                             &
!  ---  outputs:
     &      hswc,topflx,sfcflx                                          &
!! ---  optional:
     &,     HSW0,HSWB,FLXPRF,FDNCMP                                     &
     &    )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!   plyr (IMAX,NLAY) : model layer mean pressure in mb                  !
!   plvl (IMAX,NLP1) : model level pressure in mb                       !
!   tlyr (IMAX,NLAY) : model layer mean temperature in k                !
!   tlvl (IMAX,NLP1) : model level temperature in k    (not in use)     !
!   qlyr (IMAX,NLAY) : layer h2o mass mixing ratio in gm/gm  *see inside!
!   olyr (IMAX,NLAY) : layer ozone mass mixing ratio in gm/gm           !
!   gasvmr(IMAX,NLAY,:): atmospheric constent gases:                    !
!                      (check module_radiation_gases for definition)    !
!      gasvmr(:,:,1)  - co2 volume mixing ratio                         !
!      gasvmr(:,:,2)  - n2o volume mixing ratio        (not used)       !
!      gasvmr(:,:,3)  - ch4 volume mixing ratio        (not used)       !
!      gasvmr(:,:,4)  - o2  volume mixing ratio                         !
!      gasvmr(:,:,5)  - co  volume mixing ratio        (not used)       !
!      gasvmr(:,:,6)  - cfc11 volume mixing ratio      (not used)       !
!      gasvmr(:,:,7)  - cfc12 volume mixing ratio      (not used)       !
!      gasvmr(:,:,8)  - cfc22 volume mixing ratio      (not used)       !
!      gasvmr(:,:,9)  - ccl4  volume mixing ratio      (not used)       !
!      gasvmr(:,:,10) - cfc113 volume mixing ratio     (not used)       !
!   clouds(IMAX,NLAY,:): cloud profile                                  !
!                      (check module_radiation_clouds for definition)   !
!                ---  for  iflagliq > 0  ---                            !
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
!                ---  for  iflagliq = 0  ---                            !
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
!  control parameters in module "module_radsw_cntr_para":               !
!   iswrate: heating rate unit selections                               !
!            =1: output in k/day                                        !
!            =2: output in k/second                                     !
!   iaersw : flags for aerosols effect                                  !
!            =0: without aerosol effect                                 !
!            >0: include aerosol effect                                 !
!   irgassw: control flag for rare gases (o2 etc.)                      !
!            =0: do not include rare gases                              !
!            =1: include all rare gases                                 !
!                                                                       !
!                                                                       !
!  output variables:                                                    !
!   hswc  (IMAX,NLAY): total sky heating rates (k/sec or k/day)         !
!   topflx(IMAX)     : radiation fluxes at toa (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at toa                   !
!     dnflx            - total sky downward flux at toa                 !
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
!     visbm            - downward surface uv+vis direct beam flux       !
!     visdf            - downward surface uv+vis diffused flux          !
!                                                                       !
!  module parameters, control and local variables:                      !
!                                                                       !
!  =====================    end of definitions    ====================  !
!
!----------------------------------------------------------------------c
! intent out:                                                          c
!                                                                      c
! dfsw =  downward radiation at all pressure levels                    c
! fsw  =  net radiation (up-down) at all pressure levels               c
! hsw  =  radiation heating rates at all pressure layers               c
! ufsw =  upward radiation at all pressure levels                      c
!----------------------------------------------------------------------c
!
      implicit   none
 
!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1, IOVR, iflip, NPTS
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

!  ---  local variables:
      real (kind=kind_phys), dimension(NPTS,NLP1,NSOLWG) :: totco2,     &
     &       totco2str, toto2, toto2str
      real (kind=kind_phys), dimension(NPTS,NLAY,NSOLWG) :: efftauco2,  &
     &       efftauo2
      real (kind=kind_phys), dimension(NPTS,NLAY,NBANDS) :: cldext,     &
     &       cldsct, cldasymm, aerextop, aerssalb, aerasymm

      real (kind=kind_phys), dimension(NPTS,NLAY) :: rrvco2, o2mixrt,   &
     &       cliqp, cicep, crain, csnow, reliq, reice, rerain, resnow,  &
     &       cliq1, cice1, crain1, csnow1, cldfrc, press,               &
     &       wh2o, wh2ostr, wo3, tco2, to2, deltaz, deltap, delpdig,    &
     &       opdep, densitymol, temp, rh2o, qo3

      real (kind=kind_phys), dimension(NPTS,NLAY) ::                    &
     &       rlydif,  rlydif0, rlydifc, rlydir,  rlydir0, rlydirc,      &
     &       tlydif,  tlydif0, tlydifc, tlydir,  tlydir0, tlydirc,      &
     &       tlyde,   tlyde0,  tlydec,  taustr0, omgstr0, ssalb0,       &
     &       gg0, ggc, ff0, ffc,        taustrc, omgstrc, ssalbc,       &
     &       gstr0, gstrc, extopdep0, sctopdep0, extopdepc, sctopdepc,  &
     &       aeroextopdep,  aerosctopdep,  aeroasymfac,  rayopdep,      &
     &       cloudextopdep, cloudsctopdep, cloudasymfac, gasopdep

      real (kind=kind_phys), dimension(NPTS,NLP1) :: pflux, reflec,     &
     &       transm, alphaco2, alphaco2str, alphao2, alphao2str,        &
     &       sumre, sumtr, scale, reflec0, transm0,                     &
     &       sumre0, sumtr0, dfsw, dfsw0, ufsw, ufsw0, fsw, fsw0

      real (kind=kind_phys), dimension(NPTS) :: cosz1, wtfac,           &
     &       albdir, albdif

      real (kind=kind_phys) :: tem0

!! ---  for optional outputs:
      real (kind=kind_phys), dimension(NPTS) :: sdnuvb, sdnuvd,         &
     &       sdnirb, sdnird, sumuvb, sumuvd, sumirb, sumird,            &
     &       strnbm, strndf
      real (kind=kind_phys), dimension(NPTS) :: suvbf0, suvbfc
      real (kind=kind_phys), dimension(NPTS,NLP1,NBANDS) :: fswb

      logical, dimension(NPTS,NLAY) :: lcloud
      integer  :: i, j1, k, k1, IK, nb, nf, ng, np, ns

!
!===> ...  begin here
!

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )

!  ---  initial output arrays

      hswc(:,:) = _zero
      topflx = topfsw_type ( _zero,_zero,_zero )
      sfcflx = sfcfsw_type ( _zero,_zero,_zero,_zero )

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

!  --- ... initialize local variables.                                        

      IK = NPTS * NLAY

      do k = 1, NLP1
        do i = 1, IMAX
          ufsw (i,k) = _zero
          dfsw (i,k) = _zero
          ufsw0(i,k) = _zero
          dfsw0(i,k) = _zero
        enddo
      enddo

      do i = 1, NPTS
        alphaco2   (i,1) = _zero
        alphaco2str(i,1) = _zero
        alphao2    (i,1) = _zero
        alphao2str (i,1) = _zero
      enddo

      if ( lfdncmp ) then
        do i = 1, NPTS
!! ---  optional uv-b surface downward flues
          suvbf0(i) = _zero
          suvbfc(i) = _zero

!! ---  optional surface downward beam and diff fluxes
          sdnuvb(i) = _zero
          sdnuvd(i) = _zero
          sdnirb(i) = _zero
          sdnird(i) = _zero
        enddo
      endif

!  ---  the internal array is always from top to surface

      if (iflip == 0) then        ! input from toa to sfc

        do i = 1, NPTS
          j1 = idxday(i)

          do k = 1, NLAY
            press(i,k) = 1.0e2 * plyr(j1,k)  ! convert press from mb to mks unit
            pflux(i,k) = 1.0e2 * plvl(j1,k)  ! convert press from mb to mks unit
!           pflux(i,k+1) = 1.0e2 * 0.5 * (plyr(j1,k)+plyr(j1,k+1))
            temp (i,k) = tlyr(j1,k)
!test use
!           rh2o (i,k) = max(_zero, qlyr(j1,k))                   ! input mass mixing ratio
!ncep model use
            rh2o (i,k) = max(_zero, qlyr(j1,k)/(_one-qlyr(j1,k))) ! input specific humidity
            qo3  (i,k) = max(_zero, olyr(j1,k))                   ! o3 mass mixing ratio
            rrvco2 (i,k) = max(_zero, gasvmr(j1,k,1))             ! co2 vol mixing ratio
            o2mixrt(i,k) = max(_zero, gasvmr(j1,k,4))             ! o2 vol mixing ratio

            cldfrc (i,k) = clouds(j1,k,1)
            cliqp  (i,k) = clouds(j1,k,2)
            reliq  (i,k) = clouds(j1,k,3)
            cicep  (i,k) = clouds(j1,k,4)
            reice  (i,k) = clouds(j1,k,5)
            crain  (i,k) = clouds(j1,k,6)
            rerain (i,k) = clouds(j1,k,7)
            csnow  (i,k) = clouds(j1,k,8)
            resnow (i,k) = clouds(j1,k,9)
          enddo

          cosz1(i)      = cosz(j1)
          pflux(i,1)    = _zero
          pflux(i,NLP1) = 1.0e2 * plvl(i,NLP1)
        enddo

      else                        ! input data from sfc to toa

        do i = 1, NPTS
          j1 = idxday(i)

          do k = 1, NLAY
            k1 = NLP1 - k
            press(i,k)   = 1.0e2 * plyr(j1,k1)  ! convert press from mb to mks unit
            pflux(i,k+1) = 1.0e2 * plvl(j1,k1)  ! convert press from mb to mks unit
!           pflux(i,k+1) = 1.0e2 * 0.5 * (plyr(j1,k)+plyr(j1,k+1))
            temp (i,k) = tlyr(j1,k1)
!test use
!           rh2o (i,k) = max(_zero, qlyr(j1,k1))                  ! input mass mixing ratio
!ncep model use
            rh2o (i,k) = max(_zero, qlyr(j1,k1)/(_one-qlyr(j1,k)))! input specific humidity
            qo3  (i,k) = max(_zero, olyr(j1,k1))                  ! o3 mass mixing ratio
            rrvco2 (i,k) = max(_zero, gasvmr(j1,k1,1))            ! co2 vol mixing ratio
            o2mixrt(i,k) = max(_zero, gasvmr(j1,k1,4))            ! o2 vol mixing ratio

            cldfrc (i,k) = clouds(j1,k1,1)
            cliqp  (i,k) = clouds(j1,k1,2)
            reliq  (i,k) = clouds(j1,k1,3)
            cicep  (i,k) = clouds(j1,k1,4)
            reice  (i,k) = clouds(j1,k1,5)
            crain  (i,k) = clouds(j1,k1,6)
            rerain (i,k) = clouds(j1,k1,7)
            csnow  (i,k) = clouds(j1,k1,8)
            resnow (i,k) = clouds(j1,k1,9)
          enddo

          cosz1(i)      = cosz(j1)
          pflux(i,1)    = _zero
          pflux(i,NLP1) = 1.0e2 * plvl(i,1)
        enddo

      endif                       ! if_iflip

!  --- ... define deltaz in meters.
 
      do i = 1, NPTS
        deltaz(i,1) = 2.0 * con_rd * temp(i,1) / con_g
      enddo

      do k = 2, NLAY
        do i = 1, NPTS
          deltaz(i,k) = alog( pflux(i,k+1)/pflux(i,k) )                 &
     &                * con_rd * temp(i,k) / con_g
        enddo
      enddo

!  --- ... define a cloud existence variable.      

      do k = 1, NLAY
        do i = 1, NPTS
          lcloud(i,k) = cldfrc(i,k) > _zero
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, NPTS
          cliq1 (i,k) = cliqp(i,k) / deltaz(i,k)   ! convert to g/m**3
          cice1 (i,k) = cicep(i,k) / deltaz(i,k)
          crain1(i,k) = crain(i,k) / deltaz(i,k)
          csnow1(i,k) = csnow(i,k) / deltaz(i,k)
        enddo
      enddo

!  --- ... obtain cloud properties from cloudrad_package

      call cloudpar                                                     &
!  ---  inputs:
     &     ( reliq,reice,rerain,cliq1,cice1,crain1,csnow1,              &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       cldext, cldsct, cldasymm                                   &
     &     )

!  --- ... obtain aerosol properties

      if ( iaersw == 0 ) then

        do nb = 1, NBANDS
          do k = 1, NLAY
            do i = 1, NPTS
              aerextop(i,k,nb) = _zero
              aerssalb(i,k,nb) = _zero
              aerasymm(i,k,nb) = _zero
            enddo
          enddo
        enddo

      else

        if (iflip == 0) then        ! input from toa to sfc

          do nb = 1, NBANDS
            do i = 1, NPTS
              j1 = idxday(i)

              do k = 1, NLAY
                aerextop(i,k,nb) = aerosols(j1,k,nb,1)
                aerssalb(i,k,nb) = aerosols(j1,k,nb,2)
                aerasymm(i,k,nb) = aerosols(j1,k,nb,3)
              enddo
            enddo
          enddo

        else                        ! input data from sfc to toa

          do nb = 1, NBANDS
            do i = 1, NPTS
              j1 = idxday(i)

              do k = 1, NLAY
                k1 = NLP1 - k
                aerextop(i,k,nb) = aerosols(j1,k1,nb,1)
                aerssalb(i,k,nb) = aerosols(j1,k1,nb,2)
                aerasymm(i,k,nb) = aerosols(j1,k1,nb,3)
              enddo
            enddo
          enddo

        endif                       ! end if_iflip_block

      endif  !  end if_iaersw_block

!  --- ... define pressure related quantities, pressure is in mks units. 

      do k = 2, NLP1
        do i = 1, NPTS
          deltap (i,k-1) = pflux(i,k) - pflux(i,k-1)
          delpdig(i,k-1) = deltap(i,k-1)/ con_g
          scale  (i,k)   = pflux(i,k)*pflux(i,k)/ con_p0
        enddo
      enddo

!  --- ... define the scaled and unscaled co2 and o2 pathlengths in cm-
!          atm, and the unscaled h2o and o3 amounts in kgrams/meter**2. 
!          cm-atm needed as units because of c2co2 having those units.
 
      do ng = 1, NSOLWG
        do k = 2, NLP1
          do i = 1, NPTS
            totco2   (i,k,ng) = 50.0 * rrvco2(i,k-1) * scale(i,k)       &
     &                        / (con_g * rhoair * cosz1(i))
            totco2str(i,k,ng) = 100.0 * rrvco2(i,k-1) * pflux(i,k)      &
     &                        / (con_g * rhoair * cosz1(i))
            toto2    (i,k,ng) = 50.0 * o2mixrt(i,k-1) * scale(i,k)      &
     &                        / (con_g * rhoair * cosz1(i))
            toto2str (i,k,ng) = 100.0 * o2mixrt(i,k-1) * pflux(i,k)     &
     &                        / (con_g * rhoair * cosz1(i))
          enddo
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, NPTS
          wh2ostr(i,k) = rh2o(i,k) * delpdig(i,k)
          wo3(i,k)     = qo3 (i,k) * delpdig(i,k)
        enddo
      enddo

!  --- ... define the molecular density for use in calculating the
!          rayleigh optical depth.
 
      do k = 1, NLAY
        do i = 1, NPTS
          densitymol(i,k) = refquanray * press(i,k) / temp(i,k)
        enddo
      enddo

!  --- ... begin band loop, np is a counter for the pseudo-monochromatic
!          frequency point number
 
      np = 0
      do nb = 1, NBANDS
 
!  --- ... define the surface albedo (infrared value for infrared bands,
!          visible value for the remaining bands). 
 
        if ( nb <= NIRBANDS ) then     ! for nir bands
          do i = 1, NPTS
            j1 = idxday(i)
            albdir(i) = sfcalb(j1,1)
            albdif(i) = sfcalb(j1,2)
          enddo
        else                           ! for uv+vis bands
          do i = 1, NPTS
            j1 = idxday(i)
            albdir(i) = sfcalb(j1,3)
            albdif(i) = sfcalb(j1,4)
          enddo
        endif
 
!  --- ... define the local variables for the band values of aerosol and
!          cloud single scattering parameters.
!  note: the unit for the aerosol extinction is kilometer**(-1). 
 
        do k = 1, NLAY
          do i = 1, NPTS
            aeroextopdep(i,k) = aerextop(i,k,nb)
            aerosctopdep(i,k) = aerssalb(i,k,nb)*aerextop(i,k,nb)
            aeroasymfac (i,k) = aerasymm(i,k,nb)
          enddo
        enddo
 
        do k = 1, NLAY
          do i = 1, NPTS
            cloudextopdep(i,k)=cldext(i,k,nb)*deltaz(i,k)*1.0e-3
            cloudsctopdep(i,k)=cldsct(i,k,nb)*deltaz(i,k)*1.0e-3
            cloudasymfac (i,k)=cldasymm(i,k,nb)
          enddo
        enddo

!  --- ... define the rayleigh optical depths.
 
        do k = 1, NLAY
          do i = 1, NPTS
            rayopdep(i,k) = bdens(nb)*densitymol(i,k)*deltaz(i,k)
          enddo
        enddo
 
        if ( nb <= NH2OBANDS ) then
 
!  --- ... define the h2o scaled gas amounts in kgrams/meter**2
 
          do k = 1, NLAY
            do i = 1, NPTS
              wh2o(i,k) = rh2o(i,k) * delpdig(i,k)                      &
     &                  * exp( powph2o(nb)*alog(press(i,k)/p0h2o(nb)) )
            enddo
          enddo

!  --- ... calculate the "effective" co2 and o2 gas optical depths for
!          the appropriate absorbing bands.
!  note: a correction is applied to the determined transmissions for 
!          co2 and o2 if t < 0, which can occur for large zenith angles.
 
          if ( c1co2(nb) /= 1.0e-99 ) then
            do ng = 1, NSOLWG
              do k = 2, NLP1
                do i = 1, NPTS
                  alphaco2(i,k) = -c4co2(nb) + c1co2(nb)                &
     &                    * exp( c3co2(nb)*alog(totco2(i,k,ng)          &
     &                    + c2co2(nb)) )
                  alphaco2str(i,k) = -c4co2str(nb) + c1co2str(nb)       &
     &                    * exp( c3co2str(nb)*alog(totco2str(i,k,ng)    &
     &                    + c2co2str(nb)) )
                  tco2(i,k-1) =                                         &
     &                (_one - alphaco2(i,k))*(_one - alphaco2str(i,k))  &
     &              / ((_one-alphaco2(i,k-1))*(_one-alphaco2str(i,k-1)))

                  if (tco2(i,k-1) .le. _zero) tco2(i,k-1) = 1.0e-9
                  efftauco2(i,k-1,ng) = -cosz1(i)*alog(tco2(i,k-1))
                enddo
              enddo
            enddo
          else
            do ng = 1, NSOLWG
              do k = 1, NLAY
                do i = 1, NPTS
                  efftauco2(i,k,ng) = _zero
                enddo
              enddo
            enddo
	  endif
 
          if ( c1o2(nb) /= 1.0e-99 ) then
            do ng = 1, NSOLWG
              do k = 2, NLP1
                do i = 1, NPTS
                  alphao2(i,k) = -c4o2(nb) + c1o2(nb)                   &
     &                     * exp( c3o2(nb)*alog(toto2(i,k,ng)           &
     &                     + c2o2(nb)) )
                  alphao2str(i,k) = -c4o2str(nb) + c1o2str(nb)          &
     &                     * exp( c3o2str(nb)*alog(toto2str(i,k,ng)     &
     &                     + c2o2str(nb)) )
                  to2(i,k-1) =                                          &
     &                (_one - alphao2(i,k))*(_one - alphao2str(i,k))    &
     &              / ((_one-alphao2(i,k-1))*(_one-alphao2str(i,k-1)))

                  if ( to2(i,k-1) .le. _zero ) to2(i,k-1) = 1.0e-9
                  efftauo2(i,k-1,ng) = -cosz1(i)*alog(to2(i,k-1))
                enddo
              enddo
            enddo
          else
            do ng = 1, NSOLWG
              do k = 1, NLAY
                do i = 1, NPTS
                  efftauo2(i,k,ng) = _zero
                enddo
              enddo
            enddo
          endif
        endif
 
!  --- ... calculate the "effective" o2 gas optical depths for the 
!          Schuman-Runge band.

        if ( nb == NBANDS ) then
          do ng = 1, NSOLWG
            do k = 2, NLP1
              do i = 1, NPTS
                alphao2str(i,k) = -c4o2strschrun + c1o2strschrun        &
     &                   * exp( c3o2strschrun*alog(toto2str(i,k,ng)     &
     &                   + c2o2strschrun) )
                to2(i,k-1) = (_one - alphao2str(i,k))                   &
     &                     / (_one - alphao2str(i,k-1))
                if ( to2(i,k-1) .le. _zero ) to2(i,k-1) = 1.0e-9
                efftauo2(i,k-1,ng) = -cosz1(i)*alog( to2(i,k-1) )
              enddo
            enddo
          enddo
        endif

!  --- ... initialize summing arrays

        do k = 1, NLP1
          do i = 1, NPTS
            sumtr0(i,k) = _zero
            sumre0(i,k) = _zero
            sumtr (i,k) = _zero
            sumre (i,k) = _zero
          enddo
        enddo

!! --- ... optional surface fluxes arrays

        if ( lfdncmp ) then
          do i = 1, NPTS
            sumuvb(i) = _zero
            sumuvd(i) = _zero
            sumirb(i) = _zero
            sumird(i) = _zero
          enddo
        endif

!  --- ... define clear sky arrays

        if (nb >= FIRSTRAYBAND) then
          do k = 1, NLAY
            do i = 1, NPTS
              sctopdep0(i,k) = rayopdep(i,k) + aerosctopdep(i,k)
              gg0 (i,k) = aeroasymfac(i,k) * aerosctopdep(i,k)          &
     &                  / sctopdep0(i,k)
              ff0 (i,k) = aeroasymfac(i,k)**2 * aerosctopdep(i,k)       &
     &                  / sctopdep0(i,k)
              gstr0(i,k) = (gg0(i,k)-ff0(i,k)) / (_one-ff0(i,k))
            enddo
          enddo

        endif

!  --- ... define cloudy sky arrays

        do k = 1, NLAY
          do i = 1, NPTS
            if ( lcloud(i,k) ) then
              sctopdepc(i,k) = rayopdep(i,k) + aerosctopdep(i,k)        &
     &                       + cloudsctopdep(i,k) 
              ggc (i,k) = ( cloudasymfac(i,k)*cloudsctopdep(i,k)        &
     &                  + aeroasymfac(i,k)*aerosctopdep(i,k) )          &
     &                  / sctopdepc(i,k)
              ffc (i,k) = ( cloudasymfac(i,k)**2*cloudsctopdep(i,k)     &
     &                  + aeroasymfac(i,k)**2*aerosctopdep(i,k) )       &
     &                  / sctopdepc(i,k)
              gstrc(i,k) = (ggc(i,k)-ffc(i,k))/(_one-ffc(i,k))
            endif
          enddo
        enddo

!  --- ... begin frequency points in the band loop

        do nf = 1, nfreqpts(nb)
          np = np + 1
 
!  --- ... define the h2o + o3 gas optical depths.
 
          if ( strterm(np) ) then
            do k = 1, NLAY
              do i = 1, NPTS
                opdep(i,k) = kh2o(np)*wh2ostr(i,k) + ko3(np)*wo3(i,k)
              enddo
            enddo				  
          else
            do k = 1, NLAY
              do i = 1, NPTS
                opdep(i,k) = kh2o(np)*wh2o(i,k) + ko3(np)*wo3(i,k)
              enddo
            enddo				  
          end if

!  --- ... begin gaussian angle loop
 
          do ng = 1, NSOLWG

            do k = 1, NLAY
              do i = 1, NPTS
                gasopdep(i,k) = opdep(i,k) + efftauco2(i,k,ng)          &
     &                        + efftauo2(i,k,ng) 
              enddo
            enddo				  

!  --- ... clear sky mode
!  note: in this mode, the delta-eddington method is performed for all 
!          spatial points.

!  --- ... calculate the scaled single-scattering quantities for use in
!          the delta-eddington routine.

            if (nb >= FIRSTRAYBAND)  then
              do k = 1, NLAY
                do i = 1, NPTS
                  extopdep0(i,k) = gasopdep(i,k) + rayopdep(i,k)        &
     &                           + aeroextopdep(i,k)
                  ssalb0   (i,k) = sctopdep0(i,k)/extopdep0(i,k)
                  taustr0  (i,k) = extopdep0(i,k)                       &
     &                             * (_one - ssalb0(i,k)*ff0(i,k))
                  omgstr0  (i,k) = ssalb0(i,k)*(_one - ff0(i,k))        &
     &                             / (_one - ssalb0(i,k)*ff0(i,k))
                enddo
              enddo

!----------------------------------------------------------------------c
! calculate the reflection and transmission in the scattering layers   c
! using the delta-eddington method.                                    c
!----------------------------------------------------------------------c

              call deledd                                               &
!  ---  inputs:
     &    ( taustr0 ,omgstr0 ,gstr0 ,cosz1, ng, lcloud,                 &
     &      NPTS, NLAY, IK, .false.,                                    &
!  ---  outputs:
     &      rlydir0, tlydir0, rlydif0, tlydif0, tlyde0                  &
     &    )

              if (ng /= 1) then
                do k = 1, NLAY
                  do i = 1, NPTS
                    tlydif0(i,k) = _zero
                    rlydif0(i,k) = _zero
                  enddo
                enddo

                do ns = 1, NSTREAMS
                  do k = 1, NLAY
                    do i = 1, NPTS
                      tlydif0(i,k) = tlydif0(i,k)                        &
     &                          + exp( -gasopdep(i,k)/coszstr(ns) )      &
     &                          * wtstr(ns)*coszstr(ns)
                    enddo
                  enddo
                enddo
              endif

!  --- ... initialize the layer reflection and transmission arrays with
!          the non-scattering case.

            else
              tlydif0(:,:) = _zero

              do ns = 1, NSTREAMS
                do k = 1, NLAY
                  do i = 1, NPTS
                    tlydif0(i,k) = tlydif0(i,k)                         &
     &                             + exp( -gasopdep(i,k)/coszstr(ns) )  &
     &                             * wtstr(ns)*coszstr(ns)
                  enddo
                enddo
              enddo

	      do k = 1, NLAY
                do i = 1, NPTS
                  tlydir0(i,k) =  exp( -gasopdep(i,k)/cosz1(i) )
                  tlyde0 (i,k) = tlydir0(i,k)
                  rlydir0(i,k) = _zero
                  rlydif0(i,k) = _zero
                enddo
	      enddo
            endif

!  --- ... overcast sky mode 
!  note: in this mode, the delta-eddington method is performed only for
!          spatial points containing a cloud.

!  --- ... calculate the scaled single-scattering quantities for use in
!          the delta-eddington routine.

            do k = 1, NLAY
              do i = 1, NPTS
                if ( lcloud(i,k) ) then
                  extopdepc(i,k) = gasopdep(i,k) + rayopdep(i,k)        &
     &                      + aeroextopdep(i,k) + cloudextopdep(i,k)
                  ssalbc   (i,k) = sctopdepc(i,k)/extopdepc(i,k)
                  taustrc  (i,k) = extopdepc(i,k)                       &
     &                           * (_one - ssalbc(i,k)*ffc(i,k))
                  omgstrc  (i,k) = ssalbc(i,k)*(_one - ffc(i,k))        &
     &                           / (_one - ssalbc(i,k)*ffc(i,k))
                endif
              enddo
            enddo

!  --- ... calculate the reflection and transmission in the scattering
!          layers using the delta-eddington method.

            call deledd                                                 &
!  ---  inputs:
     &    ( taustrc ,omgstrc ,gstrc ,cosz1, ng, lcloud,                 &
     &      NPTS, NLAY, IK, .true.,                                     &
!  ---  outputs:
     &      rlydirc, tlydirc, rlydifc, tlydifc, tlydec                  &
     &    )

!  --- ... can reduce to only points with clouds

            if (ng /= 1) then
              do k = 1, NLAY
                do i = 1, NPTS
                  tlydifc(i,k) = tlydif0(i,k)
                  rlydifc(i,k) = rlydif0(i,k)
                enddo
              enddo
            endif
 
!  --- ... weight the reflection and transmission arrays for clear and 
!          overcast sky conditions by the cloud fraction, to calculate 
!          the resultant values.
 
            do k = 1, NLAY
              do i = 1, NPTS
                if ( .not. lcloud(i,k) ) then
                  rlydir(i,k) = rlydir0(i,k)
                  tlydir(i,k) = tlydir0(i,k)
                  rlydif(i,k) = rlydif0(i,k)
                  tlydif(i,k) = tlydif0(i,k)
                  tlyde (i,k) = tlyde0 (i,k)
                else
                  rlydir(i,k) = cldfrc(i,k) * rlydirc(i,k)              &
     &                        + (_one - cldfrc(i,k)) * rlydir0(i,k)
                  rlydif(i,k) = cldfrc(i,k) * rlydifc(i,k)              &
     &                        + (_one - cldfrc(i,k)) * rlydif0(i,k)
                  tlydir(i,k) = cldfrc(i,k) * tlydirc(i,k)              &
     &                        + (_one - cldfrc(i,k)) * tlydir0(i,k)
                  tlydif(i,k) = cldfrc(i,k) * tlydifc(i,k)              &
     &                        + (_one - cldfrc(i,k)) * tlydif0(i,k)
                  tlyde (i,k) = cldfrc(i,k) * tlydec (i,k)              &
     &                        + (_one - cldfrc(i,k)) * tlyde0 (i,k)
                endif
              enddo
            enddo
 
!  --- ... calculate the reflection and transmission at flux levels from
!          the direct and diffuse values of reflection and transmission
!          in the corresponding layers using the adding method. 
!  ---  clear-sky
 
            call Adding                                                 &
!  ---  inputs
     &    ( rlydir0,tlydir0,rlydif0,tlydif0,tlyde0,albdir,albdif,       &
     &      NPTS, NLAY, NLP1,                                           &
!  ---  outputs:
     &      reflec0,transm0                                             &
     &    )
 
!  --- total-sky

            call Adding                                                 &
!  ---  inputs
     &    ( rlydir,tlydir,rlydif,tlydif,tlyde,albdir,albdif,            &
     &      NPTS, NLAY, NLP1,                                           &
!  ---  outputs:
     &      reflec,transm                                               &
!! ---  optional outputs:
     &,     strnbm,strndf                                               &
     &    )
 
!  --- ... weight and sum the reflectance and transmittance to calculate
!          the band values.
 
            do i = 1, NPTS
              wtfac(i) = wtfreq(np)*gausswt(ng)*cosz1(i)
            enddo
 
            do k = 1, NLP1
              do i = 1, NPTS
                sumtr0(i,k) = sumtr0(i,k) + transm0(i,k)*wtfac(i)
                sumre0(i,k) = sumre0(i,k) + reflec0(i,k)*wtfac(i)

                sumtr(i,k) = sumtr(i,k) + transm(i,k)*wtfac(i)
                sumre(i,k) = sumre(i,k) + reflec(i,k)*wtfac(i)
              enddo 
            enddo

!! --- ...  optional beam and diff transmttance at surface
            if ( lfdncmp ) then
              if ( nb <= NIRBANDS ) then     ! for nir bands
                do i = 1, NPTS
                  sumirb(i) = sumirb(i) + strnbm(i)*wtfac(i)
                  sumird(i) = sumird(i) + strndf(i)*wtfac(i)
                enddo
              else                           ! for uv+vis bands
                do i = 1, NPTS
                  sumuvb(i) = sumuvb(i) + strnbm(i)*wtfac(i)
                  sumuvd(i) = sumuvd(i) + strndf(i)*wtfac(i)
                enddo
              endif
            endif

          enddo    ! end of gaussian loop (ng-loop)
        enddo  ! end of frequency points (nf-loop) in the band loop
 
!  --- ... normalize the solar flux in the band to the appropriate value
!          for the given total solar insolation. 
 
        tem0 = solcon*solflxband(nb) / solflxtotal
 
!  --- ... sum the band fluxes to calculate the total spectral values.
 
        do k = 1, NLP1
          do i = 1, NPTS
            dfsw0(i,k) = dfsw0(i,k) + sumtr0(i,k) * tem0
            ufsw0(i,k) = ufsw0(i,k) + sumre0(i,k) * tem0

            dfsw(i,k) = dfsw(i,k) + sumtr(i,k) * tem0
            ufsw(i,k) = ufsw(i,k) + sumre(i,k) * tem0
          enddo
        enddo

!! --- ...  optional spectral band net flux
        if ( lhswb ) then
          do k = 1, NLP1
            do i = 1, NPTS
              fswb(i,k,nb) = (sumtr(i,k) - sumre(i,k) ) * tem0
            enddo
          enddo
        endif

        if ( lfdncmp ) then
!! --- ...  optional surface downward fluxes
          if ( nb <= NIRBANDS ) then         ! for nir bands
            do i = 1, NPTS
              sdnirb(i) = sdnirb(i) + sumirb(i) * tem0
              sdnird(i) = sdnird(i) + sumird(i) * tem0
            enddo
          else                               ! for uv+vis bands
            do i = 1, NPTS
              sdnuvb(i) = sdnuvb(i) + sumuvb(i) * tem0
              sdnuvd(i) = sdnuvd(i) + sumuvd(i) * tem0
            enddo
          endif

!! --- ...  optional uv-b surface fluxes
          if ( nb >= NUVBSTR .and. nb <= NUVBEND ) then
            do i = 1, NPTS
              suvbf0(i) = suvbf0(i) + sumtr0(i,NLP1) * tem0
              suvbfc(i) = suvbfc(i) + sumtr (i,NLP1) * tem0
            enddo
          endif
        endif

      enddo      ! end of band loop (nb-loop)
 
!  ---  toa and sfc fluxes

      do i = 1, NPTS
        j1 = idxday(i)
        topflx(j1)%upfxc = ufsw (i,1)
        topflx(j1)%dnfxc = dfsw (i,1)
        topflx(j1)%upfx0 = ufsw0(i,1)

        sfcflx(j1)%upfxc = ufsw (i,NLP1)
        sfcflx(j1)%dnfxc = dfsw (i,NLP1)
        sfcflx(j1)%upfx0 = ufsw0(i,NLP1)
        sfcflx(j1)%dnfx0 = dfsw0(i,NLP1)
      enddo

      if ( lfdncmp ) then
        do i = 1, NPTS
          j1 = idxday(i)
!! ---  optional uv-b surface downward flux
          fdncmp(j1)%uvbf0 = suvbf0(i)
          fdncmp(j1)%uvbfc = suvbfc(i)

!! ---  optional beam and diffuse sfc fluxes
          fdncmp(j1)%nirbm = sdnirb(i)
          fdncmp(j1)%nirdf = sdnird(i)
          fdncmp(j1)%visbm = sdnuvb(i)
          fdncmp(j1)%visdf = sdnuvd(i)
        enddo
      endif

      if (iflip == 0) then        ! input data from toa to sfc

        fsw(:,:) = dfsw(:,:) - ufsw(:,:) 

!  ---  compute heating rates
        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            tem0 = radcon / deltap(i,k)
            hswc(j1,k) = tem0 * (fsw (i,k) - fsw (i,k+1))
          enddo
        enddo

!! ---  optional clear sky heating
        if ( lhsw0 ) then
          fsw0(:,:) = dfsw0(:,:) - ufsw0(:,:) 
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              tem0 = radcon / deltap(i,k)
              hsw0(j1,k) = tem0 * (fsw0(i,k) - fsw0(i,k+1))
            enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhswb ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do nb = 1, NBANDS
            do k = 1, NLAY
              tem0 = radcon / deltap(i,k)
              hswb(j1,k,nb) = tem0 * (fswb(i,k,nb) - fswb(i,k+1,nb))
            enddo
            enddo
          enddo
        endif

!! ---  optional output fluxes
        if ( lflxprf ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              flxprf(j1,k)%upfxc = ufsw (i,k)
              flxprf(j1,k)%dnfxc = dfsw (i,k)
              flxprf(j1,k)%upfx0 = ufsw0(i,k)
              flxprf(j1,k)%dnfx0 = dfsw0(i,k)
            enddo
          enddo
        endif

      else                        ! input data from sfc to toa

        fsw(:,:) = dfsw(:,:) - ufsw(:,:) 

!  ---  compute heating rates
        do i = 1, NPTS
          j1 = idxday(i)
          do k = 1, NLAY
            k1 = NLP1 - k
            tem0 = radcon / deltap(i,k)
            hswc(j1,k1) = tem0 * (fsw (i,k) - fsw (i,k+1))
          enddo
        enddo

!! ---  optional clear sky heating
        if ( lhsw0 ) then
          fsw0(:,:) = dfsw0(:,:) - ufsw0(:,:) 
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLAY
              k1 = NLP1 - k
              tem0 = radcon / deltap(i,k)
              hsw0(j1,k1) = tem0 * (fsw0(i,k) - fsw0(i,k+1))
            enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhswb ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do nb = 1, NBANDS
            do k = 1, NLAY
              k1 = NLP1 - k
              tem0 = radcon / deltap(i,k)
              hswb(j1,k1,nb) = tem0 * (fswb(i,k,nb) - fswb(i,k+1,nb))
            enddo
            enddo
          enddo
        endif

!! ---  optional output fluxes
        if ( lflxprf ) then
          do i = 1, NPTS
            j1 = idxday(i)
            do k = 1, NLP1
              k1 = NLP1 - k + 1
              flxprf(j1,k1)%upfxc = ufsw (i,k)
              flxprf(j1,k1)%dnfxc = dfsw (i,k)
              flxprf(j1,k1)%upfx0 = ufsw0(i,k)
              flxprf(j1,k1)%dnfx0 = dfsw0(i,k)
            enddo
          enddo
        endif

      endif                       ! if_flip

!
      return
!...................................
      end  subroutine swrad
!-----------------------------------



!-----------------------------------
      subroutine rswinit                                                &
!...................................

!  ---  inputs:
     &     ( icwp, me, NLAY )
!  ---  outputs: (none)
 
!  *******************************************************************  !
!                                                                       !
!    subroutine rswinit defines the time-independent quantities         !
!    associated with the incoming shortwave radiation in the multiple-  !
!    band solar radiation parameterization.                             !
!                                                                       !
!  -------------------------------------------------------------------  !
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
!  control flags in module "module_radsw_cntr_para":                    !
!     iswrate - heating rate unit selections                            !
!               =1: output in k/day                                     !
!               =2: output in k/second                                  !
!     iaersw  - flags for aerosols effect (not yet!)                    !
!               =0: without aerosol effect                              !
!               >0: include aerosol effect                              !
!                                                                       !
!  -------------------------------------------------------------------  !
!                                                                       !
!  module variables used:                                               !
!                                                                       !
!   solflxband   = the solar flux in each parameterization band         !
!   endwvnbands  = the wavenumber value for the band limits             !
!   nwvnsolar    = the number of wavenumbers in each region where the   !
!                  solar flux is constant                               !
!   solint       = the solar flux in watts per meter**2 in each         !
!                  wavenumber region where it is constant               !
!   wtfreq       = the weight associated with each exponential term     !
!                                                                       !
! define parameters for use in determining the absorption due to co2,   !
! h2o, o2 and o3 in each band.                                          !
!                                                                       !
!   nfreqpts     = the number of pseudo-monochromatic frequencies       !
!   firstrayband = the first band number where the contribution by      !
!                  rayleigh scattering is included in the solar         !
!                  calculations                                         !
!   nirbands     = the number of bands in the near-infrared (used in    !
!                  assigning the value of the surface albedo for the    !
!                  near-infrared, and the visible and ultraviolet       !
!                  regions, separately)                                 !
!   powph2o      = the scaling factor used in the fit of the h2o        !
!                  transmission function                                !
!   p0h2o        = the reference pressure (mb) used in the fit of the   !
!                  h2o transmission function                            !
!   c(n)co2(str) = coefficients for the absorptivity expression for co2 !
!                  for the pressure-scaled and non-scaled, respectively,!
!                  portions of the fit (=1.0e-99 if no absorption)      !
!   c(n)o2(str)  = coefficients for the absorptivity expression for o2  !
!                  for the pressure-scaled and non-scaled, respectively,!
!                  portions of the fit (=1.0e-99 if no absorption)      !
!   ""(schrun)   = coefficients for the absorptivity expression for the !
!                  Schuman-Runge o2 band (non-scaled only)              !
!   kh2o         =  the psuedo-absorption coefficients in cm2/gm for h2o!
!   ko3          = the absorption coefficients in cm2/gm for o3         !
!   strterm      = logical flag to indicate whether or not a h2o pseudo-!
!                  absorption coefficient is assigned a non-scaled      !
!                  (true) or pressure-scaled (false) gas amount         !
!   ptstr        = gaussian points and weights for evaluation of the    !
!                  diffuse beam.                                        !
!                                                                       !
! define the solar weights and interval counters that are used to       !
! determine the single-scattering properties for the parameterization   !
! band spectral intervals, from the specified spectral intervals for    !
! cloud water drops, ice particles, rain, and snow.                     !
!                                                                       !
!   nv1liq   = interval number for liquid droplets single-scattering    !
!              properties corresponding to the first psuedo-monochrom-  !
!              atic frequency in a given parameterization band.         !
!   nv2liq   = interval number for liquid droplet single-scattering     !
!              properties corresponding to the last psuedo-monochromatic!
!              frequency in a given parameterization band.              !
!   nv1ice   = interval number for ice particle single-scattering       !
!              properties corresponding to the first psuedo-monochrom-  !
!              atic frequency in a given parameterization band.         !
!   nv2ice   = interval number for ice particle single-scattering       !
!              properties corresponding to the last psuedo-monochromatic!
!              frequency in a given parameterization band.              !
!   nv1rain  = interval number for rain drop single-scattering          !
!              properties corresponding to the first psuedo-monochrom-  !
!              atic frequency in a given parameterization band.         !
!   nv2rain  = interval number for rain drop single-scattering          !
!              properties corresponding to the last psuedo-monochromatic!
!              frequency in a given parameterization band.              !
!   nv1snow  = interval number for snow particle single-scattering      !
!              properties corresponding to the first psuedo-monochrom-  !
!              atic frequency in a given parameterization band.         !
!   nv2isnow = interval number for snow particle single-scattering      !
!              properties corresponding to the last psuedo-monochromatic!
!              frequency in a given parameterization band.              !
!   solvliq  = solar flux in the scattering spectrum for liquid droplet !
!   solvice  = solar flux in the scattering spectrum for ice droplet    !
!   solvrain = solar flux in the scattering spectrum for rain drop      !
!   solvsnow = solar flux in the scattering spectrum for snow particle  !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: (none)

!  ---  locals:
      real (kind=kind_phys) :: freqnu(NBANDS), corrfac, pinteg,         &
     &       twopiesq, densmolrefsqt3, wavelength, freqsq, ri, sumsol

      integer :: nb, ni, nw, nw1, nw2
!
!===> ...  begin here
!

      if (me == 0) then
        print *,' - Using GFDL Shortwave Radiation, Version: ',VTAGSW

        if (iaersw == 0) then
          print *,'   --- Aerosol effect is NOT included in SW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        else
          print *,'   --- Using input aerosol parameters for SW'
        endif

        print *,'   --- Number of Gaussian quadrature weights=',NSOLWG
      endif

!  --- ...  setup constant factor radcon for heating rate

      if (iswrate .eq. 1) then
        radcon = con_g * 86400.0 / con_cp           !    (in k/day)
!       radcon = 844.10328                          !    (in k/day)
      else
        radcon = con_g / con_cp                     !    (in k/second)
!       radcon = 844.10328 / 86400.0                !    (in k/second)
      endif

!  --- ...  define the one wavenumber solar fluxes
 
      do nb = 1, NINTSOLAR
	if ( nb .eq. 1 ) then
	  nw1 = 1
        else
          nw1 = nw1 + nwvnsolar(nb-1)
        end if
	nw2 = nw1 + nwvnsolar(nb) - 1
 
	do nw = nw1, nw2
	  solarfluxtoa(nw) = solint(nb)
        end do
      end do

!  --- ...  set up gaussian weight values

      if (NSOLWG == 1) then
        gausswt = _one
      elseif (NSOLWG == 2) then
        do ni = 1, NSOLWG
          gausswt(ni) = gwt2(ni)
        enddo
      elseif (NSOLWG == 4) then
        do ni = 1, NSOLWG
          gausswt(ni) = gwt4(ni)
        enddo
      elseif (NSOLWG == 8) then
        do ni = 1, NSOLWG
          gausswt(ni) = gwt8(ni)
        enddo
      else
        print *,'  ***  ERROR in parameter NSOLWG setting !!! ***'
        print *,'  ***  valid values are 1, 2, 4, or 8 only  ***'
        stop
      endif

!  --- ...  set up co2 and o2 coefficients
 
      do ni = 1, NH2OBANDS
        c4co2   (ni) = c1co2   (ni) * c2co2   (ni)**c3co2   (ni)
        c4co2str(ni) = c1co2str(ni) * c2co2str(ni)**c3co2str(ni)
        c4o2    (ni) = c1o2    (ni) * c2o2    (ni)**c3o2    (ni)
        c4o2str (ni) = c1o2str (ni) * c2o2str (ni)**c3o2str (ni)
      enddo

      c4o2strschrun = c1o2strschrun * c2o2strschrun**c3o2strschrun
 
      solflxtotal = _zero
      do nb = 1, NBANDS
	solflxtotal = solflxtotal + solflxband(nb)
      enddo
 
!  --- ...  define the wavenumbers to evaluate rayleigh optical depth
 
      do nb = 1, NBANDS
	freqnu(nb) = 0.5 * ( endwvnbands(nb-1) + endwvnbands(nb) )
      enddo
 
! -------------------------------------------------------------------- !
! define quantities used to determine the rayleigh optical depth.      !
!                                                                      !
! notes: refquanray is the quantity which multiplies pressure /        !
!        temperature to yield the molecular density.                   !
!                                                                      !
!        bdens is the quantity which multiples the                     !
!        molecular density to yield the rayleigh scattering            !
!        coefficient.                                                  !
!                                                                      !
!        depfac = 1.39E-02 is the depolorization factor.               !
! -------------------------------------------------------------------- !
 
      corrfac = (6.0 + 3.0*depfac )/(6.0 - 7.0* depfac )

! -------------------------------------------------------------------- !
!     gamma   = depfac / (2.0 - depfac)
!     f1 = 0.75 / (1.0 + 2.0*gamma)
!     f2 = 1.0 + 3.0*gamma 
!     f3 = 1.0 - gamma
!     pinteg = 2.0*con_pi * ( 2.0*f1*f2 * ( 1.0 + f3/f2/ 3.0 ) )
! -------------------------------------------------------------------- !

      pinteg         = 4.0 * con_pi
      twopiesq       = 2.0 * con_pi**2 
      densmolrefsqt3 = 3.0 * densmolref**2
 
      do nb = 1, NBANDS
        wavelength = 1.0e+4 / freqnu(nb)
        freqsq = _one / wavelength**2
        ri = _one + 1.0e-8*(6.4328e+3 + 2.94981e+6 / (146. - freqsq)    &
     &                                + 2.55400e+4 / (41.0 - freqsq))
        bdens(nb) = twopiesq * pinteg * convfac * corrfac               &
     &          * (ri**2 - _one)**2 / (densmolrefsqt3*wavelength**4)
      enddo
 
!  --- ...  define the gaussian angles for evaluation of the diffuse beam
 
      do ni = 1, NSTREAMS
        coszstr(ni) = (ptstr(ni) + _one) * 0.5
      enddo

!  --- ...  define the solar weights and interval counters of band
!           spectral intervals for liquid cloud

      ni = 1
      nb = 1
      sumsol = _zero
      solvliq(:,:) = _zero
      nv1liq(1) = 1

      do nw = 1, endwvnbands(NBANDS)
        sumsol = sumsol + solarfluxtoa(nw)

        if ( nw == endliqwvn(ni) ) then
          solvliq(nb,ni) = sumsol
          sumsol = _zero
        endif

        if ( nw == endwvnbands(nb) ) then
          if ( nw /= endliqwvn(ni) ) then
            solvliq(nb,ni) = sumsol
            sumsol = _zero
          endif

          nv2liq(nb) = ni
          nb = nb + 1

          if ( nb <= NBANDS ) then
            if ( nw == endliqwvn(ni) ) then
              nv1liq(nb) = ni + 1
            else
              nv1liq(nb) = ni
            endif
          endif
        endif

        if ( nw .eq. endliqwvn(ni) ) ni = ni + 1
      enddo

!  --- ...  define the solar weights and interval counters of band
!           spectral intervals for ice cloud

      ni = 1
      nb = 1
      sumsol = _zero
      solvice(:,:) = _zero
      nv1ice(1) = 1

      do nw = 1, endwvnbands(NBANDS)
        sumsol = sumsol + solarfluxtoa(nw)

        if ( nw == endicewvn(ni) ) then
          solvice(nb,ni) = sumsol
          sumsol = _zero
        endif

        if ( nw == endwvnbands(nb) ) then
          if ( nw /= endicewvn(ni) ) then
            solvice(nb,ni) = sumsol
            sumsol = _zero
          endif

          nv2ice(nb) = ni
          nb = nb + 1

          if ( nb <= NBANDS ) then
            if ( nw == endicewvn(ni) ) then
              nv1ice(nb) = ni + 1
            else
              nv1ice(nb) = ni
            endif
          endif
        endif

        if ( nw == endicewvn(ni) ) ni = ni + 1
      enddo

!  --- ...  define the solar weights and interval counters of band
!           spectral intervals for rain drop

      ni = 1
      nb = 1
      sumsol = _zero
      solvrain(:,:) = _zero
      nv1rain(1) = 1

      do nw = 1, endwvnbands(NBANDS)
        sumsol = sumsol + solarfluxtoa(nw)

        if ( nw == endrainwvn(ni) ) then
          solvrain(nb,ni) = sumsol
          sumsol = _zero
        endif

        if ( nw == endwvnbands(nb) ) then
          if ( nw /= endrainwvn(ni) ) then
            solvrain(nb,ni) = sumsol
            sumsol = _zero
          endif

          nv2rain(nb) = ni
          nb = nb + 1

          if ( nb <= NBANDS ) then
            if ( nw == endrainwvn(ni) ) then
              nv1rain(nb) = ni + 1
            else
              nv1rain(nb) = ni
            endif
          endif
        endif

        if ( nw == endrainwvn(ni) ) ni = ni + 1
      enddo

!  --- ...  define the solar weights and interval counters of band
!           spectral intervals for snow particle

      ni = 1
      nb = 1
      sumsol = _zero
      solvsnow(:,:) = _zero
      nv1snow(1) = 1

      do nw = 1, endwvnbands(NBANDS)
        sumsol = sumsol + solarfluxtoa(nw)

        if ( nw == endsnowwvn(ni) ) then
          solvsnow(nb,ni) = sumsol
          sumsol = _zero
        endif

        if ( nw == endwvnbands(nb) ) then
          if ( nw /= endsnowwvn(ni) ) then
            solvsnow(nb,ni) = sumsol
            sumsol = _zero
          endif

          nv2snow(nb) = ni
          nb = nb + 1

          if ( nb <= NBANDS ) then
            if ( nw == endsnowwvn(ni) ) then
              nv1snow(nb) = ni + 1
            else
              nv1snow(nb) = ni
            endif
          endif
        endif

        if ( nw == endsnowwvn(ni) ) ni = ni + 1
      enddo

!
      return
!...................................
      end subroutine rswinit
!-----------------------------------



!-----------------------------------
      subroutine adding                                                 &
!...................................
!  ---  inputs
     &    ( rlydir,tlydir,rlydif,tlydif,tlyde,albdir,albdif,            &
     &      NPTS, NLAY, NLP1,                                           &
!  ---  outputs:
     &      reflec,transm                                               &
!  ---  optional outputs:
     &,     strnbm,strndf                                               &
     &    )
 
!  *******************************************************************  !
!                                                                       !
! calculate the reflection and transmission at flux levels from the     !
! direct and diffuse values of reflection and transmission in the       !
! corresponding layers using the adding method.                         !
!                                                                       !
! references:                                                           !
!                                                                       !
! bowen, m.m., and v. ramaswamy, effects of changes in radiatively      !
!        active species upon the lower stratospheric temperatures.,     !
!        j. geophys. res., 18909-18921, 1994.                           !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! rlydir    = the layer reflectivity to a direct incident beam          !
! tlydir    = the layer transmissivity to a direct incident beam        !
! rlydif    = the layer reflectivity to a diffuse incident beam         !
! tlydif    = the layer transmissivity to a diffuse incident beam       !
! tlyde     = the layer transmissivity (non-scattered) to the direct    !
!             incident beam                                             !
! albdir    = surface albedo for direct beam                            !
! albdif    = surface albedo for diffused radiation                     !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! reflec        = the reflectance of the scattered radiation at a       !
!                 level                                                 !
! transm        = the transmittance of the scattered radiation at a     !
!                 level                                                 !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1
 
      real (kind=kind_phys), dimension(:,:), intent(in) :: tlyde,       &
     &       rlydir, rlydif, tlydir, tlydif

      real (kind=kind_phys), dimension(:),  intent(in) :: albdir, albdif
 
!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),intent(out):: reflec, transm

!  ---  optional outputs:
      real (kind=kind_phys), dimension(:), optional, intent(out) ::     &
     &       strnbm, strndf

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: dm1tl, dm2tl,      &
     &       rdm2tl, rdndif, tdndir
      real (kind=kind_phys), dimension(NPTS,NLP1) :: rupdif, rupdir,    &
     &       tlevel, alpp, dm3, dm3r, dm3r1p
 
      integer ::  i, k
!
!===> ...  begin here
!
      do i = 1, NPTS

!  --- ...  initialization for the surface layer
 
        rupdif(i,NLP1) = albdif(i)
        rupdir(i,NLP1) = albdir(i)
 
!  --- ...  initialization for the first scattering layer
 
        tdndir(i,1) = tlydir(i,1)
        rdndif(i,1) = rlydif(i,1)
 
        tlevel(i,1) = _one
        tlevel(i,2) = tlyde(i,1)
      enddo
 
!  --- ...  define the direct transmittance

      do k = 2, NLAY
        do i = 1, NPTS
          tlevel(i,k+1) = tlevel(i,k) * tlyde(i,k) 
        enddo
      enddo
 
!  --- ...  add the inhomogeneous layers downward from the second layer
!           to the surface.  radiation incident from below for diffuse
!           beam, transmission of direct beam and conversion to diffuse.
 
      do k = 2, NLAY
        do i = 1, NPTS
          dm1tl(i,k) = tlydif(i,k) / (_one - rlydif(i,k)*rdndif(i,k-1))
          rdndif(i,k) = rlydif(i,k)                                     &
     &                + rdndif(i,k-1)*tlydif(i,k)*dm1tl(i,k)
          tdndir(i,k) = (tdndir(i,k-1) - tlevel(i,k)) *dm1tl(i,k)       &
     &                + tlevel(i,k) * ( tlydir(i,k) + rlydir(i,k)       &
     &                * rdndif(i,k-1)*dm1tl(i,k) )
        enddo
      enddo
 
!  --- ...  add the inhomogeneous layers upward from the surface to the
!           top of the atmosphere.  radiation incident from above for 
!           diffuse beam, reflection of direct beam and conversion to diffuse.
 
      do k = NLAY, 1, -1
        do i = 1, NPTS
          dm2tl(i,k) = tlydif(i,k) / (_one - rlydif(i,k)*rupdif(i,k+1))
          rdm2tl(i,k) = dm2tl(i,k) * rupdif(i,k+1)
          rupdif(i,k) = rlydif(i,k) + tlydif(i,k)* rdm2tl(i,k) 
        enddo
      enddo
 
      do k = NLAY, 1, -1
        do i = 1, NPTS
          rupdir(i,k) = rlydir(i,k)+tlyde(i,k)*rupdir(i,k+1)*dm2tl(i,k) &
     &                + (tlydir(i,k)-tlyde(i,k))*rdm2tl(i,k)
        enddo
      enddo
 
!  --- ...  add downward to calculate the resultant reflectances and 
!           transmittances at flux levels.
 
      do k = 2, NLP1
        do i = 1, NPTS
          dm3   (i,k) = _one / (_one - rupdif(i,k)*rdndif(i,k-1))
          dm3r  (i,k) = dm3(i,k)*rdndif(i,k-1)
          dm3r1p(i,k) = _one + rupdif(i,k)*dm3r(i,k)
          alpp  (i,k) = (tdndir(i,k-1) - tlevel(i,k) )*dm3(i,k)
          transm(i,k) = alpp(i,k) + tlevel(i,k)                         &
     &                * (_one + rupdir(i,k)*dm3r(i,k))
          reflec(i,k) = alpp(i,k)*rupdif(i,k)                           &
     &                + tlevel(i,k)*rupdir(i,k)*dm3r1p(i,k)
        enddo
      enddo
 

      do i = 1, NPTS
        reflec(i,1) = rupdir(i,1)
        transm(i,1) = _one
      enddo

!! --- ...  optional surface downward beam and diff fluxes

      if ( present( strnbm ) ) then
        do i = 1, NPTS
          strnbm(i) = tlevel(i,NLP1)
        enddo
      endif

      if ( present( strndf ) ) then
        do i = 1, NPTS
          strndf(i) = transm(i,NLP1) - tlevel(i,NLP1)
        enddo
      endif

!
      return
!...................................
      end subroutine adding 
!-----------------------------------



!-----------------------------------
      subroutine deledd                                                 &
!...................................

!  ---  inputs:
     &    ( taustr,omgstr,gstr,cosang,ng,lcloud                         &
     &,     NPTS, NLAY, IK, sw_with_clouds                              &
!  ---  outputs:
     &,     rlydir,tlydir,rlydif,tlydif,tlyde                           &
     &    )
 
!  *******************************************************************  !
!                                                                       !
! calculate the reflection and transmission in the scattering layers    !
! using the delta-eddington method.                                     !
!                                                                       !
! references:                                                           !
!                                                                       !
! joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington      !
!      approximation for radiative flux transfer.,j. atmos. sci.,33,    !
!      2452-2459, 1976.                                                 !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! taustr   = the scaled extinction optical depth                        !
! omgstr   = the scaled single-scattering albedo                        !
! gstr     = the scaled asymmetry factor                                !
! cosang   = the cosine of the solar zenith angle                       !
! ng       = the number of gaussian angles to compute the diurnally     !
!            averaged solar radiation                                   !
! lcloud   = flag for existence of a cloud (used only in 'ovc' mode)    !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! rlydir    = the layer reflectivity to a direct incident beam          !
! tlydir    = the layer transmissivity to a direct incident beam        !
! rlydif    = the layer reflectivity to a diffuse incident beam         !
! tlydif    = the layer transmissivity to a diffuse incident beam       !
! tlyde     = the layer transmissivity (non-scattered) to the direct    !
!             incident beam                                             !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, IK, ng

      real (kind=kind_phys), dimension(:,:), intent(in) :: taustr,      &
     &       omgstr, gstr
      real (kind=kind_phys), dimension(:),   intent(in) :: cosang

      logical, intent(in) :: sw_with_clouds, lcloud(:,:)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: tlyde,      &
     &       rlydir, rlydif, tlydir, tlydif
 
!  ---  local constants
      real (kind=kind_phys), parameter :: onedi3=1.0/3.0
      real (kind=kind_phys), parameter :: twodi3=2.0/3.0

!  ---  locals:
      real (kind=kind_phys), dimension(IK,15)       :: qq
      real (kind=kind_phys), dimension(IK,NSTREAMS) :: rr, tt
      real (kind=kind_phys), dimension(IK)          :: gstr2, taustr2,  &
     &       omgstr2, cosang2, rlydir2, tlyde2, tlydir2, sumr, sumt

      integer :: i, k, ns, nn, ntot, idxcld(NPTS,NLAY)

!
!===> ...  begin here
!

!  --- ...  overcast sky mode
!     note: in this mode, the delta-eddington method is performed only
!           for spatial points containing a cloud.

      if ( sw_with_clouds ) then
        nn = 0

        do k = 1, NLAY
          do i = 1, NPTS
            if ( lcloud(i,k) ) then
              nn = nn + 1
              idxcld(i,k) = nn
              gstr2(nn) = gstr(i,k)
              taustr2(nn) = taustr(i,k)
              omgstr2(nn) = omgstr(i,k)
              cosang2(nn) = cosang(i)
            endif
          enddo
        enddo
        ntot = nn

!  --- ...  clear sky mode
!     note: in this mode, the delta-eddington method is performed for
!           all spatial points.

      else
        nn = 0

        do k = 1, NLAY
          do i = 1, NPTS
            nn = nn + 1
            gstr2(nn) = gstr(i,k)
            taustr2(nn) = taustr(i,k)
            omgstr2(nn) = omgstr(i,k)
            cosang2(nn) = cosang(i)
          enddo
        enddo

        ntot = nn 
      endif

!     note: the following are done to avoid the conservative scattering
!           case, and to eliminate floating point errors in the
!           exponential calculations, respectively.

      do nn = 1, ntot      
        if ( omgstr2(nn) >= _one )   omgstr2(nn) = 9.9999999e-01
        if ( taustr2(nn) >= 1.0e+2 ) taustr2(nn) = 1.0e+2
      enddo
 
!  --- ...  direct quantities
 
      do nn = 1, ntot      
        qq(nn,1)  = 3.0 * ( _one - omgstr2(nn) )
        qq(nn,2)  = _one - omgstr2(nn) * gstr2(nn)
        qq(nn,3)  = sqrt( qq(nn,1)*qq(nn,2) )

        qq(nn,4)  = _one + twodi3 * sqrt( qq(nn,1)/qq(nn,2) ) 
        qq(nn,5)  = _one - twodi3 * sqrt( qq(nn,1)/qq(nn,2) )
        qq(nn,6)  = exp( -taustr2(nn) * qq(nn,3) )

        qq(nn,7)  = qq(nn, 4)/qq(nn,6 )                                 &
     &            - qq(nn, 5)*qq(nn, 5)*qq(nn, 6)/qq(nn, 4)
        qq(nn,8)  = 0.75 * omgstr2(nn)                                  &
     &            / (_one - (qq(nn,3)*cosang2(nn))**2)

        qq(nn,9)  = qq(nn,8)*cosang2(nn)                                &
     &            * (_one + gstr2(nn)*qq(nn,1)*onedi3)
        qq(nn,10) = qq(nn,8)*(_one + gstr2(nn)*qq(nn,1)*cosang2(nn)**2)

        qq(nn,11) = qq(nn,9) - twodi3*qq(nn,10)     
        qq(nn,12) = qq(nn,9) + twodi3*qq(nn,10)     
        qq(nn,13) = exp( -taustr2(nn)/cosang2(nn) )

        qq(nn,14) = (qq(nn,11)*qq(nn,13)                                &
     &            - qq(nn,12)*qq(nn,6)*qq(nn,5)/qq(nn,4)) / qq(nn,7)

        qq(nn,15) = (qq(nn,12) - qq(nn,5)*qq(nn,14) ) / qq(nn,4)

        rlydir2(nn) = qq(nn,5)*qq(nn,15) + qq(nn,4)*qq(nn,14)           &
     &              - qq(nn,11)
        tlydir2(nn) = qq(nn,13) + (qq(nn,6)*qq(nn,4)*qq(nn,15)          &
     &              + qq(nn,5)*qq(nn,14)/qq(nn,6)                       &
     &              - qq(nn,12)*qq(nn,13) )
        tlyde2(nn) = qq(nn,13)
      enddo

!  --- ...  diffuse quantities
!     notes: the number of streams for the diffuse beam is fixed at 4.
!            this calculation is done only for NSOLWG=1.
 
      if ( ng == 1 ) then   
        do ns = 1, NSTREAMS
          do nn = 1, ntot
            qq(nn,8)  = 0.75 * omgstr2(nn)                              &
     &                / ( _one - (qq(nn,3)*coszstr(ns))**2 )

            qq(nn,9)  = qq(nn,8) * coszstr(ns)                          &
     &                * ( _one + gstr2(nn)*qq(nn,1)*onedi3 )
            qq(nn,10) = qq(nn,8)                                        &
     &                * ( _one + gstr2(nn)*qq(nn,1)*coszstr(ns)**2 )

            qq(nn,11) = qq(nn,9) - twodi3*qq(nn,10)
            qq(nn,12) = qq(nn,9) + twodi3*qq(nn,10)
            qq(nn,13) = exp( -taustr2(nn)/coszstr(ns) )

            qq(nn,14) = ( qq(nn,11)*qq(nn,13) - qq(nn,12)*qq(nn,6)      &
     &                * qq(nn,5)/qq(nn,4) ) / qq(nn,7)

            qq(nn,15) = ( qq(nn,12) - qq(nn,5)*qq(nn,14) ) / qq(nn,4)
            rr(nn,ns) = qq(nn,5)*qq(nn,15) + qq(nn,4)*qq(nn,14)         &
     &                - qq(nn,11)
            tt(nn,ns) = qq(nn,13) + qq(nn,6)*qq(nn,4)*qq(nn,15)         &
     &                + qq(nn,5)*qq(nn,14)/qq(nn,6)-qq(nn,12)*qq(nn,13)
          enddo
        enddo

        do nn = 1, IK
          sumr(nn) = _zero
          sumt(nn) = _zero
        enddo
 
        do ns = 1, NSTREAMS
          do nn = 1, ntot
            sumr(nn) = sumr(nn) + rr(nn,ns) * wtstr(ns) * coszstr(ns)
            sumt(nn) = sumt(nn) + tt(nn,ns) * wtstr(ns) * coszstr(ns)
          enddo
        enddo
      endif

!  --- ...  return results in proper locations in (i,k) arrays

      if ( sw_with_clouds ) then
        nn = 0

        do k = 1, NLAY
          do i = 1, NPTS
            if ( lcloud(i,k) ) then
              nn = nn + 1

              rlydir(i,k) = rlydir2(idxcld(i,k))
              tlydir(i,k) = tlydir2(idxcld(i,k))
              tlyde (i,k) = tlyde2 (idxcld(i,k))

              if ( ng == 1 ) then
                rlydif(i,k) = sumr(idxcld(i,k))
                tlydif(i,k) = sumt(idxcld(i,k))
              endif
            endif
          enddo
        enddo

      else
        nn = 0

        do k = 1, NLAY
          do i = 1, NPTS
            nn = nn + 1

            rlydir(i,k) = rlydir2(nn)
            tlydir(i,k) = tlydir2(nn)
            tlyde (i,k) = tlyde2 (nn)

            if ( ng == 1 ) then
              rlydif(i,k) = sumr(nn)
              tlydif(i,k) = sumt(nn)
            endif
          enddo
        enddo

      endif

      if  ( sw_with_clouds .and. nn /= ntot) then
        print *,' Error in deledd',                                     &
     &          'final number of deledd points differs from original'
        stop
      endif

!
      return
!...................................
      end subroutine deledd
!-----------------------------------



!-----------------------------------
      subroutine cloudpar                                               &
!...................................

!  ---  inputs:
     &    ( reliq,reice,rerain,cliq,cice,crain,csnow                    &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     cldext,cldsct,cldasy                                        &
     &    )

!  *******************************************************************  !
!                                                                       !
!  determine the parameterization band values of the single scattering  !
!  parameters (extinction coefficient, scattering coefficient and       !
!  asymmetry factor) for clouds from the size and/or concentration of   !
!  each constituent (cloud drops, rain drops, ice crystals and snow)    !
!  present.                                                             !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! reliq     = the cloud drop effective radius in microns                !
! reice     = the ice crystal effective radius in microns               !
! rerain    = the rain drop effective radius in microns                 !
! cliq      = the cloud drop liquid water concentration in grams / m**3 !
! cice     = the cloud ice water concentration in grams / m**3          !
! crain     = the rain drop water concentration in grams / m**3         !
! csnow     = the snow concentration in grams / m**3                    !
!                                                                       !
!---------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! cldext    = the parameterization band values of the cloud             !
!             extinction coefficient in kilometer**(-1)                 !
! cldsct    = the parameterization band values of the cloud             !
!             scattering coefficient in kilometer**(-1)                 !
! cldasy    = the parameterization band values of the asymmetry         !
!             factor                                                    !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: reliq,       &
     &       reice, rerain, cliq, cice, crain, csnow

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: cldasy,   &
     &       cldext, cldsct

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY,NLIQCLDV)  ::          &
     &       cextvliq, cssavliq, casyvliq
      real (kind=kind_phys), dimension(NPTS,NLAY,NICECLDV)  ::          &
     &       cextvice, cssavice, casyvice
      real (kind=kind_phys), dimension(NPTS,NLAY,NRAINCLDV) ::          &
     &       cextvrain, cssavrain, casyvrain
      real (kind=kind_phys), dimension(NPTS,NLAY,NSNOWCLDV) ::          &
     &       cextvsnow, cssavsnow, casyvsnow

      real (kind=kind_phys), dimension(NPTS,NLAY,NBANDS)    ::          &
     &       cextbliq, cssabliq, casybliq, cextbice, cssabice, casybice,&
     &       cextbrain,cssabrain,casybrain,cextbsnow,cssabsnow,casybsnow

      real (kind=kind_phys), dimension(NPTS,NLAY)        ::             &
     &       cscatliq, cscatice, cscatrain, cscatsnow

      integer :: nb, i, k
!
!===> ...  begin here
!

!  --- ...  define the single scattering parameters for cloud drops

      call slingo                                                       &
!  ---  inputs:
     &    ( cliq, reliq, NPTS, NLAY                                     &
!  ---  outputs:
     &,     cextvliq, cssavliq, casyvliq                                &
     &    )

!  --- ...  define the single scattering parameters for ice crystals

      call fu                                                           &
!  ---  inputs:
     &    ( cice, reice, NPTS, NLAY                                     &
!  ---  outputs:
     &,     cextvice, cssavice, casyvice                                &
     &    )

!  --- ...  define the single scattering parameters for rain drops

      call savijarvi                                                    &
!  ---  inputs:
     &    ( crain, rerain, NPTS, NLAY                                   &
!  ---  outputs:
     &,     cextvrain, cssavrain, casyvrain                             &
     &    )

!  --- ...  define the single scattering parameters for snow

      call snowsw                                                       &
!  ---  inputs:
     &    ( csnow, NPTS, NLAY                                           &
!  ---  outputs:
     &,     cextvsnow, cssavsnow, casyvsnow                             &
     &    )

!  --- ...  use the thick-averaging technique to define the single-
!           scattering properties of the parameterization band spectral
!           intervals from the specified spectral intervals for cloud
!           drops, ice crystals, rain drops, and snow, respectively.

      call thickavg                                                     &
!  ---  inputs:
     &    ( nv1liq,nv2liq,solvliq                                       &
     &,     NLIQCLDV,cextvliq,cssavliq,casyvliq                         &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     cextbliq, cssabliq, casybliq                                &
     &    )

      call thickavg                                                     &
!  ---  inputs:
     &    ( nv1ice,nv2ice,solvice                                       &
     &,     NICECLDV,cextvice,cssavice,casyvice                         &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     cextbice, cssabice, casybice                                &
     &    )

      call thickavg                                                     &
!  ---  inputs:
     &    ( nv1rain,nv2rain,solvrain                                    &
     &,     NRAINCLDV,cextvrain,cssavrain,casyvrain                     &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     cextbrain, cssabrain, casybrain                             &
     &    )

      call thickavg                                                     &
!  ---  inputs:
     &    ( nv1snow,nv2snow,solvsnow                                    &
     &,     NSNOWCLDV,cextvsnow,cssavsnow,casyvsnow                     &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     cextbsnow, cssabsnow, casybsnow                             &
     &    )

!  --- ...  combine the single-scattering properties for all the constituents
!           to define the corresponding values for the clouds.

      do nb = 1,NBANDS
        cscatliq (:,:) = cssabliq (:,:,nb) * cextbliq (:,:,nb)
        cscatice (:,:) = cssabice (:,:,nb) * cextbice (:,:,nb)
        cscatrain(:,:) = cssabrain(:,:,nb) * cextbrain(:,:,nb)
        cscatsnow(:,:) = cssabsnow(:,:,nb) * cextbsnow(:,:,nb)

        cldext(:,:,nb) = cextbliq(:,:,nb) + cextbrain(:,:,nb)           &
     &                 + cextbice(:,:,nb) + cextbsnow(:,:,nb)
        cldsct(:,:,nb) = cscatliq(:,:)    + cscatrain(:,:)              &
     &                 + cscatice(:,:)    + cscatsnow(:,:)
        cldasy(:,:,nb) = ( casybliq (:,:,nb) * cscatliq (:,:)           &
     &                 +   casybice (:,:,nb) * cscatice (:,:)           &
     &                 +   casybrain(:,:,nb) * cscatrain(:,:)           &
     &                 +   casybsnow(:,:,nb) * cscatsnow(:,:) )         &
     &                 / ( cldsct(:,:,nb) + 1.0e-10 )
      enddo

!
      return
!...................................
      end subroutine cloudpar
!-----------------------------------



!-----------------------------------
      subroutine slingo                                                 &
!...................................

!  ---  inputs:
     &    ( cliq, reliq, NPTS, NLAY                                     &
!  ---  outputs:
     &,     cextvliq, cssavliq, casyvliq                                &
     &    )

!  *******************************************************************  !
!                                                                       !
!  define the single scattering parameters for cloud drops using the    !
!  Slingo parameterization for his spectral intervals.                  !
!                                                                       !
!  references:                                                          !
!                                                                       !
!  slingo, a., a gcm parameterization of the shortwave properties of    !
!      water clouds., j. atmos. sci.,46, 1419-1427, 1989.               !
!                                                                       !
!  notes: the cloud drop effective radius can only be 4.2 <= re <=      !
!        16.6 microns.                                                  !
!                                                                       !
!        the single scattering properties for wavenumbers < 2500 cm-1   !
!        are assigned the values in the first interval, since the       !
!        formulation is not defined for those wavenumbers.              !
!                                                                       !
!        the extinction coefficient is converted to kilometer**(-1)     !
!        the unit utilized by the shortwave routine Swresf.             !
!                                                                       !
!        a value of 1.0e-10 is added to the size so that no division    !
!        by zero occurs when the size is zero, in defining the          !
!        extinction coefficient.                                        !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! cliq     = the cloud drop liquid water concentration in grams / m**3  !
! reliq    = the cloud drop effective radius in microns                 !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! cextvliq = the specified spectral values of the extinction            !
!                  coefficient in kilometer**(-1) for drops             !
! cssavliq = the specified spectral values of the single-               !
!                  scattering albedo for drops                          !
! casyvliq = the specified spectral values of the asymmetry             !
!                  factor for drops                                     !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: cliq, reliq

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cextvliq, cssavliq, casyvliq

!  ---  locals:
      real (kind=kind_phys), dimension(NLIQCLDV) ::  a, b, c, d, e, f

      integer i, j, k, n, ni

      data a /-1.023e+01, 1.950e+01, 1.579e+01, 1.850e+01, 1.970e+01,   &
     &         2.237e+01, 2.463e+01, 2.551e+01, 2.589e+01, 2.632e+01,   &
     &         2.497e+01, 2.622e+01, 2.650e+01, 3.115e+01, 2.895e+01,   &
     &         2.831e+01, 2.838e+01, 2.672e+01, 2.698e+01, 2.668e+01,   &
     &         2.801e+01, 3.308e+01, 2.944e+01, 3.094e+01 /
      data b / 1.933e+03, 1.540e+03, 1.611e+03, 1.556e+03, 1.501e+03,   &
     &         1.452e+03, 1.420e+03, 1.401e+03, 1.385e+03, 1.365e+03,   &
     &         1.376e+03, 1.362e+03, 1.349e+03, 1.244e+03, 1.315e+03,   &
     &         1.317e+03, 1.300e+03, 1.320e+03, 1.315e+03, 1.307e+03,   &
     &         1.293e+03, 1.246e+03, 1.270e+03, 1.252e+03 /
      data c / 2.500e-02, 4.490e-01, 1.230e-01, 1.900e-04, 1.200e-03,   &
     &         1.200e-04, 2.400e-04, 6.200e-05,-2.800e-05,-4.600e-05,   &
     &         9.800e-06, 3.300e-06, 2.300e-06,-2.700e-07,-1.200e-07,   &
     &        -1.200e-06, 0.000e+00, 0.000e+00, 1.000e-06, 0.000e+00,   &
     &         1.000e-06,-3.000e-07,-6.500e-07, 7.900e-07 /
      data d / 1.220e-02, 1.540e-03, 9.350e-03, 2.540e-03, 2.160e-03,   &
     &         6.670e-04, 8.560e-04, 2.600e-04, 8.000e-05, 5.000e-05,   &
     &         2.100e-05, 2.800e-06, 1.700e-06, 1.400e-06, 4.400e-07,   &
     &         4.000e-07, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,   &
     &         0.000e+00, 2.360e-07, 4.330e-07, 3.690e-07 /
      data e / 7.260e-01, 8.310e-01, 8.510e-01, 7.690e-01, 7.400e-01,   &
     &         7.490e-01, 7.540e-01, 7.730e-01, 7.800e-01, 7.840e-01,   &
     &         7.830e-01, 8.060e-01, 8.090e-01, 8.040e-01, 8.180e-01,   &
     &         8.280e-01, 8.250e-01, 8.280e-01, 8.200e-01, 8.400e-01,   &
     &         8.360e-01, 8.390e-01, 8.410e-01, 8.440e-01 /
      data f / 6.652e-03, 6.102e-03, 2.814e-03, 5.171e-03, 7.469e-03,   &
     &         6.931e-03, 6.555e-03, 5.405e-03, 4.989e-03, 4.745e-03,   &
     &         5.035e-03, 3.355e-03, 3.387e-03, 3.520e-03, 2.989e-03,   &
     &         2.492e-03, 2.776e-03, 2.467e-03, 3.004e-03, 1.881e-03,   &
     &         2.153e-03, 1.946e-03, 1.680e-03, 1.558e-03 /

!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS
          if ( cliq(i,k) > _zero .and.                                  &
     &         (reliq(i,k)<4.2 .or. reliq(i,k)>16.6) ) then
            stop 'Slingo: cloud drop size out of range'
          end if
        enddo
      enddo

      do j = 1, NLIQCLDV
        cextvliq(:,:,j) = cliq(:,:) * (a(j) + b(j)/(reliq(:,:)+1.0e-10))
        cssavliq(:,:,j) = min( _one, _one-(c(j) + d(j)*reliq(:,:)) )
        casyvliq(:,:,j) = e(j) + f(j) * reliq(:,:)
      enddo

!
      return
!...................................
      end subroutine slingo
!-----------------------------------



!-----------------------------------
      subroutine fu                                                     &
!...................................

!  ---  inputs:
     &    ( cice, reice, NPTS, NLAY                                     &
!  ---  outputs:
     &,     cextvice, cssavice, casyvice                                &
     &    )

!  *******************************************************************  !
!                                                                       !
! define the single scattering parameters for ice crystals using the    !
! Fu parameterization for his spectral intervals.                       !
!                                                                       !
! references:                                                           !
!                                                                       !
! fu, q., an accurate parameterization of the solar radiative           !
!      properties of cirrus clouds for climate models., j. climate,     !
!      9, 2058-2082, 1996.                                              !
!                                                                       !
! notes: the ice crystal effective size (D^sub^ge in his paper) is      !
!        the input effective radius times 2.                            !
!                                                                       !
!        the ice crystal effective size (D^sub^ge in his paper) can     !
!        only be 18.6 <= D^sub^ge <= 130.2 microns.                     !
!                                                                       !
!        the single scattering properties for wavenumbers < 2000 cm-1   !
!        are assigned the values in the first interval, since the       !
!        formulation is not defined for those wavenumbers.              !
!                                                                       !
!        the extinction coefficient is converted to kilometer**(-1)     !
!        the unit utilized by the shortwave routine Swresf.             !
!                                                                       !
!        a value of 1.0e-10 is added to the size so that no division    !
!        by zero occurs when the size is zero, in defining the          !
!        extinction coefficient.                                        !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! cice    = the ice water concentation in grams / meter**3              !
! reice   = the ice crystal effective radius in microns                 !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! cextvice = the specified spectral values of the extinction            !
!                  coefficient for ice particles in kilometers**(-1)    !
! cssavice = the specified spectral values of the single-               !
!                  scattering albedo for ice particles                  !
! casyvice = the specified spectral values of the asymmetry             !
!                  factor for ice particles                             !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: cice, reice

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cextvice, cssavice, casyvice

!  ---  locals:
      real (kind=kind_phys), dimension(NICECLDV) :: a0fu, a1fu, b0fu,   &
     &       b1fu, b2fu, b3fu, c0fu, c1fu, c2fu, c3fu
      real (kind=kind_phys), dimension(NPTS,NLAY):: esize

      integer :: i, k, n

      data a0fu / -2.54823e-01, 1.87598e-01, 2.97295e-01, 2.34245e-01,  &
     &             4.89477e-01,-8.37325e-02, 6.44675e-01,-8.05155e-01,  &
     &             6.51659e-02, 4.13595e-01,-6.14288e-01, 7.31638e-02,  &
     &             8.10443e-02, 2.26539e-01,-3.04991e-01, 1.61983e-01,  &
     &             9.82244e-02,-3.03108e-02,-9.45458e-02, 1.29121e-01,  &
     &            -1.06451e-01,-2.58858e-01,-2.93599e-01,-2.66955e-01,  &
     &            -2.36447e-01 /
      data a1fu /  2.52909e+03, 2.51396e+03, 2.48895e+03, 2.48573e+03,  &
     &             2.48776e+03, 2.52504e+03, 2.47060e+03, 2.57600e+03,  &
     &             2.51660e+03, 2.48783e+03, 2.56520e+03, 2.51051e+03,  &
     &             2.51619e+03, 2.49909e+03, 2.54412e+03, 2.50746e+03,  &
     &             2.50875e+03, 2.51805e+03, 2.52061e+03, 2.50410e+03,  &
     &             2.52684e+03, 2.53815e+03, 2.54540e+03, 2.54179e+03,  &
     &             2.53817e+03 /
      data b0fu /  2.60155e-01, 1.96793e-01, 4.64416e-01, 9.05631e-02,  &
     &             5.83469e-04, 2.53234e-03, 2.01931e-03,-2.85518e-05,  &
     &            -1.48012e-07, 6.47675e-06,-9.38455e-06,-2.32733e-07,  &
     &            -1.57963e-07,-2.75031e-07, 3.12168e-07,-7.78001e-08,  &
     &            -8.93276e-08, 9.89368e-08, 5.08447e-07, 7.10418e-07,  &
     &             3.25057e-08,-1.98529e-07, 1.82299e-07,-1.00570e-07,  &
     &            -2.69916e-07 /
      data b1fu/   5.45547e-03, 5.75235e-03, 2.04716e-05, 2.93035e-03,  &
     &             1.18127e-03, 1.75078e-03, 1.83364e-03, 1.71993e-03,  &
     &             9.02355e-05, 2.18111e-05, 1.77414e-05, 6.41602e-06,  &
     &             1.72475e-06, 9.72285e-07, 4.93304e-07, 2.53360e-07,  &
     &             1.14916e-07, 5.44286e-08, 2.73206e-08, 1.42205e-08,  &
     &             5.43665e-08, 9.39480e-08, 1.12454e-07, 1.60441e-07,  &
     &             2.12909e-07 /
      data b2fu / -5.58760e-05,-5.29220e-05,-4.60375e-07,-1.89176e-05,  &
     &            -3.40011e-06,-8.00994e-06,-7.00232e-06,-7.43697e-06,  &
     &            -1.98190e-08, 1.83054e-09,-1.13004e-09, 1.97733e-10,  &
     &             9.02156e-11,-2.23685e-10, 1.79019e-10,-1.15489e-10,  &
     &            -1.62990e-10,-1.00877e-10, 4.96553e-11, 1.99874e-10,  &
     &            -9.24925e-11,-2.54540e-10,-1.08031e-10,-2.05663e-10,  &
     &            -2.65397e-10 /
      data b3fu /  1.97086e-07, 1.76618e-07, 2.03198e-09, 5.93361e-08,  &
     &             8.78549e-09, 2.31309e-08, 1.84287e-08, 2.09647e-08,  &
     &             4.01914e-11,-8.28710e-12, 2.37196e-12,-6.96836e-13,  &
     &            -3.79423e-13, 5.75512e-13,-7.31058e-13, 4.65084e-13,  &
     &             6.53291e-13, 4.56410e-13,-1.86001e-13,-7.81101e-13,  &
     &             4.53386e-13, 1.10876e-12, 4.99801e-13, 8.88595e-13,  &
     &             1.12983e-12 /
      data c0fu /  7.99084e-01, 7.59183e-01, 9.19599e-01, 8.29283e-01,  &
     &             7.75916e-01, 7.58748e-01, 7.51497e-01, 7.52528e-01,  &
     &             7.51277e-01, 7.52292e-01, 7.52048e-01, 7.51715e-01,  &
     &             7.52318e-01, 7.51779e-01, 7.53393e-01, 7.49693e-01,  &
     &             7.52131e-01, 7.51135e-01, 7.49856e-01, 7.48613e-01,  &
     &             7.47054e-01, 7.43546e-01, 7.40926e-01, 7.37809e-01,  &
     &             7.33260e-01 /
      data c1fu /  4.81706e-03, 4.93765e-03, 5.03025e-04, 2.06865e-03,  &
     &             1.74517e-03, 2.02709e-03, 2.05963e-03, 1.95748e-03,  &
     &             1.29824e-03, 1.14395e-03, 1.12044e-03, 1.10166e-03,  &
     &             1.04224e-03, 1.03341e-03, 9.61630e-04, 1.05446e-03,  &
     &             9.37763e-04, 9.09208e-04, 8.89161e-04, 8.90545e-04,  &
     &             8.86508e-04, 9.08674e-04, 8.90216e-04, 8.97515e-04,  &
     &             9.18317e-04 /
      data c2fu / -5.13220e-05,-4.84059e-05,-5.74771e-06,-1.59247e-05,  &
     &            -9.21314e-06,-1.17029e-05,-1.12135e-05,-1.02495e-05,  &
     &            -4.99075e-06,-3.27944e-06,-3.11826e-06,-2.91300e-06,  &
     &            -2.26618e-06,-2.13121e-06,-1.32519e-06,-2.32576e-06,  &
     &            -9.72292e-07,-6.34939e-07,-3.49578e-07,-3.44038e-07,  &
     &            -2.59305e-07,-4.65326e-07,-1.87919e-07,-2.17099e-07,  &
     &            -4.22974e-07 /
      data c3fu /  1.84420e-07, 1.65801e-07, 2.01731e-08, 5.01791e-08,  &
     &             2.15003e-08, 2.95195e-08, 2.73998e-08, 2.35479e-08,  &
     &             6.33757e-09,-2.42583e-10,-5.70868e-10,-1.37242e-09,  &
     &            -3.68283e-09,-4.24308e-09,-7.17071e-09,-3.58307e-09,  &
     &            -8.62063e-09,-9.84390e-09,-1.09913e-08,-1.10117e-08,  &
     &            -1.13305e-08,-1.05786e-08,-1.16760e-08,-1.16090e-08,  &
     &            -1.07976e-08 /

!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS
          esize(i,k) = 2.0 * reice(i,k)
          if ( cice(i,k) > _zero .and.                                  &
     &        (esize(i,k)<18.6 .or. esize(i,k)>130.2) ) then
            stop 'Fu: ice crystal size out of range'
          endif
        enddo
      enddo

      do n = 1, NICECLDV
        cextvice(:,:,n) = cice(:,:)                                     &
     &                  * ( a0fu(n) + a1fu(n)/(esize(:,:) + 1.0e-10) )
        cssavice(:,:,n) = min( _one, _one-(b0fu(n) + b1fu(n)*esize(:,:) &
     &                  + b2fu(n)*esize(:,:)**2+b3fu(n)*esize(:,:)**3) )
        casyvice(:,:,n) = c0fu(n) + c1fu(n)*esize(:,:)                  &
     &                  + c2fu(n)*esize(:,:)**2 + c3fu(n)*esize(:,:)**3
      enddo

!
      return
!...................................
      end subroutine fu
!-----------------------------------



!-----------------------------------
      subroutine savijarvi                                              &
!...................................

!  ---  inputs:
     &    ( crain, rerain, NPTS, NLAY                                   &
!  ---  outputs:
     &,     cextvrain, cssavrain, casyvrain                             &
     &    )

!  *******************************************************************  !
!                                                                       !
! define the single scattering parameters for rain drops using the      !
! Savijarvi parameterization for his spectral intervals.                !
!                                                                       !
! references:                                                           !
!                                                                       !
! savijarvi, h., shortwave optical properties of rain., tellus, 49a,    !
!      177-181, 1997.                                                   !
!                                                                       !
! notes: the rain drop effective radius can only be 16.6 < re < 5000    !
!        microns.                                                       !
!                                                                       !
!        the single scattering properties for wavenumbers < 2500 cm-1   !
!        are assigned the values in the first interval, since the       !
!        formulation is not defined for those wavenumbers.              !
!                                                                       !
!        the extinction coefficient is converted to kilometer**(-1)     !
!        the unit utilized by the shortwave routine Swresf.             !
!                                                                       !
!        a value of 1.0e-10 is added to the size so that no division    !
!        by zero occurs when the size is zero, in defining the          !
!        extinction coefficient.                                        !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! crain   = the rain drop water concentration in grams / meter**3       !
! rerain  = the rain drop effective radius in microns                   !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! cextvrain   = the specified spectral values of the extinction         !
!                   coefficient for rain in kilometers**(-1)            !
! cssavrain = the specified spectral values of the single-              !
!                   scattering albedo for rain                          !
! casyvrain = the specified spectral values of the asymmetry            !
!                   factor for rain                                     !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: crain,rerain

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cextvrain, cssavrain, casyvrain

!  ---  locals:
      real (kind=kind_phys), dimension (NRAINCLDV) ::  a, asymm, b
      real (kind=kind_phys), dimension (NPTS,NLAY) ::  rcap

      integer :: i, k, n

      data a     / 4.65e-1, 2.64e-1, 1.05e-2, 8.00e-5 /
      data b     / 1.00e-3, 9.00e-2, 2.20e-1, 2.30e-1 /
      data asymm / 9.70e-1, 9.40e-1, 8.90e-1, 8.80e-1 /
!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS
          if ( crain(i,k) > _zero .and.                                 &
     &        (rerain(i,k)<16.6 .or. rerain(i,k)>5.0e+3) ) then
            stop 'Savisolar: rain drop size out of range'
          endif
        enddo
      enddo

      rcap(:,:) = ( rerain(:,:) / 5.0e+2 ) ** 4.348e+0

      do n = 1, NRAINCLDV
        cextvrain(:,:,n) = 1.505e+3 * crain(:,:)/(rerain(:,:) + 1.0e-10)
        cssavrain(:,:,n) = min( _one, _one-(a(n)*(rcap(:,:)**b(n))) )
        casyvrain(:,:,n) = asymm(n)
      enddo

!
      return
!...................................
      end subroutine savijarvi
!-----------------------------------



!-----------------------------------
      subroutine snowsw                                                 &
!...................................

!  ---  inputs:
     &    ( csnow, NPTS, NLAY                                           &
!  ---  outputs:
     &,     cextvsnow, cssavsnow, casyvsnow                             &
     &    )

!  *******************************************************************  !
!                                                                       !
! define the single scattering parameters for snow using the Fu         !
! parameterization for his spectral intervals.                          !
!                                                                       !
! author: leo donner, gfdl, 11 Sept 98                                  !
!                                                                       !
! references:                                                           !
!                                                                       !
! fu, q., et al., (See notes from Kuo-Nan Liou, 1 Sept 98). (SNOW)      !
!                                                                       !
! notes: the single scattering properties for wavenumbers < 2500 cm-1   !
!        are assigned the values in the first interval, since the       !
!        formulation is not defined for those wavenumbers.              !
!                                                                       !
!        the extinction coefficient is in units of kilometer**(-1)      !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! csnow     = the snow concentration in grams / meter**3                !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! cextvsnow = the specified spectral values of the extinction           !
!                   coefficient for snow in kilometers**(-1)            !
! cssavsnow = the specified spectral values of the single-              !
!                   scattering albedo for snow                          !
! casyvsnow = the specified spectral values of the asymmetry            !
!                   factor for snow                                     !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: csnow

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       cextvsnow, cssavsnow, casyvsnow

!  ---  locals:
      real (kind=kind_phys), dimension(NSNOWCLDV) ::  asymm, ext, ssalb
      real (kind=kind_phys) :: cref

      data asymm / .96373, .98141, .97816, .96820, .89940, .89218 /
      data ext   / .83951, .83946, .83941, .83940, .83940, .83939 /
      data ssalb / .53846, .52579, .53156, .56192, .97115, .99991 /
      data cref  / 0.5 /

      integer :: i, k, n

      do n = 1, NSNOWCLDV
        cextvsnow(:,:,n) = ext(n) * csnow(:,:) / cref
        cssavsnow(:,:,n) = ssalb(n)
        casyvsnow(:,:,n) = asymm(n)
      enddo

!
      return
!...................................
      end subroutine snowsw
!-----------------------------------



!-----------------------------------
      subroutine thickavg                                               &
!...................................

!  ---  inputs:
     &    ( nv1, nv2, solflxv                                           &
     &,     nvs, extv, ssalbv, asymmv                                   &
     &,     NPTS, NLAY                                                  &
!  ---  outputs:
     &,     extb, ssalbb, asymmb                                        &
     &    )

!  *******************************************************************  !
!                                                                       !
! use the thick-averaging technique to define the single-scattering     !
! properties of the parameterization band spectral intervals from the   !
! specified spectral intervals of the particular scatterer.             !
!                                                                       !
! references:                                                           !
!                                                                       !
! edwards,j.m. and a. slingo, studies with a flexible new radiation     !
!      code I: choosing a configuration for a large-scale model.,       !
!      q.j.r. meteorological society, 122, 689-719, 1996.               !
!                                                                       !
! note: the 1.0e-10 factor to calculate asymmband is to prevent         !
!       division by zero.                                               !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent in:                                                            !
!                                                                       !
! nv1      = interval number for the specified single-scattering        !
!              properties corresponding to the first psuedo-            !
!              monochromatic frequency in a given parameterization      !
!              band                                                     !
! nv2      = interval number for the specified single-scattering        !
!              properties corresponding to the last psuedo-             !
!              monochromatic frequency in a given parameterization      !
!              band                                                     !
! nvs      = number of specified scattering spectral intervals          !
! extv     = the specified spectral values of the extinction            !
!              coefficient                                              !
! ssalbv   = the specified spectral values of the single-               !
!              scattering albedo                                        !
! asymmv   = the specified spectral values of the asymmetry             !
!              factor                                                   !
! solflxv  = the solar flux in each specified scattering spectral       !
!              interval                                                 !
! solflxband = the solar flux in each parameterization band             !
!                                                                       !
! --------------------------------------------------------------------- !
!                                                                       !
! intent out:                                                           !
!                                                                       !
! extb   =  the parameterization band values of the extinction          !
!              coefficient                                              !
! ssalbb =  the parameterization band values of the single-             !
!              scattering albedo                                        !
! asymmb =  the parameterization band values of the asymmetry           !
!              factor                                                   !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: nvs, NPTS, NLAY

      real (kind=kind_phys), dimension(:,:,:), intent(in) ::            &
     &       extv, ssalbv, asymmv
      real (kind=kind_phys), dimension(:,:),   intent(in) :: solflxv

      integer, dimension(:),  intent(in) :: nv1, nv2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) ::           &
     &       extb, ssalbb, asymmb

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,NLAY) :: refb, refthick,    &
     &       sp, sumk, sumomegak, sumomegakg, sumrefthick

      integer :: i, k, nb, ni
!
!===> ...  begin here
!

      do nb = 1, NBANDS

        do k = 1, NLAY
          do i = 1, NPTS
            sumk       (i,k) = _zero
            sumomegak  (i,k) = _zero
            sumomegakg (i,k) = _zero
            sumrefthick(i,k) = _zero
          enddo
        enddo

        do ni = nv1(nb), nv2(nb)
          do k = 1, NLAY
            do i = 1, NPTS
              sp(i,k) = sqrt( (_one - ssalbv(i,k,ni))                   &
     &                      / (_one - ssalbv(i,k,ni)*asymmv(i,k,ni)) )
              refthick(i,k) = (_one - sp(i,k)) / (_one + sp(i,k))
              sumrefthick(i,k) = sumrefthick(i,k)                       &
     &                         + refthick(i,k)*solflxv(nb,ni)
              sumk(i,k) = sumk(i,k) + extv(i,k,ni)*solflxv(nb,ni)
              sumomegak(i,k) = sumomegak(i,k)                           &
     &                     + ssalbv(i,k,ni)*extv(i,k,ni)*solflxv(nb,ni)
              sumomegakg(i,k) = sumomegakg(i,k) + ssalbv(i,k,ni)        &
     &                     * extv(i,k,ni)*asymmv(i,k,ni)*solflxv(nb,ni)
            enddo
          enddo
        enddo

        do k = 1, NLAY
          do i = 1, NPTS
            extb(i,k,nb) = sumk(i,k) / solflxband(nb)
            asymmb(i,k,nb) = sumomegakg(i,k)/(sumomegak(i,k)+1.0e-10)
            refb(i,k) = sumrefthick(i,k) / solflxband(nb)
            ssalbb(i,k,nb) = 4.0*refb(i,k)/( (_one+refb(i,k))**2        &
     &                     - asymmb(i,k,nb)*(_one-refb(i,k))**2 )
          enddo
        enddo

      enddo    ! end nb-loop

!
      return
!...................................
      end subroutine thickavg
!-----------------------------------

!
!........................................!
      end module module_radsw_main       !
!========================================!

