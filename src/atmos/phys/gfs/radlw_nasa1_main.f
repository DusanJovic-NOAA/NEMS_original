!!!!!  ==========================================================  !!!!!
!!!!!              nasa1 radiation package description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the nasa1 package includes these parts:                           !
!                                                                      !
!       'radlw_nasa1_param.f'                                          !
!       'radlw_nasa1_datatb.f'                                         !
!       'radlw_nasa1_main.f'                                           !
!                                                                      !
!    the 'radlw_nasa1_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_nasa1_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_tables'      -- data for each lw spectral band   !
!                                                                      !
!    the 'radlw_nasa1_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main nasa1 lw radiation routine                 !
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
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_nasa1_param.f'                                          !
!       'radlw_nasa1_datatb.f'                                         !
!       'radlw_nasa1_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!                                                                      !
!                                                                      !
!    author:                                                           !
!      ming-dah chou, code 913, nasa/goddard space flight center,      !
!      greenbelt, md 20771.                                            !
!      phone: 301-614-6192, fax: 301-614-6307,                         !
!      e-mail: chou@climate.gsfc.nasa.gov                              !
!                                                                      !
!    reference:                                                        !
!      chow,m.-d., m.j.suarez, x.-z.liang, and m.m.-h.yan, 2001: a     !
!      thermal infrared radiation parameterization for atmospheric     !
!      studies, nasa/tm-2001-104606, vol. 19                           !
!                                                                      !
!    the original program description:                                 !
!                                                                      !
!    the equation numbers noted in this code follows the latest        !
!    version (july 2002) of the nasa tech. memo. (2001), which can     !
!    be accessed at http://climate.gsfc.nasa.gov/~chou/clirad_lw       !
!                                                                      !
!    change in july 2002:                                              !
!      the effective Planck functions of a layer are separately        !
!      computed for the upward and downward emission (bu and bd).      !
!      for a optically thick cloud layer, the upward emitting          !
!      temperature will be close the cloud top temperature, and the    !
!      downward emitting temperature will be close the cloud base      !
!      temperature.                                                    !
!                                                                      !
!    recent changes:                                                   !
!      subroutines for planck functions                                !
!      subroutines for cloud overlapping                               !
!      eliminate "rflx" and "rflc". fold the flux calculations in      !
!        band 10 to that of the other bands.                           !
!      return the calculations when ibn=10 and itrace=0.               !
!      the number of aerosol types is allowed to be more than one.     !
!      include sub-grid surface variability and vegetation canopy.     !
!      include the ckd continuum absorption ceofficient as an option.  !
!                                                                      !
!    ice and liquid cloud particles are allowed to co-exist in any of  !
!      the np layers.                                                  !
!                                                                      !
!    if no information is available for the effective cloud particle   !
!      size, reff, default values of 10 micron for liquid water and 75 !
!      micron for ice can be used.                                     !
!                                                                      !
!    the maximum-random assumption is applied for cloud overlapping.   !
!      clouds are grouped into high, middle, and low clouds separated  !
!      by the level indices ict and icb.  within each of the three     !
!      groups, clouds are assumed maximally overlapped.  clouds among  !
!      the three groups are assumed randomly overlapped. the indices   !
!      ict and icb correspond approximately to the 400 mb and 700 mb   !
!      levels.                                                         !
!                                                                      !
!    various types of aerosols are allowed to be in any of the np      !
!      layers. aerosol optical properties can be specified as functions!
!      of height and spectral band. (moved out off this package)       !
!                                                                      !
!    the surface can be divided into a number of sub-regions either    !
!      with or without vegetation cover. reflectivity and emissivity   !
!      can be specified for each sub-region. (moved out off this pkg)  !
!                                                                      !
!    the ir spectrum is divided into nine bands:                       !
!       band     wavenumber (/cm)   absorber                           !
!        1           0 - 340           h2o                             !
!        2         340 - 540           h2o                             !
!        3         540 - 800           h2o,cont,co2                    !
!        4         800 - 980           h2o,cont,co2,f11,f12,f22        !
!        5         980 - 1100          h2o,cont,o3,co2,f11             !
!        6        1100 - 1215          h2o,cont,n2o,ch4,f12,f22        !
!        7        1215 - 1380          h2o,cont,n2o,ch4                !
!        8        1380 - 1900          h2o                             !
!        9        1900 - 3000          h2o                             !
!    in addition, a narrow band in the 17 micrometer region (band 10)  !
!    is added to compute flux reduction due to n2o                     !
!       10         540 - 620           h2o,cont,co2,n2o                !
!                                                                      !
!    band 3 (540-800/cm) is further divided into 3 sub-bands :         !
!      subband   wavenumber (/cm)                                      !
!        3a        540 - 620                                           !
!        3b        620 - 720                                           !
!        3c        720 - 800                                           !
!                                                                      !
!    notes:                                                            !
!      (1) scattering is parameterized for clouds and aerosols.        !
!      (2) diffuse cloud and aerosol transmissions are computed        !
!          from exp(-1.66*tau).                                        !
!      (3) If there are no clouds, flxntc=flxnt0.                      !
!      (4) plevel(1) is the pressure at the top of the model           !
!          atmosphere, and plevel(np+1) is the surface pressure.       !
!      (5) downward flux is positive and upward flux is negative.      !
!      (6) sfcem and dfdts are negative because upward flux is defined !
!          as negative.                                                !
!      (7) for questions and coding errors, please contact             !
!          ming-dah chou, code 913, nasa/goddard space flight center,  !
!          greenbelt, md 20771.                                        !
!          phone: 301-614-6192, fax: 301-614-6307,                     !
!          e-mail: chou@climate.gsfc.nasa.gov                          !
!                                                                      !
!                                                                      !
!                                                                      !
!    ncep modifications history log:                                   !
!       --- 2003,  yu-tai hou   -- received original code from nasa    !
!       jul 2006,  yu-tai hou                                          !
!                  recoded to fit in ncep unified radiation package    !
!       mar 2007,  yu-tai hou                                          !
!                  add unified lw aerosol support                      !
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
      use machine,           only : kind_phys
!     use physcons,          only : con_g, con_cp, con_gasv, con_amd,   &
!    &                              con_amo3, con_rd, con_p0, con_ttp

      use module_radlw_parameters
      use module_radlw_cntr_para
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
!     character(24), parameter :: VTAGLW='NASA-LW vjul02  jul 2006'
      character(24), parameter :: VTAGLW='NASA-LW vjul02  Apr 2007'

!  ---  constant values
      real (kind=kind_phys), parameter :: _zero = 0.0
      real (kind=kind_phys), parameter :: _one  = 1.0
      real (kind=kind_phys), parameter :: fdiff = 1.66

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

!  ---  those data will be set up only once by "rlwinit"

!  ...  heatfac is the factor for heating rates (in k/day, or k/sec 
!       set by subroutine 'rlwinit')

      real (kind=kind_phys) :: heatfac

      public lwrad, rlwinit


! =================
      contains
! =================


!-----------------------------------
      subroutine lwrad                                                  &
!...................................

!  ---  inputs:
     &     ( plyr, plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                     &
     &       clouds,iovr,aerosols,sfcemis,                              &
     &       NPTS, NLAY, NLP1, iflip, lprnt,                            &
!  ---  outputs:
     &       hlwc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HLW0,HLWB,FLXPRF                                           &
     &     )

!  ====================  defination of variables  ===================  !
!                                                                      !
!  input variables:                                                    !
!     plyr   (NPTS,NLAY)    - layer pressures (mb)                     !
!     plvl   (NPTS,NLP1)    - interface pressures (mb)                 !
!     tlyr   (NPTS,NLAY)    - layer temperature (k)                    !
!     tlvl   (NPTS,NLP1)    - interface temperatures (k)               !
!     qlyr   (NPTS,NLAY)    - layer h2o mixing ratio (gm/gm)*see inside!
!     olyr   (NPTS,NLAY)    - layer o3 mixing ratio (gm/gm) *see inside!
!     gasvmr (NPTS,NLAY,:)  - atmospheric gases amount:                !
!                       (check module_radiation_gases for definition)  !
!       gasvmr(:,:,1)   -      co2 volume mixing ratio                 !
!       gasvmr(:,:,2)   -      n2o volume mixing ratio                 !
!       gasvmr(:,:,3)   -      ch4 volume mixing ratio                 !
!       gasvmr(:,:,4)   -      o2  volume mixing ratio     (not used)  !
!       gasvmr(:,:,5)   -      co  volume mixing ratio     (not used)  !
!       gasvmr(:,:,6)   -      cfc11 volume mixing ratio               !
!       gasvmr(:,:,7)   -      cfc12 volume mixing ratio               !
!       gasvmr(:,:,8)   -      cfc22 volume mixing ratio               !
!       gasvmr(:,:,9)   -      cfccl4 volume mixing ratio  (not used)  !
!       gasvmr(:,:,10)  -      cfc113 volume mixing ratio  (not used)  !
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
!     aerosols(NPTS,NLAY,NBDLW,:) - aerosol optical properties         !
!                       (check module_radiation_aerosols for definition!
!        (:,:,:,1)          - optical depth                            !
!        (:,:,:,2)          - single scattering albedo                 !
!        (:,:,:,3)          - asymmetry parameter                      !
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
!     iaerlw                - control flag for aerosols                !
!                             =0: do not include aerosol effect        !
!                             >0: include aerosol effect               !
!     itrace                - control flag for trace gases             !
!                             (ch4,n2o,cfcs)                           !
!                             =0: do not include trace gases           !
!                             =1: include all trace gases and 2 minor  !
!                              co2 bands in the window region          !
!     iflagc                - control flag for cloud optical property  !
!                             =1: input cloud optical depth (diagnostic!
!                             =2: input cloud cwp/cip/crp...           !
!     iovcst                - control flag for cloud fraction assuption!
!                             =1: overcast assumption cld either 0 or 1!
!                             =2: allowing for fractional clouds       !
!     itrans                - flag for transm funct for h2o,co2,o3     !
!                             =1: transm funcs are computed using the  !
!                              k-distr with linear press scaling for   !
!                              all bands except band 5. cooling rates  !
!                              are not accurately for press < 10 mb.   !
!                              it is faster than itrans=2.             !
!                             =2: transm funcs in the co2,o3, and 3    !
!                              h2o bands with strong abs are computed  !
!                              using table look-up. cooling rates are  !
!                              accurately from the sfc up to 0.01 mb.  !
!                                                                      !
!  output variables:                                                   !
!     hlwc   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)  !
!     topflx (NPTS)         - radiation fluxes at top, component:      !
!                        (check module_radlw_paramters for definition) !
!        upfxc                 total sky upward flux at top (w/m2)     !
!        upfx0                 clear sky upward flux at top (w/m2)     !
!     sfcflx (NPTS)         - radiation fluxes at sfc, component:      !
!                        (check module_radlw_paramters for definition) !
!        upfxc                 total sky upward flux at sfc (w/m2)     !
!        dnfxc                 total sky downward flux at sfc (w/m2)   !
!        dnfx0                 clear sky downward flux at sfc (w/m2)   !
!                                                                      !
!! optional output variables:                                          !
!     hlwb(NPTS,NLAY,NBDLW) - spectral band total sky heating rates    !
!     hlw0   (NPTS,NLAY)    - total sky heating rate (k/day or k/sec)  !
!     flxprf (NPTS,NLP1)    - level radiative fluxes (w/m2), components!
!                        (check module_radlw_paramters for definition) !
!        upfxc                 total sky upward flux                   !
!        dnfxc                 total sky dnward flux                   !
!        upfx0                 clear sky upward flux                   !
!        dnfx0                 clear sky dnward flux                   !
!                                                                      !
!  module internal variables:                                          !
!    pl    : level pressure                                mb          !
!    ta    : layer temperature                              k          !
!    wa    : layer specific humidity                       g/g         !
!    oa    : layer ozone mixing ratio by mass              g/g         !
!    rco2  : co2 mixing ratio by volume                   pppv         !
!    rn2o  : n2o mixing ratio by volume                   pppv         !
!    rch4  : ch4 mixing ratio by volume                   pppv         !
!    cfc11 : cfc11 mixing ratio by volume                 pppv         !
!    cfc12 : cfc12 mixing ratio by volume                 pppv         !
!    cfc22 : cfc22 mixing ratio by volume                 pppv         !
!    tb    : surface air temperature                        k          !
!    tg    : land or ocean surface temperature              k          !
!    eg    : land or ocean surface emissivity             fraction     !
!    fcld  : cloud amount                                 fraction     !
!    cwc   : cloud water mixing ratio                      g/g         !
!      index 1-ice, 2-liquid, 3-rain                                   !
!    reff  : effective cloud-particle size                micrometer   !
!      index 1-ice, 2-liquid, 3-rain                                   !
!    taucl : cloud optical thickness                       --          !
!      index 1-ice, 2-liquid, 3-rain                                   !
!    ict   : lev index separating hi and mid clouds        --          !
!    icb   : lev index separating mid and low clouds       --          !
!    taual : aerosol optical thickness                     --          !
!    ssaal : aerosol single-scattering albedo              --          !
!    asyal : aerosol asymmetry factor                      --          !
!                                                                      !
!    flxntc: net downward flux, all-sky                   w/m**2       !
!    flxnt0: net downward flux, clear-sky                 w/m**2       !
!not dfdts : sensitivity of net dnwd flx to sfc temp.     w/m**2/k     !
!not sfcem : emission by the surface                      w/m**2       !
!                                                                      !
!  data used in table look-up for transmittance calculations:          !
!        c1 ,  c2,  c3 : for co2 (band 3)                              !
!        o1 ,  o2,  o3 : for  o3 (band 5)                              !
!        h11, h12, h13 : for h2o (band 1)                              !
!        h21, h22, h23 : for h2o (band 2)                              !
!        h81, h82, h83 : for h2o (band 8)                              !
!                                                                      !
!                                                                      !
!  subroutine lwrad is called by : grrad - ncep radiation driver       !
!                                                                      !
!  subroutines called by lwrad : planck, sfcflux, h2oexps, conexps,    !
!                                co2exps, co2exps, n2oexps, ch4exps,   !
!                                comexps, cfcexps, b10exps, tablup,    !
!                                h2okdis, co2kdis, n2okdis, ch4kdis,   !
!                                comkdis, cfckdis, b10kdis, cldovlp    !
!                                                                      !
!  ******************************************************************  !
!
      use module_radlw_tables
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, NLP1, iovr, iflip

      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(:,:),  intent(in) :: plvl, tlvl, &
     &       plyr, tlyr, qlyr, olyr
      real (kind=kind_phys), dimension(:,:,:),intent(in) :: gasvmr,     &
     &       clouds
      real (kind=kind_phys), intent(in) :: sfcemis(:), aerosols(:,:,:,:)

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: hlwc

      type (topflw_type),    dimension(:),   intent(out) :: topflx
      type (sfcflw_type),    dimension(:),   intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys),dimension(:,:,:),optional,intent(out):: hlwb
      real (kind=kind_phys),dimension(:,:),optional,intent(out):: hlw0
      type (proflw_type),   dimension(:,:),optional,intent(out):: flxprf

!  ---  locals:
      real (kind=kind_phys), dimension(NPTS,0:NLP1) :: blayer, bd, bu

      real (kind=kind_phys), dimension(NPTS,NLP1)   :: trant, transfc,  &
     &       trantcr, flxu, flxd, flcu, flcd, dfdts, blevel, fclr, pl,  &
     &       flxdnc, flxupc, flxntc, flxdn0, flxup0, flxnt0

      real (kind=kind_phys), dimension(NPTS,NLP1,NBDLW) :: flxntb

      real (kind=kind_phys), dimension(NPTS,NLAY)   :: pa, ta, wa, oa,  &
     &       dt, dp, dh2o, dcont, dco2, do3, dn2o, dch4, df11, df12,    &
     &       df22, f11exp, f12exp, f22exp, taua, ssaa, asya, taerlyr,   &
     &       tcldlyr, fcld, rco2, rn2o, rch4, cfc11, cfc12, cfc22

      real (kind=kind_phys), dimension(NPTS)        :: xlayer, tranal,  &
     &       tg, tb, tx, x1, x2, x3, tf11, tf12, tf22, bs, dbs, rflxs,  &
     &       sfcem, cldhi, cldmd, cldlw

      real (kind=kind_phys), dimension(NPTS,NK0)    :: th2o, tcom
      real (kind=kind_phys), dimension(NPTS,NK1)    :: tn2o, tch4
      real (kind=kind_phys), dimension(NPTS,NK2)    :: tcon
      real (kind=kind_phys), dimension(NPTS,NBDLW)  :: eg

      real (kind=kind_phys), dimension(NPTS,NLAY,NK0) :: h2oexp, comexp
      real (kind=kind_phys), dimension(NPTS,NLAY,NK1) :: n2oexp, ch4exp
      real (kind=kind_phys), dimension(NPTS,NLAY,NK2) :: conexp
      real (kind=kind_phys), dimension(NPTS,NLAY,3)   :: cwp,reff,taucl

      real (kind=kind_phys), dimension(NPTS,     NK0,2) :: tco2
      real (kind=kind_phys), dimension(NPTS,NLAY,NK0,2) :: co2exp

      real (kind=kind_phys) :: xx, yy, p1, dwe, dpe, a1, b1, fk1,       &
     &       a2, b2, fk2, w1, w2, w3, g1, g2, g3, ww, gg, ff, tauc,     &
     &       reff1, reff2

      integer, dimension(NPTS) :: it, im, ib, ict, icb
      integer :: i, j, k, k1, k2, i1, i2, i3, ibn, m, ne

      logical :: oznbnd, co2bnd, h2otbl, conbnd, n2obnd, ch4bnd,        &
     &           combnd, f11bnd, f12bnd, f22bnd, b10bnd, ltb

!
!===> ...  begin here
!

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

!  ---  the internal array is always from top to surface

      if ( iflip == 0 ) then      ! input from toa to sfc

!  ---  compute layer pressure (pa) and layer temperature minus 250K (dt)

        do k = 1, NLAY
          do i = 1, NPTS
            pl(i,k) = plvl(i,k)
            pa(i,k) = plyr(i,k)
!orig       pa(i,k) = 0.5 * (plvl(i,k) + plvl(i,k+1))
            ta(i,k) = tlyr(i,k)
!test use
!           wa(i,k) = max(2.0e-7, qlyr(i,k))                 ! input mass mixing ratio
!ncep model use
            wa(i,k) = max(2.0e-7, qlyr(i,k)/(_one-qlyr(i,k)))! input specific humidity
            oa(i,k) = max(_zero,  olyr(i,k))                 ! input mass mixing ratio
!           oa(i,k) = max(_zero,  olyr(i,k)*ramdo3)          ! input vol mixing ratio
            rco2(i,k) = gasvmr(i,k,1)
            rn2o(i,k) = gasvmr(i,k,2)
            rch4(i,k) = gasvmr(i,k,3)
            cfc11(i,k)= gasvmr(i,k,6)
            cfc12(i,k)= gasvmr(i,k,7)
            cfc22(i,k)= gasvmr(i,k,8)
          enddo
        enddo

        do i = 1, NPTS
          pl(i,NLP1) = plvl(i,NLP1)
        enddo

!  ---  assign surface-skin temp (tg) and surface-air temp (tb)

!ltb    if ( size(tlyr,dim=2) == NLP1 ) then
!         do i = 1, NPTS
!           tg(i) = tlvl(i,NLP1)
!           tb(i) = tlyr(i,NLP1)  !! assume tb is at NLP1 place
!         enddo

!         ltb = .true.
!       else
!  ---  note: we assume tb = tg to avoid extrapolation later

          do i = 1, NPTS
            tg(i) = tlvl(i,NLP1)
            tb(i) = tlvl(i,NLP1)
          enddo

!         ltb = .false.
!         ltb = .true.       !! use ground temp, avoid extrapolate later
!ltb    endif

!  ---  compute layer cloud water amount (gm/m**2)
!       index is 1 for ice, 2 for waterdrops and 3 for raindrops.

        if ( iflagc == 2 ) then     ! use prognostic cloud method
          do k = 1, NLAY
            do i = 1, NPTS
              fcld(i,k)  = clouds(i,k,1)
              cwp (i,k,2)= clouds(i,k,2)
              reff(i,k,2)= clouds(i,k,3)
              cwp (i,k,1)= clouds(i,k,4)
              reff(i,k,1)= clouds(i,k,5)
              cwp (i,k,3)= clouds(i,k,6)
              reff(i,k,3)= clouds(i,k,7)
            enddo
          enddo
        else                        ! use diagnostic cloud method
          do k = 1, NLAY
            do i = 1, NPTS
              fcld (i,k)  = clouds(i,k,1)
              taucl(i,k,1)= clouds(i,k,2)
            enddo
          enddo
        endif                       ! end if_iflagc

      else                        ! input from sfc to toa

!  ---  compute layer pressure (pa) and layer temperature minus 250K (dt)

        do k = 1, NLAY
          k1 = NLP1 - k
          do i = 1, NPTS
            pl(i,k+1) = plvl(i,k1)
            pa(i,k) = plyr(i,k1)
!orig       pa(i,k) = 0.5 * (plvl(i,k1) + plvl(i,k1+1))
            ta(i,k) = tlyr(i,k1)
!test use
!           wa(i,k) = max(2.0e-7, qlyr(i,k1))                  ! input mass mixing ratio
!ncep model use
            wa(i,k) = max(2.0e-7, qlyr(i,k1)/(_one-qlyr(i,k1)))! input specific humidity
            oa(i,k) = max(_zero,  olyr(i,k1))                  ! input mass mixing ratio
!           oa(i,k) = max(_zero,  olyr(i,k1)*ramdo3)           ! input vol mixing ratio
            rco2(i,k) = gasvmr(i,k1,1)
            rn2o(i,k) = gasvmr(i,k1,2)
            rch4(i,k) = gasvmr(i,k1,3)
            cfc11(i,k)= gasvmr(i,k1,6)
            cfc12(i,k)= gasvmr(i,k1,7)
            cfc22(i,k)= gasvmr(i,k1,8)
          enddo
        enddo

        do i = 1, NPTS
          pl(i,1) = plvl(i,NLP1)
        enddo

!  ---  assign surface-skin temp (tg) and surface-air temp (tb)

!ltb    if ( size(tlyr,dim=2) == NLP1 ) then
!         do i = 1, NPTS
!           tg(i) = tlvl(i,1)
!           tb(i) = tlyr(i,NLP1)  !! assume tb is at NLP1 place
!         enddo

!         ltb = .true.
!       else
          do i = 1, NPTS
            tg(i) = tlvl(i,1)
            tb(i) = tlvl(i,1)
          enddo

!         ltb = .false.
!         ltb = .true.       !! use ground temp, avoid extrapolate later
!ltb    endif

!  ---  compute layer cloud water amount (gm/m**2)
!       index is 1 for ice, 2 for waterdrops and 3 for raindrops.

        if ( iflagc == 2 ) then     ! use prognostic cloud method
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              fcld(i,k)  = clouds(i,k1,1)
              cwp (i,k,2)= clouds(i,k1,2)
              reff(i,k,2)= clouds(i,k1,3)
              cwp (i,k,1)= clouds(i,k1,4)
              reff(i,k,1)= clouds(i,k1,5)
              cwp (i,k,3)= clouds(i,k1,6)
              reff(i,k,3)= clouds(i,k1,7)
            enddo
          enddo
        else                        ! use diagnostic cloud method
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              fcld (i,k)  = clouds(i,k1,1)
              taucl(i,k,1)= clouds(i,k1,2)
            enddo
          enddo
        endif                       ! end if_iflagc

      endif                       ! if_iflip

!  ---  compute layer absorber amount
!       dh2o : water vapor amount (g/cm**2)
!       dcont: scaled water vapor amt for continuum abs (g/cm**2)
!       dco2 : co2 amount (cm-atm)stp
!       do3  : o3 amount (cm-atm)stp
!       dn2o : n2o amount (cm-atm)stp
!       dch4 : ch4 amount (cm-atm)stp
!       df11 : cfc11 amount (cm-atm)stp
!       df12 : cfc12 amount (cm-atm)stp
!       df22 : cfc22 amount (cm-atm)stp
!       the factor 1.02 is equal to 1000/980
!       factors 789 and 476 are for unit conversion
!       the factor 0.001618 is equal to 1.02/(.622*1013.25)
!       the factor 6.081 is equal to 1800/296

!     g1 = 10.0 / con_g
      g1 = 1.02
!     g2 = 1.0e4*(con_rd/con_g)*(con_amd/con_amo3)*con_ttp/con_p0
      g2 = 476.0
!     g3 = (10.0/con_g) * (con_gasv*1.0e6/ con_amd)
      g3 = 789.0
!     a1 = g1 / (con_eps*con_p0*1.0e-2)
      a1 = 0.001618

      do k = 1, NLAY
        do i = 1, NPTS
          dp(i,k) = pl(i,k+1) - pl(i,k)
          dt(i,k) = ta(i,k) - 250.0

          dh2o(i,k) = max( 1.0e-8, g1*wa  (i,k)*dp(i,k) )
          do3 (i,k) = max( 1.0e-6, g2*oa  (i,k)*dp(i,k) )
          dco2(i,k) = max( 1.0e-4, g3*rco2(i,k)*dp(i,k) )

          dch4(i,k) = g3 * rch4 (i,k) * dp(i,k)
          dn2o(i,k) = g3 * rn2o (i,k) * dp(i,k)
          df11(i,k) = g3 * cfc11(i,k) * dp(i,k)
          df12(i,k) = g3 * cfc12(i,k) * dp(i,k)
          df22(i,k) = g3 * cfc22(i,k) * dp(i,k)

!  ---  compute scaled water vapor amount for h2o continuum absorption
!       following eq. (4.21).

          xx = pa(i,k) * a1 * wa(i,k) * wa(i,k) * dp(i,k)
          dcont(i,k) = xx * exp( 1800.0/ta(i,k) - 6.081 )
        enddo
      enddo

!  ---  find boundaries of cloud domains

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

      do j = 1, NBDLW
        do i = 1, NPTS
          eg(i,j) = sfcemis(i)
        enddo
      enddo

!  ---  the surface (np+1) is treated as a layer filled with black clouds.
!       transfc is the transmittance between the surface and a pressure level.
!       trantcr is the clear-sky transmittance between the surface and a
!       pressure level.
 
      do i = 1, NPTS
!       sfcem(i)        = _zero
        transfc(i,NLP1) = _one
        trantcr(i,NLP1) = _one
      enddo

!  ---  initialize fluxes

      do k = 1, NLP1
        do i = 1, NPTS
          flxdnc(i,k) = _zero
          flxupc(i,k) = _zero
          flxdn0(i,k) = _zero
          flxup0(i,k) = _zero
!         dfdts (i,k) = _zero
        enddo
      enddo

!  ---  integration over spectral bands

      Lab_do_ibn : do ibn = 1, 10

!       if ( ibn == 10 .and. itrace == 0 ) return

!  ---  if h2otbl, compute h2o (line) transmittance using table look-up.
!       if conbnd, compute h2o (continuum) transmittance in bands 2-7.
!       if co2bnd, compute co2 transmittance in band 3.
!       if oznbnd, compute  o3 transmittance in band 5.
!       if n2obnd, compute n2o transmittance in bands 6 and 7.
!       if ch4bnd, compute ch4 transmittance in bands 6 and 7.
!       if combnd, compute co2-minor transmittance in bands 4 and 5.
!       if f11bnd, compute cfc11 transmittance in bands 4 and 5.
!       if f12bnd, compute cfc12 transmittance in bands 4 and 6.
!       if f22bnd, compute cfc22 transmittance in bands 4 and 6.
!       if b10bnd, compute flux reduction due to n2o in band 10.

        h2otbl = (itrans == 2) .and. (ibn==1 .or. ibn==2 .or. ibn==8)
        conbnd = (ibn >= 2) .and. (ibn <= 7)
        co2bnd = (ibn == 3)
        oznbnd = (ibn == 5)
        n2obnd = (ibn == 6) .or. (ibn == 7)
        ch4bnd = (ibn == 6) .or. (ibn == 7)
        combnd = (ibn == 4) .or. (ibn == 5)
        f11bnd = (ibn == 4) .or. (ibn == 5)
        f12bnd = (ibn == 4) .or. (ibn == 6)
        f22bnd = (ibn == 4) .or. (ibn == 6)
        b10bnd = (ibn == 10)

!  ---  blayer is the spectrally integrated planck flux of the mean layer
!       temperature derived from eq. (3.11)
!       The fitting for the planck flux is valid for the range 160-345 K.

        do k = 1, NLAY
          do i = 1, NPTS
            tx(i) = ta(i,k)
          enddo

          call planck                                                   &
!  ---  inputs:
     &     ( ibn, NPTS, tx,                                             &
!  ---  outputs:
     &       xlayer                                                     &
     &     )

          do i = 1, NPTS
            blayer(i,k) = xlayer(i)
          enddo
        enddo

!  ---  Index "0" is the layer above the top of the atmosphere.

        do i = 1, NPTS
          blayer(i,0) = _zero
        enddo

!  ---  Surface emission and reflectivity. See Section 9.

        call sfcflux                                                    &
!  ---  inputs:
     &     ( tg, eg, ibn, NPTS,                                         &
!  ---  outputs:
     &       bs, dbs, rflxs                                             &
     &     )

        do i = 1, NPTS
          blayer(i,NLP1) = bs(i)
        enddo

!  ----  interpolate Planck function at model levels (linear in p)

        do k = 2, NLAY
          do i = 1, NPTS
            blevel(i,k) = (blayer(i,k-1)*dp(i,k)+blayer(i,k)*dp(i,k-1)) &
     &                  / (dp(i,k-1) + dp(i,k))
          enddo
        enddo

!  ---  extrapolate blevel(i,1) from blayer(i,2) and blayer(i,1)

        do i = 1, NPTS
          blevel(i,1) = blayer(i,1) + (blayer(i,1) - blayer(i,2))       &
     &                * dp(i,1) / (dp(i,1) + dp(i,2))
        enddo

!  ---  note: we assume tb = tg to avoid the unsafe extrapolation

!ltb    if ( ltb ) then

!  ---  if the surface air temperature tb is known, compute blevel(i,np+1)

          call planck                                                   &
!  ---  inputs:
     &     ( ibn, NPTS, tb,                                             &
!  ---  outputs:
     &       xlayer                                                     &
     &     )

          do i = 1, NPTS
            blevel(i,NLP1) = xlayer(i)
          enddo

!       else

!  ---  otherwise, extrapolate blevel(np+1) from blayer(np-1) and blayer(np)

!         do i = 1, NPTS
!           blevel(i,NLP1) = blayer(i,NLAY)                             &
!    &                     + (blayer(i,NLAY) - blayer(i,NLAY-1))        &
!    &                     * dp(i,NLAY) / (dp(i,NLAY) + dp(i,NLAY-1))
!         enddo

!ltb    endif

!  ---  compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
!       rain optical thickness is set to 0.00307 /(gm/m**2).
!       it is for a specific drop size distribution provided by Q. Fu.

        if ( iflagc == 1 ) then

          do k = 1, NLAY
            do i = 1, NPTS
              if ( taucl(i,k,1) > 0.02 .and. fcld(i,k) > 0.01 ) then
                tcldlyr(i,k) = exp( -fdiff*taucl(i,k,1) )
              else
                tcldlyr(i,k) = _one
              endif
            enddo
          enddo

        else

          do k = 1, NLAY
            do i = 1, NPTS
              taucl(i,k,1) = cwp(i,k,1) * (aib(1,ibn)+aib(2,ibn)        &
     &                     / reff(i,k,1)**aib(3,ibn))
              taucl(i,k,2) = cwp(i,k,2) * (awb(1,ibn)+(awb(2,ibn)       &
     &                     + (awb(3,ibn)+awb(4,ibn)*reff(i,k,2))        &
     &                     * reff(i,k,2)) * reff(i,k,2))
              taucl(i,k,3) = 0.00307 * cwp(i,k,3)
            enddo
          enddo

!  ---  compute cloud single-scattering albedo and asymmetry factor for
!       a mixture of ice particles and liquid drops following 
!       eqs. (6.5), (6.6), (6.9) and (6.10).
!       single-scattering albedo and asymmetry factor of rain are set
!       to 0.54 and 0.95, respectively, based on the information provided
!       by Prof. Qiang Fu.

          do k = 1, NLAY
          do i = 1, NPTS
            tcldlyr(i,k) = _one
            tauc = taucl(i,k,1) + taucl(i,k,2) + taucl(i,k,3)

            if ( tauc > 0.02 .and. fcld(i,k) > 0.01 ) then
              reff1 = min( reff(i,k,1), 130.0 )
              reff2 = min( reff(i,k,2),  20.0 )
              reff1 = max( reff1, 20.0 )
              reff2 = max( reff2,  4.0 )

              w1 = taucl(i,k,1)*(aiw(1,ibn)+(aiw(2,ibn)+(aiw(3,ibn)     &
     &           + aiw(4,ibn)*reff1)*reff1)*reff1)
              w2 = taucl(i,k,2)*(aww(1,ibn)+(aww(2,ibn)+(aww(3,ibn)     &
     &           + aww(4,ibn)*reff2)*reff2)*reff2)
              w3 = taucl(i,k,3)*0.54
              ww = (w1+w2+w3) / tauc

              g1 = w1 * (aig(1,ibn)+(aig(2,ibn)+(aig(3,ibn)             &
     &           + aig(4,ibn)*reff1)*reff1)*reff1)
              g2 = w2 * (awg(1,ibn)+(awg(2,ibn)+(awg(3,ibn)             &
     &           + awg(4,ibn)*reff2)*reff2)*reff2)
              g3 = w3*0.95

              gg = (g1+g2+g3) / (w1+w2+w3)

!  ---  parameterization of LW scattering following Eqs. (6.11) and (6.12). 

              ff = 0.5 + (0.3739 + (0.0076+0.1185*gg)*gg ) * gg
              tauc = (_one - ww*ff) * tauc

!  ---  compute cloud diffuse transmittance. It is approximated by using 
!       a diffusivity factor of 1.66.

              tcldlyr(i,k) = exp( -fdiff*tauc )
            endif
          enddo
          enddo

        endif   ! end if_iflagc_block

!  ---  parameterization of aerosol scattering following Eqs. (6.11) and (6.12). 
!       taerlyr is the aerosol diffuse transmittance

        if ( iaerlw > 0 ) then
          taerlyr(:,:) = _one

          if ( iflip == 0 ) then        ! from toa to sfc
            do k = 1, NLAY
              do i = 1, NPTS
                taua(i,k) = aerosols(i,k,ibn,1)
                ssaa(i,k) = aerosols(i,k,ibn,2)
                asya(i,k) = aerosols(i,k,ibn,3)
              enddo
            enddo
          else                          ! from sfc to toa
            do k = 1, NLAY
              k1 = NLP1 - k
              do i = 1, NPTS
                taua(i,k) = aerosols(i,k1,ibn,1)
                ssaa(i,k) = aerosols(i,k1,ibn,2)
                asya(i,k) = aerosols(i,k1,ibn,3)
              enddo
            enddo
          endif

          do k = 1, NLAY
            do i = 1, NPTS
              if ( taua(i,k) > 0.001 ) then
                ff = 0.5 + (0.3739 + (0.0076                            &
     &             + 0.1185*asya(i,k)) * asya(i,k)) * asya(i,k)
                taua(i,k) = taua(i,k) * (_one-ssaa(i,k)*ff)
                taerlyr(i,k) = exp( -fdiff*taua(i,k) )
              endif
            enddo

!  ---  check print
!           if ( k >= NLAY-10 ) then
!           if ( ibn == 1 ) then
!             print *,'  In RADLW, k =',k
!             print *,'  TAUA_after:',taua(:,k)
!             print *,'  SSAA:',ssaa(:,k)
!             print *,'  ASYA:',asya(:,k)
!             print *,'  TAERLYR:',taerlyr(:,k)
!           endif

          enddo

        else

          taua(:,:) = _zero
          ssaa(:,:) = _zero
          asya(:,:) = _zero
          taerlyr(:,:) = _one

        endif

!  ---  compute the exponential terms (eq. 8.21) at each layer due to
!       water vapor line absorption when k-distribution is used

        if ( .not.h2otbl .and. .not.b10bnd ) then

          call h2oexps                                                  &
!  ---  inputs
     &     ( dh2o, pa, dt, xkw, aw, bw, pm, mw,                         &
     &       NPTS, NLAY, ibn,                                           &
!  ---  outputs:
     &       h2oexp                                                     &
     &     )

        endif

!  ---  compute the exponential terms (eq. 4.24) at each layer due to
!       water vapor continuum absorption.
!       ne is the number of terms used in each band to compute water 
!       vapor continuum transmittance (table 9).

        ne = 0

        if ( conbnd ) then
          ne = 1

          if ( ibn == 3 ) ne = 3

          call conexps                                                  &
!  ---  inputs:
     &     ( dcont, xke, NPTS, NLAY, ibn,                               &
!  ---  outputs:
     &       conexp                                                     &
     &     )

        endif

!  ---  compute the exponential terms (Eq. 8.21) at each layer due to
!       co2 absorption

        if ( itrans==1 .and. co2bnd ) then

          call co2exps                                                  &
!  ---  inputs:
     &     ( dco2, pa, dt, NPTS, NLAY,                                  &
!  ---  outputs:
     &       co2exp                                                     &
     &     )

        endif

!  ---  for trace gases

        if ( itrace == 1 ) then

!  ---  compute the exponential terms at each layer due to n2o absorption

          if ( n2obnd ) then

            call n2oexps                                                &
!  ---  inputs:
     &     ( dn2o, pa, dt, NPTS, NLAY, ibn,                             &
!  ---  outputs:
     &       n2oexp                                                     &
     &     )

          endif

!  ---  compute the exponential terms at each layer due to ch4 absorption

          if ( ch4bnd ) then

            call ch4exps                                                &
!  ---  inputs:
     &     ( dch4, pa, dt, NPTS, NLAY, ibn,                             &
!  ---  outputs:
     &       ch4exp                                                     &
     &     )

          endif

!  ---  compute the exponential terms due to co2 minor absorption

          if ( combnd ) then

            call comexps                                                &
!  ---  inputs:
     &     ( dco2, dt, NPTS, NLAY, ibn,                                 &
!  ---  outputs:
     &       comexp                                                     &
     &     )

          endif

!  ---  compute the exponential terms due to cfc11 absorption.
!       the values of the parameters are given in Table 7.

          if ( f11bnd ) then
            a1  = 1.26610e-3
            b1  = 3.55940e-6
            fk1 = 1.89736e+1
            a2  = 8.19370e-4
            b2  = 4.67810e-6
            fk2 = 1.01487e+1

            call cfcexps                                                &
!  ---  inputs:
     &     ( a1, b1, fk1, a2, b2, fk2, df11, dt,                        &
     &       NPTS, NLAY, ibn,                                           &
!  ---  outputs:
     &       f11exp                                                     &
     &     )

          endif

!  ---  compute the exponential terms due to cfc12 absorption.

          if ( f12bnd ) then
            a1  = 8.77370e-4
            b1  =-5.88440e-6
            fk1 = 1.58104e+1
            a2  = 8.62000e-4
            b2  =-4.22500e-6
            fk2 = 3.70107e+1

            call cfcexps                                                &
!  ---  inputs:
     &     ( a1, b1, fk1, a2, b2, fk2, df12, dt,                        &
     &       NPTS, NLAY, ibn,                                           &
!  ---  outputs:
     &       f12exp                                                     &
     &     )

          endif

!  ---  compute the exponential terms due to cfc22 absorption.

          if ( f22bnd ) then
            a1  = 9.65130e-4
            b1  = 1.31280e-5
            fk1 = 6.18536e+0
            a2  =-3.00010e-5 
            b2  = 5.25010e-7
            fk2 = 3.27912e+1

            call cfcexps                                                &
!  ---  inputs:
     &     ( a1, b1, fk1, a2, b2, fk2, df22, dt,                        &
     &       NPTS, NLAY, ibn,                                           &
!  ---  outputs:
     &       f22exp                                                     &
     &     )

          endif

!  ---  compute the exponential terms at each layer in band 10 due to
!       h2o line and continuum, co2, and n2o absorption

          if ( b10bnd ) then

            call b10exps                                                &
!  ---  inputs:
     &     ( dh2o, dcont, dco2, dn2o, pa, dt,                           &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       h2oexp, conexp, co2exp, n2oexp                             &
     &     )

          endif
        endif

!  ---  initialize fluxes

        do k = 1, NLP1
          do i = 1, NPTS
            flxu(i,k) = _zero
            flxd(i,k) = _zero
            flcu(i,k) = _zero
            flcd(i,k) = _zero
          enddo
        enddo

!  ---  for a given level, k1, compute the transmittance between this level
!       all the levels below, trant(i,k2).
!       also, compute the upward and doward blackbody emissions of a layer,
!       bu and bd.

        do i = 1, NPTS
          bd(i,0) = _zero
          bu(i,NLP1) = blayer(i,NLP1)
        enddo

        Lab_do_k1 : do k1 = 1, NLAY

!  ---  initialization
!       it, im, and ib are the numbers of cloudy layers in the high,
!       middle, and low cloud groups between levels k1 and k2.
!       cldlw, cldmd, and cldhi are the equivalent black-cloud fractions
!       of low, middle, and high troposphere.
!       tranal is the aerosol transmission function

          do i = 1, NPTS
            it(i) = 0
            im(i) = 0
            ib(i) = 0

            cldlw(i) = _zero
            cldmd(i) = _zero
            cldhi(i) = _zero
            tranal(i)= _one
          enddo

!  ---  for h2o line transmission

          if ( .not. h2otbl ) then
            do m = 1, NK0
              do i = 1, NPTS
                th2o(i,m) = _one
              enddo
            enddo
          endif

!  ---  for h2o continuum transmission

          do m = 1, NK2
            do i = 1, NPTS
              tcon(i,m) = _one
            enddo
          enddo

!  ---  for co2 transmission using k-distribution method.
!       band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
!       are combined in computing the co2 transmittance.

          if ( itrans==1 .and. co2bnd ) then
            do j = 1, 2
              do m = 1, NK0
                do i = 1, NPTS
                  tco2(i,m,j) = _one
                enddo
              enddo
            enddo
          endif

!  ---  for trace gases

          if ( itrace == 1 ) then

!  ---  for n2o transmission using k-distribution method.

            if ( n2obnd ) then
              do m = 1, NK1
                do i = 1, NPTS
                  tn2o(i,m) = _one
                enddo
              enddo
            endif

!  ---  for ch4 transmission using k-distribution method.

            if ( ch4bnd ) then
              do m = 1, NK1
                do i = 1, NPTS
                  tch4(i,m) = _one
                enddo
              enddo
            endif

!  ---  for co2-minor transmission using k-distribution method.

            if ( combnd ) then
              do m = 1, NK0
                do i = 1, NPTS
                  tcom(i,m) = _one
                enddo
              enddo
            endif

!  ---  for cfc-11 transmission using k-distribution method.

            if ( f11bnd ) then
              do i = 1, NPTS
                tf11(i) = _one
              enddo
            endif

!  ---  for cfc-12 transmission using k-distribution method.

            if ( f12bnd ) then
              do i = 1, NPTS
                tf12(i) = _one
              enddo
            endif

!  ---  for cfc-22 transmission when using k-distribution method.

            if ( f22bnd ) then
              do i = 1, NPTS
                tf22(i) = _one
              enddo
            endif

!  ---  for the transmission in band 10 using k-distribution method.

            if ( b10bnd ) then
              do m = 1, NK0-1
                do i = 1, NPTS
                  th2o(i,m) = _one
                enddo
              enddo

              do m = 1, NK0
                do i = 1, NPTS
                  tco2(i,m,1) = _one
                enddo
              enddo

              do i = 1, NPTS
                tcon(i,1) = _one
              enddo

              do m = 1, 2
                do i = 1, NPTS
                  tn2o(i,m) = _one
                enddo
              enddo
            endif
          endif

!  ---  end trace gases

          do i = 1, NPTS
            x1(i) = _zero
            x2(i) = _zero
            x3(i) = _zero
          enddo

          do k = 1, NLP1
            do i = 1, NPTS
              fclr(i,k) = _one
            enddo
          enddo

!  ---  loop over the bottom level of the region (k2)

          Lab_do_k2 : do k2 = k1+1, NLP1

!  ---  trant is the total transmittance between levels k1 and k2.

            do i = 1, NPTS
              trant(i,k2) = _one
            enddo

            if ( h2otbl ) then

!  ---  compute water vapor transmittance using table look-up.
!       the following values are taken from table 8.

              w1  = -8.0
              p1  = -2.0
              dwe = 0.3
              dpe = 0.2

              if ( ibn == 1 ) then

                call tablup                                             &
!  ---  inputs:
     &     ( dh2o, pa, dt, w1, p1, dwe, dpe, h11, h12, h13,             &
     &       NPTS, k2, nx, nh,                                          &
!  ---  in/outputs:
     &       x1, x2, x3, trant                                          &
     &     )

              endif

              if ( ibn == 2 ) then

                call tablup                                             &
!  ---  inputs:
     &     ( dh2o, pa, dt, w1, p1, dwe, dpe, h21, h22, h23,             &
     &       NPTS, k2, nx, nh,                                          &
!  ---  in/outputs:
     &       x1, x2, x3, trant                                          &
     &     )

              endif

              if ( ibn == 8 ) then

                call tablup                                             &
!  ---  inputs:
     &     ( dh2o, pa, dt, w1, p1, dwe, dpe, h81, h82, h83,             &
     &       NPTS, k2, nx, nh,                                          &
!  ---  in/outputs:
     &       x1, x2, x3, trant                                          &
     &     )

              endif

              if ( conbnd ) then
                do i = 1, NPTS
                  tcon (i,1) = tcon (i,1) *conexp(i,k2-1,1)
                  trant(i,k2)= trant(i,k2)*tcon(i,1)
                enddo
              endif

            else

!  --- compute water vapor transmittance using k-distribution

              if ( .not.b10bnd ) then

                call h2okdis                                            &
!  ---  inputs:
     &     ( fkw, gkw, h2oexp, conexp,                                  &
     &       NPTS, k2-1, ibn, ne,                                       &
!  ---  in/outputs:
     &       th2o, tcon, trant                                          &
     &     )

              endif

            endif

            if ( co2bnd ) then

              if ( itrans == 2 ) then

!  ---  compute co2 transmittance using table look-up method.
!       the following values are taken from table 8.

                w1  = -4.0
                p1  = -2.0
                dwe = 0.3
                dpe = 0.2

                call tablup                                             &
!  ---  inputs:
     &     ( dco2, pa, dt, w1, p1, dwe, dpe, c1, c2, c3,                &
     &       NPTS, k2, nx, nc,                                          &
!  ---  in/outputs:
     &       x1, x2, x3, trant                                          &
     &     )

              else

!  ---  compute co2 transmittance using k-distribution method

                call co2kdis                                            &
!  ---  inputs:
     &     ( co2exp, NPTS, k2-1,                                        &
!  ---  in/outputs:
     &       tco2, trant                                                &
     &     )

              endif

            endif

!  ---  Always use table look-up to compute o3 transmittance.
!       the following values are taken from table 8.

            if ( oznbnd ) then
              w1  = -6.0
              p1  = -2.0
              dwe = 0.3
              dpe = 0.2

              call tablup                                               &
!  ---  inputs:
     &     ( do3, pa, dt, w1, p1, dwe, dpe, o1, o2, o3,                 &
     &       NPTS, k2, nx, no,                                          &
!  ---  in/outputs:
     &       x1, x2, x3, trant                                          &
     &     )

            endif

!  ---  for trace gases

            if ( itrace == 1 ) then

!  ---  compute n2o transmittance using k-distribution method

              if ( n2obnd ) then

                call n2okdis                                            &
!  ---  inputs:
     &     ( n2oexp, NPTS, k2-1, ibn,                                   &
!  ---  in/outputs:
     &       tn2o, trant                                                &
     &     )

              endif

!  ---  compute ch4 transmittance using k-distribution method

              if ( ch4bnd ) then

                call ch4kdis                                            &
     &     ( ch4exp, NPTS, k2-1, ibn,                                   &
!  ---  in/outputs:
     &       tch4, trant                                                &
     &     )

              endif

!  ---  compute co2-minor transmittance using k-distribution method

              if ( combnd ) then

                call comkdis                                            &
!  ---  inputs:
     &     ( comexp, NPTS, k2-1, ibn,                                   &
!  ---  in/outputs:
     &       tcom, trant                                                &
     &     )

              endif

!  ---  compute cfc11 transmittance using k-distribution method

              if ( f11bnd ) then

                call cfckdis                                            &
!  ---  inputs:
     &     ( f11exp, NPTS, k2-1,                                        &
!  ---  in/outputs:
     &       tf11, trant                                                &
     &     )

              endif

!  ---  compute cfc12 transmittance using k-distribution method

              if ( f12bnd ) then

                call cfckdis                                            &
!  ---  inputs:
     &     ( f12exp, NPTS, k2-1,                                        &
!  ---  in/outputs:
     &       tf12, trant                                                &
     &     )

              endif

!  ---  compute cfc22 transmittance using k-distribution method

              if ( f22bnd ) then

                call cfckdis                                            &
!  ---  inputs:
     &     ( f22exp, NPTS, k2-1,                                        &
!  ---  in/outputs:
     &       tf22, trant                                                &
     &     )
              endif

!  ---  compute transmittance in band 10 using k-distribution method.
!       for band 10, trant is the change in transmittance due to n2o 
!       absorption.

              if ( b10bnd ) then

                call b10kdis                                            &
!  ---  inputs:
     &     ( h2oexp, conexp, co2exp, n2oexp,                            &
     &       NPTS, k2-1,                                                &
!  ---  outputs:
     &       th2o, tcon, tco2, tn2o, trant                              &
     &     )

              endif
            endif

!  ---  end trace gases  *****

!  ---  include aerosol effect

            if ( iaerlw > 0 ) then
              do i = 1, NPTS
                tranal(i) = tranal(i) * taerlyr(i,k2-1)
                trant(i,k2) = trant(i,k2) * tranal(i)
              enddo
            endif

!  ---  cloud overlapping *****

            if ( iovcst == 2 ) then

              call cldovlp                                              &
!  ---  inputs:
     &      ( fcld, tcldlyr, ict, icb,                                  &
     &        NPTS, NLAY, k2,                                           &
!  ---  in/outputs:
     &        cldhi, cldmd, cldlw, it, im, ib,                          &
     &        fclr                                                      &
     &      )

            else

              do i = 1, NPTS
                fclr(i,k2) = fclr(i,k2) * tcldlyr(i,k2-1)
              enddo

            endif

!  ---  compute upward and downward blackbody emission of a layer

            if ( k1 == 1 ) then
              do i = 1, NPTS

                xx = (blayer(i,k2-1) - blevel(i,k2-1))                  &
     &             * (blayer(i,k2-1) - blevel(i,k2))

                if ( xx > _zero ) then

!  ---  if xx>0, there is a local temperature minimum or maximum.
!       computations of bd and bu follow eq. (8.20).

                  bd(i,k2-1) = 0.5 * blayer(i,k2-1)                     &
     &                       + .25 * (blevel(i,k2-1) + blevel(i,k2))
                  bu(i,k2-1) = bd(i,k2-1)

                else

!  ---  computations of bd and bu following eqs.(8.17) and (8.18).
!       the effect of clouds on the transmission of a layer is taken
!       into account, following eq. (8.19).

                  xx = (fcld(i,k2-1)*tcldlyr(i,k2-1)                    &
     &               + (_one-fcld(i,k2-1))) * trant(i,k2)
                  yy = min( 0.9999, xx )
                  yy = max( 0.00001, yy )
                  xx = (blevel(i,k2-1) - blevel(i,k2)) / alog(yy)
                  bd(i,k2-1) = (blevel(i,k2) - blevel(i,k2-1)*yy)       &
     &                       / (_one - yy) - xx
                  bu(i,k2-1) = (blevel(i,k2-1) + blevel(i,k2))          &
     &                       - bd(i,k2-1)

                endif

              enddo
            endif

          enddo  Lab_do_k2

!  ---  upward and downward flux calculations.

          do 4000 k2 = k1+1, NLP1

            if ( k2 == k1+1 .and. ibn /= 10 ) then

!  ---  the first terms on the rhs of eqs. (8.15) and (8.16)

              do i = 1, NPTS
                flcu(i,k1) = flcu(i,k1) - bu(i,k1)
                flcd(i,k2) = flcd(i,k2) + bd(i,k1)
                flxu(i,k1) = flxu(i,k1) - bu(i,k1)
                flxd(i,k2) = flxd(i,k2) + bd(i,k1)
              enddo

            endif

!  ---  the summation terms on the rhs of eqs. (8.15) and (8.16).
!       also see eqs. (5.4) and (5.5) for band 10.

            do i = 1, NPTS
              xx = trant(i,k2) * (bu(i,k2-1) - bu(i,k2))
              flcu(i,k1) = flcu(i,k1) + xx
              flxu(i,k1) = flxu(i,k1) + xx*fclr(i,k2)

              xx = trant(i,k2) * (bd(i,k1-1) - bd(i,k1))
              flcd(i,k2) = flcd(i,k2) + xx
              flxd(i,k2) = flxd(i,k2) + xx*fclr(i,k2)
            enddo

 4000     continue

!  ---  here, fclr and trant are, respectively, the clear line-of-sight 
!       and the transmittance between k1 and the surface.

          do i = 1, NPTS
            trantcr(i,k1) = trant(i,NLP1)
            transfc(i,k1) = trant(i,NLP1)*fclr(i,NLP1)
          enddo

!  ---  compute the partial derivative of fluxes with respect to
!       surface temperature (eq. 3.12). 
!       note: upward flux is negative, and so is dfdts.

!         do i = 1, NPTS
!           dfdts(i,k1) = dfdts(i,k1) - dbs(i)*transfc(i,k1)
!         enddo

        enddo  Lab_do_k1

        if ( .not. b10bnd ) then

!  ---  for surface emission.
!       note: blayer(i,NLP1) and dbs include the surface emissivity effect.
!       both dfdts and sfcem are negative quanties.

          do i = 1, NPTS
            flcu(i,NLP1) = -blayer(i,NLP1)
            flxu(i,NLP1) = -blayer(i,NLP1)
!           sfcem(i)     = sfcem(i) - blayer(i,NLP1)
!           dfdts(i,NLP1) = dfdts(i,NLP1) - dbs(i)
          enddo

!  ---  add the flux reflected by the surface. (second term on the
!       rhs of eq. 8.16)

          do k = 1, NLP1
            do i = 1, NPTS
              flcu(i,k) = flcu(i,k)                                     &
     &                  - flcd(i,NLP1)*trantcr(i,k)*rflxs(i)
              flxu(i,k) = flxu(i,k)                                     &
     &                  - flxd(i,NLP1)*transfc(i,k)*rflxs(i)   
            enddo
          enddo

        endif

        if ( ibn == 10 .and. itrace == 0 ) then

          do k = 1, NLP1
            do i = 1, NPTS
              flcu(i,k) = _zero
              flcd(i,k) = _zero
              flxu(i,k) = _zero
              flxd(i,k) = _zero
            enddo
          enddo

        endif

!  ---  summation of fluxes over spectral bands

        do k = 1, NLP1
          do i = 1, NPTS
            flxdn0(i,k) = flxdn0(i,k) + flcd(i,k)
            flxup0(i,k) = flxup0(i,k) - flcu(i,k)

            flxdnc(i,k) = flxdnc(i,k) + flxd(i,k)
            flxupc(i,k) = flxupc(i,k) - flxu(i,k)
          enddo
        enddo

!! ---  optional spectral band heating

        if ( lhlwb ) then
          do k = 1, NLP1
          do i = 1, NPTS
            flxntb(i,k,ibn) = - flxu(i,k) - flxd(i,k)
          enddo
          enddo
        endif

      enddo  Lab_do_ibn

!  ---  prepare for final outputs

      do k = 1, NLP1
        do i = 1, NPTS
          flxntc(i,k) = flxupc(i,k) - flxdnc(i,k)
        enddo
      enddo

!! ---  optional clear sky net flux for heating rate calculations
      if ( lhlw0 ) then
        do k = 1, NLP1
        do i = 1, NPTS
          flxnt0(i,k) = flxup0(i,k) - flxdn0(i,k)
        enddo
        enddo
      endif

!  ---  output total-sky and clear-sky fluxesa and heating rates

      do i = 1, NPTS
        topflx(i)%upfxc = flxupc(i,1)
        topflx(i)%upfx0 = flxup0(i,1)

        sfcflx(i)%upfxc = flxupc(i,NLP1)
        sfcflx(i)%dnfxc = flxdnc(i,NLP1)
        sfcflx(i)%dnfx0 = flxdn0(i,NLP1)
      enddo

      if (iflip == 0) then        ! output from toa to sfc

        do k = 1, NLAY
          do i = 1, NPTS
            hlwc(i,k) = (flxntc(i,k+1) - flxntc(i,k))*heatfac           &
     &                / dp(i,k)
          enddo
        enddo

!! ---  optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
          do i = 1, NPTS
            hlw0(i,k) = (flxnt0(i,k+1) - flxnt0(i,k))*heatfac           &
     &                / dp(i,k)
          enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhlwb ) then
          do ibn = 1, NBDLW
          do k = 1, NLAY
          do i = 1, NPTS
            hlwb(i,k,ibn) = (flxntb(i,k+1,ibn) - flxntb(i,k,ibn))       &
     &                    * heatfac / dp(i,k)
          enddo
          enddo
          enddo
        endif

!! ---  optional fluxes
        if ( lflxprf ) then
          do k = 1, NLP1
          do i = 1, NPTS
            flxprf(i,k)%upfxc = flxupc(i,k)
            flxprf(i,k)%dnfxc = flxdnc(i,k)
            flxprf(i,k)%upfx0 = flxup0(i,k)
            flxprf(i,k)%dnfx0 = flxdn0(i,k)
          enddo
          enddo
        endif

      else                        ! output from sfc to toa

        do k = 1, NLAY
          k1 = NLP1 - k
          do i = 1, NPTS
            hlwc(i,k1) = (flxntc(i,k+1) - flxntc(i,k))*heatfac          &
     &                 / dp(i,k)
          enddo
        enddo

!! ---  optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              hlw0(i,k1) = (flxnt0(i,k+1) - flxnt0(i,k))*heatfac        &
     &                   / dp(i,k)
            enddo
          enddo
        endif

!! ---  optional spectral band heating
        if ( lhlwb ) then
          do ibn = 1, NBDLW
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, NPTS
              hlwb(i,k1,ibn) = (flxntb(i,k+1,ibn) - flxntb(i,k,ibn))    &
     &                       * heatfac / dp(i,k)
            enddo
          enddo
          enddo
        endif

!! ---  optional fluxes
        if ( lflxprf ) then
          do k = 1, NLP1
            k1 = NLP1 - k + 1
            do i = 1, NPTS
              flxprf(i,k1)%upfxc = flxupc(i,k)
              flxprf(i,k1)%dnfxc = flxdnc(i,k)
              flxprf(i,k1)%upfx0 = flxup0(i,k)
              flxprf(i,k1)%dnfx0 = flxdn0(i,k)
            enddo
          enddo
        endif

      endif                       ! end if_iflip_block
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
!  ---  outputs: (none)

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
!     iaerlw  - control flag for aerosols (not yet)                     !
!               =0: do not include aerosol effect                       !
!               >0: include aerosol effect                              !
!     itrace  - control flag for trace gases (ch4,n2o,cfcs, etc.)       !
!               =0: do not include trace gases                          !
!               =1: include all trace gases                             !
!     iflagc  - cloud optical properties contrl flag                    !
!               =1: input cloud opt depth from diagnostic scheme        !
!               =2: input cwp,cip, and other cloud content parameters   !
!     iovcst  - cloud fraction assumption contrl flag                   !
!               =1: overcast assumption, clouds either 0 or 1           !
!               =2: allowing for fractional clouds                      !
!                                                                       !
!  *******************************************************************  !
!
      implicit none
!
!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: none

!  ---  locals:  none

!
!===> ... begin here
!

      if (me == 0) then
        print *,' - Using NASA Longwave Radiation, Version: ', VTAGLW

        if (iaerlw > 0) then
          print *,'   --- Using input aerosol parameters for LW'
        else
          print *,'   --- Aerosol effect is NOT included in LW, all'    &
     &           ,' internal aerosol parameters are reset to zeros'
        endif

        if (itrace == 1) then
          print *,'   --- Include trace gases N2O, CH4, CFCs,',         &
     &            ' absorptions in LW'
        else
          print *,'   --- Trace gases effect is NOT included in LW'
        endif
      endif

!  --- ...  check cloud flags for consistency

      if ((icwp == 0 .and. iflagc /= 1) .or.                            &
     &    (icwp == 1 .and. iflagc == 1)) then
        print *, ' *** Model cloud scheme inconsistent with LW',        &
     &           ' radiation cloud radiative property setup !!'
        stop
      endif

!  --- ...  setup constant factor and heating rate
!           the 1.0e-2 is to convert pressure from mb to N/m**2

      if (ilwrate == 1) then
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
!       heatfac = con_g * 864.0 / con_cp            !   (in k/day)
        heatfac = 8.441874                          !   (in k/day)
      else
!       heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
        heatfac = 8.441874 / 86400.0                !   (in k/second)
      endif

!
      return
!...................................
      end subroutine rlwinit
!-----------------------------------



!-----------------------------------
      subroutine planck                                                 &
!...................................

!  ---  inputs:
     &     ( ibn, NPTS, t,                                              &
!  ---  outputs:
     &       xlayer                                                     &
     &     )

!-----Compute spectrally integrated Planck flux
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: ibn                     ! spectral band index
      integer, intent(in) :: NPTS                    ! no of points
      real (kind=kind_phys), intent(in) :: t(:)      ! temperature (K)

!  ---  outputs:
      real (kind=kind_phys), intent(out):: xlayer(:) ! planck flux (w/m2)

!  ---  locals:
      integer :: i

      real (kind=kind_phys) :: cb(6,NBDLW)

!  ---  the following coefficients are given in table 2 for computing  
!       spectrally integrated planck fluxes using eq. (3.11)

      data   cb   /                                                     &
     &  5.3443e+0,-2.0617e-1,2.5333e-3,-6.8633e-6,1.0115e-8,-6.2672e-12,&
     &  2.7148e+1,-5.4038e-1,2.9501e-3,2.7228e-7,-9.3384e-9, 9.9677e-12,&
     & -3.4860e+1,1.1132e+0,-1.3006e-2,6.4955e-5,-1.1815e-7, 8.0424e-11,&
     & -6.0513e+1,1.4087e+0,-1.2077e-2,4.4050e-5,-5.6735e-8, 2.5660e-11,&
     & -2.6689e+1,5.2828e-1,-3.4453e-3,6.0715e-6, 1.2523e-8,-2.1550e-11,&
     & -6.7274e+0,4.2256e-2, 1.0441e-3,-1.2917e-5,4.7396e-8,-4.4855e-11,&
     &  1.8786e+1,-5.8359e-1,6.9674e-3,-3.9391e-5,1.0120e-7,-8.2301e-11,&
     &  1.0344e+2,-2.5134e+0,2.3748e-2,-1.0692e-4,2.1841e-7,-1.3704e-10,&
     & -1.0482e+1,3.8213e-1,-5.2267e-3,3.4412e-5,-1.1075e-7, 1.4092e-10,&
     &  1.6769e+0,6.5397e-2,-1.8125e-3,1.2912e-5,-2.6715e-8, 1.9792e-11/
!
!===> ...  begin here
!
      do i = 1, NPTS
        xlayer(i) = t(i) * (t(i) * (t(i) * (t(i) * (t(i)*cb(6,ibn)      &
     &            + cb(5,ibn)) + cb(4,ibn)) + cb(3,ibn)) + cb(2,ibn))   &
     &            + cb(1,ibn)
      enddo
!
      return
!...................................
      end subroutine planck
!-----------------------------------



!-----------------------------------
      subroutine plancd                                                 &
!...................................

!  ---  inputs:
     &     ( ibn, NPTS, t,                                              &
!  ---  outputs:
     &       dbdt                                                       &
     &     )

!  ---  compute the derivative of Planck flux wrt temperature
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: ibn                    ! spectral band index
      integer, intent(in) :: NPTS                   ! no of points
      real (kind=kind_phys), intent(in) :: t(:)     ! temperature (K)

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: dbdt(:) ! derivative of Planck flux
                                                    ! wrt temperature
!  ---  locals:
      integer i

      real (kind=kind_phys) :: dcb(5,NBDLW)

!  ---  coefficients for computing the derivative of planck function
!       with respect to temperature (eq. 3.12).
!       dcb(1)=1*cb(2), dcb(2)=2*cb(3), dcb(3)=3*cb(4) ...  etc

      data   dcb   /                                                    &
     &  -2.0617E-01, 5.0666E-03,-2.0590E-05, 4.0460E-08,-3.1336E-11,    &
     &  -5.4038E-01, 5.9002E-03, 8.1684E-07,-3.7354E-08, 4.9839E-11,    &
     &   1.1132E+00,-2.6012E-02, 1.9486E-04,-4.7260E-07, 4.0212E-10,    &
     &   1.4087E+00,-2.4154E-02, 1.3215E-04,-2.2694E-07, 1.2830E-10,    &
     &   5.2828E-01,-6.8906E-03, 1.8215E-05, 5.0092E-08,-1.0775E-10,    &
     &   4.2256E-02, 2.0882E-03,-3.8751E-05, 1.8958E-07,-2.2428E-10,    &
     &  -5.8359E-01, 1.3935E-02,-1.1817E-04, 4.0480E-07,-4.1150E-10,    &
     &  -2.5134E+00, 4.7496E-02,-3.2076E-04, 8.7364E-07,-6.8520E-10,    &
     &   3.8213E-01,-1.0453E-02, 1.0324E-04,-4.4300E-07, 7.0460E-10,    &
     &   6.5397E-02,-3.6250E-03, 3.8736E-05,-1.0686E-07, 9.8960E-11 /
!
!===> ...  begin here
!
      do i = 1, NPTS
        dbdt(i) = t(i) * (t(i) * (t(i) * (t(i)*dcb(5,ibn)+dcb(4,ibn))   &
     &          + dcb(3,ibn)) + dcb(2,ibn)) + dcb(1,ibn)
      enddo
!
      return
!...................................
      end subroutine plancd
!-----------------------------------



!-----------------------------------
      subroutine h2oexps                                                &
!...................................

!  ---  inputs
     &     ( dh2o, pa, dt, xkw, aw, bw, pm, mw,                         &
     &       NPTS, NLAY, ib,                                            &
!  ---  outputs:
     &       h2oexp                                                     &
     &     )

!**********************************************************************
!  compute exponentials for water vapor line absorption
!  in individual layers using eqs. (8.21) and (8.22).
!
!  --- input parameters
!  layer water vapor amount for line absorption (dh2o) 
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  absorption coefficients for the first k-distribution
!     function due to h2o line absorption (xkw)
!  coefficients for the temperature and pressure scaling (aw,bw,pm)
!  ratios between neighboring absorption coefficients for
!     h2o line absorption (mw)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  6 exponentials for each layer  (h2oexp)
!
!  note: that the 3 sub-bands in band 3 use the same set of xkw, aw,
!  and bw,  therefore, h2oexp for these sub-bands are identical.
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib, mw(:)

      real (kind=kind_phys), dimension(:,:), intent(in) :: dh2o, pa, dt
      real (kind=kind_phys), dimension(:),   intent(in) :: xkw,aw,bw,pm

!  --- outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: h2oexp

!  ---  locals:
      integer :: i, k, ik

      real (kind=kind_phys) :: xh
!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS

!  ---  xh is the scaled water vapor amount for line absorption
!       computed from eq. (4.4).

          xh = dh2o(i,k) * (pa(i,k)/500.0)**pm(ib)                      &
     &       * ( _one + (aw(ib) + bw(ib)*dt(i,k)) * dt(i,k) )

!  ---  h2oexp is the water vapor transmittance of the layer k
!       due to line absorption

          h2oexp(i,k,1) = exp( -xh*xkw(ib) )

        enddo
      enddo

!  ---  compute transmittances from eq. (8.22)

      do ik = 2, NK0

        if ( mw(ib) == 6 ) then

          do k = 1, NLAY
            do i = 1, NPTS
              xh = h2oexp(i,k,ik-1) * h2oexp(i,k,ik-1)
              h2oexp(i,k,ik) = xh*xh*xh
            enddo
          enddo

        elseif ( mw(ib) == 8 ) then

          do k = 1, NLAY
            do i = 1, NPTS
              xh = h2oexp(i,k,ik-1) * h2oexp(i,k,ik-1)
              xh = xh*xh
              h2oexp(i,k,ik) = xh*xh
            enddo
          enddo

        elseif ( mw(ib) == 9 ) then

          do k = 1, NLAY
            do i = 1, NPTS
              xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
              h2oexp(i,k,ik) = xh*xh*xh
            enddo
          enddo

        else

          do k = 1, NLAY
            do i = 1, NPTS
              xh = h2oexp(i,k,ik-1) * h2oexp(i,k,ik-1)
              xh = xh*xh
              xh = xh*xh
              h2oexp(i,k,ik) = xh*xh
            enddo
          enddo

        endif
       enddo
!
      return
!...................................
      end subroutine h2oexps
!-----------------------------------



!-----------------------------------
      subroutine conexps                                                &
!...................................

!  ---  inputs:
     &     ( dcont, xke, NPTS, NLAY, ib,                                &
!  ---  outputs:
     &       conexp                                                     &
     &     )

!**********************************************************************
!  compute exponentials for continuum absorption in individual layers.
!
!  --- input parameters
!  layer scaled water vapor amount for continuum absorption (dcont) 
!  absorption coefficients for the first k-distribution function
!     due to water vapor continuum absorption (xke)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  1 or 3 exponentials for each layer (conexp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib

      real (kind=kind_phys), intent(in) :: dcont(:,:), xke(:)

!  --- outputs:
      real (kind=kind_phys), intent(out) :: conexp(:,:,:)

!  ---  locals:
      integer :: i, k

!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS
          conexp(i,k,1) = exp( -dcont(i,k)*xke(ib) )
        enddo
      enddo

      if ( ib == 3 ) then

!  ---  the absorption coefficients for sub-bands 3b and 3a are, respectively,
!       two and four times the absorption coefficient for sub-band 3c (table 9).
!       note that conexp(i,k,3) is for sub-band 3a. 

        do k = 1, NLAY
          do i = 1, NPTS
            conexp(i,k,2) = conexp(i,k,1) * conexp(i,k,1)
            conexp(i,k,3) = conexp(i,k,2) * conexp(i,k,2)
          enddo
        enddo

      endif
!
      return
!...................................
      end subroutine conexps
!-----------------------------------



!-----------------------------------
      subroutine co2exps                                                &
!...................................

!  ---  inputs:
     &     ( dco2, pa, dt, NPTS, NLAY,                                  &
!  ---  outputs:
     &       co2exp                                                     &
     &     )

!**********************************************************************
!  compute co2 exponentials for individual layers.
!
!  --- input parameters
!  layer co2 amount (dco2)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!
!  --- output parameters
!  6 exponentials for each layer (co2exp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: dco2, pa, dt

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) :: co2exp

!  ---  locals:
      integer :: i, k

      real (kind=kind_phys) :: xc

!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS

!  ---  the scaling parameters are given in table 3, and values of
!       the absorption coefficient are given in table 10.
!       scaled co2 amount for band-wings (sub-bands 3a and 3c)

          xc = dco2(i,k) * (pa(i,k)/300.0)**0.5                         &
     &       * (_one + (0.0182 + 1.07e-4*dt(i,k)) * dt(i,k))

!  ---  six exponentials by powers of 8 (See eqs. 8.21, 8.22 and table 10).

          co2exp(i,k,1,1) = exp( -xc*2.656e-5 )

          xc = co2exp(i,k,1,1) * co2exp(i,k,1,1)
          xc = xc*xc
          co2exp(i,k,2,1) = xc*xc

          xc = co2exp(i,k,2,1) * co2exp(i,k,2,1)
          xc = xc*xc
          co2exp(i,k,3,1) = xc*xc

          xc = co2exp(i,k,3,1) * co2exp(i,k,3,1)
          xc = xc*xc
          co2exp(i,k,4,1) = xc*xc

          xc = co2exp(i,k,4,1) * co2exp(i,k,4,1)
          xc = xc*xc
          co2exp(i,k,5,1) = xc*xc

          xc = co2exp(i,k,5,1) * co2exp(i,k,5,1)
          xc = xc*xc
          co2exp(i,k,6,1) = xc*xc

!  ---  for band-center region (sub-band 3b)

          xc = dco2(i,k) * (pa(i,k)/30.0)**0.85                         &
     &       * (_one + (0.0042 + 2.00e-5*dt(i,k)) * dt(i,k))

          co2exp(i,k,1,2) = exp( -xc*2.656e-3 )

          xc = co2exp(i,k,1,2) * co2exp(i,k,1,2)
          xc = xc*xc
          co2exp(i,k,2,2) = xc*xc

          xc = co2exp(i,k,2,2) * co2exp(i,k,2,2)
          xc = xc*xc
          co2exp(i,k,3,2) = xc*xc

          xc = co2exp(i,k,3,2) * co2exp(i,k,3,2)
          xc = xc*xc
          co2exp(i,k,4,2) = xc*xc

          xc = co2exp(i,k,4,2) * co2exp(i,k,4,2)
          xc = xc*xc
          co2exp(i,k,5,2) = xc*xc

          xc = co2exp(i,k,5,2) * co2exp(i,k,5,2)
          xc = xc*xc
          co2exp(i,k,6,2) = xc*xc

        enddo
      enddo
!
      return
!...................................
      end subroutine co2exps
!-----------------------------------



!-----------------------------------
      subroutine n2oexps                                                &
!...................................

!  ---  inputs:
     &     ( dn2o, pa, dt, NPTS, NLAY, ib,                              &
!  ---  outputs:
     &       n2oexp                                                     &
     &     )

!**********************************************************************
!   compute n2o exponentials for individual layers 
!
!  --- input parameters
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  2 or 4 exponentials for each layer (n2oexp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib

      real (kind=kind_phys), dimension(:,:), intent(in) :: dn2o, pa, dt

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: n2oexp

!  ---  locals:
      integer :: i, k

      real (kind=kind_phys) :: xc, xc1, xc2
!
!===> ...  begin here
!
!  ---  scaling and absorption data are given in table 5.
!       transmittances are computed using eqs. (8.21) and (8.22).

      do k = 1, NLAY
        do i = 1, NPTS

!  ---  four exponential by powers of 21 for band 6.

          if ( ib == 6 ) then

            xc = dn2o(i,k) * (_one + (1.9297e-3                         &
     &         + 4.3750e-6*dt(i,k)) * dt(i,k))
            n2oexp(i,k,1) = exp( -xc*6.31582e-2 )

            xc = n2oexp(i,k,1) * n2oexp(i,k,1) * n2oexp(i,k,1)
            xc1 = xc*xc
            xc2 = xc1*xc1
            n2oexp(i,k,2) = xc*xc1*xc2

!  ---  four exponential by powers of 8 for band 7

          else

            xc = dn2o(i,k) * (pa(i,k)/500.0)**0.48                      &
     &         * (_one + (1.3804e-3 + 7.4838e-6*dt(i,k)) * dt(i,k))
            n2oexp(i,k,1) = exp( -xc*5.35779e-2 )

            xc = n2oexp(i,k,1) * n2oexp(i,k,1)
            xc = xc*xc
            n2oexp(i,k,2) = xc*xc

            xc = n2oexp(i,k,2) * n2oexp(i,k,2)
            xc = xc*xc
            n2oexp(i,k,3) = xc*xc

            xc = n2oexp(i,k,3) * n2oexp(i,k,3)
            xc = xc*xc
            n2oexp(i,k,4) = xc*xc

          endif

        enddo
       enddo
!
      return
!...................................
      end subroutine n2oexps
!-----------------------------------



!-----------------------------------
      subroutine ch4exps                                                &
!...................................

!  ---  inputs:
     &     ( dch4, pa, dt, NPTS, NLAY, ib,                              &
!  ---  outputs:
     &       ch4exp                                                     &
     &     )

!**********************************************************************
!  compute ch4 exponentials for individual layers
!
!  --- input parameters
!  layer ch4 amount (dch4)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  1 or 4 exponentials for each layer (ch4exp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib

      real (kind=kind_phys), dimension(:,:), intent(in) :: dch4, pa, dt

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: ch4exp

!  ---  locals:
      integer :: i, k

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  scaling and absorpton data are given in table 5

      do k = 1, NLAY
        do i = 1, NPTS

!  ---  four exponentials for band 6

          if ( ib == 6 ) then

            xc = dch4(i,k) * (_one + (1.7007e-2                         &
     &         + 1.5826e-4*dt(i,k)) * dt(i,k))
            ch4exp(i,k,1) = exp( -xc*5.80708e-3 )

!  ---  four exponentials by powers of 12 for band 7

          else

            xc = dch4(i,k) * (pa(i,k)/500.0)**0.65                      &
     &         * (_one + (5.9590e-4 - 2.2931e-6*dt(i,k)) * dt(i,k))
           ch4exp(i,k,1) = exp( -xc*6.29247e-2 )

            xc = ch4exp(i,k,1) * ch4exp(i,k,1) * ch4exp(i,k,1)
            xc = xc*xc
            ch4exp(i,k,2) = xc*xc

            xc = ch4exp(i,k,2) * ch4exp(i,k,2) * ch4exp(i,k,2)
            xc = xc*xc
            ch4exp(i,k,3) = xc*xc

            xc = ch4exp(i,k,3) * ch4exp(i,k,3) * ch4exp(i,k,3)
            xc = xc*xc
            ch4exp(i,k,4) = xc*xc

          endif

        enddo
      enddo

      return
!...................................
      end subroutine ch4exps
!-----------------------------------



!-----------------------------------
      subroutine comexps                                                &
!...................................

!  ---  inputs:
     &     ( dcom, dt, NPTS, NLAY, ib,                                  &
!  ---  outputs:
     &       comexp                                                     &
     &     )

!**********************************************************************
!   compute co2-minor exponentials for individual layers using 
!   eqs. (8.21) and (8.22).
!
!  --- input parameters
!  layer co2 amount (dcom)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  6 exponentials for each layer (comexp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib

      real (kind=kind_phys), dimension(:,:), intent(in) :: dcom, dt

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: comexp

!  ---  locals:
      integer :: i, k, ik

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  scaling and absorpton data are given in table 6

      do k = 1, NLAY
        do i = 1, NPTS

          if ( ib == 4 ) then
            xc = dcom(i,k) * (_one + (3.5775e-2                         &
     &         + 4.0447e-4*dt(i,k)) * dt(i,k))
          endif

          if ( ib == 5 ) then
            xc = dcom(i,k) * (_one + (3.4268e-2                         &
     &         + 3.7401e-4*dt(i,k)) * dt(i,k))
          endif

          comexp(i,k,1) = exp( -xc*1.922e-7 )

          do ik = 2, NK0
            xc = comexp(i,k,ik-1) * comexp(i,k,ik-1)
            xc = xc*xc
            comexp(i,k,ik) = xc * comexp(i,k,ik-1)
          enddo

        enddo
      enddo
!
      return
!...................................
      end subroutine comexps
!-----------------------------------



!-----------------------------------
      subroutine cfcexps                                                &
!...................................

!  ---  inputs:
     &     ( a1, b1, fk1, a2, b2, fk2, dcfc, dt,                        &
     &       NPTS, NLAY, ib,                                            &
!  ---  outputs:
     &       cfcexp                                                     &
     &     )

!**********************************************************************
!   compute cfc(-11, -12, -22) exponentials for individual layers.
!
!  --- input parameters
!  parameters for computing the scaled cfc amounts
!     for temperature scaling (a1,b1,a2,b2)
!  the absorption coefficients for the
!     first k-distribution function due to cfcs (fk1,fk2)
!  layer cfc amounts (dcfc)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!  spectral band (ib)
!
!  --- output parameters
!  1 exponential for each layer (cfcexp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, ib

      real (kind=kind_phys), dimension(:,:), intent(in) :: dcfc, dt
      real (kind=kind_phys), intent(in) :: a1, b1, fk1, a2, b2, fk2

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: cfcexp

!  ---  locals:
      integer :: i, k

      real (kind=kind_phys) :: xf
!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS

!  ---  compute the scaled cfc amount (xf) and exponential (cfcexp)

          if ( ib == 4 ) then
            xf = dcfc(i,k) * (_one + (a1 + b1*dt(i,k)) * dt(i,k))
            cfcexp(i,k) = exp( -xf*fk1 )
          else
            xf = dcfc(i,k) * (_one + (a2 + b2*dt(i,k)) * dt(i,k))
            cfcexp(i,k) = exp( -xf*fk2 )
          endif

        enddo
      enddo
!
      return
!...................................
      end subroutine cfcexps
!-----------------------------------



!-----------------------------------
      subroutine b10exps                                                &
!...................................

!  ---  inputs:
     &     ( dh2o, dcont, dco2, dn2o, pa, dt,                           &
     &       NPTS, NLAY,                                                &
!  ---  outputs:
     &       h2oexp, conexp, co2exp, n2oexp                             &
     &     )

!**********************************************************************
!  compute band3a exponentials for individual layers
!
!  --- input parameters
!  layer h2o amount for line absorption (dh2o)
!  layer h2o amount for continuum absorption (dcont)
!  layer co2 amount (dco2)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!  number of grid intervals (NPTS)
!  number of layers (NLAY)
!
!  --- output parameters
!  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: dh2o, dcont, &
     &        dn2o, dco2, pa, dt

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: h2oexp,   &
     &        conexp, n2oexp
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) :: co2exp

!  ---  locals:
      integer :: i, k

      real (kind=kind_phys) :: xx, xx1, xx2, xx3
!
!===> ...  begin here
!
      do k = 1, NLAY
        do i = 1, NPTS

!  ---  compute scaled h2o-line amount for band 10 (eq. 4.4 and table 3).

          xx = dh2o(i,k) * (pa(i,k)/500.0)                              &
     &       * (_one + (0.0149 + 6.20e-5*dt(i,k)) * dt(i,k))

!  ---  six exponentials by powers of 8

          h2oexp(i,k,1) = exp( -xx*0.10624 )

          xx = h2oexp(i,k,1) * h2oexp(i,k,1)
          xx = xx*xx
          h2oexp(i,k,2) = xx*xx

          xx = h2oexp(i,k,2) * h2oexp(i,k,2)
          xx = xx*xx
          h2oexp(i,k,3) = xx*xx

          xx = h2oexp(i,k,3) * h2oexp(i,k,3)
          xx = xx*xx
          h2oexp(i,k,4) = xx*xx

          xx = h2oexp(i,k,4) * h2oexp(i,k,4)
          xx = xx*xx
          h2oexp(i,k,5) = xx*xx

!  ---  one exponential of h2o continuum for sub-band 3a (table 9).

          conexp(i,k,1) = exp( -dcont(i,k)*109.0 )

!  ---  scaled co2 amount for the band 10 (eq. 4.4, tables 3 and 6).

          xx = dco2(i,k) * (pa(i,k)/300.0)**0.5                         &
     &       * (_one + (0.0179 + 1.02e-4*dt(i,k)) * dt(i,k))

!  ---  six exponentials by powers of 8

          co2exp(i,k,1,1) = exp( -xx*2.656e-5 )

          xx = co2exp(i,k,1,1) * co2exp(i,k,1,1)
          xx = xx*xx
          co2exp(i,k,2,1) = xx*xx

          xx = co2exp(i,k,2,1) * co2exp(i,k,2,1)
          xx = xx*xx
          co2exp(i,k,3,1) = xx*xx

          xx = co2exp(i,k,3,1) * co2exp(i,k,3,1)
          xx = xx*xx
          co2exp(i,k,4,1) = xx*xx

          xx = co2exp(i,k,4,1) * co2exp(i,k,4,1)
          xx = xx*xx
          co2exp(i,k,5,1) = xx*xx

          xx = co2exp(i,k,5,1) * co2exp(i,k,5,1)
          xx = xx*xx
          co2exp(i,k,6,1) = xx*xx

!  ---  compute the scaled n2o amount for band 10 (table 5).

          xx = dn2o(i,k) * (_one + (1.4476e-3 + 3.6656e-6*dt(i,k))      &
     &       * dt(i,k))

!  ---  two exponentials by powers of 58

          n2oexp(i,k,1) = exp( -xx*0.25238 )

          xx  = n2oexp(i,k,1) * n2oexp(i,k,1)
          xx1 = xx*xx
          xx1 = xx1*xx1
          xx2 = xx1*xx1
          xx3 = xx2*xx2
          n2oexp(i,k,2) = xx*xx1*xx2*xx3

        enddo
      enddo
!
      return
!...................................
      end subroutine b10exps
!-----------------------------------



!-----------------------------------
      subroutine tablup                                                 &
!...................................

!  ---  inputs:
     &     ( dw, p, dt, w1, p1, dwe, dpe, coef1, coef2, coef3,          &
     &       NPTS, k2, nx, nh,                                          &
!  ---  in/outputs:
     &       s1, s2, s3, tran                                           &
     &     )

!**********************************************************************
!   compute water vapor, co2 and o3 transmittances between level
!   k1 and and level k2 for m soundings, using table look-up.
!
!   calculations follow eq. (4.16).
!
!  --- input ---------------------
!  layer absorber amount (dw)
!  layer pressure (p)
!  deviation of layer temperature from 250K (dt)
!  first value of absorber amount (log10) in the table (w1) 
!  first value of pressure (log10) in the table (p1) 
!  size of the interval of absorber amount (log10) in the table (dwe)
!  size of the interval of pressure (log10) in the table (dpe)
!  pre-computed coefficients (coef1, coef2, and coef3)
!  number of grid intervals (NPTS)
!  index for level (k2)
!  number of pressure intervals in the table (nx)
!  number of absorber amount intervals in the table (nh)
!
!  --- updated ---------------------
!
!  column integrated absorber amount (s1)
!  absorber-weighted column pressure (s2)
!  absorber-weighted column temperature (s3)
!  transmittance (tran)
!
!  note: units of s1 are g/cm**2 for water vapor and
!       (cm-atm)stp for co2 and o3.
!   
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k2, nx, nh

      real (kind=kind_phys), dimension(:,:), intent(in) :: dw, p, dt,   &
     &       coef1, coef2, coef3
      real (kind=kind_phys), intent(in) :: w1, p1, dwe, dpe

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:),   intent(inout) :: s1,s2,s3
      real (kind=kind_phys), dimension(:,:), intent(inout) :: tran

!  ---  locals:
      integer :: i, iw, ip

      real (kind=kind_phys) :: we, pe, fw, fp, pa, pb, pc, ax, ba, bb,  &
     &       t1, ca, cb, t2, x1, x2, x3
!
!===> ...  begin here
!
!  ---  compute effective pressure (x2) and temperature (x3) following 
!       eqs. (8.28) and (8.29)

      do i = 1, NPTS

        s1(i) = s1(i) + dw(i,k2-1)
        s2(i) = s2(i) + p (i,k2-1) * dw(i,k2-1)
        s3(i) = s3(i) + dt(i,k2-1) * dw(i,k2-1)

        x1 = s1(i)
        x2 = s2(i) / s1(i)
        x3 = s3(i) / s1(i)

!  ---  normalize we and pe

        we = (log10(x1) - w1) / dwe
        pe = (log10(x2) - p1) / dpe

!  ---  restrict the magnitudes of the normalized we and pe.

        fw = nh - 1
        fp = nx - 1
        we = min( we, fw )
        pe = min( pe, fp )

!  ---  assign iw and ip and compute the distance of we and pe 
!       from iw and ip.

        iw = int(we + _one)
        iw = max(2, min(iw, nh-1))
        fw = we - float(iw-1)

        ip = int(pe + _one)
        ip = max(1, min(ip, nx-1))
        fp = pe - float(ip-1)

!  ---  linear interpolation in pressure

        pa = coef1(ip,iw-1)*(_one-fp) + coef1(ip+1,iw-1)*fp
        pb = coef1(ip,  iw)*(_one-fp) + coef1(ip+1,  iw)*fp
        pc = coef1(ip,iw+1)*(_one-fp) + coef1(ip+1,iw+1)*fp

!  ---  quadratic interpolation in absorber amount for coef1

        ax = (-pa*(_one-fw) + pc*(_one+fw))*fw*0.5 + pb*(_one-fw*fw)

!  ---  linear interpolation in absorber amount for coef2 and coef3

        ba = coef2(ip,  iw)*(_one-fp) + coef2(ip+1,  iw)*fp
        bb = coef2(ip,iw+1)*(_one-fp) + coef2(ip+1,iw+1)*fp
        t1 = ba*(_one-fw) + bb*fw

        ca = coef3(ip,  iw)*(_one-fp) + coef3(ip+1,  iw)*fp
        cb = coef3(ip,iw+1)*(_one-fp) + coef3(ip+1,iw+1)*fp
        t2 = ca*(_one-fw) + cb*fw

!  ---  update the total transmittance between levels k1 and k2

        tran(i,k2) = (ax + (t1 + t2*x3)*x3) * tran(i,k2)
!       tran(i,k2) = min(tran(i,k2), 0.9999999)
!       tran(i,k2) = max(tran(i,k2), 0.0000001)

      enddo
!
      return
!...................................
      end subroutine tablup
!-----------------------------------



!-----------------------------------
      subroutine h2okdis                                                &
!...................................

!  ---  inputs:
     &     ( fkw, gkw, h2oexp, conexp,                                  &
     &       NPTS, k, ib, ne,                                           &
!  ---  in/outputs:
     &       th2o, tcon, tran                                           &
     &     )

!**********************************************************************
!   compute water vapor transmittance between levels k1 and k2 for
!   m soundings, using the k-distribution method.
!
!  --- input parameters
!  planck-weighted k-distribution function due to
!    h2o line absorption (fkw)
!  planck-weighted k-distribution function due to
!    h2o continuum absorption (gkw)
!  exponentials for line absorption (h2oexp) 
!  exponentials for continuum absorption (conexp) 
!  number of grid intervals (NPTS)
!  current level (k)
!  spectral band (ib)
!  number of terms used in each band to compute water vapor
!    continuum transmittance (ne)
!
!---- updated parameters
!  transmittance between levels k1 and k2 due to
!    water vapor line absorption (th2o)
!  transmittance between levels k1 and k2 due to
!    water vapor continuum absorption (tcon)
!  total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k, ib, ne

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: conexp,    &
     &       h2oexp
      real (kind=kind_phys), dimension(:,:),   intent(in) :: fkw, gkw

!  ---  in/outputs:
      real (kind=kind_phys), dimension(:,:), intent(inout) :: th2o,     &
     &       tcon, tran

!  ---  locals:
      integer :: i, j

      real (kind=kind_phys) :: trnth2o
!
!===> ...  begin here
!
!  ---  tco2 are the six exp factors between levels k1 and k2 
!       tran is the updated total transmittance between levels k1 and k2
!       th2o is the 6 exp factors between levels k1 and k2 due to h2o
!            line absorption. 
!       tcon is the 3 exp factors between levels k1 and k2 due to h2o
!            continuum absorption.
!       trnth2o is the total transmittance between levels k1 and k2 due
!            to both line and continuum absorption.

!  ---  compute th2o following eq. (8.23).

      do j = 1, NK0
        do i = 1, NPTS
          th2o(i,j) = th2o(i,j) * h2oexp(i,k,j)
        enddo
      enddo


      if ( ne == 0 ) then

!  ---  compute trnh2o following eq. (8.25). fkw is given in table 4.

        do i = 1, NPTS
          trnth2o = (fkw(1,ib)*th2o(i,1) + fkw(2,ib)*th2o(i,2)          &
     &            +  fkw(3,ib)*th2o(i,3) + fkw(4,ib)*th2o(i,4)          &
     &            +  fkw(5,ib)*th2o(i,5) + fkw(6,ib)*th2o(i,6))

          tran(i,k+1) = tran(i,k+1) * trnth2o
        enddo

      elseif ( ne == 1 ) then

!  ---  compute trnh2o following eqs. (8.25) and (4.27).

        do i = 1, NPTS
          tcon(i,1) = tcon(i,1)*conexp(i,k,1)

          trnth2o = (fkw(1,ib)*th2o(i,1) + fkw(2,ib)*th2o(i,2)          &
     &            +  fkw(3,ib)*th2o(i,3) + fkw(4,ib)*th2o(i,4)          &
     &            +  fkw(5,ib)*th2o(i,5) + fkw(6,ib)*th2o(i,6))         &
     &            * tcon(i,1)

          tran(i,k+1) = tran(i,k+1) * trnth2o
        enddo

      else

!  ---  for band 3. this band is divided into 3 subbands.

        do i = 1, NPTS
          tcon(i,1) = tcon(i,1) * conexp(i,k,1)
          tcon(i,2) = tcon(i,2) * conexp(i,k,2)
          tcon(i,3) = tcon(i,3) * conexp(i,k,3)

!  ---  compute trnh2o following eqs. (4.29) and (8.25).

          trnth2o = (gkw(1,1)*th2o(i,1) + gkw(2,1)*th2o(i,2)            &
     &            +  gkw(3,1)*th2o(i,3) + gkw(4,1)*th2o(i,4)            &
     &            +  gkw(5,1)*th2o(i,5) + gkw(6,1)*th2o(i,6))*tcon(i,1) &
     &            + (gkw(1,2)*th2o(i,1) + gkw(2,2)*th2o(i,2)            &
     &            +  gkw(3,2)*th2o(i,3) + gkw(4,2)*th2o(i,4)            &
     &            +  gkw(5,2)*th2o(i,5) + gkw(6,2)*th2o(i,6))*tcon(i,2) &
     &            + (gkw(1,3)*th2o(i,1) + gkw(2,3)*th2o(i,2)            &
     &            +  gkw(3,3)*th2o(i,3) + gkw(4,3)*th2o(i,4)            &
     &            +  gkw(5,3)*th2o(i,5) + gkw(6,3)*th2o(i,6))*tcon(i,3)

          tran(i,k+1) = tran(i,k+1) * trnth2o
        enddo

      endif
!
      return
!...................................
      end subroutine h2okdis
!-----------------------------------



!-----------------------------------
      subroutine co2kdis                                                &
!...................................

!  ---  inputs:
     &     ( co2exp, NPTS, k,                                           &
!  ---  in/outputs:
     &       tco2, tran                                                 &
     &     )

!**********************************************************************
!   compute co2 transmittances between levels k1 and k2 for
!   m soundings, using the k-distribution method with linear
!   pressure scaling.
!
!  --- input parameters
!   exponentials for co2 absorption (co2exp)
!   number of grid intervals (NPTS)
!   current level (k)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k

      real (kind=kind_phys), intent(in) :: co2exp(:,:,:,:)

!  ---  in/outputs:
      real (kind=kind_phys), intent(inout) :: tco2(:,:,:), tran(:,:)

!  ---  locals:
      integer :: i

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  tco2 is the 6 exp factors between levels k1 and k2 computed
!       from eqs. (8.23) and (8.25). also see eq. (4.30).
!       the k-distribution functions are given in table 10.

      do i = 1, NPTS

!  ---  band-wings

        tco2(i,1,1) = tco2(i,1,1) * co2exp(i,k,1,1)
        xc =      0.1395*tco2(i,1,1)

        tco2(i,2,1) = tco2(i,2,1) * co2exp(i,k,2,1)
        xc = xc + 0.1407*tco2(i,2,1)

        tco2(i,3,1) = tco2(i,3,1) * co2exp(i,k,3,1)
        xc = xc + 0.1549*tco2(i,3,1)

        tco2(i,4,1) = tco2(i,4,1) * co2exp(i,k,4,1)
        xc = xc + 0.1357*tco2(i,4,1)

        tco2(i,5,1) = tco2(i,5,1) * co2exp(i,k,5,1)
        xc = xc + 0.0182*tco2(i,5,1)

        tco2(i,6,1) = tco2(i,6,1) * co2exp(i,k,6,1)
        xc = xc + 0.0220*tco2(i,6,1)

!  ---  band-center region

        tco2(i,1,2) = tco2(i,1,2) * co2exp(i,k,1,2)
        xc = xc + 0.0766*tco2(i,1,2)

        tco2(i,2,2) = tco2(i,2,2) * co2exp(i,k,2,2)
        xc = xc + 0.1372*tco2(i,2,2)

        tco2(i,3,2) = tco2(i,3,2) * co2exp(i,k,3,2)
        xc = xc + 0.1189*tco2(i,3,2)

        tco2(i,4,2) = tco2(i,4,2) * co2exp(i,k,4,2)
        xc = xc + 0.0335*tco2(i,4,2)

        tco2(i,5,2) = tco2(i,5,2) * co2exp(i,k,5,2)
        xc = xc + 0.0169*tco2(i,5,2)

        tco2(i,6,2) = tco2(i,6,2) * co2exp(i,k,6,2)
        xc = xc + 0.0059*tco2(i,6,2)

        tran(i,k+1) = tran(i,k+1) * xc

      enddo
!
      return
!...................................
      end subroutine co2kdis
!-----------------------------------


!-----------------------------------
      subroutine n2okdis                                                &
!...................................

!  ---  inputs:
     &     ( n2oexp, NPTS, k, ib,                                       &
!  ---  in/outputs:
     &       tn2o, tran                                                 &
     &     )

!**********************************************************************
!   compute n2o transmittances between levels k1 and k2 for
!   m soundings, using the k-distribution method with linear
!   pressure scaling.
!
!  --- input parameters
!   exponentials for n2o absorption (n2oexp)
!   number of grid intervals (NPTS)
!   current level (k)
!   spectral band (ib)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k, ib

      real (kind=kind_phys), intent(in) :: n2oexp(:,:,:)

!  ---  in/outputs:
      real (kind=kind_phys), intent(inout) :: tn2o(:,:), tran(:,:)

!  ---  locals:
      integer :: i

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  tn2o is computed from eq. (8.23). 
!       xc is the total n2o transmittance computed from (8.25)
!       the k-distribution functions are given in table 5.

      do i = 1, NPTS

!  ---  band 6

        if ( ib == 6 ) then

          tn2o(i,1) = tn2o(i,1) * n2oexp(i,k,1)
          xc =      0.940414*tn2o(i,1)

          tn2o(i,2) = tn2o(i,2) * n2oexp(i,k,2)
          xc = xc + 0.059586*tn2o(i,2)

!  ---  band 7

        else

          tn2o(i,1) = tn2o(i,1) * n2oexp(i,k,1)
          xc =      0.561961*tn2o(i,1)

          tn2o(i,2) = tn2o(i,2) * n2oexp(i,k,2)
          xc = xc + 0.138707*tn2o(i,2)

          tn2o(i,3) = tn2o(i,3) * n2oexp(i,k,3)
          xc = xc + 0.240670*tn2o(i,3)

          tn2o(i,4) = tn2o(i,4) * n2oexp(i,k,4)
          xc = xc + 0.058662*tn2o(i,4)

        endif

        tran(i,k+1) = tran(i,k+1) * xc

      enddo
!
      return
!...................................
      end subroutine n2okdis
!-----------------------------------


!-----------------------------------
      subroutine ch4kdis                                                &
!...................................

!  ---  inputs:
     &     ( ch4exp, NPTS, k, ib,                                       &
!  ---  in/outputs:
     &       tch4, tran                                                 &
     &     )

!**********************************************************************
!   compute ch4 transmittances between levels k1 and k2 for
!    m soundings, using the k-distribution method with
!    linear pressure scaling.
!
!  --- input parameters
!   exponentials for ch4 absorption (ch4exp)
!   number of grid intervals (NPTS)
!   current level (k)
!   spectral band (ib)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to ch4 absorption
!     for the various values of the absorption coefficient (tch4)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k, ib

      real (kind=kind_phys), intent(in) :: ch4exp(:,:,:)

!  ---  in/outputs:
      real (kind=kind_phys), intent(inout) :: tch4(:,:), tran(:,:)

!  ---  locals:
      integer :: i

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  tch4 is computed from eq. (8.23). 
!       xc is the total ch4 transmittance computed from (8.25)
!       the k-distribution functions are given in table 5.

      do i = 1, NPTS

!  ---  band 6

        if ( ib == 6 ) then

          tch4(i,1) = tch4(i,1) * ch4exp(i,k,1)
          xc = tch4(i,1)

!  ---  band 7

        else

          tch4(i,1)=tch4(i,1) * ch4exp(i,k,1)
          xc =      0.610650*tch4(i,1)

          tch4(i,2)=tch4(i,2) * ch4exp(i,k,2)
          xc = xc + 0.280212*tch4(i,2)

          tch4(i,3)=tch4(i,3) * ch4exp(i,k,3)
          xc = xc + 0.107349*tch4(i,3)

          tch4(i,4)=tch4(i,4) * ch4exp(i,k,4)
          xc = xc + 0.001789*tch4(i,4)

        endif

        tran(i,k+1) = tran(i,k+1) * xc

      enddo
!
      return
!...................................
      end subroutine ch4kdis
!-----------------------------------


!-----------------------------------
      subroutine comkdis                                                &
!...................................

!  ---  inputs:
     &     ( comexp, NPTS, k, ib,                                       &
!  ---  in/outputs:
     &       tcom, tran                                                 &
     &     )

!**********************************************************************
!  compute co2-minor transmittances between levels k1 and k2
!   for m soundings, using the k-distribution method
!   with linear pressure scaling.
!
!  --- input parameters
!   exponentials for co2-minor absorption (comexp)
!   number of grid intervals (NPTS)
!   current level (k)
!   spectral band (ib)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to co2-minor absorption
!     for the various values of the absorption coefficient (tcom)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k, ib

      real (kind=kind_phys), intent(in) :: comexp(:,:,:)

!  ---  in/outputs:
      real (kind=kind_phys), intent(inout) :: tcom(:,:), tran(:,:)

!  ---  locals:
      integer :: i

      real (kind=kind_phys) :: xc
!
!===> ...  begin here
!
!  ---  tcom is computed from eq. (8.23). 
!       xc is the total co2 transmittance computed from (8.25)
!       the k-distribution functions are given in table 6.

      do i = 1, NPTS

!  ---  band 4

        if ( ib == 4 ) then

          tcom(i,1) = tcom(i,1) * comexp(i,k,1)
          xc =      0.12159*tcom(i,1)

          tcom(i,2) = tcom(i,2) * comexp(i,k,2)
          xc = xc + 0.24359*tcom(i,2)

          tcom(i,3) = tcom(i,3) * comexp(i,k,3)
          xc = xc + 0.24981*tcom(i,3)

          tcom(i,4) = tcom(i,4) * comexp(i,k,4)
          xc = xc + 0.26427*tcom(i,4)

          tcom(i,5) = tcom(i,5) * comexp(i,k,5)
          xc = xc + 0.07807*tcom(i,5)

          tcom(i,6) = tcom(i,6) * comexp(i,k,6)
          xc = xc + 0.04267*tcom(i,6)

!  ---  band 5

        else

          tcom(i,1) = tcom(i,1) * comexp(i,k,1)
          xc =      0.06869*tcom(i,1)

          tcom(i,2) = tcom(i,2) * comexp(i,k,2)
          xc = xc + 0.14795*tcom(i,2)

          tcom(i,3) = tcom(i,3) * comexp(i,k,3)
          xc = xc + 0.19512*tcom(i,3)

          tcom(i,4) = tcom(i,4) * comexp(i,k,4)
          xc = xc + 0.33446*tcom(i,4)

          tcom(i,5) = tcom(i,5) * comexp(i,k,5)
          xc = xc + 0.17199*tcom(i,5)

          tcom(i,6) = tcom(i,6) * comexp(i,k,6)
          xc = xc + 0.08179*tcom(i,6)

        endif

        tran(i,k+1) = tran(i,k+1) * xc

      enddo
!
      return
!...................................
      end subroutine comkdis
!-----------------------------------


!-----------------------------------
      subroutine cfckdis                                                &
!...................................

!  ---  inputs:
     &     ( cfcexp, NPTS, k,                                           &
!  ---  in/outputs:
     &       tcfc, tran                                                 &
     &     )

!**********************************************************************
!  compute cfc-(11,12,22) transmittances between levels k1 and k2
!   for m soundings, using the k-distribution method with
!   linear pressure scaling.
!
!  --- input parameters
!   exponentials for cfc absorption (cfcexp)
!   number of grid intervals (NPTS)
!   current level (k)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to cfc absorption
!     for the various values of the absorption coefficient (tcfc)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k

      real (kind=kind_phys), intent(in) :: cfcexp(:,:)

!  ---  in/outputs:
      real (kind=kind_phys), intent(inout) :: tcfc(:), tran(:,:)

!  ---  locals:
      integer :: i
!
!===> ...  begin here
!
!  ---  tcfc is the exp factors between levels k1 and k2. 

      do i = 1, NPTS
        tcfc(i) = tcfc(i) * cfcexp(i,k)
        tran(i,k+1) = tran(i,k+1) * tcfc(i)
      enddo
!
      return
!...................................
      end subroutine cfckdis
!-----------------------------------


!-----------------------------------
      subroutine b10kdis                                                &
!...................................

!  ---  inputs:
     &     ( h2oexp, conexp, co2exp, n2oexp,                            &
     &       NPTS, k,                                                   &
!  ---  outputs:
     &       th2o, tcon, tco2, tn2o, tran                               &
     &     )

!**********************************************************************
!
!   compute h2o (line and continuum),co2,n2o transmittances between
!   levels k1 and k2 for m soundings, using the k-distribution
!   method with linear pressure scaling.
!
!  --- input parameters
!   number of grid intervals (NPTS)
!   current level (k)
!   exponentials for h2o line absorption (h2oexp)
!   exponentials for h2o continuum absorption (conexp)
!   exponentials for co2 absorption (co2exp)
!   exponentials for n2o absorption (n2oexp)
!
!  --- updated parameters
!   transmittance between levels k1 and k2 due to h2o line absorption
!     for the various values of the absorption coefficient (th2o)
!   transmittance between levels k1 and k2 due to h2o continuum
!     absorption for the various values of the absorption
!     coefficient (tcon)
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!**********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, k

      real (kind=kind_phys), dimension(:,:,:),  intent(in) :: h2oexp,   &
     &       conexp, n2oexp
      real (kind=kind_phys), dimension(:,:,:,:),intent(in) :: co2exp

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),   intent(out) :: th2o,     &
     &       tcon, tn2o, tran
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: tco2

!  ---  locals:
      integer :: i
      real (kind=kind_phys) :: xx
!
!===> ...  begin here
!

!  ---  for h2o line. the k-distribution functions are given in table 4.

      do i = 1, NPTS
        th2o(i,1) = th2o(i,1) * h2oexp(i,k,1)
        xx =      0.3153 * th2o(i,1)

        th2o(i,2) = th2o(i,2) * h2oexp(i,k,2)
        xx = xx + 0.4604 * th2o(i,2)

        th2o(i,3) = th2o(i,3) * h2oexp(i,k,3)
        xx = xx + 0.1326 * th2o(i,3)

        th2o(i,4) = th2o(i,4) * h2oexp(i,k,4)
        xx = xx + 0.0798 * th2o(i,4)

        th2o(i,5) = th2o(i,5) * h2oexp(i,k,5)
        xx = xx + 0.0119 * th2o(i,5)

        tran(i,k+1) = xx
      enddo

!  ---  for h2o continuum. note that conexp(i,k,3) is for subband 3a.

      do i = 1, NPTS
        tcon(i,  1) = tcon(i,  1) * conexp(i,k,3)
        tran(i,k+1) = tran(i,k+1) * tcon  (i,1)
      enddo
 
!  ---  for co2 (table 6)

      do i = 1, NPTS
        tco2(i,1,1) = tco2(i,1,1) * co2exp(i,k,1,1)
        xx =      0.2673 * tco2(i,1,1)

        tco2(i,2,1) = tco2(i,2,1) * co2exp(i,k,2,1)
        xx = xx + 0.2201 * tco2(i,2,1)

        tco2(i,3,1) = tco2(i,3,1) * co2exp(i,k,3,1)
        xx = xx + 0.2106 * tco2(i,3,1)

        tco2(i,4,1) = tco2(i,4,1) * co2exp(i,k,4,1)
        xx = xx + 0.2409 * tco2(i,4,1)

        tco2(i,5,1) = tco2(i,5,1) * co2exp(i,k,5,1)
        xx = xx + 0.0196 * tco2(i,5,1)

        tco2(i,6,1) = tco2(i,6,1) * co2exp(i,k,6,1)
        xx = xx + 0.0415 * tco2(i,6,1)

        tran(i,k+1) = tran(i,k+1) * xx
      enddo

!  ---  for n2o (table 5)

      do i = 1, NPTS
        tn2o(i,1) = tn2o(i,1) * n2oexp(i,k,1)
        xx =      0.970831 * tn2o(i,1)

        tn2o(i,2) = tn2o(i,2) * n2oexp(i,k,2)
        xx = xx + 0.029169 * tn2o(i,2)

        tran(i,k+1) = tran(i,k+1) * (xx - _one)
      enddo
!
      return
!...................................
      end subroutine b10kdis
!-----------------------------------


!-----------------------------------
      subroutine cldovlp                                                &
!...................................

!  ---  inputs:
     &      ( fcld, tcldlyr, ict, icb,                                  &
     &        NPTS, NLAY, k2,                                           &
!  ---  in/outputs:
     &        cldhi, cldmd, cldlw, it, im, ib,                          &
     &        fclr                                                      &
     &      )

!***********************************************************************
!     compute the fractional clear line-of-sight between levels k1
!     and k2 following Eqs.(6.18)-(6.21).
!
!  ---  input parameters
!  fcld:    fractional cloud cover of a layer
!  tcldlyr: transmittance of a cloud layer
!  ict:     the level separating high and middle clouds
!  icb:     the level separating middle and low clouds
!  it:      number of cloudy layers in the high-cloud group
!  im:      number of cloudy layers in the middle-cloud group
!  ib:      number of cloudy layers in the low-cloud group
!  NPTS:    number of soundings
!  NLAY:    number of layers
!  k2:      index for the level
!  
!  ---  in/output parameter
!  cldhi:
!  cldmd:
!  cldlw:
!  fclr:    clear line-of-sight between levels k1 and k2
!
!***********************************************************************
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NPTS, NLAY, k2, ict(:), icb(:)

      real (kind=kind_phys), intent(in) :: fcld(:,:), tcldlyr(:,:)

!  ---  in/outputs:
      integer, intent(inout) :: it(:), im(:), ib(:)

      real (kind=kind_phys), intent(inout) :: cldhi(:),cldmd(:),cldlw(:)
      real (kind=kind_phys), intent(inout) :: fclr(:,:)

!  ---  locals:
      integer :: i, j, k, ii
      integer, dimension(NPTS,NLAY) :: itx, imx, ibx
!
!===> ...  begin here
!
      do i = 1, NPTS

!  ---  for high clouds "it" is the number of high-cloud layers

        if ( k2 <= ict(i) ) then
          if ( fcld(i,k2-1) > 0.001 ) then

            it(i) = it(i) + 1
            ii    = it(i)
            itx(i,ii) = k2 - 1

            if ( ii == 1 ) go to 11

!  ---  rearrange the order of cloud layers with increasing cloud amount

            do k = 1, ii-1
              j = itx(i,k)

              if ( fcld(i,j) > fcld(i,k2-1) ) then
                do j = ii-1, k, -1
                  itx(i,j+1) = itx(i,j)
                enddo

                itx(i,k) = k2 - 1
                go to 11
              endif
            enddo

   11       continue

!  ---  compute equivalent black-body high cloud amount

            cldhi(i) = _zero
            do k = 1, ii
              j = itx(i,k)
              cldhi(i) = fcld(i,j) - tcldlyr(i,j)*(fcld(i,j) - cldhi(i))
            enddo

          endif
        endif

!  ---  for middle clouds "im" is the number of middle-cloud layers

        if ( k2 > ict(i) .and. k2 <= icb(i) ) then
          if ( fcld(i,k2-1) > 0.001 ) then

            im(i) = im(i) + 1
            ii    = im(i)
            imx(i,ii) = k2 - 1

            if ( ii == 1 ) go to 21

!  ---  rearrange the order of cloud layers with increasing cloud amount

            do k = 1, ii-1
              j = imx(i,k)

              if ( fcld(i,j) > fcld(i,k2-1) ) then
                do j = ii-1, k, -1
                  imx(i,j+1) = imx(i,j)
                enddo

                imx(i,k) = k2 - 1
                go to 21
              endif
            enddo

   21       continue

!  ---  compute equivalent black-body middle cloud amount

            cldmd(i) = _zero
            do k = 1, ii
              j = imx(i,k)
              cldmd(i) = fcld(i,j) - tcldlyr(i,j)*(fcld(i,j) - cldmd(i))
            enddo

          endif
        endif

!  ---  for low clouds "ib" is the number of low-cloud layers

        if ( k2 > icb(i) ) then
          if ( fcld(i,k2-1) > 0.001 ) then

            ib(i) = ib(i) + 1
            ii    = ib(i)
            ibx(i,ii) = k2 - 1

            if ( ii == 1 ) go to 31

!  ---  rearrange the order of cloud layers with increasing cloud amount

            do k = 1, ii-1
              j = ibx(i,k)

              if ( fcld(i,j) > fcld(i,k2-1) ) then
                do j = ii-1, k, -1
                  ibx(i,j+1) = ibx(i,j)
                enddo

                ibx(i,k) = k2 - 1
                go to 31
              endif
            enddo

   31       continue

!  ---  compute equivalent black-body low cloud amount

            cldlw(i) = _zero
            do k = 1, ii
              j = ibx(i,k)
              cldlw(i) = fcld(i,j) - tcldlyr(i,j)*(fcld(i,j) - cldlw(i))
            enddo

          endif
        endif

!  ---  fclr is the equivalent clear fraction between levels k1 and k2
!       assuming the three cloud groups are randomly overlapped.
!       it follows eqs. (6.20) and (6.21).

        fclr(i,k2) = (_one - cldhi(i)) * (_one - cldmd(i))              &
     &             * (_one - cldlw(i))   

      enddo
!
      return
!...................................
      end subroutine cldovlp
!-----------------------------------


!-----------------------------------
      subroutine sfcflux                                                &
!...................................

!  ---  inputs:
     &     ( tg, eg, ibn, NPTS,                                         &
!  ---  outputs:
     &       bs, dbs, rflxs                                             &
     &     )

!***********************************************************************
!  compute emission and reflection by an inhomogeneous surface
!  with vegetation cover.
!
!  ---  Input parameters
!  index for the spectral band (IBN)
!  number of grid box (NPTS)
!  ground temperature (tg)
!  ground emissivity (eg)
!
!  ---  Output parameters
!  Emission by the surface (bs)
!  Derivative of bs rwt temperature (dbs)
!  Reflection by the surface (rflxs)
!
!**********************************************************************
!
       implicit none

!  ---  inputs:
       integer, intent(in) :: ibn, NPTS

       real (kind=kind_phys), intent(in) :: tg(:), eg(:,:)

!  ---  outputs:
       real (kind=kind_phys), dimension(:), intent(out) :: bs,dbs,rflxs

!  ---  locals:
       integer :: i

!
!===> ...  begin here
!

      call planck                                                       &
!  ---  inputs:
     &     ( ibn, NPTS, tg,                                             &
!  ---  outputs:
     &       bs                                                         &
     &     )

      call plancd                                                       &
!  ---  inputs:
     &     ( ibn, NPTS, tg,                                             &
!  ---  outputs:
     &       dbs                                                        &
     &     )

      do i = 1, NPTS
        rflxs(i) = _one - eg(i,IBN)
      enddo
!
      return
!...................................
      end subroutine sfcflux
!-----------------------------------

!
!........................................!
      end module module_radlw_main       !
!========================================!

