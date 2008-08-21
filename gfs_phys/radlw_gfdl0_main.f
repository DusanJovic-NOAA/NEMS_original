!!!!!  ==========================================================  !!!!!
!!!!!            gfdl0 radiation package description               !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the gfdl0 package includes these parts:                           !
!                                                                      !
!       'radlw_gfdl0_param.f'                                          !
!       'radlw_gfdl0_datatb.f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    the 'radlw_gfdl0_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_gfdl0_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_banddata'    -- data for spectral bands          !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!                                                                      !
!    the 'radlw_gfdl0_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main rrtm lw radiation routine                  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                     !
!            clouds,iovr,aerosols,sfemis,                              !
!            IMAX, NLAY, NLP1, iflip, lprnt,                           !
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
!       'radlw_gfdl0_param,f'                                          !
!       'radlw_gfdl0_datatb,f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!                                                                      !
!                                                                      !
!    original program descriptions:                                    !
!                                                                      !
!       - subroutine fst88 is the main computation module of the long- !
!       wave radiation code. in it all "emissivity" calculations,      !
!       including calls to table lookup subroutines. also after        !
!       calling subroutine "spa88", final combined heating rates and   !
!       ground flux are obtained.                                      !
!       - subroutine spa88 is called to obtain exact cts for water,    !
!       co2 and o3, and approxinate cts co2 and o3 calculations.       !
!       - subroutine e1e290 computes the exchange terms in the flux    !
!       equation for logwave radiation for all terms except the exchange
!       with the top of the atmosphere. the method is a table lookup   !
!       on a precomputed e2 function.                                  !
!       - subroutine e290 computes the exchange terms in the flux      !
!       equation for logwave radiation for all terms except the exchange
!       with the top of the atmosphere. the method is a table lookup   !
!       on a precomputed e2 function.                                  !
!       - subroutine lwr88 computes temperature-corrected co2          !
!       transmission functions and also computes the pressure grid and !
!       layer optical paths.                                           !
!                                                                      !
!                                                                      !
!    references:                                                       !
!      1) schwarzkopf, m. d., and s. b. fels, "the simplified exchange !
!         method revisited: an accurate, rapid method for computation  !
!         of infrared cooling rates and fluxes," jgr, 96 (1991).       !
!      2) schwarzkopf, m. d., and s. b. fels, "improvements to the     !
!         algorithm for computing co2 transmissivities and cooling     !
!         rates," jgr, 90 (1985) 10541-10550.                          !
!      3) fels, s.b., "simple strategies for inclusion of voigt        !
!         effects in infrared cooling calculations," appl. opt., 18    !
!         (1979), 2634-2637.                                           !
!      4) fels, s. b., and m. d. schwarzkopf, "the simplified exchange !
!         approximation: a new method for radiative transfer           !
!         calculations," jas, 32 (1975), 1475-1488.                    !
!                                                                      !
!                                                                      !
!    ncep modifications history log:                                   !
!                                                                      !
!       --- 1989,  original code from gfdl                             !
!                                                                      !
!       jul 2004,  yu-tai hou                                          !
!                  modified to confirm fortran 90/95 standard and      !
!                  written in modular form                             !
!       apr 2005,  yu-tai hou                                          !
!                  minor modifications on module structures            !
!       apr 2007,  yu-tai hou                                          !
!                  standardize optional outputs                        !
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
      use physcons,         only : con_g, con_cp, con_sbc, con_p0,      &
     &                             con_amd, con_amo3
      use module_iounitdef, only : NICO2TR

      use module_radlw_parameters
      use module_radlw_cntr_para
      use module_radlw_banddata
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
!
!     character(24), parameter :: VTAGLW='GFDL-LW v88.2   aug 2004'
      character(24), parameter :: VTAGLW='GFDL-LW v88.2   Apr 2007'

!  ---  constant values (in cgs unit)
      real (kind=kind_phys), parameter :: ginv       = 1.0e-2 / con_g
      real (kind=kind_phys), parameter :: cgs_p0     = 10.0 * con_p0
      real (kind=kind_phys), parameter :: p0inv      = 1.0 / cgs_p0
      real (kind=kind_phys), parameter :: gp0inv     = ginv * p0inv
      real (kind=kind_phys), parameter :: diffctr    = 1.66
      real (kind=kind_phys), parameter :: ramdo3     = con_amo3/con_amd

!  ---  index no of the 40 wide bands used in combined wide band cals. telling
!       which bands between 160-560 cm-1 are included in the first 8 combined bands
      integer, dimension(40) :: iband

      data iband    / 2, 1, 2, 2, 1, 2, 1, 3, 2, 2, 3, 2, 2, 4, 2, 4,   &
     &                2, 3, 3, 2, 4, 3, 4, 3, 7, 5, 6, 7, 6, 5, 7, 6,   &
     &                7, 8, 6, 6, 8, 8, 8, 8 /

!  ---  those data will be set up only once by "rlwinit"

      real (kind=kind_phys), dimension(28,NBLY) :: source, dsrce

      real (kind=kind_phys), dimension(N5040) :: em1v,em1vw,em3v, t2,t4

      real (kind=kind_phys) :: t1(N5040+1), radcon, rco2tr

      real (kind=kind_phys), allocatable, dimension(:) :: stemp, gtemp, &
     &       cdt31,  co231,  c2d31,  cdt38,  co238,  c2d38,             &
     &       cdt71,  co271,  c2d71,  cdt78,  co278,  c2d78,             &
     &       cdtm51, co2m51, c2dm51, cdtm58, co2m58, c2dm58

      real (kind=kind_phys), allocatable, dimension(:,:) ::             &
     &       cdt51,  co251,  c2d51,  cdt58,  co258,  c2d58

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

      public lwrad, rlwinit


! =================
      contains
! =================

!-----------------------------------
      subroutine lwrad                                                  &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,iovr,aerosols,sfemis,                               &
     &       IMAX, NLAY, NLP1, iflip, lprnt,                            &
!  ---  outputs:
     &       hlwc,topflx,sfcflx                                         &
!! ---  optional:
     &,      HLW0,HLWB,FLXPRF                                           &
     &     )

!  *******************************************************************  !
!                                                                       !
!     subroutine lwrad computes temperature-corrected co2 transmission  !
!     functions and also computes the pressure grid and layer optical   !
!     paths.                                                            !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!     plyr   (IMAX,NLAY)    - layer pressures (mb)                      !
!     plvl   (IMAX,NLP1)    - interface pressures (mb)                  !
!     tlyr   (IMAX,NLAY)    - layer temperature (k)                     !
!     tlvl   (IMAX,NLP1)    - interface temperatures (k)                !
!     qlyr   (IMAX,NLAY)    - layer h2o mixing ratio (gm/gm)*see inside !
!     olyr   (IMAX,NLAY)    - layer o3 mixing ratio (gm/gm) *see inside !
!     gasvmr (IMAX,NLAY,:)  - atmospheric gases amount, (not used!)     !
!                       (not used! co2 is fixed for trans func table)   !
!                       (check module_radiation_gases for definition)   !
!     clouds (IMAX,NLAY,:)  - cloud profiles                            !
!                       (check module_radiation_clouds for definition)  !
!                ---  for  iflagcld > 0  ---                            !
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
!                ---  for  iflagliq = 0  ---                            !
!        clouds(:,:,1)  -  layer total cloud fraction                   !
!        clouds(:,:,2)  -  layer cloud optical depth                    !
!        clouds(:,:,3)  -  layer cloud single scattering albedo         !
!        clouds(:,:,4)  -  layer cloud asymmetry factor                 !
!     iovr                  - control flag for cloud overlapping        !
!                             =0: random overlapping clouds             !
!                             =1: max/ran overlapping clouds            !
!     aerosols(IMAX,NLAY,NBDLW,:) - aerosol optical prop. **not used!!**!
!                       (check module_radiation_aerosols for definition)!
!        (:,:,:,1)      - optical depth                                 !
!        (:,:,:,2)      - single scattering albedo                      !
!        (:,:,:,3)      - asymmetry parameter                           !
!     sfemis (IMAX)         - surface emissivity         **not used!!** !
!     IMAX                  - total number of horizontal points         !
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
!not  iaerlw                - control flag for aerosols  **not used!!** !
!                             =0: do not include aerosol effect         !
!                             >0: include aerosol effect                !
!     iflagcld              - cloud optical properties control flag     !
!                             =0: input cld opt dep                     !
!                             =1: input cwp,cip,rew,rei                 !
!not  ico2tran              - control flag for co2 data sources         !
!                             =0: use given pre-calc trans table        !
!                             =1: compute co2 trans data (not yet) ***  !
!                                                                       !
!  output variables:                                                    !
!     hlwc   (IMAX,NLAY)    - total sky heating rate (k/day or k/sec)   !
!     topflx (IMAX)         - radiation fluxes at top, component:       !
!                       (check module_radlw_paramters for definition)   !
!        upfxc                 total sky upward flux at top (w/m2)      !
!        upfx0                 clear sky upward flux at top (w/m2)      !
!     sfcflx (IMAX)         - radiation fluxes at sfc, component:       !
!                       (check module_radlw_paramters for definition)   !
!        upfxc                 total sky upward flux at sfc (w/m2)      !
!        dnfxc                 total sky downward flux at sfc (w/m2)    !
!        dnfx0                 clear sky downward flux at sfc (w/m2)    !
!                                                                       !
!! optional output variables:                                           !
!     hlwb(IMAX,NLAY,NBDLW) - spectral band total sky heating rates     !
!     hlw0   (IMAX,NLAY)    - total sky heating rate (k/day or k/sec)   !
!     flxprf (IMAX,NLP1)    - level radiative fluxes (w/m2), components !
!                       (check module_radlw_paramters for definition)   !
!        netfc                 total sky net lw flux                    !
!        netf0                 clear sky net lw flux                    !
!                                                                       !
!                                                                       !
!                                                                       !
!  subroutine lwrad is called by : grrad                                !
!                                                                       !
!  subroutines called by lwrad : cldprp, fst88                          !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!

      implicit none

!  ---  constant parameters:
!       b0,b1,b2,b3 are coefficients used to correct for the use of 250k
!       in the planck function used in evaluating planck-weighted co2
!       transmission functions. (see ref. 4)

      real (kind=kind_phys), parameter :: b0 = -0.51926410e-4
      real (kind=kind_phys), parameter :: b1 = -0.18113332e-3
      real (kind=kind_phys), parameter :: b2 = -0.10680132e-5
      real (kind=kind_phys), parameter :: b3 = -0.67303519e-7

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1, iovr, iflip

      logical, intent(in) :: lprnt

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
      real (kind=kind_phys), dimension(IMAX,NLP1,NLP1) :: cldfac, co21
      real (kind=kind_phys), dimension(IMAX,2*NLAY+1)  :: empl

      real (kind=kind_phys), dimension(IMAX,NLP1) :: plv, tlv, ply, tly,&
     &       toto3, tphio3, totphi, totvo2, co2sp1, co2sp2, co2r1,      &
     &       co2r2, dco2d1, dco2d2, d2cd21, d2cd22, cntval, tdav,       &
     &       tstdav, vsum3, dift, tlsqu

      real (kind=kind_phys), dimension(IMAX,NLAY) ::  qo3, qh2o, delp,  &
     &       delp2, vv, co2nbl, co2mr, co2md, co2m2d, var1, var2,       &
     &       var3, var4, cldfrac, cwp, cip, rew, rei, cdat1

      real (kind=kind_phys), dimension(IMAX) :: a1, a2, emx1, emx2,     &
     &       topfup, topfu0, sfcfup, grnflx, grnfl0

      real (kind=kind_phys) :: texpsl, tem1, tem2, vsum2, co2r, dco2dt, &
     &       d2cdt2

      real (kind=kind_phys), dimension(IMAX,NLP1) :: flxnet
      real (kind=kind_phys), dimension(IMAX,NLAY) :: heatra

!! ---  for optional clear sky heating and flux profile
      real (kind=kind_phys), dimension(IMAX,NLAY) :: heatr0
      real (kind=kind_phys), dimension(IMAX,NLP1) :: flxne0

      integer :: i, k, k1, kk, kp, LLP1

!
!===> ... begin here
!

      LLP1 = 2 * NLAY + 1

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

!  ---  the internal array is always from top to surface
!       gfdl lw use cgs unit, need a factor of 1.0e3 for pressure

      if (iflip == 0) then        ! input from toa to sfc

        do k = 1, NLAY
          do i = 1, IMAX
            ply (i,k) = plyr(i,k) * 1.0e3
            tly (i,k) = tlyr(i,k)
            plv (i,k) = plvl(i,k) * 1.0e3
            tlv (i,k) = tlvl(i,k)
!test use
!           qh2o(i,k) = max(0.0, qlyr(i,k))                     ! input mass mixing ratio
!ncep model use
            qh2o(i,k) = max(0.0, qlyr(i,k)/(1.0-qlyr(i,k)))     ! input specific humidity
            qo3 (i,k) = max(0.0, olyr(i,k))                     ! input mass mixing ratio
!           qo3 (i,k) = max(0.0, olyr(i,k)*ramdo3)              ! input vol mixing ratio
          enddo
        enddo

        do i = 1, IMAX
          plv(i,1)    = 0.0
          tlv(i,1)    = tly (i,1)
          plv(i,NLP1) = plvl(i,NLP1) * 1.0e3
          tlv(i,NLP1) = tlvl(i,NLP1)
          ply(i,NLP1) = plv (i,NLP1)
          tly(i,NLP1) = tlvl(i,NLP1)
        enddo

        if (iflagcld > 0) then   ! use prognostic cloud method
          do k = 1, NLAY
            do i = 1, IMAX
              cldfrac(i,k)= clouds(i,k,1)
              cwp  (i,k)  = clouds(i,k,2)
              rew  (i,k)  = clouds(i,k,3)
              cip  (i,k)  = clouds(i,k,4)
              rei  (i,k)  = clouds(i,k,5)
            enddo
          enddo
        else                     ! use diagnostic cloud method
          do k = 1, NLAY
            do i = 1, IMAX
              cldfrac(i,k)= clouds(i,k,1)
              cdat1(i,k)  = clouds(i,k,2)
            enddo
          enddo
        endif                    ! end if_iflagcld

      else                        ! input data from sfc to toa

        do k = 1, NLAY
          k1 = NLP1 - k

          do i = 1, IMAX
            ply (i,k) = plyr(i,k1) * 1.0e3
            tly (i,k) = tlyr(i,k1)
            plv (i,k) = plvl(i,k1+1) * 1.0e3
            tlv (i,k) = tlvl(i,k1+1)
!test use
!           qh2o(i,k) = max(0.0, qlyr(i,k1))                    ! input mass mixing ratio
!ncep model use
            qh2o(i,k) = max(0.0, qlyr(i,k1)/(1.0-qlyr(i,k1)))   ! input specific humidity
            qo3 (i,k) = max(0.0, olyr(i,k1))                    ! input mass mixing ratio
!           qo3 (i,k) = max(0.0, olyr(i,k1)*ramdo3)             ! input vol mixing ratio
          enddo
        enddo

        do i = 1, IMAX
          plv(i,1)    = 0.0
          tlv(i,1)    = tlvl(i,NLP1)
          plv(i,NLP1) = plvl(i,1) * 1.0e3
          tlv(i,NLP1) = tlvl(i,1)
          ply(i,NLP1) = plv (i,NLP1)
          tly(i,NLP1) = tlvl(i,1)
        enddo

        if (iflagcld > 0) then   ! use prognostic cloud method
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, IMAX
              cldfrac(i,k)= clouds(i,k1,1)
              cwp  (i,k)  = clouds(i,k1,2)
              rew  (i,k)  = clouds(i,k1,3)
              cip  (i,k)  = clouds(i,k1,4)
              rei  (i,k)  = clouds(i,k1,5)
            enddo
          enddo
        else                     ! use diagnostic cloud method
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, IMAX
              cldfrac(i,k)= clouds(i,k1,1)
              cdat1(i,k)  = clouds(i,k1,2)
            enddo
          enddo
        endif                    ! end if_iflagcld

      endif                       ! end if_iflip

!  ---  compute flux pressures (plv) and differences (delp2,delp)
!       compute flux level temperatures (tlv) and continuum temperature
!       corrections (texpsl)

!old code
      do k = 2, NLAY
        do i = 1, IMAX
          plv(i,k) = 0.5*(ply(i,k-1) + ply(i,k))
          tlv(i,k) = 0.5*(tly(i,k-1) + tly(i,k))
        enddo
      enddo
!new code will use inputs directly

      do k = 1, NLAY
        do i = 1, IMAX
          delp2(i,k) = plv(i,k+1) - plv(i,k)
          delp (i,k) = 1.0 / delp2(i,k)
        enddo
      enddo

!  ---  compute cloud transmission function matrix

      call cldprp                                                       &
!  ---  inputs:
     &     ( cldfrac,cwp,cip,rew,rei,cdat1,                             &
     &       IMAX, NLAY, NLP1, iovr,                                    &
!  ---  outputs:
     &       cldfac                                                     &
     &     )

!  ---  compute optical paths for h2o and o3, using the diffusivity
!       approximation for the angular integration (1.66). obtain the
!       unweighted values(var1,var3) and the weighted values(var2,var4).
!       the quantities h3m4(.0003) and h3m3(.003) appearing in the var2 and
!       var4 expressions are the approximate voigt corrections for h2o and
!       o3,respectively.

      tem1 = diffctr * ginv

      do k = 1, NLAY
        do i = 1, IMAX

!  ---  vv is the layer-mean pressure (in atm),which is not the same as
!       the data level pressure

!old code
          vv  (i,k) = 0.5 * (plv(i,k+1) + plv(i,k)) * p0inv
!new code
!         vv  (i,k) = ply(i,k) * p0inv
          var1(i,k) = delp2(i,k) * qh2o(i,k) * tem1
          var3(i,k) = delp2(i,k) * qo3 (i,k) * tem1
          var2(i,k) = var1(i,k)  * (vv(i,k) + 0.0003)
          var4(i,k) = var3(i,k)  * (vv(i,k) + 0.003)

!  ---   compute optical path for the h2o continuum, using roberts coeffs.
!        (betinw),and temp. correction (texpsl). the diffusivity factor
!        (which cancels out in this expression) is assumed to be 1.66. the
!        use of the diffusivity factor has been shown to be a significant
!        source of error in the continuum calcs.,but the time penalty of
!        an angular integration is severe.
!        temp. correction coeff texpsl = 1800.(1./tly-1./296.)

          texpsl  = exp(1800.0/tly(i,k) - 6.081081081)
          cntval(i,k) = texpsl * qh2o(i,k) * var2(i,k) * betinw         &
     &                / (qh2o(i,k) + 0.622)
        enddo
      enddo

!  ---  compute summed optical paths for h2o,o3 and continuum

      do i = 1, IMAX
        totphi(i,1) = 0.0
        toto3 (i,1) = 0.0
        tphio3(i,1) = 0.0
        totvo2(i,1) = 0.0
      enddo

      do k = 2, NLP1
        do i = 1, IMAX
          totphi(i,k) = totphi(i,k-1) + var2(i,k-1)
          toto3 (i,k) = toto3(i,k-1)  + var3(i,k-1)
          tphio3(i,k) = tphio3(i,k-1) + var4(i,k-1)
          totvo2(i,k) = totvo2(i,k-1) + cntval(i,k-1)
        enddo
      enddo

!  ---  emx1 is the additional pressure-scaled mass from ply(NLAY) to
!       plv(NLAY). it is used in nearby layer and emiss calculations.
!       emx2 is the additional pressure-scaled mass from ply(NLAY) to
!       plv(NLP1). it is used in calculations between flux levels NLAY and NLP1.

      tem1 = gp0inv * diffctr

      do i = 1, IMAX
        tem2    = qh2o(i,NLAY) * ply(i,NLAY) * tem1
        emx1(i) = tem2 * (ply(i,NLAY) - plv(i,NLAY))
        emx2(i) = tem2 * (plv(i,NLP1) - ply(i,NLAY))
      enddo

!  ---  empl is the pressure scaled mass from plv(k) to ply(k) or to
!       ply(k+1)

      do k = 1, NLAY
        do i = 1, IMAX
          empl(i,k+1) = qh2o(i,k) * plv(i,k+1) * tem1                   &
     &                * (plv(i,k+1) - ply(i,k))
        enddo
      enddo

      do k = 1, NLAY-1
        kk = k + NLP1
        k1 = k + 1

        do i = 1, IMAX
          empl(i,kk) = qh2o(i,k1) * plv(i,k1) * tem1                    &
     &               * (ply(i,k1) - plv(i,k1))
        enddo
      enddo

      do i = 1, IMAX
        empl(i,1)    = var2(i,NLAY)
        empl(i,LLP1) = empl(i,LLP1-1)
      enddo

!  ---  compute weighted temperature (tdav) and pressure (tstdav) integrals
!       for use in obtaining temp. difference bet. sounding and std.
!       temp. sounding (dift)

      do i = 1, IMAX
        tstdav(i,1) = 0.0
        tdav(i,1)   = 0.0
      enddo

      do k = 1, NLP1
        do i = 1, IMAX
          vsum3(i,k) = tly(i,k) - stemp(k)
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IMAX
          vsum2         = gtemp(k)    * delp2(i,k)
          tstdav(i,k+1) = tstdav(i,k) + vsum2
          tdav(i,k+1)   = tdav(i,k)   + vsum2 * vsum3(i,k)
        enddo
      enddo

!  ---  evaluate coefficients for co2 pressure interpolation (a1,a2)

      tem2 = 1.0 / 202649.902
      do i = 1, IMAX
        a1(i) = (ply(i,NLP1) - 810600.098) * tem2
        a2(i) = (cgs_p0 - ply(i,NLP1)) * tem2
      enddo

!  ---  perform co2 pressure interpolation on all inputted transmission
!       functions and temp. derivatives
!       successively computing co2r,dco2dt and d2cdt2 is done to save
!       storage (at a slight loss in computation time)

      do k = 1, NLP1
        do i = 1, IMAX
          co2r1 (i,k) =           a1(i)*co231(k) + a2(i)*co238(k)
          d2cd21(i,k) = 1.0e-3 * (a1(i)*c2d31(k) + a2(i)*c2d38(k))
          dco2d1(i,k) = 1.0e-2 * (a1(i)*cdt31(k) + a2(i)*cdt38(k))
          co2r2 (i,k) =           a1(i)*co271(k) + a2(i)*co278(k)
          d2cd22(i,k) = 1.0e-3 * (a1(i)*c2d71(k) + a2(i)*c2d78(k))
          dco2d2(i,k) = 1.0e-2 * (a1(i)*cdt71(k) + a2(i)*cdt78(k))
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IMAX
          co2mr (i,k) =           a1(i)*co2m51(k) + a2(i)*co2m58(k)
          co2md (i,k) = 1.0e-2 * (a1(i)*cdtm51(k) + a2(i)*cdtm58(k))
          co2m2d(i,k) = 1.0e-3 * (a1(i)*c2dm51(k) + a2(i)*c2dm58(k))
        enddo
      enddo

!  ---  compute co2 temperature interpolations for all bands,using dift
!       the case where k=1 is handled first. we are now replacing
!       3-dimensional arrays by 2-d arrays, to save space. thus this
!       calculation is for (i,kp,1)

      do kp = 2, NLP1
        do i = 1, IMAX
          dift(i,kp) = tdav(i,kp) / tstdav(i,kp)
        enddo
      enddo

      do i = 1, IMAX
        co21(i,1,1) = 1.0
        co2sp1(i,1) = 1.0
        co2sp2(i,1) = 1.0
      enddo

      do kp = 2, NLP1
        do i = 1, IMAX

!  ---  calculations for kp>1 for k=1

          co2r         =         a1(i)*co251(kp,1) + a2(i)*co258(kp,1)
          dco2dt       = 1.0e-2*(a1(i)*cdt51(kp,1) + a2(i)*cdt58(kp,1))
          d2cdt2       = 1.0e-3*(a1(i)*c2d51(kp,1) + a2(i)*c2d58(kp,1))
          co21(i,kp,1) = co2r+dift(i,kp)*(dco2dt+0.5*dift(i,kp)*d2cdt2)

!  ---  calculations for (effectively) kp=1,k>kp. these use the
!       same value of dift due to symmetry

          co2r         =         a1(i)*co251(1,kp) + a2(i)*co258(1,kp)
          dco2dt       = 1.0e-2*(a1(i)*cdt51(1,kp) + a2(i)*cdt58(1,kp))
          d2cdt2       = 1.0e-3*(a1(i)*c2d51(1,kp) + a2(i)*c2d58(1,kp))
          co21(i,1,kp) = co2r+dift(i,kp)*(dco2dt+0.5*dift(i,kp)*d2cdt2)
        enddo
      enddo

!  ---  the transmission functions used in spa88 may be computed now.
!       (in the 250 loop,dift really should be (i,1,k), but dift is
!       invariant with respect to k,kp,and so (i,1,k)=(i,k,1))

      do k = 2, NLP1
        do i = 1, IMAX
          co2sp1(i,k) = co2r1(i,k) + dift(i,k)*( dco2d1(i,k)            &
     &                + 0.5*dift(i,k)*d2cd21(i,k) )
          co2sp2(i,k) = co2r2(i,k) + dift(i,k)*( dco2d2(i,k)            &
     &                + 0.5*dift(i,k)*d2cd22(i,k) )
        enddo
      enddo

!  ---  next the case when k=2...l

      do k = 2, NLAY
        do kp = k+1, NLP1
          do i = 1, IMAX
           dift(i,kp) = (tdav(i,kp) - tdav(i,k))                        &
     &                / (tstdav(i,kp) - tstdav(i,k))

           co2r         =         a1(i)*co251(kp,k) + a2(i)*co258(kp,k)
           dco2dt       = 1.0e-2*(a1(i)*cdt51(kp,k) + a2(i)*cdt58(kp,k))
           d2cdt2       = 1.0e-3*(a1(i)*c2d51(kp,k) + a2(i)*c2d58(kp,k))
           co21(i,kp,k) = co2r+dift(i,kp)*(dco2dt+0.5*dift(i,kp)*d2cdt2)

           co2r         =         a1(i)*co251(k,kp) + a2(i)*co258(k,kp)
           dco2dt       = 1.0e-2*(a1(i)*cdt51(k,kp) + a2(i)*cdt58(k,kp))
           d2cdt2       = 1.0e-3*(a1(i)*c2d51(k,kp) + a2(i)*c2d58(k,kp))
           co21(i,k,kp) = co2r+dift(i,kp)*(dco2dt+0.5*dift(i,kp)*d2cdt2)
          enddo
        enddo
      enddo

!  ---  finally the case when k=kp,k=2..nlp1

      do k = 2, NLP1
        do i = 1, IMAX
          dift(i,k)   = 0.5*(vsum3(i,k) + vsum3(i,k-1))

          co2r        =         a1(i)*co251(k,k) + a2(i)*co258(k,k)
          dco2dt      = 1.0e-2*(a1(i)*cdt51(k,k) + a2(i)*cdt58(k,k))
          d2cdt2      = 1.0e-3*(a1(i)*c2d51(k,k) + a2(i)*c2d58(k,k))
          co21(i,k,k) = co2r+dift(i,k)*(dco2dt+0.5*dift(i,k)*d2cdt2)
        enddo
      enddo

!  ---  we aren't doing nbl tfs on the 100 cm-1 bands .

      do k = 1, NLAY
        do i = 1, IMAX
          co2nbl(i,k) = co2mr(i,k) + vsum3(i,k)*( co2md(i,k)            &
     &                + 0.5*vsum3(i,k)*co2m2d(i,k) )
        enddo
      enddo

!  ---  compute temp. coefficient based on tlv(k) (see ref.2)

      do k = 1, NLP1
        do i = 1, IMAX
          if (tlv(i,k) <= 250.0) then
            tem1       = tlv(i,k) - 250.0
            tlsqu(i,k) = b0 + tem1 * (b1 + tem1 * (b2 + b3*tem1))
          else
            tlsqu(i,k) = b0
          endif
        enddo
      enddo

!  ---  apply to all co2 tfs

      do k = 1, NLP1
        do kp = 1, NLP1
          do i = 1, IMAX
            co21(i,kp,k) = co21(i,kp,k)*(1.0-tlsqu(i,kp)) + tlsqu(i,kp)
          enddo
        enddo

        do i = 1, IMAX
          co2sp1(i,k) = co2sp1(i,k)*(1.0-tlsqu(i,1)) + tlsqu(i,1)
          co2sp2(i,k) = co2sp2(i,k)*(1.0-tlsqu(i,1)) + tlsqu(i,1)
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IMAX
          co2nbl(i,k) = co2nbl(i,k)*(1.0-tlsqu(i,k)) + tlsqu(i,k)
        enddo
      enddo

      call fst88                                                        &
!  ---  inputs:
     &     ( plv,ply,tlv,tly,delp,delp2,cldfac,                         &
     &       co21,co2nbl,co2sp1,co2sp2,var1,var2,var3,var4,cntval,      &
     &       toto3,tphio3,totphi,totvo2,emx1,emx2,empl,                 &
     &       IMAX,NLAY,NLP1, LLP1,                                      &
!  ---  output:
     &       heatra,grnflx,topfup,grnfl0,topfu0,flxnet                  &
!! ---  for optional clear sky output:
     &,      heatr0,flxne0                                              &
     &     )

!  ---  prepare for output

      do i = 1, IMAX
        topflx(i)%upfxc = topfup(i)
        topflx(i)%upfx0 = topfu0(i)

        sfcfup(i)       = con_sbc * tlv(i,NLP1)**4
        sfcflx(i)%upfxc = sfcfup(i)
        sfcflx(i)%dnfxc = sfcfup(i) - grnflx(i)
        sfcflx(i)%dnfx0 = sfcfup(i) - grnfl0(i)
      enddo

      if (iflip == 0) then        ! output from toa to sfc

        do k = 1, NLAY
          do i = 1, IMAX
            hlwc(i,k) = heatra(i,k)
          enddo
        enddo

!! ---  for optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
          do i = 1, IMAX
            hlw0(i,k) = heatr0(i,k)
          enddo
          enddo
        endif

!! ---  optional spectral heating rate
!!      *** currently only one broad band
        if ( lhlwb ) then
          do k = 1, NLAY
          do i = 1, IMAX
            hlwb(i,k,1) = heatra(i,k)
          enddo
          enddo
        endif

!! ---  for optional net flux profiles
        if ( lflxprf ) then
          do k = 1, NLP1
          do i = 1, IMAX
            flxprf(i,k)%netfc = flxnet(i,k)
            flxprf(i,k)%netf0 = flxne0(i,k)
          enddo
          enddo
        endif

      else                        ! output from sfc to toa

        do k = 1, NLAY
          k1 = NLP1 - k
          do i = 1, IMAX
            hlwc(i,k) = heatra(i,k1)
          enddo
        enddo

!! ---  for optional clear sky heating rate
        if ( lhlw0 ) then
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, IMAX
              hlw0(i,k) = heatr0(i,k1)
            enddo
          enddo
        endif

!! ---  optional spectral heating rate
!!      *** currently only one broad band
        if ( lhlwb ) then
          do k = 1, NLAY
            k1 = NLP1 - k
            do i = 1, IMAX
              hlwb(i,k,1) = heatra(i,k1)
            enddo
          enddo
        endif

!! ---  for optional net flux profiles
        if ( lflxprf ) then
          do k = 1, NLP1
            k1 = NLP1 - k + 1
            do i = 1, IMAX
              flxprf(i,k)%netfc = flxnet(i,k1)
              flxprf(i,k)%netf0 = flxne0(i,k1)
            enddo
          enddo
        endif

      endif                       ! if_iflip

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

!  *******************************************************************  !
!                                                                       !
!   subroutine rlwinit reads co2 transmission data from unit 'NICO2TR'  !
!   those data are pre-calculated for the model's vertical structure.   !
!   then it initialize some variables for lw radiation calculation.     !
!                                                                       !
!  *******************************************************************  !
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
!  subroutines called by rlwinit : lwtable                              !
!                                                                       !
!  modification history log:                                            !
!                                                                       !
!      mar     1990: k. a. campana    --- original program conrad       !
!                                                                       !
!      aug  9, 2004: yu-tai hou                                         !
!                    modified to be as part of the radiation module     !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: icwp, me, NLAY

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_io4), dimension(NLAY+1) :: stmp4, gtmp4

      real (kind=kind_io4), dimension(NLAY+1,NLAY+1,6) :: co22d4
      real (kind=kind_io4), dimension(NLAY+1,6)   :: co21d34, co21d74
      real (kind=kind_io4), dimension(NLAY,  6)   :: co21d4

      real (kind=kind_io4) :: rco24

      integer   :: i, j, k, kk, NLP1

!
!                 co2 data tables for users vertical coordinate
!
!   the following module blocks contain pretabulated co2 transmission
!       functions, evaluated using the methods of fels and
!       schwarzkopf (1981) and schwarzkopf and fels (1985),
!-----  the 2-dimensional arrays are
!                    co2 transmission functions and their derivatives
!        from 109-level line-by-line calculations made using the 1982
!        mcclatchy tape (12511 lines),consolidated,interpolated
!        to the nmc mrf vertical coordinatte,and re-consolidated to a
!        200 cm-1 bandwidth. the interpolation method is described in
!        schwarzkopf and fels (j.g.r.,1985).
!-----  the 1-dim arrays are
!                  co2 transmission functions and their derivatives
!          for tau(i,i+1),i=1,l,
!            where the values are not obtained by quadrature,but are the
!            actual transmissivities,etc,between a pair of pressures.
!          these used only for nearby layer calculations including qh2o.
!-----  the weighting function gtemp=plv(k)**0.2*(1.+plv(k)/30000.)**0.8/
!         1013250.,where plv(k)=pressure,nmc mrf(new)  l18 data levels for
!         pstar=1013250.
!-----  stemp is us standard atmospheres,1976,at data pressure levels
!        using nmc mrf sigmas,where pstar=1013.25 mb (ptz program)

!
!====> ... begin here
!

      NLP1 = NLAY + 1

      if (me == 0) then
        print *,' - Using GFDL Longwave Radiation, Version: ', VTAGLW
        print *,'   --- Clouds are randomly overlapped in this version'
!       if ( iaerlw > 0 ) then
          print *,'   This version of GFDL-LW does not include ',       &
     &            'aerosols effect.'
!       else
!         print *,'   --- Aerosol effect is NOT included in LW, all'    &
!    &           ,' internal aerosol parameters are reset to zeros'
!       endif
        print *,'   --- Rare gases effect is NOT included in LW'
        print *,'   --- Included absorbers are CO2, O3, and H2O in LW'
        print *,'   --- No up/down, only net flux profiles produced'
      endif

      if (ilwrate == 1) then
!       radcon = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
!       radcon = con_g * 864.0 / con_cp            !   (in k/day)
        radcon = 8.427
      else
!       radcon = con_g * 1.0e-2 / con_cp           !   (in k/second)
        radcon = 8.427 / 86400.0
      endif

!  ---  allocate module arrays

      if (.not. allocated(stemp) ) allocate (stemp (NLP1))
      if (.not. allocated(gtemp) ) allocate (gtemp (NLP1))

      if (.not. allocated(cdt31) ) allocate (cdt31 (NLP1))
      if (.not. allocated(co231) ) allocate (co231 (NLP1))
      if (.not. allocated(c2d31) ) allocate (c2d31 (NLP1))
      if (.not. allocated(cdt38) ) allocate (cdt38 (NLP1))
      if (.not. allocated(co238) ) allocate (co238 (NLP1))
      if (.not. allocated(c2d38) ) allocate (c2d38 (NLP1))

      if (.not. allocated(cdt71) ) allocate (cdt71 (NLP1))
      if (.not. allocated(co271) ) allocate (co271 (NLP1))
      if (.not. allocated(c2d71) ) allocate (c2d71 (NLP1))
      if (.not. allocated(cdt78) ) allocate (cdt78 (NLP1))
      if (.not. allocated(co278) ) allocate (co278 (NLP1))
      if (.not. allocated(c2d78) ) allocate (c2d78 (NLP1))

      if (.not. allocated(cdtm51)) allocate (cdtm51(NLP1))
      if (.not. allocated(co2m51)) allocate (co2m51(NLP1))
      if (.not. allocated(c2dm51)) allocate (c2dm51(NLP1))
      if (.not. allocated(cdtm58)) allocate (cdtm58(NLP1))
      if (.not. allocated(co2m58)) allocate (co2m58(NLP1))
      if (.not. allocated(c2dm58)) allocate (c2dm58(NLP1))

      if (.not. allocated(cdt51) ) allocate (cdt51 (NLP1,NLP1))
      if (.not. allocated(co251) ) allocate (co251 (NLP1,NLP1))
      if (.not. allocated(c2d51) ) allocate (c2d51 (NLP1,NLP1))
      if (.not. allocated(cdt58) ) allocate (cdt58 (NLP1,NLP1))
      if (.not. allocated(co258) ) allocate (co258 (NLP1,NLP1))
      if (.not. allocated(c2d58) ) allocate (c2d58 (NLP1,NLP1))

!  ---  read in pre-computed co2 transmission data

      rewind  NICO2TR

      read(NICO2TR) (stmp4(i),i=1,NLP1)
      read(NICO2TR) (gtmp4(i),i=1,NLP1)

      do k = 1, 6
        read(NICO2TR) (co21d4(i,k),i=1,NLAY)
      enddo

      do k = 1, 6
        read(NICO2TR) ((co22d4(i,j,k),i=1,NLP1),j=1,NLP1)
      enddo

      do k = 1, 6
        read(NICO2TR) (co21d34(i,k),i=1,NLP1)
      enddo

      do k = 1, 6
        read(NICO2TR) (co21d74(i,k),i=1,NLP1)
      enddo

!  ---  read co2 concentration in ppmv

      read(NICO2TR,end=31) rco24
   31 continue

      stemp  = stmp4
      gtemp  = gtmp4

      do i = 1, NLAY
        cdtm51(i) = co21d4(i,1)
        co2m51(i) = co21d4(i,2)
        c2dm51(i) = co21d4(i,3)
        cdtm58(i) = co21d4(i,4)
        co2m58(i) = co21d4(i,5)
        c2dm58(i) = co21d4(i,6)
      enddo

      do j = 1, NLP1
        do i = 1, NLP1
          cdt51(i,j) = co22d4(i,j,1)
          co251(i,j) = co22d4(i,j,2)
          c2d51(i,j) = co22d4(i,j,3)
          cdt58(i,j) = co22d4(i,j,4)
          co258(i,j) = co22d4(i,j,5)
          c2d58(i,j) = co22d4(i,j,6)
        enddo
      enddo

      do i = 1, NLP1
        cdt31(i) = co21d34(i,1)
        co231(i) = co21d34(i,2)
        c2d31(i) = co21d34(i,3)
        cdt38(i) = co21d34(i,4)
        co238(i) = co21d34(i,5)
        c2d38(i) = co21d34(i,6)

        cdt71(i) = co21d74(i,1)
        co271(i) = co21d74(i,2)
        c2d71(i) = co21d74(i,3)
        cdt78(i) = co21d74(i,4)
        co278(i) = co21d74(i,5)
        c2d78(i) = co21d74(i,6)
      enddo

      rco2tr   = rco24 * 1.0e-6

      rewind NICO2TR

!     print 66, NICO2TR
!  66 format(1h ,' --- Read co2 transmission functions from unit ',i2)

!  ---  define tables for lw radiation

      call lwtable
!  ---  inputs: ( none )
!  ---  outputs:( none )

      return
!...................................
      end subroutine rlwinit
!-----------------------------------


!-----------------------------------
      subroutine lwtable
!...................................

!  ---  inputs: ( none )
!  ---  output: ( none )

!  *******************************************************************  !
!                                                                       !
!    purpose:  compute table entries used in the longwave radiation     !
!              program. also calculated are indices used in strip-      !
!              mining and for some pre-computable functions.            !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!                                                                       !
!                                                                       !
!                                                                       !
!  outputs:                                                             !
!      to the module variables source, dsrce                            !
!                                                                       !
!  subroutine lwtable is called by : lwrinit                            !
!                                                                       !
!  subroutines called by lwtable : none                                 !
!                                                                       !
!  *******************************************************************  !
!

      implicit none

!  ---  inputs: ( none )

!  ---  outputs:( none )

!  ---  locals:
      real (kind=kind_phys), dimension(28,180) :: sum,sum3,pertsm,sumwde

      real (kind=kind_phys), dimension(28,NBLW):: src1nb, dbdtnb
      real (kind=kind_phys), dimension(28,NBLX):: srcwd

      real (kind=kind_phys), dimension(NBLW) :: alfanb, arotnb, centnb, &
     &       delnb
      real (kind=kind_phys), dimension(181)  :: zmass, zroot
      real (kind=kind_phys), dimension(180)  :: expo, fac, x2
      real (kind=kind_phys), dimension(28)   :: sc, dsc, xtemv, x, x1,  &
     &       tfour, fortcu, srcs, sum4, sum6, sum7, sum8, sum4wd, s2,   &
     &       r1, r2, t3, r1wd
      real (kind=kind_phys), dimension(30)   :: cnusb, dnusb

      real (kind=kind_phys) :: cent, del, bdhi, bdlo, anu, c1

      integer :: i, i1, ia, j, jp, n, nsubds, nsb


!  ---  compute local quantities for narrow-bands

      do n = 1, NBLW
        centnb(n) = 0.5 * (bandlo(n) + bandhi(n))
        delnb(n)  = bandhi(n) - bandlo(n)
      enddo

!  ---  begin table computations here  ---
!       compute temps, masses for table entries
!       note: the dimensioning and initialization of xtemv and other arrays
!       with dimension of 28 imply a restriction of model temperatures from
!       100k to 370k.
!       the dimensioning of zmass,zroot and other arrays with dimension of
!       180 imply a restriction of model h2o amounts such that optical paths
!       are between 10**-16 and 10**2, in cgs units.

      zmass(1) = 1.0e-16
      do j = 1, 180
        jp        = j + 1
        zroot(j)  = sqrt( zmass(j) )
        zmass(jp) = zmass(j) * 1.258925411
      enddo

      do i = 1, 28
        xtemv(i)  = 90.0 + 10.0*i
        tfour(i)  = 1.0 / (xtemv(i)*xtemv(i)*xtemv(i)*xtemv(i))
        fortcu(i) = 1.0 / (4.0*xtemv(i)*xtemv(i)*xtemv(i))
      enddo

!  ---  the computation of source,dsrce is  needed only for the combined
!       wide-band case.to obtain them,the source must be computed for
!       each of the (nblx) wide bands(=srcwd) then combined (using iband)
!       into source.

      do n = 1, NBLY
        do i = 1, 28
          source(i,n) = 0.0
        enddo
      enddo

      do n = 1, NBLX
        do i = 1, 28
          srcwd(i,n) = 0.0
        enddo
      enddo

!  ---  begin freq. loop (on n)

      do n = 1, NBLX

!  ---  the 160-1200 band cases
        if (n <= 46) then
          cent = centnb(n+16)
          del  = delnb(n+16)
          bdlo = bandlo(n+16)
          bdhi = bandhi(n+16)
        endif

!  ---  the 2270-2380 band case
        if (n == NBLX) then
          cent = centnb(NBLW)
          del  = delnb (NBLW)
          bdlo = bandlo(NBLW)
          bdhi = bandhi(NBLW)
        endif

!  ---  for purposes of accuracy, all evaluations of planck fctns are made
!       on 10 cm-1 intervals, then summed into the (nblx) wide bands.

        nsubds = (del - 1.0e-3) / 10 + 1
        do nsb = 1, nsubds
          if (nsb /= nsubds) then
            cnusb(nsb) = 10.0*(nsb - 1) + bdlo + 5.0
            dnusb(nsb) = 10.0
          else
            cnusb(nsb) = 0.5 * (10.0*(nsb - 1) +bdlo + bdhi)
            dnusb(nsb) = bdhi - (10.0*(nsb - 1) + bdlo)
          endif

          c1 = 3.7412e-5 * cnusb(nsb)**3

!  ---  begin temp. loop (on i)

          do i = 1, 28
            x(i)       = 1.4387 * cnusb(nsb) / xtemv(i)
            x1(i)      = exp( x(i) )
            srcs(i)    = c1 / (x1(i) - 1.0)
            srcwd(i,n) = srcwd(i,n) + srcs(i)*dnusb(nsb)
          enddo
        enddo
      enddo                          ! end of n loop!

!  ---  the following loops create the combined wide band quantities
!       source and dsrce

      do n = 1, 40
        do i = 1, 28
          source(i,iband(n)) = source(i,iband(n)) + srcwd(i,n)
        enddo
      enddo

      do n = 9, NBLY
        do i = 1, 28
          source(i,n) = srcwd(i,n+32)
        enddo
      enddo

      do n = 1, NBLY
        do i = 1, 27
          dsrce(i,n) = (source(i+1,n) - source(i,n)) * 0.1
        enddo
      enddo

      do n = 1, NBLW
        alfanb(n) = brndm(n) * arndm(n)
        arotnb(n) = sqrt( alfanb(n) )
      enddo

!  ---  first compute planck fctns (src1nb) and derivatives (dbdtnb) for
!       use in table evaluations. these are different from source,dsrce
!       because different frequency pts are used in evaluation, the freq.
!       ranges are different, and the derivative algorithm is different.

      do n = 1, NBLW
        cent = centnb(n)
        del  = delnb(n)

!  ---  note: at present, the ia loop is only used for ia=2. the loop struct
!       is kept so that in the future, we may use a quadrature scheme for
!       the planck fctn evaluation, rather than use the mid-band frequency.

        do ia = 1, 3
          anu = cent + 0.5*(ia - 2)*del
          c1  = 3.7412e-5*anu*anu*anu + 1.0e-20

!  ---  temperature loop

          do i = 1, 28
            x(i)   = 1.4387 * anu / xtemv(i)
            x1(i)  = exp( x(i) )
            sc(i)  = c1 / ((x1(i) - 1.0) + 1.0e-20)
            dsc(i) = sc(i)*sc(i)*x(i)*x1(i) / (xtemv(i)*c1)
          enddo

          if (ia == 2) then
            do i = 1, 28
              src1nb(i,n) = del * sc(i)
              dbdtnb(i,n) = del * dsc(i)
            enddo
          endif
        enddo
      enddo

!  ---  next compute r1,r2,s2,and t3- coefficients used for e3 function
!       when the optical path is less than 10-4. in this case, we assume a
!       different dependence on (zmass).
!       also obtain r1wd, which is r1 summed over the 160-560 cm-1 range

      do i = 1, 28
        sum4(i)   = 0.0
        sum6(i)   = 0.0
        sum7(i)   = 0.0
        sum8(i)   = 0.0
        sum4wd(i) = 0.0
      enddo

      do n = 1, NBLW
        cent = centnb(n)

!  ---  perform summations for freq. ranges of 0-560,1200-2200 cm-1 for
!       sub4, sum6, sum7, sum8

        if (cent < 560.0 .or. cent > 1200.0 .and. cent <= 2200.0) then
          do i = 1, 28
            sum4(i) = sum4(i) + src1nb(i,n)
            sum6(i) = sum6(i) + dbdtnb(i,n)
            sum7(i) = sum7(i) + dbdtnb(i,n) * arotnb(n)
            sum8(i) = sum8(i) + dbdtnb(i,n) * alfanb(n)
          enddo
        endif

!  ---  perform summations over 160-560 cm-1 freq range for e1 calcs (sum4wd

        if (cent > 160.0 .and. cent < 560.0) then
          do i = 1, 28
            sum4wd(i) = sum4wd(i) + src1nb(i,n)
          enddo
        endif
      enddo

      do i = 1, 28
        r1(i)   = sum4(i)   * tfour(i)
        r2(i)   = sum6(i)   * fortcu(i)
        s2(i)   = sum7(i)   * fortcu(i)
        t3(i)   = sum8(i)   * fortcu(i)
        r1wd(i) = sum4wd(i) * tfour(i)
      enddo

      do j = 1, 180
        do i = 1, 28
          sum   (i,j) = 0.0
          pertsm(i,j) = 0.0
          sum3  (i,j) = 0.0
          sumwde(i,j) = 0.0
        enddo
      enddo

!  ---  frequency loop begins

      do n = 1, NBLW
        cent = centnb(n)

!  ---  perform calculations for freq. ranges of 0-560,1200-2200 cm-1

        if (cent < 560.0 .or. cent > 1200.0 .and. cent <= 2200.0) then
          do j = 1, 180
            x2(j)   = arotnb(n) * zroot(j)
            expo(j) = exp( -x2(j) )
          enddo

          do j = 1, 180
            if (x2(j) >= 100.0) then
              expo(j) = 0.0
            endif
          enddo

          do j = 121, 180
            fac(j) = zmass(j)*(1.0-(1.0+x2(j))*expo(j))/(x2(j)*x2(j))
          enddo

          do j = 1, 180
            do i = 1, 28
              sum   (i,j) = sum   (i,j) + src1nb(i,n)*expo(j)
              pertsm(i,j) = pertsm(i,j) + dbdtnb(i,n)*expo(j)
            enddo
          enddo

          do j = 121, 180
            do i = 1, 28
              sum3(i,j) = sum3(i,j) + dbdtnb(i,n)*fac(j)
            enddo
          enddo
        endif

!  ---  compute sum over 160-560 cm-1 range for use in e1 calcs (sumwde)

        if (cent > 160.0 .and. cent < 560.0) then
          do j = 1, 180
            do i = 1, 28
              sumwde(i,j) = sumwde(i,j) + src1nb(i,n)*expo(j)
            enddo
          enddo
        endif
      enddo

      i1 = 0
      do j = 1, 180
        do i = 1, 28
          i1 = i1 + 1
          em1v(i1) = sum(i,j)    * tfour(i)
          t1  (i1) = pertsm(i,j) * fortcu(i)
        enddo
      enddo

      i1 = 120 * 28
      do j = 121, 180
        do i = 1, 28
          i1 = i1 + 1
          em3v(i1) = sum3(i,j) * fortcu(i)
        enddo
      enddo

      i1 = 0
      do j = 1, 179
        do i = 1, 28
          i1 = i1 + 1
          t2(i1) = (t1(i1+28) - t1(i1)) * 10.0
        enddo
      enddo

      i1 = 0
      do j = 1, 180
        do i = 1, 27
          t4(i1+i) = (t1(i1+i+1) - t1(i1+i)) * 0.1
        enddo

        i1 = i1 + 28
        t4(i1) = 0.0
      enddo

      i1 = 179 * 28
      do i = 1, 28
        i1 = i1 + 1
        t2(i1) = 0.0
      enddo

      i1 = 0
      do j = 1, 2
        do i = 1, 28
          i1 = i1 + 1
          em1v(i1) = r1(i)
        enddo
      enddo

      i1 = 0
      do j = 1, 120
        do i = 1, 28
          i1 = i1 + 1
          em3v(i1) = r2(i)/2.0 - s2(i)*sqrt(zmass(j))/3.0               &
     &             + t3(i)*zmass(j)/8.0
        enddo
      enddo

      i1 = 120 * 28
      do j = 121, 180
        do i = 1, 28
          i1 = i1 + 1
          em3v(i1) = em3v(i1) / zmass(j)
        enddo
      enddo

!  ---  now compute e1 tables for 160-560 cm-1 bands only. we use r1wd
!       and sumwde obtained above.

      i1 = 0
      do j = 1, 180
        do i = 1, 28
          i1 = i1 + 1
          em1vw(i1) = sumwde(i,j) * tfour(i)
        enddo
      enddo

      i1 = 0
      do j = 1, 2
        do i = 1, 28
          i1 = i1 + 1
          em1vw(i1) = r1wd(i)
        enddo
      enddo

      return
!...................................
      end subroutine lwtable
!-----------------------------------



!-----------------------------------
      subroutine cldprp                                                 &
!...................................

!  ---  inputs:
     &     ( cldfrac,cwp,cip,rew,rei,cdat1,                             &
     &       IMAX, NLAY, NLP1, iovr1,                                   &
!  ---  output:
     &       cldfac                                                     &
     &     )

!  *******************************************************************  !
!                                                                       !
!    subroutine cldprp first puts input layered cloud into gfdl cloud-  !
!    radiation package required cloud structure, then calculate the     !
!    effective cloud amount.  next it computes cloud transmission       !
!    functions for the longwave radiation.                              !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!     cldfrac(IMAX,NLAY) - layer cloud fraction                         !
!        - - -  for iflagcld > 0  (prognostic cloud sckeme)  - - -      !
!     cwp  (IMAX,NLAY)   - layer cloud liquid water path  (g/m**2)      !
!     rew  (IMAX,NLAY)   - effective radius for water cloud (micron)    !
!     cip  (IMAX,NLAY)   - layer cloud ice water path  (g/m**2)         !
!     rei  (IMAX,NLAY)   - effective radius for ice cloud (micron)      !
!     cdat1(IMAX,NLAY)   - not used                                     !
!        - - -  for iflagcld = 0  (diagnostic cloud sckeme)  - - -      !
!     cdat1(IMAX,NLAY)   - input cloud optical depth                    !
!     cwp  (IMAX,NLAY)   - not used                                     !
!     rew  (IMAX,NLAY)   - not used                                     !
!     cip  (IMAX,NLAY)   - not used                                     !
!     rei  (IMAX,NLAY)   - not used                                     !
!                                                                       !
!     IMAX               - number of horizontal points                  !
!     NLAY/NLP1          - number of vertical layers/levels             !
!     iovr               - cloud overlapping control variable           !
!                          =0: random overlapping clouds                !
!                          =1: max/ran overlapping clouds               !
!                                                                       !
!  control parameter in module "module_radlw_cntr_para":                !
!    iflagcld            - =0: diagnostic cloud scheme, cloud optical   !
!                              depth is input                           !
!                          =1: prognostic cloud scheme, cwp/cip and     !
!                              rew/rei are input                        !
!                                                                       !
!  outputs:                                                             !
!    cldfac(IMAX,NLP1,NLP1)                                             !
!                        - cloud transmission function matrix           !
!                                                                       !
!  subroutine cldprp is called by : lwrad                               !
!                                                                       !
!  subroutines called by cldprp : none                                  !
!                                                                       !
!                                                                       !
!  modification history log:                                            !
!                                                                       !
!     dec 1988 - original clo88/clo89 code written by bert katz and     !
!                modified by dan schwarzkopf.                           !
!                                                                       !
!     feb 1993 - yu-tai hou                                             !
!                created the original version of cldprp for mrf model   !
!                to compute cloud radiative properties after davis      !
!                (1982) and harshvardhan et al. (1987).                 !
!     nov 1995 - yu-tai hou                                             !
!                modified to provide mixed cloud overlapping scheme.    !
!     aug 1998 - yu-tai hou                                             !
!                modified cloud emissivity for lw radiation calculation !
!                based on kiehl et al. 1998,j.clim.                     !
!                                                                       !
!     aug 2004 - recode to combine clo89 with cldprp and in fortran     !
!                90/95 standard                                         !
!                                                                       !
!  *******************************************************************  !
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1, iovr1

      real (kind=kind_phys), dimension(:,:), intent(in) :: cldfrac,     &
     &       cwp, cip, rew, rei, cdat1

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: cldfac

!  ---  locals:
      real (kind=kind_phys), dimension(IMAX,NLP1) :: camt
      real (kind=kind_phys), dimension(NLP1,NLP1) :: cldipt
      real (kind=kind_phys), dimension(IMAX) :: xamt, tauc
      real (kind=kind_phys) :: cldrow(NLP1), xcld, awlw, ailw, tau0

      integer, dimension(IMAX,NLP1) :: ktop, kbtm
      integer, dimension(IMAX)      :: nclds
      integer, dimension(IMAX)      :: kcld, mbtm, mthk
      integer :: i, ir, j, k, kb, k1, k2, kp, kt, nc, iovr

      logical, dimension(IMAX) :: firstc, topgap
!
!===> ... begin here
!
      iovr = 0           ! currently only random overlapping is used
!     iovr = iovr1

!  --- ...  put layered cloud fraction into gfdl's cloud-radiation
!           form. where ktop and kbtm values are model layers for
!           cloud tops and bases (from toa to surface). the k-index
!           of ktop, kbtm, and camt is arranged as k=1 is for the
!           surface, k=2 is for the lowest cloud, k=3 is for the 2nd
!           lowest cloud, and so on.  nclds indecates how many
!           individual clouds in the atmospheric column.

      do i = 1, IMAX
        kcld(i)   = 2
        ktop(i,1) = NLP1
        kbtm(i,1) = NLP1
        camt(i,1) = 1.0
      enddo

      do k = 2, NLP1
        do i = 1, IMAX
          ktop(i,k) = 1
          kbtm(i,k) = 1
          camt(i,k) = 0.0
        enddo
      enddo

!  ---  loop over mdl layers (bottom up)

      if (iovr == 0) then             ! random overlapping clouds

        do k = 1, NLAY
          k2 = NLP1 - k

          do i = 1, IMAX
            if (cldfrac(i,k2) > 0.0) then
              k1 = kcld(i)
              ktop(i,k1)  = k2
              kbtm(i,k1)  = k2
              kcld(i)     = k1 + 1

              if (iflagcld == 0) then       ! use input cloud optical depth
                camt(i,k1) = cldfrac(i,k2)                              &
     &                     * (1.0 - exp( -0.75*cdat1(i,k2) ))
              else                          ! compute cloud optical depth
                awlw    = 0.078 * cwp(i,k2)
                ailw    = (0.0029 + 1.0/rei(i,k2)) * cip(i,k2)
                camt(i,k1) = cldfrac(i,k2) * max( 0.0, min( 1.0,        &
     &                       1.0 - exp( -1.66*(awlw+ailw) ) ))
              endif
            endif
          enddo   ! end i_loop

        enddo     ! end k_loop

      elseif (iovr == 1) then         ! max/ran overlapping clouds

        do i = 1, IMAX
          firstc(i) = .true.
          topgap(i) = .false.
          mbtm(i) = 1
          mthk(i) = 0
          tauc(i) = 0.0
        enddo

        do k = 1, NLAY
          k2 = NLP1 - k

          do i = 1, IMAX

            if (cldfrac(i,k2) > 0.0) then
              k1 = kcld(i)

              if (iflagcld == 0) then       ! use input cloud optical depth
                tau0 = cldfrac(i,k2)*cdat1(i,k2)
              else                          ! compute cloud optical depth
                awlw = 0.078 * cwp(i,k2)
                ailw = (0.0029 + 1.0/rei(i,k2)) * cip(i,k2)
                tau0 = cldfrac(i,k2)*(awlw + ailw)
              endif                         ! end if_iflagcld block

              if (k < NLAY) then
                topgap(i) = cldfrac(i,k2-1) < 0.001
              else
                topgap(i) = .true.
              endif

              if (firstc(i)) then
                firstc(i) = .false.
                xamt(i) = cldfrac(i,k2)
                mbtm(i) = k2
                mthk(i) = 1
                tauc(i) = tau0
              else
                xamt(i) = max( xamt(i), cldfrac(i,k2) )
                xamt(i) = xamt(i) + cldfrac(i,k2)
                tauc(i) = tauc(i) + tau0
                mthk(i) = mthk(i) + 1
              endif

              if (topgap(i)) then
                xcld = 1.0 / float(mthk(i))
                if (iflagcld == 0) then
                  camt(i,k1) = xamt(i)*xcld                             &
     &                       * (1.0 - exp( -0.75*tauc(i)*xcld ))
                else
                  camt(i,k1) = xamt(i)*xcld                             &
     &                       * (1.0 - exp( -1.66*tauc(i)*xcld ))
                endif

                ktop(i,k1) = k2
                kbtm(i,k1) = mbtm(i)
                kcld(i)    = k1 + 1
                xamt(i)    = 0.0
                tauc(i)    = 0.0
                mbtm(i)    = 1
                mthk(i)    = 0
                firstc(i)  = .true.
              endif

            endif   ! end if_cldfrac block

          enddo   ! end do_i loop
        enddo     ! end do_k loop

      endif       ! end if_iovr block

!  ---  record number of cloud layers

      do i = 1, IMAX
        nclds(i) = kcld(i) - 2
      enddo

      lab_do_ir : do ir = 1, IMAX

        lab_if_nclds0 : if (nclds(ir) == 0) then

          cldipt(:,:) = 1.0

        else  lab_if_nclds0

          lab_if_nclds1 : if (nclds(ir) >= 1) then

            xcld = 1.0 - camt(ir,2)
            k1   = ktop(ir,2) + 1
            k2   = kbtm(ir,2)

            do j = 1, NLP1
              if (j <= k2) then
                cldrow(j) = xcld
              else
                cldrow(j) = 1.0
              endif
            enddo

            kb = max(k1, k2+1)
            do k = kb, NLP1
              do kp = 1, NLP1
                cldipt(kp,k) = cldrow(kp)
              enddo
            enddo

            do j = 1, NLP1
              if (j < k1) then
                cldrow(j) = 1.0
              else
                cldrow(j) = xcld
              endif
            enddo

            kt = min(k1-1, k2)
            do k = 1, kt
              do kp = 1, NLP1
                cldipt(kp,k) = cldrow(kp)
              enddo
            enddo

            if (k2+1 <= k1-1) then
              do j = k2+1, k1-1
                do i = 1, NLP1
                  cldipt(i,j) = 1.0
                enddo
              enddo
            else if(k1 <= k2) then
              do j = k1, k2
                do i = 1, NLP1
                  cldipt(i,j) = xcld
                enddo
              enddo
            endif

          endif  lab_if_nclds1

          lab_if_nclds2 : if (nclds(ir) >= 2) then

            lab_do_nc : do nc = 2, nclds(ir)
              xcld = 1.0 - camt(ir,nc+1)
              k1   = ktop(ir,nc+1) + 1
              k2   = kbtm(ir,nc+1)

              do j = 1, NLP1
                if (j <= k2) then
                  cldrow(j) = xcld
                else
                  cldrow(j) = 1.0
                endif
              enddo

              kb = max(k1, k2+1)
              do k = kb, NLP1
                do kp = 1, NLP1
                  cldipt(kp,k) = cldipt(kp,k) * cldrow(kp)
                enddo
              enddo

              do j = 1, NLP1
                if (j < k1) then
                  cldrow(j) = 1.0
                else
                  cldrow(j) = xcld
                endif
              enddo

              kt = min(k1-1, k2)
              do k = 1, kt
                do kp = 1, NLP1
                  cldipt(kp,k) = cldipt(kp,k) * cldrow(kp)
                enddo
              enddo

              if (k1 <= k2) then
                do j = k1, k2
                  do i = 1, NLP1
                    cldipt(i,j) = cldipt(i,j) * xcld
                  enddo
                enddo
              endif
            enddo  lab_do_nc

          endif  lab_if_nclds2

        endif lab_if_nclds0

        do j = 1, NLP1
          do i = 1, NLP1
            cldfac(ir,i,j) = cldipt(i,j)
          enddo
        enddo

      enddo  lab_do_ir

      return
!...................................
      end subroutine cldprp
!-----------------------------------



!-----------------------------------
      subroutine fst88                                                  &
!...................................

!  ---  inputs:
     &     ( plv,ply,tlv,tly,delp,delp2,cldfac,                         &
     &       co21,co2nbl,co2sp1,co2sp2,var1,var2,var3,var4,cntval,      &
     &       toto3,tphio3,totphi,totvo2,emx1,emx2,empl,                 &
     &       IMAX, L, LP1, LLP1,                                        &
!  ---  output:
     &       heatra,grnflx,topflx,grnfl0,topfl0,flxnet                  &
!! ---  for optional clear sky output
     &,      heatr0,flxne0                                              &
     &     )

!  *******************************************************************  !
!                                                                       !
!    purpose:  subroutine fst88 is the main computation module of the   !
!              longwave radiation code. in it all "emissivity"          !
!              calculations, including calls to table lookup subroutines!
!              also, after calling subroutine "spa88", final combined   !
!              heating rates and ground flux are obtained.              !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  subroutine fst88 is called by : lwrad                                !
!                                                                       !
!  subroutines called by fst88 : e1e290, e290, spa88                    !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!   passed variables:                                                   !
!         in e3v88:                                                     !
!     emd     =  e3 function for h2o lines (0-560,1200-2200 cm-1)       !
!                computed in e3v88                                      !
!     tpl     =  temperature input for e3 calculation in e3v88          !
!     empl    =  h2o amount,input for e3 calculation in e3v88           !
!                (computed in lwrad; stored in kdacom.h)                !
!         in e1e288:                                                    !
!     e1cts1  =  e1 function for the (i+1)th level using the            !
!                temperature of the ith data level,computed over        !
!                the frequency range 0-560,1200-2200 cm-1. (e1cts1-     !
!                e1ctw1) is used in obtaining the flux at the top       !
!                in the 0-160,1200-2200 cm-1 range (flx1e1).            !
!     e1cts2  =  e1 function for the ith level, using the temp. of      !
!                the ith data level,computed over the frequency range   !
!                0-560,1200-2200 cm-1. (e1cts2-e1ctw2) is also used     !
!                in obtaining the flux at the top in the 0-160,.        !
!                1200-2200 cm-1 range.                                  !
!     e1flx   =  e1 fctn. for the ith level,using the temperature at    !
!                the top of the atmosphere. computed over the freq.     !
!                range 0-560,1200-2200 cm-1. used for q(approx) term.   !
!                (in module block tfcom)                                !
!     e1ctw1  =  like e1cts1,but computed over the 160-560 cm-1 range   !
!                and used for q(approx,cts) calculation                 !
!     e1ctw2  =  like e1cts2,but computed over the 160-560 cm-1 range   !
!                and used for q(approx,cts) calculation                 !
!     fxo     =  temperature index used for e1 function and also        !
!                used for source function calc. in fst88.               !
!     dt      =  temp. diff.between model temps. and temps. at          !
!                tabular values of e1 and source fctns. used in         !
!                fst88 and in e1 function calc.                         !
!     fxoe2   =  temperature index used for e2 function                 !
!     dte2    =  temp. diff. between model temp. and temps. at          !
!                tabular values of e2 function.                         !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, L, LP1, LLP1       

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: cldfac

      real (kind=kind_phys), dimension(:,:),   intent(in) ::            &
     &       plv, ply, tlv, tly, delp, delp2, var1, var2, var3, var4,   &
     &       co2nbl, co2sp1, co2sp2, cntval, toto3, tphio3, totphi,     &
     &       totvo2, empl

      real (kind=kind_phys), dimension(:), intent(in) :: emx1, emx2

!  ---  input/output:
      real (kind=kind_phys), dimension(:,:,:), intent(inout) :: co21

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:),intent(out):: heatra,flxnet

      real (kind=kind_phys), dimension(:),  intent(out):: topflx,grnflx,&
     &       topfl0, grnfl0

!! ---  for optional clear sky outputs:
      real (kind=kind_phys), dimension(:,:),intent(out) :: heatr0

!  ---  locals:
      real (kind=kind_phys), dimension(IMAX,LP1) :: avephi, emiss,      &
     &       emissb, e1flx, co2sp, to3sp, oss, css, ss1, ss2, tc,       &
     &       dtc, csour, avvo2, over1d, to31d, cont1d, avmo3, avpho3,   &
     &       vtmp3, delpr1, delpr2, emisdg, contdg, to3dg, vsum1,       &
     &       e1cts1, e1ctw1, fxo, fxoe2, dt, dte2, flx, totevv, cnttau, &
     &       vtmp30, vsum10, flx0, flxne0

      real (kind=kind_phys), dimension(IMAX,L)   :: excts, ctso3, cts,  &
     &       to3spc, e1cts2, e1ctw2, rlog, excts0, ctso30, cts0

      real (kind=kind_phys), dimension(IMAX,LLP1):: c, c2, alp, tpl

      real (kind=kind_phys), dimension(IMAX,2)   :: emspec, fxosp, dtsp

      real (kind=kind_phys), dimension(IMAX)     :: gxcts, flx1e1,      &
     &       gxcts0, flx1e10

      real (kind=kind_phys), dimension(IMAX,LP1,NBLY) :: sorc

      real (kind=kind_phys) :: vtmp, fac1, tem, tmp3, du, fyo, dt3,     &
     &       ww1, ww2, fxo3, csub2

      integer :: i, k, k1, kk, kp, kk1, kkk, klen, LL, LM1, LLM1,       &
     &           it, ival, item

      logical, parameter :: lthick = .false.
!
!===> ... begin here
!
      LM1  = L - 1
      LL   = LLP1 - 1
      LLM1 = LL - 1


!  ---  first section is table lookup for source function and derivative
!       (b and db/dt).  below, decrementing the index by 9 accounts for
!       the tables beginning at t=100k.
!
!                                 ******* e1 source *******
      do k = 1, LP1
        do i = 1, IMAX
          vtmp     = aint( tly(i,k)*0.1 )
          fxo(i,k) = vtmp  - 9.0
          dt (i,k) = tly(i,k) - 10.0*vtmp

          item         = fxo(i,k)

!       source function for 14 combined bands
!         band 9  - (560-670 cm-1)  band 10 - (670-800 cm-1)
!         band 11 - (800-900 cm-1)  band 12 - (900-990 cm-1)
!         band 13 - (990-1070 cm-1) band 14 - (1070-1200 cm-1)

          sorc(i,k,1)  = source(item,1)  + dt(i,k)*dsrce(item,1)  ! band 1
          sorc(i,k,2)  = source(item,2)  + dt(i,k)*dsrce(item,2)  ! band 2
          sorc(i,k,3)  = source(item,3)  + dt(i,k)*dsrce(item,3)  ! band 3
          sorc(i,k,4)  = source(item,4)  + dt(i,k)*dsrce(item,4)  ! band 4
          sorc(i,k,5)  = source(item,5)  + dt(i,k)*dsrce(item,5)  ! band 5
          sorc(i,k,6)  = source(item,6)  + dt(i,k)*dsrce(item,6)  ! band 6
          sorc(i,k,7)  = source(item,7)  + dt(i,k)*dsrce(item,7)  ! band 7
          sorc(i,k,8)  = source(item,8)  + dt(i,k)*dsrce(item,8)  ! band 8
          sorc(i,k,9)  = source(item,9)  + dt(i,k)*dsrce(item,9)  ! band 9
          sorc(i,k,10) = source(item,10) + dt(i,k)*dsrce(item,10) ! band 10
          sorc(i,k,11) = source(item,11) + dt(i,k)*dsrce(item,11) ! band 11
          sorc(i,k,12) = source(item,12) + dt(i,k)*dsrce(item,12) ! band 12
          sorc(i,k,13) = source(item,13) + dt(i,k)*dsrce(item,13) ! band 13
          sorc(i,k,14) = source(item,14) + dt(i,k)*dsrce(item,14) ! band 14
        enddo
      enddo

!  ---  temp. indices for e2 (kp=1 layer not used in flux calculations)

      do k = 1, L
        do i = 1, IMAX
          vtmp       = aint( tlv(i,k+1)*0.1 )
          fxoe2(i,k) = vtmp - 9.0
          dte2 (i,k) = tlv(i,k+1) - 10.0*vtmp
        enddo
      enddo

!  ---  special case to handle kp=LP1 layer and special e2 calcs.

      do i = 1, IMAX
        fxoe2(i,LP1) = fxo(i,L)
        dte2 (i,LP1) = dt(i,L)
        fxosp(i,1)   = fxoe2(i,LM1)
        fxosp(i,2)   = fxo(i,LM1)
        dtsp (i,1)   = dte2(i,LM1)
        dtsp (i,2)   = dt(i,LM1)
      enddo

!  ---  obtain special source functions for the 15 um band (csour) and
!       the window region (ss1).  also compute tly**4 (tc) and vertical
!       temperature differences (oss,css,ss2,dtc). all these will be
!       used later in flux computations.

      do k = 1, LP1
        do i = 1, IMAX
          ss1  (i,k) = sorc(i,k,11) + sorc(i,k,12) + sorc(i,k,14)
          csour(i,k) = sorc(i,k,9) + sorc(i,k,10)
          vtmp       = tly(i,k) * tly(i,k)
          tc   (i,k) = vtmp * vtmp
        enddo
      enddo

      do k = 1, L
        k1 = k + 1

        do i = 1, IMAX
          oss(i,k1) = sorc(i,k1,13) - sorc(i,k,13)
          css(i,k1) = csour(i,k1)   - csour(i,k)
          dtc(i,k1) = tc(i,k1)      - tc(i,k)
          ss2(i,k1) = ss1(i,k1)     - ss1(i,k)
        enddo
      enddo

!  ---  the followimg is a drastic rewrite of the radiation code to
!       (largely) eliminate three-dimensional arrays. the code works
!       on the following principles:
!
!       let k = fixed flux level, kp = varying flux level then
!           flux(k)=sum over kp : (deltab(kp)*tau(kp,k))
!              over all kp's, from 1 to LP1.
!
!       we can break down the calculations for all k's as follows:
!       for all k's k=1 to LP1:
!           flux(k)=sum over kp : (deltab(kp)*tau(kp,k))      (1)
!              over all kp's, from k+1 to LP1
!       and for kp from k+1 to LP1:
!           flux(kp) = deltab(k)*tau(k,kp)                    (2)
!
!       now if tau(k,kp)=tau(kp,k) (symmetrical arrays) we can compute
!       a 1-dimensional array tau1d(kp) from k+1 to LP1, each time k
!       is incremented. equations (1) and (2) then become:
!
!           tau1d(kp) = (values for tau(kp,k) at the particular k)
!           flux(k) = sum over kp : (deltab(kp)*tau1d(kp))    (3)
!           flux(kp) = deltab(k)*tau1d(kp)                    (4)
!
!       the terms for tau (k,k) and other special terms (for nearby
!       layers) must, of course, be handled separately, and with care.
!
!       compute "upper triangle" transmission functions for the 9.6 um
!       band (to3sp) and the 15 um band (over1d). also, the stage 1...
!       compute o3 ,over transmission fctns and avephi

      do k = 1, L
        do i = 1, IMAX
          avephi(i,k) = totphi(i,k+1)
        enddo
      enddo

!  ---  in order to properly evaluate emiss integrated over the (LP1)
!       layer, a special evaluation of emiss is done. this requires
!       a special computation of avephi, and it is stored in the
!       (otherwise vacant) LP1'th position

      do i = 1, IMAX
        avephi(i,LP1) = avephi(i,LM1) + emx1(i)
      enddo

!  ---  compute fluxes for k=1

      call e1e290                                                       &
!  ---  inputs:
     &     ( fxo,fxoe2,dt,dte2,avephi,em1v,em1vw,t1,t2,t4,              &
     &       IMAX, L, LP1,                                              &
!  ---  outputs:
     &       e1cts1,e1cts2,e1flx,e1ctw1,e1ctw2,emiss                    &
     &     )


      do k = 1, L
        k1 = k + 1

        do i = 1, IMAX
          fac1        = bo3rnd * tphio3(i,k1) / toto3(i,k1)
          to3spc(i,k) = 0.5 * (fac1*(sqrt(1.0                           &
     &                + (4.0*ao3rnd*toto3(i,k1))/fac1) - 1.0))

!  ---  for k=1, to3sp is used instead of to31d (they are equal in this
!       case); to3sp is passed to spa90, while to31d is a work-array.

          to3sp(i,k)  = exp( -1.0*( to3spc(i,k) + sko3r*totvo2(i,k1) ))
          over1d(i,k) = exp( -1.0*( sqrt(ab15wd*totphi(i,k1))           &
     &                            + skc1r*totvo2(i,k1) ))

!  ---  because all continuum transmissivities are obtained from the 2-d
!       quantity cnttau (and its reciprocal totevv) we store both of
!       these here. for k=1, cont1d equals cnttau

          cnttau(i,k) = exp( -1.0*totvo2(i,k1) )
          totevv(i,k) = 1.0 / max(cnttau(i,k), 1.0e-25)
        enddo
      enddo

      do k = 1, L
        do i = 1, IMAX
          co2sp(i,k+1) = over1d(i,k)*co21(i,1,k+1)
        enddo
      enddo

      do k = 1, L
        k1 = k + 1

        do i = 1, IMAX
          co21(i,k1,1) = co21(i,k1,1)*over1d(i,k)
        enddo
      enddo

!  ---  rlog is the nbl amount for the 15 um band calculation

      do i = 1, IMAX
        rlog(i,1) = over1d(i,1)*co2nbl(i,1)
      enddo

!  ---  the terms when kp=1 for all k are the photon exchange with the
!       top of the atmosphere, and are obtained differently than the
!       other calculations

      do k = 2, LP1
        do i = 1, IMAX
          tem  = tc(i,1)*e1flx(i,k)        + ss1(i,1)*cnttau(i,k-1)     &
     &         + sorc(i,1,13)*to3sp(i,k-1) + csour(i,1)*co2sp(i,k)
          flx (i,k) = tem * cldfac(i,1,k)
          flx0(i,k) = tem
        enddo
      enddo

      do i = 1, IMAX
        flx (i,1) = tc(i,1)*e1flx(i,1)+ss1(i,1)+sorc(i,1,13)+csour(i,1)
        flx0(i,1) = flx(i,1)
      enddo

!  ---  the kp terms for k=1...

      do k = 2, LP1
        do i = 1, IMAX
          tem      = oss(i,k)*to3sp(i,k-1) + ss2(i,k)*cnttau(i,k-1)     &
     &             + css(i,k)*co21(i,k,1)  + dtc(i,k)*emiss(i,k-1)
          flx (i,1) = flx (i,1) + tem*cldfac(i,k,1)
          flx0(i,1) = flx0(i,1) + tem
        enddo
      enddo

!  ---  subroutine spa88 is called to obtain exact cts for water
!       co2 and o3, and approximate cts co2 and o3 calculations.

      call spa88                                                        &
!  ---  inputs:
     &     ( plv,ply,tly,var1,var2,delp,cldfac,                         &
     &       totvo2,to3sp,to3spc,co2sp,co2sp1,co2sp2,sorc,csour,        &
     &       IMAX, L, LP1,                                              &
!  ---  output:
     &       excts,ctso3,gxcts,excts0,ctso30,gxcts0                     &
     &     )

!  ---  this section computes the emissivity cts heating rates for 2
!       emissivity bands: the 0-160,1200-2200 cm-1 band and the 800-
!       990,1070-1200 cm-1 band. the remaining cts comtributions are
!       contained in ctso3, computed in spa88.

      do i = 1, IMAX
        vtmp3 (i,1) = 1.0
        vtmp30(i,1) = 1.0
      enddo

      do k = 1, L
        do i = 1, IMAX
          vtmp3 (i,k+1) = cnttau(i,k)*cldfac(i,k+1,1)
          vtmp30(i,k+1) = cnttau(i,k)
        enddo
      enddo

      do k = 1, L
        do i = 1, IMAX
          cts (i,k) = tc(i,k)*( e1ctw2(i,k)*cldfac(i,k+1,1)             &
     &                       - e1ctw1(i,k)*cldfac(i,k,1) )              &
     &                       + ss1(i,k)*(vtmp3 (i,k+1)-vtmp3 (i,k))
          cts0(i,k) = tc(i,k)*( e1ctw2(i,k) - e1ctw1(i,k) )             &
     &                       + ss1(i,k)*(vtmp30(i,k+1)-vtmp30(i,k))
        enddo
      enddo

      do k = 1, L
        do i = 1, IMAX
          vtmp3 (i,k)=tc(i,k)*(cldfac(i,k,1)*(e1cts1(i,k)-e1ctw1(i,k))  &
     &                      - cldfac(i,k+1,1)*(e1cts2(i,k)-e1ctw2(i,k)))
          vtmp30(i,k)=tc(i,k)*((e1cts1(i,k) - e1ctw1(i,k))              &
     &                      -  (e1cts2(i,k) - e1ctw2(i,k)))
        enddo
      enddo

      do i = 1, IMAX
        tem = tc(i,LP1) * (e1cts1(i,LP1)-e1ctw1(i,LP1))
        flx1e1 (i) = tem * cldfac(i,LP1,1)
        flx1e10(i) = tem
      enddo

      do k = 1, L
        do i = 1, IMAX
          flx1e1 (i) = flx1e1 (i) + vtmp3 (i,k)
          flx1e10(i) = flx1e10(i) + vtmp30(i,k)
        enddo
      enddo

!  ---  now repeat flux calculations for the k=2..LM1  cases.
!       calculations for flux level L and LP1 are done separately, as all
!       emissivity and co2 calculations are special cases or nearby layers.

      do k = 2, LM1
        KLEN = k

        do kk = 1, LP1-k
          do i = 1, IMAX
            avephi(i,kk+k-1) = totphi(i,kk+k) - totphi(i,k)
          enddo
        enddo

        do i = 1, IMAX
          avephi(i,LP1) = avephi(i,LM1) + emx1(i)
        enddo

!  ---  compute emissivity fluxes (e2) for this case. note that we have
!       omitted the nearby later case (emiss(i,k,k)) as well as all cases
!       with k=L or LP1. but these cases have always been handled as
!       special cases, so we may as well compute their fluxes separastely.

        call e290                                                       &
!  ---  inputs:
     &     ( avephi, fxoe2, dte2, t1, t2, t4,                           &
     &       IMAX, L, LP1, KLEN,                                        &
!  ---  outputs:
     &       emiss, emissb                                              &
     &     )

        do kk = 1, LP1-k
          kkk = kk + k
          kk1 = kkk - 1

          do i = 1, IMAX
            avmo3 (i,kk1) = toto3 (i,kkk) - toto3 (i,k)
            avpho3(i,kk1) = tphio3(i,kkk) - tphio3(i,k)
            avvo2 (i,kk1) = totvo2(i,kkk) - totvo2(i,k)
            cont1d(i,kk1) = cnttau(i,kk1) * totevv(i,k-1)
          enddo
        enddo

        do kk = 1, LP1-k
          kkk = kk  + k
          kk1 = kkk - 1

          do i = 1, IMAX
            fac1 = bo3rnd * avpho3(i,kk1) / avmo3(i,kk1)
            vtmp = 0.5*( fac1*( sqrt(1.0                                &
     &           + (4.0*ao3rnd*avmo3(i,kk1))/fac1 ) - 1.0 ))
            to31d(i,kk1)  = exp( -1.0*( vtmp+sko3r*avvo2(i,kk1) ))
            over1d(i,kk1) = exp( -1.0*( sqrt(ab15wd*avephi(i,kk1))      &
     &                               +  skc1r*avvo2(i,kk1) ))
            co21(i,kkk,k) = over1d(i,kk1)*co21(i,kkk,k)
          enddo
        enddo

        do kp = k+1, LP1
          do i = 1, IMAX
            co21(i,k,kp) = over1d(i,kp-1)*co21(i,k,kp)
          enddo
        enddo

!  ---  rlog is the NBL amount for the 15 um band calculation

        do i = 1, IMAX
          rlog(i,k) = over1d(i,k)*co2nbl(i,k)
        enddo

!  ---  the kp terms for arbirrary k..

        do kp = k+1, LP1
          do i = 1, IMAX
            tem = oss(i,kp)*to31d(i,kp-1) + ss2(i,kp)*cont1d(i,kp-1)    &
     &          + css(i,kp)*co21 (i,kp,k) + dtc(i,kp)*emiss(i,kp-1)
            flx (i,k) = flx (i,k) + tem*cldfac(i,kp,k)
            flx0(i,k) = flx0(i,k) + tem
          enddo
        enddo

        do kp = k+1, LP1
          do i = 1, IMAX
            tem = oss(i,k)*to31d(i,kp-1) + ss2(i,k)*cont1d(i,kp-1)      &
     &          + css(i,k)*co21 (i,k,kp) + dtc(i,k)*emissb(i,kp-1)
            flx (i,kp) = flx (i,kp) + tem*cldfac(i,k,kp)
            flx0(i,kp) = flx0(i,kp) + tem
          enddo
        enddo

      enddo

!  ---  now do k=L case. since the kp loop is length 1, many simplifications
!       occur. also, the co2 quantities (as well as the emiss quantities)
!       are computed in the nbl sedction; therefore, we want only over, to3
!       and cont1d (over(i,L),to31d(i,L) and cont1d(i,L) according to the
!       notation. thus no call is made to the e290 subroutine.
!
!       the third section calculates boundary layer and nearby layer
!       corrections to the transmission functions obtained above. methods
!       are given in ref. (4).
!
!       the remaining calculations are for :
!          1) the (k,k) terms, k=2,LM1;        2) the (L,L) term;
!          3) the (L,LP1) term;                4) the (LP1,L) term;
!          5) the (LP1,LP1) term.
!
!       each is uniquely handled; different flux terms are computed
!       differently
!
!       fourth section obtains water transmission functions used in q(approx)
!       calculations and also makes nbl corrections:
!          1) emiss (i,j) is the transmission function matrix
!          2) "nearby layer" corrections (emiss(i,i)) are obtained using
!             subroutine e3v88 (inline coded)
!          3) special values at the surface (emiss(L,LP1),emiss(LP1,L),
!             emiss(LP1,LP1)) are calculated.

      do i = 1, IMAX
        tpl(i,1)    = tly(i,L)
        tpl(i,LP1)  = 0.5*(tlv(i,LP1) + tly(i,L))
        tpl(i,LLP1) = 0.5*(tlv(i,L)   + tly(i,L))

!  ---  e2 functions are required in the nbl calculations for 2 cases,
!       denoted (in old code) as (L,LP1) and (LP1,LP1)

        avephi(i,1) = var2(i,L)
        avephi(i,2) = var2(i,L) + empl(i,L)
      enddo

      do k = 2, L
        do i = 1, IMAX
          tpl(i,k)   = tlv(i,k)
          tpl(i,k+L) = tlv(i,k)
        enddo
      enddo

!  ---  subroutine e2spec computes the exchange terms in the flux equation
!       for longwave radiation for 2 terms used for nearby layer compu-
!       tations. the method is a table lookup on a pre-computed e2 function
!       (defined in ref. (4)).

      do  k = 1, 2
        do  i = 1, IMAX
          tmp3       = log10(avephi(i,k)) + 16.0
          fyo        = aint(tmp3*10.0)
          du         = tmp3 - 0.1*fyo
          ival       = 28.0*fyo + fxosp(i,k)
          emiss(i,k) = t1(ival)  + du*t2(ival) + dtsp(i,k)*t4(ival)
        enddo
      enddo

!  ---  compute nearby (NBL) layer transmissivities for h2o using a table
!       lookup of the pre-computed e3 function ( described in ref. (4)).

      do k = 1, LLP1
        do i = 1, IMAX
          fxo3 = aint(tpl(i,k)*0.1)
          tmp3 = log10(empl(i,k)) + 16.0
          dt3  = tpl(i,k) - 10.0*fxo3
          fyo  = aint(tmp3*10.0)
          du   = tmp3 - 0.1*fyo

!  ---  obtain index for table lookup; this value will have to be
!       decremented by 9 to account for table temps starting at 100k.

          it   = fxo3 + fyo*28.0
          ww1  = 10.0 - dt3
          ww2  = 0.1 - du
          tpl(i,k) = ww2 * (ww1*em3v(it-9)  + dt3*em3v(it-8))           &
     &             + du  * (ww1*em3v(it+19) + dt3*em3v(it+20))
        enddo
      enddo

!  ---  compute nearby layer and special-case transmissivities for emiss
!       using methods for h2o given in ref. (4)

      do k = 2, L
        do i = 1, IMAX
          emisdg(i,k) = tpl(i,k+L) + tpl(i,k)
        enddo
      enddo

!  ---  note that emx1/2 (pressure scaled paths) are now computed in lwrad

      do i = 1, IMAX
        emspec(i,1) = (tpl(i,1)*empl(i,1)-tpl(i,LP1)*empl(i,LP1))       &
     &              /  emx1(i) + 0.25*(emiss(i,1)+emiss(i,2))
        emisdg(i,LP1) = 2.0*tpl(i,LP1)
        emspec(i,2) = 2.0*(tpl(i,1)*empl(i,1)-tpl(i,LLP1)*empl(i,LLP1)) &
     &              / emx2(i)
      enddo

      do i = 1, IMAX
        fac1 = bo3rnd * var4(i,L) / var3(i,L)
        vtmp = 0.5*( fac1*( sqrt(1.0                                    &
     &       + (4.0*ao3rnd*var3(i,L))/fac1) - 1.0 ))
        to31d(i,L)  = exp( -1.0*( vtmp+sko3r*cntval(i,L) ))
        over1d(i,L) = exp( -1.0*( sqrt(ab15wd*var2(i,L))                &
     &                         +  skc1r*cntval(i,L) ))
        cont1d(i,L) = cnttau(i,L)*totevv(i,LM1)
        rlog(i,L)   = over1d(i,L)*co2nbl(i,L)
      enddo

      do k = 1, L
        k1 = k + 1

        do i = 1, IMAX
          rlog(i,k)    = log(rlog(i,k))
          delpr2(i,k1) = delp(i,k)*(plv(i,k1) - ply(i,k))
          tpl(i,k)     = -sqrt(delpr2(i,k1)) * rlog(i,k)
        enddo
      enddo

      do k = 1, LM1
        k1 = k + 1

        do i = 1, IMAX
          delpr1(i,k1) = delp(i,k1)*(ply(i,k1) - plv(i,k1))
          tpl(i,k+L)   = -sqrt(delpr1(i,k1))*rlog(i,k1)
        enddo
      enddo

      do i = 1, IMAX
        tpl(i,LL)   = -rlog(i,L)
        tpl(i,LLP1) = -rlog(i,L)                                        &
     &              * sqrt( delp(i,L)*(plv(i,LP1) - ply(i,LM1)) )
      enddo

!  ---  the first computation is for the 15 um band,with the for the
!       combined h2o and co2 transmission function.

      do k = 1, LLP1
        do i = 1, IMAX
          c(i,k)=tpl(i,k)*(-0.66667+tpl(i,k)*(0.25-0.066667*tpl(i,k)))
        enddo
      enddo

      do i = 1, IMAX
        co21(i,LP1,LP1) = 1.0 + c(i,L)
        co21(i,LP1,L)   = 1.0 + ( delp2(i,L)*c(i,LL)                    &
     &                  - (ply(i,L) - plv(i,L))*c(i,LLM1) )             &
     &                  / (plv(i,LP1) - ply(i,L))
        co21(i,L,LP1)   = 1.0 + ( (plv(i,LP1) - ply(i,LM1))*c(i,LLP1)   &
     &                  - (plv(i,LP1) - ply(i,L))*c(i,L) )              &
     &                  / (ply(i,L) - ply(i,LM1))
      enddo

      do k = 2, L
        do i = 1, IMAX
          co21(i,k,k) = 1.0 + 0.5*(c(i,LM1+k) + c(i,k-1))
        enddo
      enddo

!  ---  compute nearby-layer transmissivities for the o3 band and for the
!       one-band continuum band (to3 and emiss2). the sf2 function is
!       used. the method is the same as described for co2 in ref (4).

      do k = 1, LM1
        k1 = k + 1

        do i = 1, IMAX
          tpl(i,k1)  = cntval(i,k1) * delpr1(i,k1)
          tpl(i,k+L) = cntval(i,k)  * delpr2(i,k1)
        enddo
      enddo

!  ---  the sf2 function in prev. versions is now explicitly evaluated

      do k = 1, LLM1-1
        do i = 1, IMAX
          tem       = tpl(i,k+1)
          csub2     = sko3r*tem
          c(i,k+1)  = tem  *(-0.5 + tem  *(0.166666 - tem*0.0416666))
          c2(i,k+1) = csub2*(-0.5 + csub2*(0.166666 - csub2*0.0416666))
        enddo
      enddo

      do i = 1, IMAX
        contdg(i,LP1) = 1.0 + c (i,LLM1)
        to3dg (i,LP1) = 1.0 + c2(i,LLM1)
      enddo

      do k = 2, L
        do i = 1, IMAX
          contdg(i,k) = 1.0 + 0.5*(c (i,k) + c (i,LM1+k))
          to3dg (i,k) = 1.0 + 0.5*(c2(i,k) + c2(i,LM1+k))
        enddo
      enddo

!  ---  now obtain fluxes for the diagonal terms...

      do k = 2, LP1
        do i = 1, IMAX
          tem = dtc(i,k)*emisdg(i,k) + ss2(i,k)*contdg(i,k)             &
     &        + oss(i,k)*to3dg(i,k)  + css(i,k)*co21(i,k,k)
          flx (i,k) = flx (i,k) + tem*cldfac(i,k,k)
          flx0(i,k) = flx0(i,k) + tem
        enddo
      enddo

!  ---  for the two off-diagonal terms...

      do i = 1, IMAX
        tem = css(i,LP1)*co21(i,LP1,L) + dtc(i,LP1)*emspec(i,2)         &
     &      + oss(i,LP1)*to31d(i,L)    + ss2(i,LP1)*cont1d(i,L)
        flx (i,L) = flx (i,L) + tem*cldfac(i,LP1,L)
        flx0(i,L) = flx0(i,L) + tem

        tem = css(i,L)*co21(i,L,LP1) + oss(i,L)*to31d(i,L)              &
     &      + ss2(i,L)*cont1d(i,L)   + dtc(i,L)*emspec(i,1)
        flx (i,LP1) = flx (i,LP1) + tem*cldfac(i,L,LP1)
        flx0(i,LP1) = flx0(i,LP1) + tem
      enddo

!  ---  final section obtains emissivity heating rates, total heating
!       rates and the flux at the ground

      do k = 1, L
        do i = 1, IMAX
          vsum1 (i,k) = flx (i,k+1) - flx (i,k) - cts (i,k)             &
     &                - ctso3 (i,k) + excts (i,k)
          vsum10(i,k) = flx0(i,k+1) - flx0(i,k) - cts0(i,k)             &
     &                - ctso30(i,k) + excts0(i,k)

          heatra(i,k) = vsum1(i,k) * radcon * delp(i,k)

!        print *,' heatra=',heatra(i,k),' flx=', flx(i,k+1),flx(i,k)
!    &,' cts=',cts(i,k),' ctso3=',ctso3(i,k)
!    &,' excts=',excts(i,k),' k=',k,' tem=',radcon*delp(i,k)

        enddo
      enddo

!! ---  for optional clear sky
      if ( lhlw0 ) then
        do k = 1, L
        do i = 1, IMAX
          heatr0(i,k) = vsum10(i,k) * radcon * delp(i,k)
        enddo
        enddo
      endif

!  ---  calculate the flux at each flux level using the flux at the top
!       (flx1e1+gxcts) and the integral of the heating rates (vsum1)
!       a factor of 1.0e-3 is used to convert cgs unit to mks unit

      do i = 1, IMAX
        topflx(i)   = (flx1e1 (i) + gxcts (i)) * 1.0e-3
        flxnet(i,1) = topflx(i)
        topfl0(i)   = (flx1e10(i) + gxcts0(i)) * 1.0e-3
        flxne0(i,1) = topfl0(i)
      enddo

!  ---  only the surface value of flux (grnflx) is needed unless the
!       thick cloud section is invoked.

      do k = 1, L
        do i = 1, IMAX
          flxnet(i,k+1) = flxnet(i,k) + vsum1(i,k) *1.0e-3
          flxne0(i,k+1) = flxne0(i,k) + vsum10(i,k)*1.0e-3
        enddo
      enddo

      do i = 1, IMAX
        grnflx(i) = flxnet(i,LP1)
        grnfl0(i) = flxne0(i,LP1)
      enddo

!  ***************************************************************
!  *   thick cloud section no longer used ....k.a.c. sep96       *
!  ***************************************************************

!     if (lthick) then

!  ---  this is the thick cloud section.optionally,if thick cloud
!       fluxes are to be "convectively adjusted",ie,df/dp is constant,
!       for cloudy part of grid point, the following code is executed.

!  ---  first,count the number of clouds along the lat. row. skip the
!       entire thick cloud computation if there are no clouds.

!       icnt = 0
!       do i = 1, IMAX
!         icnt = icnt + nclds(i)
!       enddo

!       if (icnt /= 0) then

!  ---  find the maximum number of clouds in the latitude row

!         kclds = nclds(1)

!         do i = 2, IMAX
!           kclds = max(nclds(i), kclds)
!         enddo

!  ---  obtain the pressures and fluxes of the top and bottom of the
!       nc'th cloud (it is assumed that all ktop and kbtm's have
!       been defined!).

!         do kk = 1, kclds
!           kmin = LP1
!           kmax = 0

!           lab_do_i1 : do i = 1, IMAX
!             j1 = ktop(i,kk+1)
!             if (j1 == 1) exit lab_do_i1
!             j3 = kbtm(i,kk+1)
!             if (j3 > j1) then
!               ptop(i) = plv(i,j1)
!               pbot(i) = plv(i,j3+1)
!               ftop(i) = flxnet(i,j1)
!               fbot(i) = flxnet(i,j3+1)

!  ---  obtain the "flux derivative" df/dp (delptc)

!               delptc(i) = (ftop(i) - fbot(i))/(ptop(i) - pbot(i))
!               kmin = min(kmin, j1)
!               kmax = max(kmax, j3)
!             endif
!           enddo  lab_do_i1

!           kmin = kmin + 1

!  ---  calculate the tot. flux chg. from the top of the cloud, for
!       all levels.

!           do k = kmin, kmax
!             lab_do_i2 : do i = 1, IMAX
!               if (ktop(i,kk+1) == 1) exit lab_do_i2
!               if (ktop(i,kk+1) < k .and. k <= kbtm(i,kk+1)) then
!                 z1(i,k) = (plv(i,k) - ptop(i))*delptc(i) + ftop(i)
!orig             flxnet(i,k) = flxnet(i,k)*(1.0 - camt(i,kk+1))        &
!orig&                        + z1(i,k)*camt(i,kk+1)
!                 flxnet(i,k) = z1(i,k)
!               endif
!             enddo  lab_do_i2
!           enddo
!         enddo  ! end_do_kk

!  ---  using this flux chg. in the cloudy part of the grid box, obtain
!       the new fluxes, weighting the clear and cloudy fluxes:again, only
!       the fluxes in thick-cloud levels will eventually be used.

!         do k = 1, LP1
!           do i = 1, IMAX
!             flxnet(i,k) = flxnet(i,k)*(1.0 - camt(i,nc))              &
!    &                    + z1(i,k)*camt(i,nc)
!           enddo
!         enddo

!  ---  merge flxthk into flxnet for appropriate levels.

!         do k = 1, LP1
!           do i = 1, IMAX
!             if (k > itop(i) .and. k <= ibot(i)                        &
!    &                        .and. (nc-1) <= nclds(i))  then
!               flxnet(i,k) = flxthk(i,k)
!             endif
!           enddo
!         enddo

!  ---  end of cloud loop  ---

!       endif  ! end_if_icnt

!  ---  the final step is to recompute the heating rates based on the
!       revised fluxes:

!       do 6101 k = 1, L
!         do 6101 i = 1, IMAX
!           heatra(i,k) = radcon*(flxnet(i,k+1)-flxnet(i,k))*delp(i,k)
!         enddo
!       enddo

!     endif    ! end_if_lthick

!  ***************************************************************
!  *   thick cloud section no longer used ....k.a.c. sep96       *
!  ***************************************************************

      return
!...................................
      end subroutine fst88
!-----------------------------------



!-----------------------------------
      subroutine e290                                                   &
!...................................

!  ---  inputs:
     &     ( avephi, fxoe2, dte2, t1, t2, t4,                           &
     &       IMAX, L, LP1, KLEN,                                        &
!  ---  outputs:
     &       emiss, emissb                                              &
     &     )

!  *******************************************************************  !
!                                                                       !
!    purpose:  compute the exchange terms in the flux equation for      !
!              longwave radiation for all terms except the exchange     !
!              with the top of the atmosphere. the method is a table    !
!              lookup on a pre-computed e2 function (defined in         !
!              ref. (4)).                                               !
!                                                                       !
!    calculations are done in the frequency range:                      !
!              0-560,1200-2200 cm-1   for q(approx)                     !
!                                                                       !
!    motivation for these calculations is in references (1) and (4).    !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!                                                                       !
!                                                                       !
!                                                                       !
!  outputs:                                                             !
!     taucloud - cloud optical depth                        NBANDS*L    !
!                                                                       !
!                                                                       !
!  subroutine e290 is called by : fst88                                 !
!                                                                       !
!  subroutines called by e290 : none                                    !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, L, LP1, KLEN

      real (kind=kind_phys), dimension(:,:), intent(in) :: avephi,      &
     &       fxoe2, dte2

      real (kind=kind_phys), dimension(:),   intent(in) :: t1, t2, t4

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out):: emiss,emissb

!  ---  locals:
      real (kind=kind_phys), dimension(IMAX,LP1) :: dt, du, fyo

      real (kind=kind_phys)                      :: tmp3

      integer, dimension(IMAX,LP1) :: ival

      integer                      :: i, k, k1, item

!
!===> ... begin here
!

!  ---  first we obtain the emissivities as a function of water amount
!       (index fyo). this part of the code thus generates the e2 function.

      do k = 1, L-KLEN+2
        k1 = k + KLEN - 1

        do i = 1, IMAX
          tmp3       = log10( avephi(i,k1) ) + 16.0
          fyo  (i,k) = aint( tmp3*10.0 )
          du   (i,k) = tmp3  - 0.1*fyo(i,k)
          fyo  (i,k) = 28.0 * fyo(i,k)
          item       = fyo(i,k) + fxoe2(i,k1)
          emiss(i,k1)= t1(item) + du(i,k)*t2(item) + dte2(i,k1)*t4(item)
        enddo
      enddo

!  ---  the special case emiss(i,l) (layer kp) is obtained now
!       by averaging the values for l and lp1:
!       note that emiss(i,lp1) is not useful after this point.

      do i = 1, IMAX
        emiss(i,L) = 0.5 * (emiss(i,L) + emiss(i,LP1))
      enddo

!  ---  calculations for kp=KLEN and varying k; results are in emissb.
!       in this case, the temperature index is unchanged, always being
!       fxo(i,KLEN-1); the water index changes, but is symmetrical with
!       that for the varying kp case.note that the special case is not
!       involved here. (fixed level) k varies from (KLEN+1) to LP1;
!       results are in emissb(i,(KLEN) to l)

      do k = 1, LP1-KLEN
        do i = 1, IMAX
          dt  (i,k) = dte2(i,KLEN-1)
          ival(i,k) = fyo(i,k) + fxoe2(i,KLEN-1)
        enddo
      enddo

      do k = 1, LP1-KLEN
        k1 = k + KLEN - 1

        do i = 1, IMAX
          item = ival(i,k)
          emissb(i,k1) = t1(item) + du(i,k)*t2(item) + dt(i,k)*t4(item)
        enddo
      enddo

      return
!...................................
      end subroutine e290
!-----------------------------------



!-----------------------------------
      subroutine e1e290                                                 &
!...................................

!  ---  inputs:
     &     ( fxoe1,fxoe2,dte1,dte2,avephi,em1v,em1vw,t1,t2,t4,          &
     &       IMAX, L, LP1,                                              &
!  ---  outputs:
     &       g1,g2,g3,g4,g5,emiss                                       &
     &     )

!  *******************************************************************  !
!                                                                       !
!    purpose:  compute the exchange terms in the flux equation for all  !
!              terms except the exchange with the top of the atmosphere.!
!                                                                       !
!    the method is a table lookup on a pre-computed e2 function (defined!
!    in ref. (4)).                                                      !
!                                                                       !
!    the e1 function  calculations (formerly done in subroutine e1v88   !
!    compute the flux resulting from the exchange of photons between a  !
!    layer and the top of the atmosphere.  the method is a table lookup !
!    on a pre-computed e1 function.                                     !
!                                                                       !
!    calculations are done in two frequency ranges:                     !
!       1) 0-560,1200-2200 cm-1   for q(approx)                         !
!       2) 160-560 cm-1           for q(approx,cts)                     !
!                                                                       !
!    motivation for these calculations is in references (1) and (4).    !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!                                                                       !
!                                                                       !
!                                                                       !
!                                                                       !
!                                                                       !
!  outputs:                                                             !
!                                                                       !
!  subroutine e1e290 is called by : fst88                               !
!                                                                       !
!  subroutines called by e1e290 : none                                  !
!                                                                       !
!  *******************************************************************  !
!

      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, L, LP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: fxoe1,       &
     &       fxoe2, dte1, dte2, avephi

      real (kind=kind_phys), dimension(:),   intent(in) :: em1v,        &
     &       em1vw, t1, t2, t4

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: g1, g2,     &
     &       g3, g4, g5, emiss

!  ---  locals:
      real (kind=kind_phys), dimension(IMAX,LP1) :: fyo, du, ww1, ww2

      real (kind=kind_phys) :: tmp3, tem1, tem2, tem3, tem4

      integer, dimension(IMAX,3*L+2) :: it1

      integer :: i, k, k1, k2, ll, llp1, ival, item

!
!===> ... begin here
!
      ll   = L + L
      llp1 = ll + 1

!  ---  first we obtain the emissivities as a function of temperature
!       (index fxo) and water amount (index fyo). this part of the code
!       thus generates the e2 function. the fxo indices have been
!       obtained in fst88, for convenience.
!
!  ---  this subroutine evaluates the k=1 case only--

      do k = 1, LP1
        do i = 1, IMAX
          tmp3       = log10( avephi(i,k) ) + 16.0
          fyo  (i,k) = aint( tmp3*10.0 )
          du   (i,k) = tmp3 - 0.1*fyo(i,k)
          fyo  (i,k) = 28.0 * fyo(i,k)
          ival       = fyo(i,k) + fxoe2(i,k)
          emiss(i,k) = t1(ival) + du(i,k)*t2(ival) + dte2(i,k)*t4(ival)
        enddo
      enddo

!  ---  the special case emiss(i,L) is obtained now by averaging the
!       values for L and LP1:

      do i = 1, IMAX
        emiss(i,L) = 0.5*(emiss(i,L) + emiss(i,LP1))
      enddo

!  ---  calculations for the k=1 layer are not performed, as the
!       radiation code assumes that the top flux layer (above the
!       top data level) is isothermal, and hence contributes nothing
!       to the fluxes at other levels.
!
!  ---  the following is the calculation for the e1 function, formerly
!       done in subroutine e1v88. the move to e1e288 is due to the
!       savings in obtaining index values (the temp. indices have been
!       obtained in fst88, while the u-indices are obtained in the e2
!       calcs.,with k=1).
!
!  ---  for terms involving top layer, du is not known; in fact, we
!       use index 2 to repersent index 1 in prev. code. this means that
!       the it1 index 1 and llp1 has to be calculated separately. the
!       index llp2 gives the same value as 1; it can be omitted.

      do i = 1, IMAX
        it1(i,1) = fxoe1(i,1)
        ww1(i,1) = 10.0 - dte1(i,1)
        ww2(i,1) = 0.1
      enddo

      do k = 1, L
        k1 = k + 1
        k2 = k + LP1

        do i = 1, IMAX
          it1(i,k1) = fyo(i,k) + fxoe1(i,k1)
          it1(i,k2) = fyo(i,k) + fxoe1(i,k)
          ww1(i,k1) = 10.0 - dte1(i,k1)
          ww2(i,k1) = 0.1  - du(i,k)
        enddo
      enddo

      do k = 1, L
        do i = 1, IMAX
          it1(i,k+llp1) = fyo(i,k) + fxoe1(i,1)
        enddo
      enddo

!  ---  g3(i,1) has the same values as g1 (and did all along)

      do i = 1, IMAX
        tem1 = ww1(i,1) * ww2(i,1)
        tem2 = ww2(i,1) * dte1(i,1)
        item = it1(i,1)
        g1(i,1) = tem1 * em1v (item) + tem2 * em1v (item+1)
        g4(i,1) = tem1 * em1vw(item) + tem2 * em1vw(item+1)
        g3(i,1) = g1(i,1)
      enddo

      do k = 1, L
        k1 = k + 1

        do i = 1, IMAX
          tem1 = ww1(i,k1)  * ww2(i,k1)
          tem2 = ww2(i,k1)  * dte1(i,k1)
          tem3 = ww1(i,k1)  * du(i,k)
          tem4 = dte1(i,k1) * du(i,k)
          item = it1(i,k1)

          g1(i,k1) = tem1 * em1v(item)    + tem2 * em1v(item+1)         &
     &             + tem3 * em1v(item+28) + tem4 * em1v(item+29)
          g4(i,k1) = tem1 * em1vw(item)   + tem2 * em1vw(item+1)        &
     &             + tem3 * em1vw(item+28)+ tem4 * em1vw(item+29)

          tem1 = ww1(i,k)  * ww2(i,k1)
          tem2 = ww2(i,k1) * dte1(i,k)
          tem3 = ww1(i,k)  * du(i,k)
          tem4 = dte1(i,k) * du(i,k)
          item = it1(i,LP1+k)

          g2(i,k) = tem1 * em1v(item)    + tem2 * em1v(item+1)          &
     &            + tem3 * em1v(item+28) + tem4 * em1v(item+29)
          g5(i,k) = tem1 * em1vw(item)   + tem2 * em1vw(item+1)         &
     &            + tem3 * em1vw(item+28)+ tem4 * em1vw(item+29)
        enddo
      enddo

      do k = 2, LP1
        do i = 1, IMAX
          item = it1(i,ll+k)
          g3(i,k) = ww1(i,1) * ww2(i,k)  * em1v(item)                   &
     &            + ww2(i,k) * dte1(i,1) * em1v(item+1)                 &
     &            + ww1(i,1) * du(i,k-1) * em1v(item+28)                &
     &            + dte1(i,1)* du(i,k-1) * em1v(item+29)
        enddo
      enddo

      return
!...................................
      end subroutine e1e290
!-----------------------------------



!-----------------------------------
      subroutine spa88                                                  &
!...................................

!  ---  inputs:
     &     ( plv,ply,tly,var1,var2,delp,cldfac,                         &
     &       totvo2,to3sp,to3spc,co2sp,co2sp1,co2sp2,sorc,csour,        &
     &       IMAX, L, LP1,                                              &
!  ---  output:
     &       excts,ctso3,gxcts,excts0,ctso30,gxcts0                     &
     &     )

!  *******************************************************************  !
!                                                                       !
!    purpose:  compute layer net fluxes (heating rates) for exact cts   !
!              and 15um and 9.6um bands cts (o3+co2).                   !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  inputs:                                                              !
!                                                                       !
!                                                                       !
!                                                                       !
!                                                                       !
!  outputs:                                                             !
!                                                                       !
!  *******************************************************************  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX, L, LP1

      real (kind=kind_phys), dimension(:,:,:),intent(in) :: cldfac, sorc

      real (kind=kind_phys), dimension(:,:),  intent(in) :: plv, ply,   &
     &       tly, var1, var2, delp, totvo2, to3sp, to3spc, co2sp,       &
     &       co2sp1, co2sp2, csour

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:), intent(out) :: excts,ctso3,&
     &       excts0, ctso30

      real (kind=kind_phys), dimension(:),   intent(out) :: gxcts,gxcts0

!  ---  locals:
      real (kind=kind_phys), dimension(IMAX,LP1) :: ctmp, ctmp3, ctmp2, &
     &       ctmp0, ctmp30, ctmp20

      real (kind=kind_phys), dimension(IMAX,L)   :: phitmp, psitmp,     &
     &       tt, x, y, topm, topphi

      real (kind=kind_phys) :: f, ff, ag, agg, fac1, fac2, tem
      integer i, k, lm1, ib

!
!===> ... begin here
!
      lm1 = L - 1

!  ---  initiallization

      do i = 1, IMAX
        ctmp  (i,1) = 1.0
        ctmp2 (i,1) = 1.0
        ctmp3 (i,1) = 1.0
        gxcts (i)   = 0.0
        gxcts0(i)   = 0.0
        ctmp0 (i,1) = 1.0
        ctmp20(i,1) = 1.0
        ctmp30(i,1) = 1.0
      enddo

      do k = 1, L
        do i = 1, IMAX
          excts (i,k) =  0.0
          excts0(i,k) =  0.0

          x(i,k) = tly(i,k) - 250.0
          y(i,k) = x(i,k) * x(i,k)
        enddo
      enddo

!  --- ...  begin loop on frequency bands

      do ib = 1, 14

!  ---    obtain temperature correction (capphi,cappsi),then multiply
!         by optical path (var1,var2) to compute temperature-corrected
!         optical path and mean pressure for a layer (phitmp,psitmp)

        do k = 1, L
          do i = 1, IMAX
            f   = 0.044194 * (apcm (ib)*x(i,k) + bpcm (ib)*y(i,k))
            ff  = 0.044194 * (atpcm(ib)*x(i,k) + btpcm(ib)*y(i,k))
            ag  = (1.418191 + f )*f  + 1.0
            agg = (1.418191 + ff)*ff + 1.0

            ag  = ag * ag      !  ag ** 2
            ag  = ag * ag      !  ag ** 4
            ag  = ag * ag      !  ag ** 8
            agg = agg * agg
            agg = agg * agg
            agg = agg * agg

            phitmp(i,k) = var1(i,k) * (ag *ag )  ! ag ** 16
            psitmp(i,k) = var2(i,k) * (agg*agg)
          enddo
        enddo

!  ---  obtain optical path,mean pressure from the top to the pressure
!       plv(k) (topm,topphi)

        do i = 1, IMAX
          topm  (i,1) = phitmp(i,1)
          topphi(i,1) = psitmp(i,1)
        enddo

        do k = 2, L
          do i = 1, IMAX
            topm  (i,k) = topm  (i,k-1) + phitmp(i,k)
            topphi(i,k) = topphi(i,k-1) + psitmp(i,k)
          enddo
        enddo

!  ---  tt is the cloud-free cts transmission function

        if (ib < 5) then

          do k = 1, L
            do i = 1, IMAX
              fac1    = acomb(ib)*topm(i,k)
              fac2    = fac1*topm(i,k) / (bcomb(ib)*topphi(i,k))
              tt(i,k) = exp( -1.0*fac1 / sqrt(1.0+fac2) )
            enddo
          enddo

        elseif (ib < 9 .or. (ib > 10 .and. ib /= 13)) then

          do k = 1, L
            do i = 1, IMAX
              fac1    = acomb(ib)*topm(i,k)
              fac2    = fac1*topm(i,k) / (bcomb(ib)*topphi(i,k))
              tt(i,k) = exp( -1.0 * (fac1/sqrt(1.0+fac2)                &
     &                            +  betacm(ib)*totvo2(i,k+1)*sko2d) )
            enddo
          enddo

        elseif (ib == 9) then

          do k = 1, L
            do i = 1, IMAX
              fac1    = acomb(ib)*topm(i,k)
              fac2    = fac1*topm(i,k) / (bcomb(ib)*topphi(i,k))
              tt(i,k) = exp( -1.0 * (fac1/sqrt(1.0+fac2)                &
     &                            +  betacm(ib)*totvo2(i,k+1)*sko2d) )  &
     &                * co2sp1(i,k+1)
            enddo
          enddo

        elseif (ib == 10) then

          do k = 1, L
            do i = 1, IMAX
              fac1    = acomb(ib)*topm(i,k)
              fac2    = fac1*topm(i,k) / (bcomb(ib)*topphi(i,k))
              tt(i,k) = exp( -1.0 * (fac1/sqrt(1.0+fac2)                &
     &                            +  betacm(ib)*totvo2(i,k+1)*sko2d) )  &
     &                * co2sp2(i,k+1)
            enddo
          enddo

        elseif (ib == 13) then

          do k = 1, L
            do i = 1, IMAX
              fac1    = acomb(ib)*topm(i,k)
              fac2    = fac1*topm(i,k) / (bcomb(ib)*topphi(i,k))
              tt(i,k) = exp( -1.0 * (fac1/sqrt(1.0+fac2)                &
     &                            +  betacm(ib)*totvo2(i,k+1)*sko2d     &
     &                            +  to3spc(i,k)) )
            enddo
          enddo

        endif       ! end_if_ib

        do k = 1, L
          do i = 1, IMAX
            ctmp (i,k+1) = tt(i,k)*cldfac(i,k+1,1)
            ctmp0(i,k+1) = tt(i,k)
          enddo
        enddo

!  ---  excts is the net cts cooling flux accumulated over frequency bands

        do k = 1, L
          do i = 1, IMAX
            excts (i,k) = excts(i,k)                                    &
     &                  + sorc(i,k,ib)*(ctmp (i,k+1) - ctmp (i,k))
            excts0(i,k) = excts0(i,k)                                   &
     &                  + sorc(i,k,ib)*(ctmp0(i,k+1) - ctmp0(i,k))
          enddo
        enddo

!  ---  gxcts is the exact cts top flux accumulated over frequency bands

        do i = 1, IMAX
          tem = tt(i,L)*sorc(i,L,ib) + ( 0.5*delp(i,L)                  &
     &        * (tt(i,lm1)*(plv(i,LP1) - ply(i,L))                      &
     &        +  tt(i,L)  *(plv(i,LP1) + ply(i,L) - 2.0*plv(i,L))) )    &
     &        * (sorc(i,LP1,ib) - sorc(i,L,ib))
          gxcts (i) = gxcts (i) + tem * cldfac(i,LP1,1)
          gxcts0(i) = gxcts0(i) + tem
        enddo

      enddo      ! band loop ends here!


!  ---  obtain cts flux at the top by integration of heating rates and
!       using cts flux at the bottom (current value of gxcts). note
!       that the pressure quantities and conversion factors have not
!       been included either in excts or in gxcts. these cancel out, thus
!       reducing computations!

      do k = 1, L
        do i = 1, IMAX
          gxcts (i) = gxcts (i) - excts (i,k)
          gxcts0(i) = gxcts0(i) - excts0(i,k)
        enddo
      enddo

!  ---  compute approximate cts net flux for 15um and 9.6 um bands

      do k = 1, L
        do i = 1, IMAX
          ctmp2 (i,k+1) = co2sp(i,k+1) * cldfac(i,k+1,1)
          ctmp3 (i,k+1) = to3sp(i,k)   * cldfac(i,k+1,1)
          ctmp20(i,k+1) = co2sp(i,k+1)
          ctmp30(i,k+1) = to3sp(i,k)
        enddo
      enddo

      do k = 1, L
        do i = 1, IMAX
          ctso3 (i,k) = csour(i,k)  *(ctmp2 (i,k+1) - ctmp2 (i,k))      &
     &                + sorc(i,k,13)*(ctmp3 (i,k+1) - ctmp3 (i,k))
          ctso30(i,k) = csour(i,k)  *(ctmp20(i,k+1) - ctmp20(i,k))      &
     &                + sorc(i,k,13)*(ctmp30(i,k+1) - ctmp30(i,k))
        enddo
      enddo

      return
!...................................
      end subroutine spa88
!-----------------------------------


!
!........................................!
      end module module_radlw_main       !
!========================================!
