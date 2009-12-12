       subroutine gloopr
!*   &    ( grid_gr,
     &    ( grid_fld, g3d_fld,
     &     lats_nodes_r,global_lats_r, lonsperlar, phour,
     &     xlon,xlat,coszdg,COSZEN,
     &     SLMSK,SNWDPH,SNCOVR,SNOALB,ZORL,TSEA,HPRIME,SFALB,
     &     ALVSF,ALNSF ,ALVWF ,ALNWF,FACSF ,FACWF,CV,CVT ,
     &     CVB  ,SWH,HLW,SFCNSW,SFCDLW,
     &     FICE ,TISFC, SFCDSW,            
     &     TSFLW,FLUXR, phy_f3d,slag,sdec,cdec,NBLCK,KDT,
     &     global_times_r)

!! Code Revision:
!! Oct 11 2009       Sarah Lu, grid_gr is replaced by grid_fld
!! Oct 16 2009       Sarah Lu, grid_fld%tracers used
!! Dec 01 2009       Sarah Lu, update fcld (instant cloud cover) in addition
!!                             to cldcov (cumulative cloud cover)
!! Dec 09 2009       Sarah Lu, (1) g3d_fld added to calling argument; (2) grrad
!!                   returns instant cloud cover (cldcov_v); the accumulative 
!!                   and instant cloud cover fields are updated after grrad call
!! Dec 11 2009       Sarah Lu, ldiag3d removed from grrad calling argument
!!
cc
!#include "f_hpm.h"
!
      USE MACHINE              ,     ONLY : kind_phys,
     &                                      kind_grid,
     &                                      kind_evod
      USE FUNCPHYS             ,     ONLY : fpkap
      USE PHYSCONS, fv => con_fvirt, rerth => con_rerth 
      USE PHYSCONS, rk => con_rocp	! hmhj

      use module_radiation_driver,   only : radinit, grrad
      use module_radiation_astronomy,only : astronomy
      USE gfs_phy_tracer_config,     only : gfs_phy_tracer
!
!! ---  for optional spectral band heating outputs
!!    use module_radsw_parameters,   only : NBDSW
!!    use module_radlw_parameters,   only : NBDLW
!
      use resol_def,            ONLY: levs, levr, latr, lonr, lotgr,
     &                                g_t, g_p, g_q, g_dp, g_ps, 
     &                                ntcw, ntoz, ncld, num_p3d, 
     &                                nmtvr, ntrac, levp1, nfxr, g_dpdt,
     &                                lgocart
      use layout1,              ONLY: me, nodes, lats_node_r, 
     &                                lats_node_r_max, ipt_lats_node_r
      use gg_def,               ONLY: coslat_r, sinlat_r
      use date_def,             ONLY: idate
      use namelist_physics_def, ONLY: lsswr, iaer, lslwr, sashal, 
     &                                lssav, flgmin, ldiag3d,
     &                                iovr_lw, iovr_sw, isol, iems, 
     &                                ialb, fhlwr, fhswr, ico2, ngptc
      use d3d_def ,             ONLY: cldcov
      use gfs_physics_gridgr_mod, ONLY: Grid_Var_Data
      use gfs_physics_g3d_mod,    ONLY: G3D_Var_Data
!
      implicit none
!
      real (kind=kind_phys), parameter :: QMIN =1.0e-10
      real (kind=kind_phys), parameter :: Typical_pgr = 95.0
      real (kind=kind_phys), parameter :: cons0 = 0.0,  cons2 = 2.0
      real (kind=kind_phys), parameter :: pt01=0.01
!
!  --- ...  inputs:
      integer, intent(in) :: lats_nodes_r(nodes)
      integer, intent(in) :: global_lats_r(latr), lonsperlar(latr)

!*    real(kind=kind_grid) grid_gr(lonr*lats_node_r_max,lotgr)
      TYPE(Grid_Var_Data)       :: grid_fld 
      TYPE(G3D_Var_Data)        :: g3d_fld 

      integer, intent(in) :: NBLCK


      real (kind=kind_phys), dimension(LONR,LATS_NODE_R), intent(in) :: &
     &                       xlon, xlat, slmsk, snwdph, zorl, tsea,     &
     &                       alvsf, alnsf, alvwf, alnwf, facsf, facwf,  &
     &                       cv, cvt, cvb, FICE, tisfc, sncovr, snoalb

      real (kind=kind_phys), intent(in) ::                              &
     &                    hprime(NMTVR,LONR,LATS_NODE_R), phour,        &
     &                    phy_f3d(NGPTC,LEVS,NBLCK,LATS_NODE_R,NUM_P3D)
!

      real (kind=kind_phys), intent(inout) ::                           &
     &                    fluxr (NFXR,LONR,LATS_NODE_R)

      integer, intent(in) :: KDT
!  --- ...  outputs:
      real(kind=kind_evod), intent(out) ::                              &
     &                    global_times_r(latr,NODES)

      real (kind=kind_phys), intent(out) ::                             &
     &                    swh(NGPTC,LEVS,NBLCK,LATS_NODE_R),            &
     &                    hlw(NGPTC,LEVS,NBLCK,LATS_NODE_R)

      real (kind=kind_phys),dimension(LONR,LATS_NODE_R), intent(out) :: &
     &                    coszdg, coszen, sfcnsw, sfcdlw, tsflw,        &
     &                    sfcdsw, SFALB

      real (kind=kind_phys), intent(out) :: slag, sdec, cdec

!! --- ...  optional spectral band heating rates
!!    real (kind=kind_phys), optional, intent(out) ::                   &
!!   &                 htrswb(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R),      &
!!   &                 htrlwb(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R)

!  --- ...  locals:
      real(kind=kind_phys) :: prsl(NGPTC,LEVS),  prslk(NGPTC,LEVS),     &
     &                        prsi(NGPTC,LEVP1), prsik(NGPTC,LEVP1)

      real (kind=kind_phys) :: si_loc(LEVR+1)

      real (kind=kind_phys) ::                                          &
     &                        gt(NGPTC,LEVR),                           &
     &                        gr(NGPTC,LEVR), gr1(NGPTC,LEVR,NTRAC-1)

      real (kind=kind_phys) :: f_ice(NGPTC,LEVS), f_rain(NGPTC,LEVS),   &
     &                        r_rime(NGPTC,LEVS)

      real (kind=kind_phys) :: cldcov_v(NGPTC,LEVS), hprime_v(NGPTC),   &
     &                        fluxr_v(NGPTC,NFXR), vvel(NGPTC,LEVS)
      real (kind=kind_phys) :: flgmin_l(ngptc), work1, work2

      real (kind=kind_phys) :: rinc(5), dtsw, dtlw, solcon, raddt

      real (kind=kind_phys), save :: facoz

      integer :: njeff, lon, lan, lat, iblk, lons_lat, istrt
      integer :: idat(8), jdat(8), DAYS(13), iday, imon, midmon, id
      integer :: kk, jlonr, ilan

      integer, save :: icwp, k1oz, k2oz, midm, midp

!  ---  number of days in a month
      data DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /

!  --- ...  control parameters: 
!           (some of the them may be moved into model namelist)

!  ---  ISOL controls solar constant data source
!!    integer, parameter :: ISOL  = 0  ! use prescribed solar constant
!     integer, parameter :: ISOL  = 1  ! use varying solar const with 11-yr cycle

!  ---  ICO2 controls co2 data source for radiation
!     integer, parameter :: ICO2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ICO2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ICO2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  IALB controls surface albedo for sw radiation
!!    integer, parameter :: IALB = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: IALB = 1   ! use modis derived alb (to be developed)

!  ---  IEMS controls surface emissivity for lw radiation
!!    integer, parameter :: IEMS = 0   ! use fixed value of 1.0
!     integer, parameter :: IEMS = 1   ! use varying sfc emiss, based on sfc type
!  ---  IAER  controls aerosols scheme selections
!     integer, parameter :: IAER  = 1  ! opac climatology, without volc forcing
!     integer, parameter :: IAER  =11  ! opac climatology, with volcanic forcing
!     integer, parameter :: IAER  = 2  ! gocart prognostic, without volc forcing
!     integer, parameter :: IAER  =12  ! gocart prognostic, with volcanic forcing

!  ---  IOVR controls cloud overlapping method in radiation:
!     integer, parameter :: IOVR_SW = 0  ! sw: random overlap clouds
!!    integer, parameter :: IOVR_SW = 1  ! sw: max-random overlap clouds

!     integer, parameter :: IOVR_LW = 0  ! lw: random overlap clouds
!!    integer, parameter :: IOVR_LW = 1  ! lw: max-random overlap clouds

!  ---  iflip indicates model vertical index direction:
!     integer, parameter :: IFLIP = 0    ! virtical profile index from top to bottom
      integer, parameter :: IFLIP = 1    ! virtical profile index from bottom to top
!
!    The following parameters are from gbphys
!
      real (kind=kind_phys), parameter :: dxmax=-16.118095651,          &
     &                dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)

      integer :: ierr, dimg
      integer :: i, j, k, n

      logical :: lslag, change, lprnt
      data  lslag / .false. /,    lprnt / .false. /
      logical, save :: first
      data  first / .true. /

!  ---  timers:
      real*8 :: rtc, timer1, timer2
!
!===> *** ...  begin here
!
!!
      integer              kap,kar,kat,kau,kav,kdrlam
      integer              ksd,ksplam,kspphi,ksq,ksr,kst
      integer              ksu,ksv,ksz,node
!!
!     print *,' enter gloopr '
!
      idat = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc = 0.
! test repro
!     rinc(2) = fhour
      rinc(2) = phour
!     print *,' idate ',idate
!     print *,' rinc ',rinc
      call w3movdat(rinc, idat, jdat)
!     print *,' jdat ',jdat
!
      if (ntoz .le. 0) then                ! Climatological Ozone!
!
      if(me .eq. 0) WRITE (6,989) jdat(1),jdat(2),jdat(3),jdat(5)
  989 FORMAT(' UPDATING OZONE FOR ', I4,I3,I3,I3)
!
        IDAY   = jdat(3)
        IMON   = jdat(2)
        MIDMON = DAYS(IMON)/2 + 1
        CHANGE = FIRST .OR.
     &          ( (IDAY .EQ. MIDMON) .AND. (jdat(5).EQ.0) )
!
        IF (CHANGE) THEN
            IF (IDAY .LT. MIDMON) THEN
               K1OZ = MOD(IMON+10,12) + 1
               MIDM = DAYS(K1OZ)/2 + 1
               K2OZ = IMON
               MIDP = DAYS(K1OZ) + MIDMON
            ELSE
               K1OZ = IMON
               MIDM = MIDMON
               K2OZ = MOD(IMON,12) + 1
               MIDP = DAYS(K2OZ)/2 + 1 + DAYS(K1OZ)
            ENDIF
        ENDIF
!
        IF (IDAY .LT. MIDMON) THEN
           ID = IDAY + DAYS(K1OZ)
        ELSE
           ID = IDAY
        ENDIF
        FACOZ = real (ID-MIDM) / real (MIDP-MIDM)
      endif
!
      if (first) then
!
        si_loc(1)=1.0
        do k=1,levr-1
!*        si_loc(k+1)=si_loc(k)-grid_gr(1,g_dp+k-1)/grid_gr(1,g_ps)
          si_loc(k+1)=si_loc(k)-grid_fld%dp(1,1,k)/grid_fld%ps(1,1) 
        enddo
        si_loc(levr+1)=0.0

!  --- determin prognostic/diagnostic cloud scheme

        icwp   = 0
        if (NTCW > 0) icwp = 1
           
        first = .false.
           
      endif         ! end_if_first
!
!===> *** ...  radiation initialization
!
      dtsw  = 3600.0 * fhswr
      dtlw  = 3600.0 * fhlwr

      raddt = min(dtsw, dtlw)

      call radinit                                                      &
!  ---  input:
     &     ( si_loc, LEVR, IFLIP, NUM_P3D,                              &
     &       ISOL, ICO2, ICWP, IALB, IEMS, IAER, jdat, me )
!  ---  output: ( none )
                                                                                                            
!
!===> *** ...  astronomy for sw radiation calculation.
!
      call astronomy                                                    &
!  ---  inputs:
     &     ( lonsperlar, global_lats_r, sinlat_r, coslat_r, xlon,       &
     &       fhswr, jdat,                                               &
     &       LONR, LATS_NODE_R, LATR, IPT_LATS_NODE_R, lsswr, me,       &
!  ---  outputs:
     &       solcon, slag, sdec, cdec, coszen, coszdg                   &
     &      )

!
!===> *** ...  spectrum to grid transformation for radiation calculation.
!     -----------------------------------
cc
!     call f_hpmstart(61,"gr delnpe")
!     call f_hpmstop(61)
!     call f_hpmstart(62,"gr delnpo")
!     call f_hpmstop(62)
cc
cc
!     call f_hpmstart(63,"gr dezouv dozeuv")
!
!     call f_hpmstop(63)
cc
!     CALL countperf(0,5,0.)
!     CALL synctime()
!     CALL countperf(1,5,0.)
!!
!     CALL countperf(0,1,0.)
cc
!     call f_hpmstart(67,"gr sumfln")
!     call f_hpmstop(67)
cc
!     CALL countperf(1,1,0.)
cc
!     CALL countperf(0,1,0.)                                            ! hmhj
!
!     call f_hpmstart(68,"gr sumder2")                                  ! hmhj
!     call f_hpmstop(68)                                                ! hmhj
!
      CALL countperf(1,1,0.)                                            ! hmhj
!
!
!===> *** ...  starting latitude loop
!
      do lan=1,lats_node_r
cc
         lat = global_lats_r(ipt_lats_node_r-1+lan)
cc
         lons_lat = lonsperlar(lat)

         jlonr = (lan-1) * lonr

!!
!$omp parallel do schedule(dynamic,1) private(lon,j,k)
!$omp+private(istrt,njeff,iblk,n)
!$omp+private(vvel,gt,gr,gr1)
!$omp+private(cldcov_v,hprime_v,fluxr_v,f_ice,f_rain,r_rime)
!$omp+private(prslk,prsl,prsik,prsi)

        DO lon=1,lons_lat,NGPTC
!!
          NJEFF   = MIN(NGPTC,lons_lat-lon+1)
          ISTRT   = lon
          if (NGPTC.ne.1) then
            IBLK  = lon/NGPTC+1
          else
            IBLK  = lon
          endif
!

          do i = 1, njeff
!           ilan = jlonr + istrt+i-1
!           prsi(i,1) = grid_gr(ilan,g_ps)
            prsi(i,1) = grid_fld%ps(lon+i-1,lan)
          enddo
          do k = 1, LEVS
            do i = 1, njeff
!             ilan = jlonr + istrt+i-1
!             gt(i,k)    = grid_gr(ilan,g_t+k-1)
!             gr(i,k)    = max(qmin,grid_gr(ilan,g_q+k-1))
!             prsl(i,k)  = grid_gr(ilan,g_p+k-1)
!             vvel(i,k)  = grid_gr(ilan,g_dpdt+k-1)
!             prsi(i,k+1)= prsi(i,k)-                                      &
!     &                    grid_gr(ilan,g_dp+k-1)
              gt(i,k)    = grid_fld%t(lon+i-1,lan,k)                           
!*            gr(i,k)    = max(qmin,grid_fld%q(lon+i-1,lan,k))                 
              gr(i,k)= max(qmin,grid_fld%tracers(1)%flds(lon+i-1,lan,k))                 
              prsl(i,k)  = grid_fld%p(lon+i-1,lan,k)                         
              vvel(i,k)  = grid_fld%dpdt(lon+i-1,lan,k)                        
              prsi(i,k+1)= prsi(i,k)-                                      &  
     &                     grid_fld%dp(lon+i-1,lan,k)          
            enddo
          enddo
          do i = 1, njeff
            prsi (i,levs+1) = 0.0
            prsik(i,levs+1) = 0.0
          enddo
          do k = 1, levs
            do i = 1, njeff
              prslk(i,k) = (prsl(i,k)*pt01)**rk
              prsik(i,k) = (prsi(i,k)*pt01)**rk
            enddo
          enddo
!
!       Remaining tracers
!
         do n = 1, NTRAC-1
!           kk = g_q + n*levs
            do k = 1, LEVS
              do i = 1, njeff
!               ilan = jlonr + istrt+i-1
!               gr1(i,k,n) = grid_gr(ilan,kk+k-1)
!*              gr1(i,k,ntoz-1) = grid_fld%oz(lon+i-1,lan,k)  
!*              gr1(i,k,ntcw-1) = grid_fld%cld(lon+i-1,lan,k) 
                gr1(i,k,n) = grid_fld%tracers(n+1)%flds(lon+i-1,lan,k)
              enddo
            enddo
          enddo

!!
!.....
          if (levr .lt. levs) then
            do i=1,njeff
              prsi(i,levr+1) = prsi(i,levp1)
              prsl(i,levr)   = (prsi(i,levp1)+prsi(i,levr)) * 0.5
              prsik(i,levr+1) = prslk(i,levp1)
              prslk(i,levr)   = fpkap(prsl(i,levr)*1000.0)
            enddo
          endif
          do i=1,njeff
            hprime_v(i) = hprime(1,istrt+i-1,lan)
          enddo
!
          do k=1,nfxr
            do i=1,njeff
              fluxr_v(i,k) = fluxr(k,istrt+i-1,lan)
            enddo
          enddo
          if (NUM_P3D == 3) then
            do k = 1, LEVR
              do i = 1, njeff
                f_ice (i,k) = phy_f3d(i,k,iblk,lan,1)
                f_rain(i,k) = phy_f3d(i,k,iblk,lan,2)
                r_rime(i,k) = phy_f3d(i,k,iblk,lan,3)
              enddo
            enddo
          endif
          work1 = (log(coslat_r(lat)/(lons_lat*latr)) - dxmin) * dxinv
          work1 = max(0.0, min(1.0,work1))
          work2 = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
          do i=1,njeff
            flgmin_l(i) = work2
          enddo
 
!  *** ...  calling radiation driver
 
!
!     lprnt = me .eq. 0 .and. kdt .ge. 120
!     if (lprnt) then
!     if (kdt .gt. 85) then
!     print *,' calling grrad for me=',me,' lan=',lan,' lat=',lat
!    &,' num_p3d=',num_p3d,' snoalb=',snoalb(lon,lan),' lon=',lon
!    &,' tsea=',tsea(lon,lan),' sncovr=',sncovr(lon,lan),
!    &' snwdph=',snwdph(lon,lan)
!


          call grrad
!  ---  inputs:
     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             &
     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 &
     &       snwdph(lon,lan),sncovr(lon,lan),snoalb(lon,lan),           &
     &       zorl(lon,lan),hprime_v,                                    &
     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              &
     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              &
                                          ! fice FOR SEA-ICE XW Nov04
     &       fice(lon,lan),tisfc(lon,lan),                              &
     &       solcon,coszen(lon,lan),coszdg(lon,lan),k1oz,k2oz,facoz,    &
     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     &
     &       IOVR_SW,IOVR_LW,f_ice,f_rain,r_rime,flgmin_l,              &
     &       NUM_P3D,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,                   &
     &       dtlw,dtsw, lsswr,lslwr,lssav,sashal,                       &
     &       NGPTC,njeff,LEVR,IFLIP, me, lprnt,                         &
!  ---  outputs:
     &       swh(1,1,iblk,lan),sfcnsw(lon,lan),sfcdsw(lon,lan),         & ! sfcdsw FOR SEA-ICE XW Nov04
     &       sfalb(lon,lan),                                            & ! lu [+1L]: add sfalb
     &       hlw(1,1,iblk,lan),sfcdlw(lon,lan),tsflw(lon,lan),          & 
     &       cldcov_v,                                                  & ! return instant cloud cover
!  ---  input/output:
     &       fluxr_v                                                    & 
!! ---  optional outputs:
!!   &,      HTRSWB=htrswb(1,1,1,iblk,lan),                             &
!!   &,      HTRLWB=htrlwb(1,1,1,iblk,lan)                              &
     &     )
!
! grrad routine computes cldcov_v (instant 3D cloud cover)    -- Sarah Lu
! if ldiag3d is T, update cldcov (accumulative 3D cloud cover)
! if lgocart is T, update fcld (instant 3D cloud cover)

          if (ldiag3d) then
            do k=1,levr
              do i=1,njeff
                cldcov(k,istrt+i-1,lan) = cldcov(k,istrt+i-1,lan) +     &
     &                                    cldcov_v(i,k) * raddt
              enddo
            enddo
          endif
          if (lgocart) then
            do k=1,levr
              do i=1,njeff
                g3d_fld%fcld(istrt+i-1,lan,k) = cldcov_v(i,k) 
              enddo
            enddo
          endif

          do k=1,nfxr
            do i=1,njeff
              fluxr(k,istrt+i-1,lan) = fluxr_v(i,k)
            enddo
          enddo
          if (levr .lt. levs) then
            do k=levr+1,levs
              do i=1,njeff
                hlw(i,k,iblk,lan) = hlw(i,levr,iblk,lan)
                swh(i,k,iblk,lan) = swh(i,levr,iblk,lan)
              enddo
            enddo
          endif
 
c$$$          write(2900+lat,*) ' ilon = ',istrt
c$$$          write(2900+lat,'("swh",T16,"hlw")')
c$$$      do k=1,levs
c$$$         write(2900+lat,
c$$$     .         '(e10.3,T16,e10.3,T31,e10.3)')
c$$$     . swh(1,k,iblk,lan),hlw(1,k,iblk,lan)
c$$$       enddo
 
!!
          CALL countperf(1,12,0.)
          ENDDO
!
      enddo
cc
      call f_hpmstop(69)
!!
      CALL countperf(0,5,0.)
      CALL synctime()
      CALL countperf(1,5,0.)
!sela print*,'completed gloopr_v kdt=',kdt
!     print *,' end of gloopr '
!!
      return
      end subroutine gloopr
