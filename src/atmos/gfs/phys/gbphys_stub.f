      SUBROUTINE GBPHYS(IM,IX,levs,lsoil,lsm,ntrac,ncld &
     , ntoz,ntcw,nmtvr,lonr,latr,jcap,ras,nlons,xkt2,nrcm,pre_rad &
     , UGRS,VGRS,PGR,TGRS,QGRS,vvel &
     , GT0,GQ0,GU0,GV0,sinlat,coslat,rcs2 &
     , prsi,prsl,prslk,prsik,phii,phil,dpshc,fhour,lssav,solhr &
     , lsfwd,clstp,dtp,dtf,poz,prdout,ko3,pl_coeff &
     , HICE,FICE,TISFC,SFCDSW &                         ! FOR SEA-ICE - XW Nov04
     , TPRCP, SRFLAG &						! lu add
     , SLC   ,SNWDPH,SLOPE ,SHDMIN,SHDMAX,SNOALB,SFALB &	! lu add
     , CHH,CMM,EPI,DLWSFCI,ULWSFCI,USWSFCI,DSWSFCI,DTSFCI &
     , DQSFCI,GFLUXI,SRUNOFF,T1,Q1,U1,V1,ZLVL,EVBSA,EVCWA &
     , TRANSA,SBSNOA,SNOWCA,SOILM &
     , RAIN,RAINC,WET1   &                                      ! added for gocart
     , TSEA  ,SHELEG,SNCOVR, TG3   ,ZORL  ,CV    ,CVB   ,CVT   &
     , SLMSK ,VFRAC ,CANOPY,F10M  ,VTYPE ,STYPE ,UUSTAR,FFMM  ,FFHH   &
     , TMPMIN,TMPMAX,GESHEM,DUSFC ,DVSFC ,DTSFC ,DQSFC ,DLWSFC,ULWSFC &
     , GFLUX ,RUNOFF,EP    ,CLDWRK,DUGWD ,DVGWD ,PSMEAN,BENGSH,XLON   &
     , COSZEN,SFCNSW,XLAT & 
     , SFCDLW,TSFLW ,PSURF ,U10M  ,V10M  ,T2M   ,Q2M   &
     , HPBL  ,PWAT  ,SWH,HLW,SMC,STC,HPRIME,slag,sdec,cdec &
     , acv,acvb,acvt &
     , phy_f3d, phy_f2d, num_p3d, num_p2d, flgmin &
     , DT3DT, DQ3DT, DU3DT, DV3DT, DQDT_V                          &
     , upd_mf, dwn_mf, det_mf, LDIAG3D                             &
     , flipv, me,kdt,lat,oro &
     , crtrh, ncw, old_monin,cnvgwd,ccwf, sashal,newsas )

! Coupling insertion->
!<-Coupling insertion
!
      implicit none
!
      integer  levs,lsoil,lsm,ix,im,ntrac,ncld,ntoz,ntcw,nmtvr,lonr &
     ,         latr,jcap,nlons(im),num_p3d,num_p2d,nrcm,lat &
     ,        pl_coeff
!
      LOGICAL lssav,lsfwd, old_monin, cnvgwd
      integer levshc(im), levshcm     ! Needed for pry version
      integer ncw(2)
!     real(kind=8) dtp,dtf,FHOUR,solhr, prsshc
      real(kind=8) dtp,dtf,FHOUR,solhr, dpshc(im), crtrh(3)
      real(kind=8) flgmin(2), ccwf
      real, parameter :: fhourpr=0.0



      real(kind=8) UGRS(IX,LEVS),      VGRS(IX,LEVS) &
     ,                     TGRS(IX,LEVS),      qgrs(IX,levs,ntrac) &
     ,                     VVEL(IX,LEVS) &
!
     ,                     GT0(IX,LEVS),       GU0(IX,LEVS) &
     ,                     GV0(IX,LEVS),       gq0(IX,levs,ntrac) &
!
     ,                     DEL(IX,LEVS),       PRSI(IX,LEVS+1) &
     ,                     PRSL(IX,LEVS),      PRSLK(IX,LEVS) &
     ,                     PRSIK(IX,LEVS+1),   PHII(IX,LEVS+1) &
     ,                     PHIL(IX,LEVS) &
     ,                     PGR(IM),            XKT2(IX,nrcm) &
     ,                    ccwfac(im)

      real(kind=8) RCS2(IM), SINLAT(IM), COSLAT(IM),clstp

      real(kind=8) SMC(IX,LSOIL), STC(IX,LSOIL),    SWH(IX,LEVS) &
     ,                    HICE(IM),      FICE(IM),         SFCDSW(IM) & ! SEA-ICE
     ,                    TISFC(IM) &
     ,                    HLW(IX,LEVS),  HPRIME(IX,NMTVR), TSEA(IM) &
     ,                    SHELEG(IM),    TG3(IM),          ZORL(IM) &
     ,                    SNCOVR(IM) &
     ,                    CV(IM),        CVB(IM),          CVT(IM) &
     ,                    COSZEN(IM),    PWAT(IM),         SLMSK(IM) &
     ,                    VFRAC(IM),     CANOPY(IM),       F10M(IM) &
     ,                    VTYPE(IM),     STYPE(IM),        UUSTAR(IM) &
     ,                    FFMM(IM),      FFHH(IM),         TMPMIN(IM) &
     ,                    TMPMAX(IM),    GESHEM(IM),       DUSFC(IM) &
     ,                    DVSFC(IM),     DTSFC(IM),        DQSFC(IM) &
     ,                    DLWSFC(IM),    ULWSFC(IM),       GFLUX(IM) &
     ,                    RUNOFF(IM),    EP(IM),           CLDWRK(IM) &
     ,                    DUGWD(IM),     DVGWD(IM),        PSMEAN(IM) &
     ,                    BENGSH(IM),    XLON(IM),         SFCNSW(IM) &
     ,                    SFCDLW(IM),    TSFLW(IM),        PSURF(IM) &
     ,                    U10M(IM),      V10M(IM),         T2M(IM) &
     ,                    Q2M(IM),       HPBL(IM),         xlat(IM) &
     ,                    TPRCP(IM),     SRFLAG(IM) &
     ,                    SLC(IX,LSOIL) &
     ,                    SNWDPH(IM),    SHDMIN(IM),      SHDMAX(IM) &
     ,                    SNOALB(IM),    SLOPE(IM),       SFALB(IM) &
     , CHH(IM),CMM(IM),EPI(IM),DLWSFCI(IM),ULWSFCI(IM),USWSFCI(IM) &
     , DSWSFCI(IM),DTSFCI(IM),DQSFCI(IM),GFLUXI(IM),SRUNOFF(IM),T1(IM) &
     , Q1(IM),U1(IM),V1(IM),ZLVL(IM) &
     , EVBSA(IM),EVCWA(IM),TRANSA(IM),SBSNOA(IM),SNOWCA(IM),SOILM(IM) &
     , RAIN(IM), RAINC(IM),WET1(IM)                                   &
     ,                    phy_f3d(IX,LEVS,num_p3d), phy_f2d(IX,num_p2d) &
     ,                    acv(IM),       acvb(IM), acvt(IM) &
     ,                    oro(im)
      real(kind=8) slag,sdec,cdec
      real(kind=8) dt3dt(IX,levs,6),  dq3dt(IX,levs,5+pl_coeff) &
     ,                     du3dt(IX,levs,4),  dv3dt(IX,levs,4) &
!    ,                    cumflx(ix,levs) &
     ,                    upd_mf(ix,levs), dwn_mf(ix,levs) &
     ,                    det_mf(ix,levs)

! add local working array for total moisture tendency
      real(kind=8) dqdt_v(IX,levs)
!
      integer me, kdt
      logical RAS,LDIAG3D,pre_rad,sashal,newsas
      real(kind=8) CLW(IX,LEVS,2) &
     ,                    garea(im), dlength(im)
!
      integer KO3
      real(kind=8) poz(KO3), prdout(IX,ko3,pl_coeff)
!
      real(kind=8) RHC(IM,LEVS), SR(IM,levs) &
     ,                    xncw(IM),  rhbbot, rhbtop, rhpbl &
     ,                    dxmax,     dxmin,  dxinv &
     ,                    ud_mf(ix,levs), dd_mf(ix,levs) &
     ,                    dt_mf(ix,levs)
      logical flipv
!
      RETURN
      END

