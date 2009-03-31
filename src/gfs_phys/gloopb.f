      subroutine gloopb
     &    ( grid_gr,
     x     lats_nodes_r,global_lats_r,lonsperlar,
     &     tstep,phour,sfc_fld, flx_fld, SFALB,xlon,
     &     swh,hlw,hprime,slag,sdec,cdec,
     &     ozplin,jindx1,jindx2,ddy,
     &     phy_f3d, phy_f2d,xlat,nblck,kdt,
     &     global_times_b)
!!
! #include "f_hpm.h"
!!
      use resol_def
      use layout1
      use gg_def
      use vert_def
      use date_def
      use namelist_physics_def
      use coordinate_def                                                ! hmhj
      use module_ras , only : ras_init
      use physcons, grav => con_g , rerth => con_rerth, rk => con_rocp  ! hmhj
      use ozne_def
      use d3d_def
      use gfs_physics_sfc_flx_mod
      use mersenne_twister
      include '../../inc/mpif.h'
      implicit none
!
      real(kind=kind_grid) grid_gr(lonr*lats_node_r_max,lotgr)
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld

      integer nlons_v(ngptc)
!
      integer id,njeff,istrt,lon,iblk,kdt
!!
      integer nblck
!!
      real(kind=kind_phys)    phour
      real(kind=kind_phys)    prsl(ngptc,levs)
      real(kind=kind_phys)   prslk(ngptc,levs), dpshc(ngptc)
      real(kind=kind_phys)    prsi(ngptc,levs+1),phii(ngptc,levs+1)
      real(kind=kind_phys)   prsik(ngptc,levs+1),phil(ngptc,levs)
!!
      real (kind=kind_rad) gu(ngptc,levs),gv(ngptc,levs)
      real (kind=kind_rad) gq(ngptc),gt(ngptc,levs)
      real (kind=kind_rad) gr(ngptc,levs,ntrac)
      real (kind=kind_rad) vvel(ngptc,levs),rand
      real (kind=kind_rad) adt(ngptc,levs),adr(ngptc,levs,ntrac)
      real (kind=kind_rad) adu(ngptc,levs),adv(ngptc,levs)
!!
      real (kind=kind_rad) xlon(lonr,lats_node_r)
      real (kind=kind_rad) xlat(lonr,lats_node_r)
      real (kind=kind_rad) 
     &                     hprime(nmtvr,lonr,lats_node_r),
     &                     fluxr(nfxr,lonr,lats_node_r),
     &                     sfalb(lonr,lats_node_r)
      real (kind=kind_rad)  swh(ngptc,levs,nblck,lats_node_r)
      real (kind=kind_rad)  hlw(ngptc,levs,nblck,lats_node_r)
!!
      real  (kind=kind_phys)
     &     phy_f3d(ngptc,levs,nblck,lats_node_r,num_p3d),
     &     phy_f2d(lonr,lats_node_r,num_p2d)
!!
      real (kind=kind_phys) dtphys,dtp,dtf
      real (kind=kind_evod) tstep
!!
      integer              lats_nodes_r(nodes)
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
cc
      integer              i,j,k,kk,n
      integer              l,lan,lat,jlonr,ilan
      integer              lon_dim,lons_lat
      integer              nsphys
!
      real(kind=kind_evod) solhr,clstp
cc
!timers______________________________________________________---
 
      real*8 rtc ,timer1,timer2
      real(kind=kind_evod) global_times_b(latr,nodes)
 
!timers______________________________________________________---
cc
      real(kind=kind_phys) cons0,cons2     !constant
      logical, parameter :: flipv = .true.
      real(kind=kind_phys), parameter :: pt01=0.01
cc
c for nasa ozone production and distruction rates:(input throu fixio)
      real ozplin(latsozp,levozp,pl_coeff,timeoz)
      integer jindx1(lats_node_r),jindx2(lats_node_r)!for ozone interpolaton
      real ddy(lats_node_r)                          !for ozone interpolaton
      real ozplout(levozp,lats_node_r,pl_coeff),
     &     ozplout_v(ngptc,levozp,pl_coeff)
!!
      real(kind=kind_phys), allocatable :: acv(:,:),acvb(:,:),acvt(:,:)
      save acv,acvb,acvt
!!
      real (kind=kind_phys) rannum(lonr*latr*nrcm)
!     real (kind=kind_phys) rannum(lonr,latr,nrcm)
     &,                     xkt2(lonr,lats_node_r,nrcm)
      integer iseed, nrc, seed0
      integer nf0,nf1,ind,nt,indod,indev
      real(kind=kind_evod) fd2
      real(kind=kind_evod) wrk(1)
      logical first
      data first/.true./
!     save    krsize, first, nrnd,seed0
      save    first, seed0
!
      real(kind=kind_phys), parameter :: cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9
!
! add more

      real (kind=kind_evod) slag,sdec,cdec

!     real (kind=kind_phys) dt3dt(ngptc,levs,6,nblck,lats_node_r),
!    &  dq3dt(ngptc,levs,5+pl_coeff,nblck,lats_node_r),
!    &  du3dt(ngptc,levs,3,nblck,lats_node_r),
!    &  dv3dt(ngptc,levs,3,nblck,lats_node_r)

      real(kind=kind_evod) sinlat_v(ngptc),coslat_v(ngptc),rcs2_v(ngptc)
      real(kind=kind_evod) smc_v(ngptc,lsoil),stc_v(ngptc,lsoil)
      real(kind=kind_evod) slc_v(ngptc,lsoil)
      real(kind=kind_evod) hprime_v(ngptc,nmtvr)
      real(kind=kind_evod) phy_f3dv(ngptc,LEVS,num_p3d)
      real(kind=kind_evod) phy_f2dv(ngptc,num_p2d)
      real(kind=kind_evod) rannum_v(ngptc,nrcm)

!
      if (first) then
!
!       call random_seed(size=krsize)
!       if (me.eq.0) print *,' krsize=',krsize
!       allocate (nrnd(krsize))

        allocate (acv(lonr,lats_node_r))
        allocate (acvb(lonr,lats_node_r))
        allocate (acvt(lonr,lats_node_r))
!
        seed0 = idate(1) + idate(2) + idate(3) + idate(4)

!       nrnd  = idate(1) + idate(3) * 24
!       call random_seed(generator=2)
!       call random_seed(put=nrnd)

        call random_setseed(seed0)
        call random_number(wrk)
        seed0 = seed0 + nint(wrk(1)*1000.0)
!
        if (me .eq. 0) then
          print *,' seed0=',seed0,' idate=',idate,' wrk=',wrk
          if (num_p3d .eq. 3) print *,' USING Ferrier-MICROPHYSICS'
          if (num_p3d .eq. 4) print *,' USING ZHAO-MICROPHYSICS'
        endif
        if (fhour .eq. 0.0) then
          do j=1,lats_node_r
            do i=1,lonr
              phy_f2d(i,j,num_p2d) = 0.0
            enddo
          enddo
        endif
       
        if (ras) call ras_init(levs, me)
       
        first = .false.

      endif

cc
      dtphys=3600.
      nsphys=max(int(2*tstep/dtphys+0.9999),1)
      dtp=2*tstep/nsphys
      dtf=0.5*dtp
      if(lsfwd) dtf=dtp
c
      solhr=mod(phour+idate(1),cons_24)
c...  set switch for saving convective clouds
      if(lscca.and.lsswr) then
        clstp=1100+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
      elseif(lscca) then
        clstp=0100+min(fhswr,fhour,cons_99)  !accumulate,convert
      elseif(lsswr) then
        clstp=1100                       !initialize,accumulate
      else
        clstp=0100                       !accumulate
      endif
!
!
      iseed = mod(100.0*sqrt(fhour*3600),cons_1p0d9) + 1 + seed0

!     nrnd  = iseed
!     call random_seed(generator=2)
!     call random_seed(put=nrnd)

      call random_setseed(iseed)
      call random_number(rannum)
      do nrc=1,nrcm
        do j=1,lats_node_r
           lat=global_lats_r(ipt_lats_node_r-1+j)
           do i=1,lonr
!             xkt2(i,j,nrc) = rannum(i,lat,nrc)
              xkt2(i,j,nrc) = rannum(i+(lat-1)*lonr+(nrc-1)*latr)
           enddo
        enddo
      enddo
!
      if(.not.random_xkt2)then
        if(kdt.lt.3) print*,'random_xkt2=',random_xkt2
        xkt2 = 0.6
      endif
!     print *,' xkt2=',xkt2(1:5,1:5,1),' kdt=',kdt
!     if (kdt .eq. 1) print *,' xkt2=',xkt2(1:5,1:5,1)
!
! doing ozone i/o and latitudinal interpolation to local gauss lats
c      ifozphys=.true.
 
      if (ntoz .gt. 0) then
       call ozinterpol(me,lats_node_r,lats_node_r,idate,fhour,
     &                 jindx1,jindx2,ozplin,ozplout,ddy)
      endif

!!
c ----------------------------------------------------
c
!
      do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lon_dim = lon_dims_r(lan)
         lons_lat = lonsperlar(lat)
         jlonr = (lan-1)*lonr

!$omp parallel do  schedule(dynamic,1) private(lon)
!$omp+private(njeff,istrt,iblk)
!$omp+private(ilan, prsi, gq, gu, gv, gt, gr, prsl)
!$omp+private(prsik, prslk, kk, vvel, phil, dpshc, ozplout_v)
!$omp+private(smc_v, stc_v, slc_v, hprime_v, phy_f3dv, phy_f2dv)
!$omp+private(rannum_v, adu, adv, adt, adr)
!$omp+private(nlons_v, sinlat_v, coslat_v, rcs2_v, phii, pl_pres)

        do lon=1,lons_lat,ngptc
!!
          njeff=min(ngptc,lons_lat-lon+1)
          istrt=lon
          if (ngptc.ne.1) then
            iblk=lon/ngptc+1
          else
            iblk=lon
          endif
!!
          do i = 1, njeff
            ilan=jlonr+istrt+i-1
            prsi(i,1) = grid_gr(ilan,g_ps)
            gq(i)     = prsi(i,1)
          enddo
          do k = 1, LEVS
            do i = 1, njeff
              ilan=jlonr+istrt+i-1
              gu(i,k)    = grid_gr(ilan,g_u+k-1)
              gv(i,k)    = grid_gr(ilan,g_v+k-1)
              gt(i,k)    = grid_gr(ilan,g_t+k-1)
              prsl(i,k)  = grid_gr(ilan,g_p+k-1)
              vvel(i,k)  = grid_gr(ilan,g_dpdt+k-1)
              prsi(i,k+1)= prsi(i,k)-
     &                     grid_gr(ilan,g_dp+k-1)
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
          do n = 1, NTRAC
            kk = g_q + (n-1)*levs
            do k = 1, LEVS
              do i = 1, njeff
                ilan=jlonr+istrt+i-1
                gr(i,k,n) = grid_gr(ilan,kk+k-1)
              enddo
            enddo
          enddo

      do i=1,ngptc
        phil(i,levs)  = 0.0 ! will force calculation of geopotential in gbphys.
      enddo
      if (gen_coord_hybrid .and. thermodyn_id == 3) then
        do i=1,ngptc
          prslk(i,1) = 0.0 ! will force calculation of geopotential in gbphys.
          prsik(i,1) = 0.0 ! will force calculation of geopotential in gbphys.
        enddo
      endif
      do i=1,njeff
        dpshc(i) = 0.3 * prsi(i,1)
      enddo
      nlons_v  = lons_lat
      sinlat_v = sinlat_r(lat)
      coslat_v = coslat_r(lat)
      rcs2_v   = rcs2_r(min(lat,latr-lat+1))

      if (ntoz .gt. 0) then
        do j=1,pl_coeff
          do k=1,levozp
            do i=1,ngptc
              ozplout_v(i,k,j) = ozplout(k,lan,j)
            enddo
          enddo
        enddo
      endif

      do k=1,lsoil
        do i=1,njeff
          smc_v(i,k) = sfc_fld%smc(k,istrt+i-1,lan)
          stc_v(i,k) = sfc_fld%stc(k,istrt+i-1,lan)
          slc_v(i,k) = sfc_fld%slc(k,istrt+i-1,lan)
        enddo
      enddo
      do k=1,nmtvr
        do i=1,njeff
          hprime_v(i,k) = hprime(k,istrt+i-1,lan)
        enddo
      enddo
!!
      do j=1,num_p3d
        do k=1,levs
          do i=1,njeff
            phy_f3dv(i,k,j) = phy_f3d(i,k,iblk,lan,j)
          enddo
        enddo
      enddo
      do j=1,num_p2d
        do i=1,njeff
          phy_f2dv(i,j) = phy_f2d(istrt+i-1,lan,j)
        enddo
      enddo
      do j=1,nrcm
        do i=1,njeff
          rannum_v(i,j) = xkt2(istrt+i-1,lan,j)
        enddo
      enddo

!
            call gbphys(njeff,ngptc,levs,lsoil,lsm,ntrac,ncld,
     &      ntoz,ntcw,nmtvr,lonr,latr,jcap,ras,nlons_v,rannum_v,nrcm,
     &      pre_rad,
     &      gu,gv,gq,gt,gr,vvel,
     &      adt,adr,adu,adv,
     &      sinlat_v,coslat_v,rcs2_v,
     &      prsi,prsl,prslk,prsik,phii,phil,dpshc,fhour,lssav,solhr,
     &      lsfwd,clstp,dtp,dtf,
     &      pl_pres,ozplout_v,levozp,pl_coeff,
     &      sfc_fld%hice(istrt,lan),sfc_fld%fice(istrt,lan),
     &      sfc_fld%tisfc(istrt,lan),
     &      flx_fld%sfcdsw(istrt,lan),                      ! SEA-ICE - Nov04
     +      sfc_fld%tprcp(istrt,lan),sfc_fld%srflag(istrt,lan),
     +      slc_v,sfc_fld%snwdph(istrt,lan),sfc_fld%slope(istrt,lan),
     &      sfc_fld%shdmin(istrt,lan),
     +      sfc_fld%shdmax(istrt,lan),sfc_fld%snoalb(istrt,lan),
     &      sfalb(istrt,lan),
     +     flx_fld%chh(istrt,lan),flx_fld%cmm(istrt,lan),
     +     flx_fld%epi(istrt,lan),flx_fld%dlwsfci(istrt,lan),
     +     flx_fld%ulwsfci(istrt,lan),flx_fld%uswsfci(istrt,lan),
     +     flx_fld%dswsfci(istrt,lan),flx_fld%dtsfci(istrt,lan),
     +     flx_fld%dqsfci(istrt,lan),flx_fld%gfluxi(istrt,lan),
     &     flx_fld%srunoff(istrt,lan),
     +     flx_fld%t1(istrt,lan),flx_fld%q1(istrt,lan),
     +     flx_fld%u1(istrt,lan),flx_fld%v1(istrt,lan),
     +     flx_fld%zlvl(istrt,lan),flx_fld%evbsa(istrt,lan),
     +     flx_fld%evcwa(istrt,lan),flx_fld%transa(istrt,lan),
     +     flx_fld%sbsnoa(istrt,lan),flx_fld%snowca(istrt,lan),
     +     flx_fld%soilm(istrt,lan),
     &     sfc_fld%tsea(istrt,lan),sfc_fld%sheleg(istrt,lan),
     &     sfc_fld%sncovr(istrt,lan),sfc_fld%tg3(istrt,lan),
     &     sfc_fld%zorl(istrt,lan),sfc_fld%cv(istrt,lan),
     &     sfc_fld%cvb(istrt,lan),sfc_fld%cvt(istrt,lan),
     &     sfc_fld%slmsk(istrt,lan),sfc_fld%vfrac(istrt,lan),
     &     sfc_fld%canopy(istrt,lan),
     &     sfc_fld%f10m(istrt,lan),sfc_fld%vtype(istrt,lan),
     &     sfc_fld%stype(istrt,lan),sfc_fld%uustar(istrt,lan),
     &     sfc_fld%ffmm(istrt,lan),sfc_fld%ffhh(istrt,lan),
     &     flx_fld%tmpmin(istrt,lan),flx_fld%tmpmax(istrt,lan),
     &     flx_fld%geshem(istrt,lan),
     &     flx_fld%dusfc(istrt,lan) ,flx_fld%dvsfc(istrt,lan) ,
     &     flx_fld%dtsfc(istrt,lan) ,
     &     flx_fld%dqsfc(istrt,lan),flx_fld%dlwsfc(istrt,lan),
     &     flx_fld%ulwsfc(istrt,lan),
     &     flx_fld%gflux(istrt,lan),flx_fld%runoff(istrt,lan),
     &     flx_fld%ep(istrt,lan),
     &     flx_fld%cldwrk(istrt,lan),flx_fld%dugwd(istrt,lan),
     &     flx_fld%dvgwd(istrt,lan),flx_fld%psmean(istrt,lan),
     &     flx_fld%bengsh(istrt,lan),
     &     xlon(istrt,lan),flx_fld%coszen(istrt,lan),
     &     flx_fld%sfcnsw(istrt,lan),
     +     xlat(istrt,lan),
     &     flx_fld%sfcdlw(istrt,lan),flx_fld%tsflw(istrt,lan) ,
     &     flx_fld%psurf(istrt,lan),flx_fld%u10m(istrt,lan),
     &     flx_fld%v10m(istrt,lan),
     &     sfc_fld%t2m(istrt,lan),sfc_fld%q2m(istrt,lan),
     &     flx_fld%hpbl(istrt,lan),flx_fld%pwat(istrt,lan),
     &     swh(1,1,iblk,lan),hlw(1,1,iblk,lan),
     &     smc_v,stc_v,hprime_v,
     &     slag,sdec,cdec,
     &     acv(istrt,lan),acvb(istrt,lan),acvt(istrt,lan),
     &     phy_f3dv,phy_f2dv,num_p3d,num_p2d,flgmin,

     &     dt3dt(1,1,1,iblk,lan), dq3dt(1,1,1,iblk,lan),
     &     du3dt(1,1,1,iblk,lan), dv3dt(1,1,1,iblk,lan),
     &     upd_mf(1,1,iblk,lan), dwn_mf(1,1,iblk,lan),
     &     det_mf(1,1,iblk,lan), ldiag3d,
     &     flipv,me,kdt,lat,sfc_fld%oro(istrt,lan),
! Coupling insertion->
!    > lssav_cc,
!    > DLWSFC_cc(istrt,lan),ULWSFC_cc(istrt,lan),SWSFC_cc(istrt,lan),
!    > DUSFC_cc(istrt,lan),DVSFC_cc(istrt,lan),DTSFC_cc(istrt,lan),
!    > DQSFC_cc(istrt,lan),PRECR_cc(istrt,lan),
!<-Coupling insertion
     &     crtrh, ncw,  old_monin, cnvgwd, ccwf, sashal, newsas )     ! hmhj

!
!!
      do k=1,lsoil
        do i=1,njeff
          sfc_fld%smc(k,istrt+i-1,lan) = smc_v(i,k)
          sfc_fld%stc(k,istrt+i-1,lan) = stc_v(i,k)
          sfc_fld%slc(k,istrt+i-1,lan) = slc_v(i,k)
        enddo
      enddo
!!
      do j=1,num_p3d
        do k=1,levs
          do i=1,njeff
            phy_f3d(i,k,iblk,lan,j) = phy_f3dv(i,k,j)
          enddo
        enddo
      enddo
      do j=1,num_p2d
        do i=1,njeff
          phy_f2d(istrt+i-1,lan,j) = phy_f2dv(i,j)
        enddo
      enddo

       do k = 1, LEVS
         do i = 1, njeff
           ilan=jlonr+istrt+i-1
           grid_gr(ilan,g_u+k-1) = adu(i,k)
           grid_gr(ilan,g_v+k-1) = adv(i,k)
           grid_gr(ilan,g_t+k-1) = adt(i,k)
         enddo
       enddo
       do n = 1, NTRAC
         kk = g_q + (n-1)*levs
         do k = 1, LEVS
           do i = 1, njeff
             ilan=jlonr+istrt+i-1
             grid_gr(ilan,kk+k-1) = adr(i,k,n)
           enddo
         enddo
       enddo
!!
       enddo   !lon
!
      enddo   !lan
cc
      call countperf(0,4,0.)
      call synctime()
      call countperf(1,4,0.)
!!
      return
      end
