      subroutine gloopb
!*   &    ( grid_gr,
     &    ( grid_fld, g3d_fld,                               
     x     lats_nodes_r,global_lats_r,lonsperlar,
     &     tstep,phour,sfc_fld, flx_fld, nst_fld, SFALB,xlon,
     &     swh,hlw,hprime,slag,sdec,cdec,
     &     ozplin,jindx1,jindx2,ddy,
     &     phy_f3d, phy_f2d,xlat,nblck,kdt,
     &     global_times_b,fscav)
!!
!! Code Revision:
!! Sep    2009       Shrinivas Moorthi added nst_fld
!! Oct 11 2009       Sarah Lu, grid_gr replaced by gri_fld
!! Oct 16 2009       Sarah Lu, grid_fld%tracers used
!! Nov 18 2009       Sarah Lu, rain/rainc added to gbphys call arg
!! Dec 14 2009       Sarah Lu, add g3d_fld to calling argument,
!!                             update dqdt after gbphys returns dqdt_v
!! July   2010       Shrinivas Moorthi - Updated for new physics
!! Aug    2010       Shrinivas Moorthi - Recoded 3d diagnostic arrays so that
!                              trap will not occur on call to gbphys
!! Oct 18 2010       Shrinivas Moorthi - Added fscav
!! Dec 23 2010       Sarah Lu, add lgocart to gbphys call arg
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
!-> Coupling insertion
!     USE SURFACE_cc
!<- Coupling insertion
      use d3d_def
      use gfs_physics_sfc_flx_mod
      use gfs_physics_nst_var_mod
      use gfs_physics_gridgr_mod, ONLY: Grid_Var_Data
      use gfs_physics_g3d_mod,    ONLY: G3D_Var_Data            
      use mersenne_twister
      include 'mpif.h'
      implicit none
!
!  **********************************************************************
!      The following arrays are for coupling to MOM4, but temporarily 
!      dimensioned here to make the code work.  Need to figure out how
!      to handel these  -- Moorthi
!
       real (kind=kind_phys) DLWSFC_cc(lonr,latr), ULWSFC_cc(lonr,latr)
     &,                      DTSFC_cc(lonr,latr),  SWSFC_cc(lonr,latr)
     &,                      DUSFC_cc(lonr,latr),  DVSFC_cc(lonr,latr)
     &,                      DQSFC_cc(lonr,latr),  PRECR_cc(lonr,latr)
 
     &,                      XMU_cc(lonr,latr),    DLW_cc(lonr,latr)
     &,                      DSW_cc(lonr,latr),    SNW_cc(lonr,latr)
     &,                      LPREC_cc(lonr,latr)
       logical lssav_cc
!  **********************************************************************
!
!
!*    real(kind=kind_grid) grid_gr(lonr*lats_node_r_max,lotgr)
      TYPE(Grid_Var_Data)       :: grid_fld 
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld
      TYPE(G3D_Var_Data)        :: g3d_fld 		

!
      integer id,njeff,lon,iblk,kdt,item
!!
      integer nblck
!!
      real(kind=kind_phys)    phour
      real(kind=kind_phys)    prsl(ngptc,levs)
      real(kind=kind_phys)   prslk(ngptc,levs), dpshc(ngptc)
      real(kind=kind_phys)    prsi(ngptc,levs+1),phii(ngptc,levs+1)
      real(kind=kind_phys)   prsik(ngptc,levs+1),phil(ngptc,levs)
!!
      real (kind=kind_rad) gu(ngptc,levs), gv(ngptc,levs)
      real (kind=kind_rad) gt(ngptc,levs), pgr(ngptc)
      real (kind=kind_rad) gr(ngptc,levs,ntrac)
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
     &     phy_f2d(lonr,lats_node_r,num_p2d), fscav(ntrac-ncld-1)
!!
      real (kind=kind_phys) dtphys,dtp,dtf
      real (kind=kind_evod) tstep
!!
      integer              lats_nodes_r(nodes)
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
!
      integer              i,j,k,kk,n
      integer              l,lan,lat,ii,lonrbm,jj
!     integer              l,lan,lat,jlonr,ilan,ii,lonrb2
      integer              lon_dim,lons_lat
      integer              nsphys
!
      real(kind=kind_evod) solhr,clstp
!
!timers______________________________________________________---
 
      real*8 rtc ,timer1,timer2
      real(kind=kind_evod) global_times_b(latr,nodes)
 
!timers______________________________________________________---
!
      logical, parameter :: flipv = .true.
      real(kind=kind_phys), parameter :: pt01=0.01, pt00001=1.0e-5
     &,                                  thousnd=1000.0
!
! for nrl/nasa ozone production and distruction rates:(input through fixio)
! ---------------------------------------------------
      integer jindx1(lats_node_r),jindx2(lats_node_r)    !for ozone interpolaton
      real(kind=kind_phys) ozplin(latsozp,levozp,pl_coeff,timeoz)
     &,                    ddy(lats_node_r)              !for ozone interpolaton
     &,                    ozplout(levozp,lats_node_r,pl_coeff)
!!
      real(kind=kind_phys), allocatable :: acv(:,:),acvb(:,:),acvt(:,:)
      save acv,acvb,acvt
!!
!     integer, parameter :: maxran=5000
!     integer, parameter :: maxran=3000
      integer, parameter :: maxran=6000, maxsub=6, maxrs=maxran/maxsub
      type (random_stat) :: stat(maxrs)
      real (kind=kind_phys), allocatable, save :: rannum_tank(:,:,:)
      real (kind=kind_phys)                    :: rannum(lonr*latr)
      integer iseed, nrc, seed0, kss, ksr, indxr(nrcm), iseedl
      integer nf0,nf1,ind,nt,indod,indev
      real(kind=kind_evod) fd2, wrk(1), wrk2(nrcm)

      logical first
      data first/.true./
!     save    krsize, first, nrnd,seed0
      save    first, seed0
!
      real(kind=kind_phys), parameter :: cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9
!
      real(kind=kind_phys) slag,sdec,cdec

!!
      integer nlons_v(ngptc)
      real(kind=kind_phys) smc_v(ngptc,lsoil),stc_v(ngptc,lsoil)
     &,                    slc_v(ngptc,lsoil)
     &,                    vvel(ngptc,levs)
     &,                    hprime_v(ngptc,nmtvr)
      real(kind=kind_phys) phy_f3dv(ngptc,LEVS,num_p3d),
     &                     phy_f2dv(ngptc,num_p2d)
     &,                    rannum_v(ngptc,nrcm)
      real(kind=kind_phys) sinlat_v(ngptc),coslat_v(ngptc)
     &,                    ozplout_v(ngptc,levozp,pl_coeff)
      real(kind=kind_phys) rqtk(ngptc)
      real(kind=kind_phys) dt3dt_v(ngptc,levs,6), du3dt_v(ngptc,levs,4)
     &,                    dv3dt_v(ngptc,levs,4)
     &,                    dq3dt_v(ngptc,levs,5+pl_coeff)
      real(kind=kind_phys) upd_mfv(ngptc,levs), dwn_mfv(ngptc,levs)
     &,                    det_mfv(ngptc,levs), dkh_v(ngptc,levs)
     &,                    rnp_v(ngptc,levs)

! local working array for moisture tendency 
      real(kind=kind_phys) dqdt_v(ngptc,LEVS) 

      real(kind=kind_phys) work1, qmin, tem
      parameter (qmin=1.0e-10)

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

        call random_setseed(seed0)
        call random_number(wrk)
        seed0 = seed0 + nint(wrk(1)*thousnd)
!
        if (.not. newsas) then  ! random number needed for RAS and old SAS
          if (random_clds) then ! create random number tank
!                                 -------------------------
            if (.not. allocated(rannum_tank))
     &                allocate (rannum_tank(lonr,maxran,lats_node_r))
!           lonrb2 = lonr / 2
            lonrbm = lonr / maxsub
            if (me == 0) write(0,*)' maxran=',maxran,' maxrs=',maxrs,
     &          'maxsub=',maxsub,' lonrbm=',lonrbm
!$OMP       parallel do private(nrc,iseedl,rannum,lat,i,j,k,ii,jj,kk)
            do nrc=1,maxrs
              iseedl = seed0 + nrc - 1
              call random_setseed(iseedl,stat(nrc))
              call random_number(rannum,stat(nrc))
              do j=1,lats_node_r
                lat  = global_lats_r(ipt_lats_node_r-1+j)
                jj = (lat-1)*lonr
                do k=1,maxsub
                  kk = k - 1
                  do i=1,lonr
                    ii = kk*lonrbm + i
                    if (ii > lonr) ii = ii - lonr
                    rannum_tank(i,nrc+kk*maxrs,j) = rannum(ii+jj)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
!
        if (me  ==  0) then
          write(0,*)' seed0=',seed0,' idate=',idate,' wrk=',wrk
          if (num_p3d == 3) write(0,*)' USING Ferrier-MICROPHYSICS'
          if (num_p3d == 4) write(0,*)' USING ZHAO-MICROPHYSICS'
        endif
        if (fhour == 0.0) then
          do j=1,lats_node_r
            do i=1,lonr
              phy_f2d(i,j,num_p2d) = 0.0
            enddo
          enddo
        endif
       
        if (ras) call ras_init(levs, me)
       
        first = .false.

      endif

!
      dtphys = 3600.
      nsphys = max(int(2*tstep/dtphys+0.9999),1)
      dtp    = (tstep+tstep)/nsphys
      dtf    = 0.5*dtp
      if(lsfwd) dtf = dtp
!
      solhr = mod(phour+idate(1),cons_24)

! **************  Ken Campana Stuff  ********************************
!...  set switch for saving convective clouds
      if(lscca.and.lsswr) then
        clstp = 1100+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
      elseif(lscca) then
        clstp = 0100+min(fhswr,fhour,cons_99)  !accumulate,convert
      elseif(lsswr) then
        clstp = 1100                           !initialize,accumulate
      else
        clstp = 0100                           !accumulate
      endif
! **************  Ken Campana Stuff  ********************************
!
!
      iseed = mod(100.0*sqrt(fhour*3600),cons_1p0d9) + 1 + seed0

      if (.not. newsas) then  ! random number needed for RAS and old SAS
        call random_setseed(iseed)
        call random_number(wrk2)
        if (random_clds) then
          do nrc=1,nrcm
            indxr(nrc) = max(1, min(nint(wrk2(nrc)*maxran)+1,maxran))
          enddo
        endif
      endif
!
! do ozone i/o and latitudinal interpolation to local gaussian lats
!
      if (ntoz > 0) then
       call ozinterpol(me,lats_node_r,lats_node_r,idate,fhour,
     &                 jindx1,jindx2,ozplin,ozplout,ddy)
      endif
!
! ----------------------------------------------------
!
      do lan=1,lats_node_r
         lat      = global_lats_r(ipt_lats_node_r-1+lan)
         lon_dim  = lon_dims_r(lan)
!        pwatp    = 0.
         lons_lat = lonsperlar(lat)
!        jlonr    = (lan-1)*lonr

!     write(0,*)' lan=',lan,' lats_node_r=',lats_node_r,' lons_lat='
!    &,lons_lat,' lat=',lat,' lonsperlar=',lonsperlar(lat)

!$omp parallel do  schedule(dynamic,1) private(lon)
!$omp+private(hprime_v,stc_v,smc_v,slc_v)
!$omp+private(nlons_v,sinlat_v,coslat_v,ozplout_v,rannum_v)
!$omp+private(prslk,prsl,prsik,prsi,phil,phii,dpshc)
!$omp+private(gu,gv,gt,gr,vvel)
!$omp+private(adt,adr,adu,adv,pgr,rqtk)
!$omp+private(phy_f3dv,phy_f2dv)
!$omp+private(dt3dt_v,du3dt_v,dv3dt_v,dq3dt_v,dqdt_v)
!$omp+private(upd_mfv,dwn_mfv,det_mfv,dkh_v,rnp_v)
!$omp+private(njeff,iblk,i,j,k,n,item)
!!$omp+private(njeff,iblk,ilan,i,j,k,n,item)
!!!$omp+private(temlon,temlat,lprnt,ipt)


        do lon=1,lons_lat,ngptc
!!
          njeff = min(ngptc,lons_lat-lon+1)
          iblk  = (lon-1)/ngptc + 1
!!
          do i = 1, njeff
!           ilan      = jlonr + lon + i - 1
!*          prsi(i,1) = grid_gr(ilan,g_ps)
            prsi(i,1) = grid_fld%ps(lon+i-1,lan)
            pgr(i)    = prsi(i,1)
 
!     write(0,*)' lan=',lan,' pgr=',pgr(i),' i=',i,' njeff=',njeff
!     print *,' lan=',lan,' pgr=',pgr(i),' grid_gr=',grid_gr(ilan,g_ps)
!    &,' i=',i,' lan=',lan
          enddo
          do k = 1, LEVS
            do i = 1, njeff
              item = lon+i-1
              gu(i,k)     = grid_fld%u(item,lan,k)        
              gv(i,k)     = grid_fld%v(item,lan,k)        
              gt(i,k)     = grid_fld%t(item,lan,k)      
              prsl(i,k)   = grid_fld%p(item,lan,k)      
              vvel(i,k)   = grid_fld%dpdt(item,lan,k)   
              prsi(i,k+1) = prsi(i,k) - grid_fld%dp(item,lan,k)   
            enddo
          enddo
          do i = 1, njeff
            prsi (i,levs+1) = 0.0
            prsik(i,levs+1) = 0.0
          enddo
          do n = 1, NTRAC
            do k = 1, LEVS
              do i = 1, njeff
                gr(i,k,n)= grid_fld%tracers(n)%flds(lon+i-1,lan,k)
              enddo
            enddo
          enddo

          do i=1,njeff
            phil(i,levs) = 0.0 ! will force calculation of geopotential in gbphys.
            dpshc(i)     = 0.3 * prsi(i,1)
!
            nlons_v(i)   = lons_lat
            sinlat_v(i)  = sinlat_r(lat)
            coslat_v(i)  = coslat_r(lat)
          enddo

          if (gen_coord_hybrid .and. thermodyn_id == 3) then
            do i=1,njeff
              prslk(i,1) = 0.0 ! forces calculation of geopotential in gbphys
              prsik(i,1) = 0.0 ! forces calculation of geopotential in gbphys
            enddo
          else
            do k = 1, levs
              do i = 1, njeff
                prslk(i,k) = (prsl(i,k)*pt00001)**rk
                prsik(i,k) = (prsi(i,k)*pt00001)**rk
              enddo
            enddo
          endif

          if (ntoz .gt. 0) then
            do j=1,pl_coeff
              do k=1,levozp
                do i=1,njeff
                  ozplout_v(i,k,j) = ozplout(k,lan,j)
                enddo
              enddo
            enddo
          endif

          do k=1,lsoil
            do i=1,njeff
              item = lon+i-1
              smc_v(i,k) = sfc_fld%smc(k,item,lan)
              stc_v(i,k) = sfc_fld%stc(k,item,lan)
              slc_v(i,k) = sfc_fld%slc(k,item,lan)
            enddo
          enddo
          do k=1,nmtvr
            do i=1,njeff
              hprime_v(i,k) = hprime(k,lon+i-1,lan)
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
              phy_f2dv(i,j) = phy_f2d(lon+i-1,lan,j)
            enddo
          enddo
          if (.not. newsas) then
            if (random_clds) then
              do j=1,nrcm
                do i=1,njeff
                  rannum_v(i,j) = rannum_tank(lon+i-1,indxr(j),lan)
                enddo
              enddo
            else
              do j=1,nrcm
                do i=1,njeff
                  rannum_v(i,j) = 0.6    ! This is useful for debugging
                enddo
              enddo
            endif
          endif
          if (ldiag3d) then
            do k=1,6
              do j=1,levs
                do i=1,njeff
                  dt3dt_v(i,j,k) = dt3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
            do k=1,4
              do j=1,levs
                do i=1,njeff
                  du3dt_v(i,j,k) = du3dt(i,j,k,iblk,lan)
                  dv3dt_v(i,j,k) = dv3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
          endif
          if (ldiag3d .or. lggfs3d) then
            do k=1,5+pl_coeff
              do j=1,levs
                do i=1,njeff
                  dq3dt_v(i,j,k) = dq3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
          endif
          if (lggfs3d) then
            do j=1,levs
              do i=1,njeff
                upd_mfv(i,j) = upd_mf(i,j,iblk,lan)
                dwn_mfv(i,j) = dwn_mf(i,j,iblk,lan)
                det_mfv(i,j) = det_mf(i,j,iblk,lan)
                dkh_v(i,j)   = dkh(i,j,iblk,lan)
                rnp_v(i,j)   = rnp(i,j,iblk,lan)
              enddo
            enddo
          endif
!
!     write(0,*)' before gbphys:', njeff,ngptc,levs,lsoil,lsm,          &
!    &      ntrac,ncld,ntoz,ntcw,                                       &
!    &      nmtvr,nrcm,levozp,lonr,latr,jcap,num_p3d,num_p2d,           &
!    &      kdt,lat,me,pl_coeff,ncw,flgmin,crtrh,cdmbgwd
!    &,' ccwf=',ccwf,' dlqf=',dlqf
!     write(0,*)' tisfc=',sfc_fld%tisfc(1:20,lan),' lan=',lan,' lon=',  &
!    &            lon
!     write(0,*) ' stc_v=',stc_v(1:5,1),' xlonlat=',xlon(lon,lan),
!    &xlat(lon,lan)
!     if (lan == 2) print *,' pgr=',pgr(1:5)
!     if (lan == 2) print *,' pgr=',pgr(45:55)
!
      lssav_cc = lssav      ! temporary assighment - neede to be revisited
!
!     if (lan == 1) call mpi_quit(4444)
          call gbphys                                                   &
!  ---  inputs:
     &    ( njeff,ngptc,levs,lsoil,lsm,ntrac,ncld,ntoz,ntcw,            &
     &      nmtvr,nrcm,levozp,lonr,latr,jcap,num_p3d,num_p2d,           &
     &      kdt,lat,me,pl_coeff,nlons_v,ncw,flgmin,crtrh,cdmbgwd,       &
     &      ccwf,dlqf,ctei_rm,clstp,dtp,dtf,fhour,solhr,                &
     &      slag,sdec,cdec,sinlat_v,coslat_v,pgr,gu,gv,                 &
     &      gt,gr,vvel,prsi,prsl,prslk,prsik,phii,phil,                 &
     &      rannum_v,ozplout_v,pl_pres,dpshc,                           &
     &      hprime_v, xlon(lon,lan),xlat(lon,lan),                      &
     &      sfc_fld%slope (lon,lan),    sfc_fld%shdmin(lon,lan),        &
     &      sfc_fld%shdmax(lon,lan),    sfc_fld%snoalb(lon,lan),        &
     &      sfc_fld%tg3   (lon,lan),    sfc_fld%slmsk (lon,lan),        &
     &      sfc_fld%vfrac (lon,lan),    sfc_fld%vtype (lon,lan),        &
     &      sfc_fld%stype (lon,lan),    sfc_fld%uustar(lon,lan),        &
     &      sfc_fld%oro   (lon,lan),    flx_fld%coszen(lon,lan),        &
     &      flx_fld%sfcdsw(lon,lan),    flx_fld%sfcnsw(lon,lan),        &
     &      flx_fld%sfcdlw(lon,lan),    flx_fld%tsflw (lon,lan),        &
     &      flx_fld%sfcemis(lon,lan),   sfalb(lon,lan),                 &
     &      swh(1,1,iblk,lan),hlw(1,1,iblk,lan),                        &
!    &      ras,pre_rad,ldiag3d,lggfs3d,lssav,                          &
!    &      ras,pre_rad,ldiag3d,lggfs3d,lssav,lssav_cc,                 &
     &      ras,pre_rad,ldiag3d,lggfs3d,lgocart,lssav,lssav_cc,         &
     &      bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,psautco,prautco,evpco,  &
     &      flipv,old_monin,cnvgwd,shal_cnv,sashal,newsas,cal_pre,      &
     &      mom4ice,mstrat,trans_trac,nst_fcst,moist_adj,fscav,         &
     &      thermodyn_id, sfcpress_id, gen_coord_hybrid,                &
!  ---  input/outputs:
     &      sfc_fld%hice  (lon,lan),    sfc_fld%fice  (lon,lan),        &
     &      sfc_fld%tisfc (lon,lan),    sfc_fld%tsea  (lon,lan),        &
     &      sfc_fld%tprcp (lon,lan),    sfc_fld%cv    (lon,lan),        &
     &      sfc_fld%cvb   (lon,lan),    sfc_fld%cvt   (lon,lan),        &
     &      sfc_fld%srflag(lon,lan),    sfc_fld%snwdph(lon,lan),        &
     &      sfc_fld%sheleg(lon,lan),    sfc_fld%sncovr(lon,lan),        &
     &      sfc_fld%zorl  (lon,lan),    sfc_fld%canopy(lon,lan),        &
     &      sfc_fld%ffmm  (lon,lan),    sfc_fld%ffhh  (lon,lan),        &
     &      sfc_fld%f10m  (lon,lan),    flx_fld%srunoff(lon,lan),       &
     &      flx_fld%evbsa (lon,lan),    flx_fld%evcwa (lon,lan),        &
     &      flx_fld%snohfa(lon,lan),    flx_fld%transa(lon,lan),        &
     &      flx_fld%sbsnoa(lon,lan),    flx_fld%snowca(lon,lan),        &
     &      flx_fld%soilm (lon,lan),    flx_fld%tmpmin(lon,lan),        &
     &      flx_fld%tmpmax(lon,lan),    flx_fld%dusfc (lon,lan),        &
     &      flx_fld%dvsfc (lon,lan),    flx_fld%dtsfc (lon,lan),        &
     &      flx_fld%dqsfc (lon,lan),    flx_fld%geshem(lon,lan),        &
     &      flx_fld%gflux (lon,lan),    flx_fld%dlwsfc(lon,lan),        &
     &      flx_fld%ulwsfc(lon,lan),    flx_fld%suntim(lon,lan),        &
     &      flx_fld%runoff(lon,lan),    flx_fld%ep    (lon,lan),        &
     &      flx_fld%cldwrk(lon,lan),    flx_fld%dugwd (lon,lan),        &
     &      flx_fld%dvgwd (lon,lan),    flx_fld%psmean(lon,lan),        &
     &      flx_fld%bengsh(lon,lan),    flx_fld%spfhmin(lon,lan),       &
     &      flx_fld%spfhmax(lon,lan),                                   &
     &      flx_fld%rain(lon,lan),      flx_fld%rainc(lon,lan),         &
     &      dt3dt_v, dq3dt_v,  du3dt_v, dv3dt_v, dqdt_v,                & ! added for GOCART
     &      acv(lon,lan), acvb(lon,lan), acvt(lon,lan),                 &
     &      slc_v, smc_v, stc_v, upd_mfv, dwn_mfv, det_mfv, dkh_v,rnp_v,&
     &      phy_f3dv, phy_f2dv,                                         &
     &      DLWSFC_cc(lon,lan),  ULWSFC_cc(lon,lan),                    &
     &      DTSFC_cc(lon,lan),   SWSFC_cc(lon,lan),                     &
     &      DUSFC_cc(lon,lan),   DVSFC_cc(lon,lan),                     &
     &      DQSFC_cc(lon,lan),   PRECR_cc(lon,lan),                     &

     &      nst_fld%xt(lon,lan),        nst_fld%xs(lon,lan),            &
     &      nst_fld%xu(lon,lan),        nst_fld%xv(lon,lan),            &
     &      nst_fld%xz(lon,lan),        nst_fld%zm(lon,lan),            &
     &      nst_fld%xtts(lon,lan),      nst_fld%xzts(lon,lan),          &
     &      nst_fld%d_conv(lon,lan),    nst_fld%ifd(lon,lan),           &
     &      nst_fld%dt_cool(lon,lan),   nst_fld%Qrain(lon,lan),         &
!  ---  outputs:
     &      adt, adr, adu, adv,                                         &
     &      sfc_fld%t2m   (lon,lan),    sfc_fld%q2m   (lon,lan),        &
     &      flx_fld%u10m  (lon,lan),    flx_fld%v10m  (lon,lan),        &
     &      flx_fld%zlvl  (lon,lan),    flx_fld%psurf (lon,lan),        &
     &      flx_fld%hpbl  (lon,lan),    flx_fld%pwat  (lon,lan),        &
     &      flx_fld%t1    (lon,lan),    flx_fld%q1    (lon,lan),        &
     &      flx_fld%u1    (lon,lan),    flx_fld%v1    (lon,lan),        &
     &      flx_fld%chh   (lon,lan),    flx_fld%cmm   (lon,lan),        &
     &      flx_fld%dlwsfci(lon,lan),   flx_fld%ulwsfci(lon,lan),       &
     &      flx_fld%dswsfci(lon,lan),   flx_fld%uswsfci(lon,lan),       &
     &      flx_fld%dtsfci(lon,lan),    flx_fld%dqsfci(lon,lan),        &
     &      flx_fld%gfluxi(lon,lan),    flx_fld%epi   (lon,lan),        &
     &      flx_fld%smcwlt2(lon,lan),   flx_fld%smcref2(lon,lan),       &
     &      flx_fld%wet1(lon,lan),                                      &
!hchuang code change [+3L] 11/12/2007 : add 2D
     &     flx_fld%gsoil(lon,lan),      flx_fld%gtmp2m(lon,lan),        &
     &     flx_fld%gustar(lon,lan),     flx_fld%gpblh(lon,lan),         &
     &     flx_fld%gu10m(lon,lan),      flx_fld%gv10m(lon,lan),         &
     &     flx_fld%gzorl(lon,lan),      flx_fld%goro(lon,lan),          &

     &      XMU_cc(lon,lan), DLW_cc(lon,lan), DSW_cc(lon,lan),          &
     &      SNW_cc(lon,lan), LPREC_cc(lon,lan),                         &

     &      nst_fld%Tref(lon,lan),       nst_fld%z_c(lon,lan),          &
     &      nst_fld%c_0 (lon,lan),       nst_fld%c_d(lon,lan),          &
     &      nst_fld%w_0 (lon,lan),       nst_fld%w_d(lon,lan),          &
     &      rqtk                                                        &! rqtkD
     &      )
!         if(kdt==100) then
!      print *,'in gloopb,aft gbphys,kdt=',kdt,'lat=',lat,lon,'smcwlt=',
!     &     flx_fld%smcwlt2(lon:lon+3,lan),
!     &    'loc=',minloc(flx_fld%smcwlt2(lon:lon+njeff-1,lan))
!         endif
!
!!
          do k=1,lsoil
            do i=1,njeff
              item = lon+i-1
              sfc_fld%smc(k,item,lan) = smc_v(i,k)
              sfc_fld%stc(k,item,lan) = stc_v(i,k)
              sfc_fld%slc(k,item,lan) = slc_v(i,k)
            enddo
          enddo
          if (ldiag3d) then
            do k=1,6
              do j=1,levs
                do i=1,njeff
                  dt3dt(i,j,k,iblk,lan) = dt3dt_v(i,j,k)
                enddo
              enddo
            enddo
            do k=1,4
              do j=1,levs
                do i=1,njeff
                  du3dt(i,j,k,iblk,lan) = du3dt_v(i,j,k)
                  dv3dt(i,j,k,iblk,lan) = dv3dt_v(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (ldiag3d .or. lggfs3d) then
            do k=1,5+pl_coeff
              do j=1,levs
                do i=1,njeff
                  dq3dt(i,j,k,iblk,lan) = dq3dt_v(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (lggfs3d) then
            do j=1,levs
              do i=1,njeff
                upd_mf(i,j,iblk,lan) = upd_mfv(i,j)
                dwn_mf(i,j,iblk,lan) = dwn_mfv(i,j)
                det_mf(i,j,iblk,lan) = det_mfv(i,j)
                dkh(i,j,iblk,lan)    = dkh_v(i,j)
                rnp(i,j,iblk,lan)    = rnp_v(i,j)
              enddo
            enddo
          endif
!!
!! total moist tendency (kg/kg/s): from local to global array
!!
      if (lgocart) then
        do k=1,levs
          do i=1,njeff
            g3d_fld%dqdt(lon+i-1,lan,k) = dqdt_v(i,k) 
          enddo        
        enddo         
      endif          
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
          phy_f2d(lon+i-1,lan,j) = phy_f2dv(i,j)
        enddo
      enddo

       do k = 1, LEVS
         do i = 1, njeff
           item = lon+i-1
           grid_fld%u(item,lan,k) = adu(i,k)            
           grid_fld%v(item,lan,k) = adv(i,k)         
           grid_fld%t(item,lan,k) = adt(i,k)
         enddo
       enddo
       do n = 1, NTRAC
         do k = 1, LEVS
           do i = 1, njeff
             grid_fld%tracers(n)%flds(lon+i-1,lan,k)= adr(i,k,n)
           enddo
         enddo
       enddo
!!
!     write(0,*)' adu=',adu(1,:)
!     write(0,*)' adv=',adv(1,:)
!     write(0,*)' adt=',adt(1,:)

       enddo                                   !lon
!
      enddo                                    !lan
!
      call countperf(0,4,0.)
      call synctime()
      call countperf(1,4,0.)
!!
!      write(0,*)' returning from gloopb for kdt=',kdt
!      if (kdt >1) call mpi_quit(3333)
      return
      end
