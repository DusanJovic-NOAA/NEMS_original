      subroutine shalcnv(im,ix,km,jcap,delt,del,prsl,ps,phil,ql,
     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,
     &     dot,ncloud,hpbl)
!    &     dot,ncloud,hpbl,me)
!
      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap
     &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c
     &,             cvap => con_cvap, cliq => con_cliq
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer            im, ix,  km, jcap, ncloud,
     &                   kbot(im), ktop(im), kcnv(im)
!    &,                  me
      real(kind=kind_phys) delt
      real(kind=kind_phys) ps(im),     del(ix,km),  prsl(ix,km),
     &                     ql(ix,km,2),q1(ix,km),   t1(ix,km),
     &                     u1(ix,km),  v1(ix,km),   rcs(im),
     &                     cldwrk(im), rn(im),      slimsk(im), 
     &                     dot(ix,km), phil(ix,km), hpbl(im)
!
      integer              i, indx, jmn, k, latd, lond, km1
      integer              kpbl(im)
cj
      real(kind=kind_phys) alpha,   alphal,  alphas,
     &                     beta,    betal,   betas,
     &                     c0,      cpoel,   dellat,  delta,
     &                     desdt,   dg,
     &                     dh,      dhh,     dlnsig,  dp,
     &                     dq,      dqsdp,   dqsdt,   dt,
     &                     dt2,     dtmax,   dtmin,   dv1h,
     &                     dv1q,    dv2h,    dv2q,    dv1u,
     &                     dv1v,    dv2u,    dv2v,    dv3q,
     &                     dv3h,    dv3u,    dv3v,    alphat,
     &                     dz,      dz1,     e1, 
     &                     el2orc,  elocp,   es,      etah,
     &                     evef,    evfact,  evfactl, fact1,
     &                     fact2,   factor,  fjcap,   fkm,
     &                     g,       gamma,   pprime,
     &                     qc,      qlk,     qrch,    qs,
     &                     rfact,   shear,   tem1,    dthk,
     &                     tem2,    terr,    val,     val1,
     &                     val2,    w1,      w1l,     w1s,
     &                     w2,      w2l,     w2s,     w3,
     &                     w3l,     w3s,     w4,      w4l,
     &                     w4s,     xdby,    xpw,     xpwd,
     &                     xqc,     xqrch,   mbdt,    tem,
     &                     ptem,    ptem1,   pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im), ktcon(im),
     &                     ktcon1(im), kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     acrt(im),   acrtfct(im),
     &                     delhbar(im), delq(im),   delq2(im),
     &                     delqbar(im), delqev(im), deltbar(im),
     &                     deltv(im),   dtconv(im), edt(im),
     &                     fld(im),     cincr(im),
     &                     pbcdif(im),  pdot(im),   po(im,km),
     &                     qcond(im),   qevap(im),  hmax(im),
     &                     rntot(im),   vshear(im), xaa0(im),
     &                     xk(im),      xlamb(im),  xlamue(im),  
     &                     xlamdt(im),  xmb(im),    xmbmax(im),
     &                     delubar(im), delvbar(im)
c
      real(kind=kind_phys) cincrmax, cincrmin
cc
c  physical parameters
      parameter(g=grav)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,
     &          el2orc=hvap*hvap/(rv*cp))
      parameter(terr=0.,c0=.002,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cincrmax=180.,cincrmin=120.)
      parameter(dthk=25.)
c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km),
     &                     uo(im,km),      vo(im,km),     qeso(im,km)
c  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im),    tvo(im,km),
      real(kind=kind_phys) qlko_ktcon(im), dellal(im),
     &                     dbyo(im,km),    zo(im,km),
     &                     heo(im,km),     heso(im,km),
     &                     dellah(im,km),  dellaq(im,km),
     &                     dellau(im,km),  dellav(im,km), hcko(im,km),
     &                     ucko(im,km),    vcko(im,km),   qcko(im,km),
     &                     eta(im,km),     zi(im,km),     pwo(im,km),
     &                     tx1(im),        sumx(im),      dzmax(im)
!
      logical totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) pcrit(15), acritt(15), acrit(15)
!!    save pcrit, acritt
      data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,
     &           350.,300.,250.,200.,150./
      data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
c  gdas derived acrit
c     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
c    &            .743,.813,.886,.947,1.138,1.377,1.896/
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
cc--------------------------------------------------------------------
      real(kind=kind_phys) cons_0               !constant
      cons_0          =         0.d0            !constant
cc--------------------------------------------------------------------
!
      km1 = km - 1
c
c  initialize arrays
c
      do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i).eq.1) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        rn(i)=0.
        cldwrk(i) = 0.
        kbcon(i)=km
        ktcon(i)=1
        kb(i)=km
        dtconv(i) = 3600.
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
        dellal(i) = 0.
        edt(i)  = 0.
        acrt(i) = 0.
        acrtfct(i) = 1.
        aa1(i)  = 0.
        xaa0(i) = 0.
        vshear(i) = 0.
      enddo
c
      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo
      dt2 = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
c  model tunable parameters are all here
      mbdt    = 10.
      alphal  = .5
      alphas  = .5
      alphat  = .1
!     betal   = .15
!     betas   = .15
      betal   = .05
      betas   = .05
c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3 
      w2l     = -4.e-2
      w3l     = -5.e-3 
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
c
c  define top layer for search of the downdraft originating layer
c  and the maximum thetae for updraft
c
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!     
      do k = 1, km1
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.50) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.40) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
c
c  hydrostatic height assume zero terr
c
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
        enddo
      enddo
c
c  pbl height
c
      do i=1,im
        flg(i) = .true.
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.zi(i,k).ge.hpbl(i)) then
             kpbl(i) = k
             flg(i) = .false.
          endif
        enddo
      enddo
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c   convert surface pressure to mb from cb
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
             pfld(i,k) = prsl(i,k) * 10.0
             eta(i,k)  = 1.
             hcko(i,k) = 0.
             qcko(i,k) = 0.
             ucko(i,k) = 0.
             vcko(i,k) = 0.
             dbyo(i,k) = 0.
             pwo(i,k)  = 0.
             to(i,k)   = t1(i,k)
             qo(i,k)   = q1(i,k)
             uo(i,k)   = u1(i,k)
             vo(i,k)   = v1(i,k)
          endif
        enddo
      enddo
c
c  column variables
c  p is pressure of the layer (mb)
c  t is temperature at t-dt (k)..tn
c  q is mixing ratio at t-dt (kg/kg)..qn
c  to is temperature at t+dt (k)... this is after advection and turbulan
c  qo is mixing ratio at t+dt (kg/kg)..q1
c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
c
c  compute moist static energy
c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
c
c  determine level with largest moist static energy within pbl
c  this is the level where updraft starts
c
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i).and.k .lt. kpbl(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
c
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
c
c  look for the level of free convection as cloud base
c
      do i=1,im
        flg(i)    = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.k.le.kbm(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
c
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  determine critical convective inhibition
c  as a function of vertical velocity at cloud base.
c
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            ptem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            ptem = - (pdot(i) + w4) / (w4 - w3)
          else
            ptem = 0.
          endif
          val1    =             -1.
          ptem = max(ptem,val1)
          val2    =             1.
          ptem = min(ptem,val2)
          ptem = 1. - ptem
          ptem1= .5*(cincrmax-cincrmin)
          cincr(i) = cincrmax - ptem * ptem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(pbcdif(i).gt.cincr(i)) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  determine convective cloud top as the level of zero buoyancy
c  without considering entrainment and detrainment
c
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = 1
      enddo
      do k = 2, km
      do i=1,im
        if (flg(i).and.k .le. kbm(i)) then
          if(k.gt.kbcon(i).and.heo(i,kb(i)).lt.heso(i,k)) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
c
      do i = 1, im
        if(cnvflg(i)) then
          pbcdif(i) = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          if(pbcdif(i).le.0..or.pbcdif(i).gt.300.) cnvflg(i)=.false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c if cloud top is less than the pbl top, which will be the case
c of stratocumulus topped pbl, the shallow convection is turned
c off to prevent stratocumulus from being unrealistically destroyed.
c
      do i = 1, im
        if(cnvflg(i)) then
          if(pfld(i,ktcon(i)).ge.pfld(i,kpbl(i))) then
             cnvflg(i)=.false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  determine entrainment rate between kb and kbcon
c
      do i = 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
        endif
      enddo
      do k = 1, km1
      do i = 1, im
        if(cnvflg(i).and.k.ge.kb(i).and.k.lt.kbcon(i)) then
          dz = zi(i,k+1) - zi(i,k)
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        alpha = alphas
        if(slimsk(i).eq.1.) alpha = alphal
        if(cnvflg(i)) then
          dz  = sumx(i)/float(kbcon(i)-kb(i))
          tem = 1./float(kbcon(i)-kb(i))
          ptem= 1./alpha
          xlamb(i) = (ptem**tem - 1.) / dz
        endif
      enddo
c
c  work up updraft cloud properties
c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
c
c  cloud property in subcloud layer is modified by the entrainment process
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.kbcon(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = xlamb(i) * dz
              factor = 1. + tem
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = (hcko(i,k-1)+tem*0.5*
     &                     (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = (ucko(i,k-1)+ptem*uo(i,k)
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = (vcko(i,k-1)+ptem*vo(i,k)
     &                     +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k)-heso(i,k)
            endif
          endif
        enddo
      enddo
c
c  determine updraft mass flux for subcloud layer
c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              eta(i,k) = eta(i,k+1) / (1. + xlamb(i) * dz)
            endif
          endif
        enddo
      enddo
c
c  compute the average of the entrainment rate for the cloud layer
c  which is inversely proportional to the height
c  (~0.55/z: Siebesma et al. 2003)
c
      do i = 1, im
        if(cnvflg(i)) then
          dz = zi(i,ktcon(i))-zi(i,kbcon(i))
          xlamue(i)=0.55*log(zi(i,ktcon(i))/zi(i,kbcon(i)))/dz
        endif
      enddo
c
c  enhance the detrainment rate for the updrafts for mass flux to
c  decrease above cloud base.
c  the mass flux at cloud top is assumed to is reduced to
c  (alphat*100)% of one at cloud base.
c
      do i = 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
          dzmax(i) = 0.
        endif
      enddo
      do k = 1, km1
      do i = 1, im
        if(cnvflg(i)) then
          if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
            dz = zi(i,k+1) - zi(i,k)
            if(dz.gt.dzmax(i)) dzmax(i) = dz
            sumx(i) = sumx(i) + dz
          endif
        endif
      enddo
      enddo
c
      do i = 1, im
        if(cnvflg(i)) then
          dz  = sumx(i)/float(ktcon(i)-kbcon(i))
          tem = 1./float(ktcon(i)-kbcon(i))
          ptem= (alphat**tem - 1.) / dz
          xlamdt(i) = xlamue(i) - ptem
          xlamdt(i) = max(xlamdt(i),xlamue(i))
          xlamdt(i) = min(xlamdt(i),xlamue(i)+1./dzmax(i))
        endif
      enddo
c
c  compute updraft property above cloud base
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kbcon(i).and.k.le.ktcon(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = xlamue(i) * dz
              tem1 = 0.5 * xlamdt(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*
     &                     (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k)
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k)
     &                     +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
c
c  renew the cloud top
c
      do i=1,im
        flg(i) = cnvflg(i)
        if(flg(i)) kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.ktcon(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          ptem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(ptem.gt.dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
c
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon1(i) = ktcon(i)
      enddo
      do k = 2, km
      do i=1,im
        if (flg(i).and.k .le. ktcon(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon1(i) = k
             flg(i)    = .false.
          endif
        endif
      enddo
      enddo
c
      do i = 1, im
        if(cnvflg(i)) then
          ktcon(i) = ktcon1(i)
          if(pfld(i,ktcon(i)).ge.pfld(i,kpbl(i))) then
             cnvflg(i)=.false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  specify upper limit of mass flux at cloud base
c
      do i = 1, im
        if(cnvflg(i)) then
          xmbmax(i) = .1
        endif
      enddo
c
c  determine updraft mass flux for cloud layer
c
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k.gt.kbcon(i).and.k.le.ktcon(i)) then
              dz        = zi(i,k) - zi(i,k-1)
              ptem      = xlamue(i) - xlamdt(i)
              eta(i,k)  = eta(i,k-1) * (1 + ptem * dz)
            endif
          endif
        enddo
      enddo
c
c  compute cloud moisture property and precipitation
c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
!         rhbar(i) = 0.
        endif
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              dz1   = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
cj
              if(k.le.kbcon(i)) then
                 tem  = xlamb(i) * dz
                 tem1 = 0.
              else 
                 tem  = xlamue(i) * dz
                 tem1 = 0.5 * xlamdt(i) * dz
              endif
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*
     &                     (qo(i,k)+qo(i,k-1)))/factor
cj
              dq = eta(i,k) * qcko(i,k) - eta(i,k) * qrch
!
!             dq = qcko(i,k) - qrch
c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
c
c  below lfc check if there is excess moisture to release latent heat
c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                qlk = dq / (eta(i,k) + etah * c0 * dz)
!
!               qlk = dq / (1. + c0 * dz)
!               qlk = c0 * dz * dq
!               qlk = min(qlk,dq)
!
                aa1(i) = aa1(i) - dz1 * g * qlk
                qc = qlk + qrch
!
!               qc = qrch + dq - qlk
!               pwo(i,k) = eta(i,k-1) * c0 * dz * qlk
!               pwo(i,k) = eta(i,k-1) * qlk
!
                pwo(i,k) = etah * c0 * dz * qlk
                qcko(i,k)= qc
              endif
            endif
          endif
        enddo
      enddo
c
!     do i = 1, im
!       if(cnvflg(i)) then
!         indx = ktcon(i) - kb(i) - 1
!         rhbar(i) = rhbar(i) / float(indx)
!       endif
!     enddo
c
c  this section is ready for cloud water
c
      if(ncloud.gt.0) then
c
c  compute liquid and vapor separation at cloud top
c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i)
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k-1) - qrch
c
c  check if there is excess moisture to release latent heat
c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k-1) = qrch
          endif
        endif
      enddo
      endif
c
c  calculate cloud work function at t+dt
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +
     &                 dz1 * (g / (cp * to(i,k)))
     &                 * dbyo(i,k) / (1. + gamma)
     &                 * rfact
              val = 0.
              aa1(i)=aa1(i)+
     &                 dz1 * g * delta *
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) then
           cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c--- compute precipitation efficiency in terms of windshear
c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear=rcs(i) * sqrt((uo(i,k)-uo(i,k-1)) ** 2
     &                          + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
        endif
      enddo
c
c--- what would the change be, that a cloud with unit mass
c--- will do to the environment?
c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
c
c--- changed due to subsidence and entrainment
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = .5 * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = .5 * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)
c
              if(k.le.kbcon(i)) then
                tem  = xlamb(i)
                tem1 = 0.
              else
                tem  = xlamue(i)
                tem1 = xlamdt(i)
              endif
cj
              dellah(i,k) = dellah(i,k) +
     &     ( eta(i,k)*dv1h - eta(i,k-1)*dv3h
     &    -  tem*eta(i,k-1)*dv2h*dz
     &    +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz
     &         ) *g/dp
cj
              dellaq(i,k) = dellaq(i,k) +
     &     ( eta(i,k)*dv1q - eta(i,k-1)*dv3q
     &    -  tem*eta(i,k-1)*dv2q*dz
     &    +  tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz
     &         ) *g/dp
cj
              dellau(i,k) = dellau(i,k) +
     &     ( eta(i,k)*dv1u - eta(i,k-1)*dv3u
     &    -  tem*eta(i,k-1)*dv2u*dz
     &    +  tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz
     &    -  pgcon*eta(i,k-1)*(dv1u-dv3u)
     &         ) *g/dp
cj
              dellav(i,k) = dellav(i,k) +
     &     ( eta(i,k)*dv1v - eta(i,k-1)*dv3v
     &    -  tem*eta(i,k-1)*dv2v*dz
     &    +  tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz
     &    -  pgcon*eta(i,k-1)*(dv1v-dv3v)
     &         ) *g/dp
cj
            endif
          endif
        enddo
      enddo
c
c------- cloud top
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *
     &                     (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *
     &                     (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *
     &                     (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *
     &                     (vcko(i,indx-1) - dv1v) * g / dp
c
c  cloud water
c
          dellal(i) = eta(i,indx-1) * qlko_ktcon(i) * g / dp
        endif
      enddo
c
c------- final changed variable per unit mass flux
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i).and.k .le. kmax(i)) then
            if(k.gt.ktcon(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon(i)) then
              qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              to(i,k) = dellat * mbdt + t1(i,k)
              val   =           1.e-10
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c--- the above changed environment is now used to calulate the
c--- effect the arbitrary cloud (with unit mass flux)
c--- would have on the stability,
c--- which then is used to calculate the real mass flux,
c--- necessary to keep this change in balance with the large-scale
c--- destabilization.
c
c--- environmental conditions again, first heights
c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k)+epsm1*qeso(i,k))
            val       =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
c
c--- moist static energy
c
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (pfld(i,k+1) - pfld(i,k))
            es = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime = pfld(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1 * qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)   = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                    cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
c         heo(i,k) = min(heo(i,k),heso(i,k))
        endif
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
        endif
      enddo
c
c**************************** static control
c
c------- moisture and cloud work functions
c
      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
        endif
      enddo
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              if(k.le.kbcon(i)) then
                 tem  = xlamb(i) * dz
                 tem1 = 0.
              else
                 tem  = xlamue(i) * dz
                 tem1 = 0.5 * xlamdt(i) * dz
              endif
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*
     &                     (heo(i,k)+heo(i,k-1)))/factor
            endif
          endif
        enddo
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              xdby = hcko(i,k) - heso(i,k)
              val  =          0.
              xdby = max(xdby,val)
              xqrch = qeso(i,k)
     &              + gamma * xdby / (hvap * (1. + gamma))
cj
              if(k.le.kbcon(i)) then
                 tem  = xlamb(i) * dz
                 tem1 = 0.
              else
                 tem  = xlamue(i) * dz
                 tem1 = 0.5 * xlamdt(i) * dz
              endif
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*
     &                     (qo(i,k)+qo(i,k-1)))/factor
cj
              dq = eta(i,k) * qcko(i,k) - eta(i,k) * xqrch
!
!             dq = qcko(i,k) - xqrch
c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                qlk = dq / (eta(i,k) + etah * c0 * dz)
!
!               qlk = dq / (1. + c0 * dz)
!               qlk = c0 * dz * dq
!               qlk = min(qlk,dq)
!
                xaa0(i) = xaa0(i) - (zo(i,k+1) - zo(i,k)) * g * qlk
                xqc = qlk + xqrch
!
!               xqc = xqrch + dq - qlk
!               xpw = eta(i,k-1) * c0 * dz * qlk
!               xpw = eta(i,k-1) * qlk
!
                xpw = etah * c0 * dz * qlk
                qcko(i,k) = xqc
              endif
            endif
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma
     &                 * to(i,k) / hvap
              xdby = hcko(i,k) - heso(i,k)
              xaa0(i) = xaa0(i)
     &                + dz1 * (g / (cp * to(i,k)))
     &                * xdby / (1. + gamma)
     &                * rfact
              val=0.
              xaa0(i)=xaa0(i)+
     &                 dz1 * g * delta *
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
c
c  calculate critical cloud work function
c
      do i = 1, im
        if(cnvflg(i)) then
          if(pfld(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.-pfld(i,ktcon(i)))
     &              /(975.-pcrit(15))
          else if(pfld(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850. - pfld(i,ktcon(i)))/50.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*
     &           (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
        endif
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
c
c  modify critical cloud workfunction by cloud base vertical velocity
c
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          val1    =             -1.
          acrtfct(i) = max(acrtfct(i),val1)
          val2    =             1.
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = 1. - acrtfct(i)
c
c  modify acrtfct(i) by colume mean rh if rhbar(i) is greater than 80 percent
c
c         if(rhbar(i).ge..8) then
c           acrtfct(i) = acrtfct(i) * (.9 - min(rhbar(i),.9)) * 10.
c         endif
c
c  modify adjustment time scale by cloud base vertical velocity
c
          dtconv(i) = dt2 + max((1800. - dt2),cons_0) *     !constant
     &                (pdot(i) - w2) / (w1 - w2)
c         dtconv(i) = max(dtconv(i), dt2)
c         dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
c
        endif
      enddo
c
c--- large scale forcing
c
      do i= 1, im
        if(cnvflg(i)) then
!!        f = aa1(i) / dtconv(i)
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
!!        fld(i)=aa1(i)/dtconv(i)
!!        fld(i)=(aa1(i)-acrt(i))/dtconv(i)
          if(fld(i).le.0.) cnvflg(i) = .false.
        endif
        if(cnvflg(i)) then
c         xaa0(i) = max(xaa0(i),0.)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) cnvflg(i) = .false.
        endif
c
c--- kernel, cloud base mass flux
c
        if(cnvflg(i)) then
          xmb(i) = -fld(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!           to(i,k) = t1(i,k)
!           qo(i,k) = q1(i,k)
!           uo(i,k) = u1(i,k)
!           vo(i,k) = v1(i,k)
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c--- feedback: simply the changes from the cloud with unit mass flux
c---           multiplied by  the mass flux necessary to keep the
c---           equilibrium with the larger-scale.
c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
c
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.ktcon(i).and.k.gt.kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
c
c evaporating rain
c
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k.lt.ktcon(i).and.k.gt.kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i).and.k.le.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
c             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.
     &           delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                rn(i) = rn(i) - .001 * qevap(i) * dp / g
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
cj
!     do i = 1, im
!     if(cnvflg(i)) then
!     if(me.eq.31.and.cnvflg(i)) then
!     if(me.eq.8.and.cnvflg(i)) then
c       print *, ' kb,kbot,ktop =',
c    &            kb(i),kbcon(i),ktcon(i)
c       print *, 'pkb pbot, ptop =', po(i,kb(i)),
c    &po(i,kbcon(i)),po(i,ktcon(i))
c       print *, '   eta ='
c       print *, (eta(i,k),k=1,kmax(i))
!       print *, ' shallow delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip xmb= ', hvap*rn(i)*1000./dt2, xmb(i)
!     endif
!     enddo
cj
c
c  cloud water
c
      if (ncloud.gt.0) then
      do i = 1, im
        if (cnvflg(i)) then
           k = ktcon(i)
           tem  = dellal(i) * xmb(i) * dt2
           tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
           if (ql(i,k,2) .gt. -999.0) then
             ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
             ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
           else
             ql(i,k,1) = ql(i,k,1) + tem
           endif
           dp = 1000. * del(i,k)
           dellal(i) = dellal(i) * xmb(i) * dp / g
        endif
      enddo
      endif
c
      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i).lt.0..and..not.flg(i)) rn(i) = 0.
          if(rn(i).le.0.) then
            rn(i) = 0.
          endif
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 0
          cldwrk(i) = aa1(i)
        endif
      enddo
!!
      return
      end
