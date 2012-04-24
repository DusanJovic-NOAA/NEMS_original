      SUBROUTINE dfi_fixwr(iflag,

     & xt,xs,xu,xv,xz,zm,xtts,xzts,dt_cool,z_c,
     & c_0,c_d,w_0,w_d,d_conv,ifd,Tref,Qrain,

     & hice,fice,tisfc,                                ! FOR SEA-ICE - XW Nov04
     & tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,
     & alvsf,alvwf,alnsf,alnwf,vfrac,canopy,f10m,vtype,stype,
     & facsf,facwf,uustar,ffmm,ffhh,tprcp,srflag,
     & slc,snwdph,slope,shdmin,shdmax,snoalb,sncovr)

!
!***********************************************************************
!     PURPOSE:
!      save or retrieve fixed fields in digifilt
!
!     REVISION HISOTRY:
!     2011-12-05  J.Wang    Adopted from fixwr in /nwprod/sorc/global_fcst.fd
!
!***********************************************************************
!
      use resol_def
      use layout1
      implicit none
      integer,intent(in) ::  iflag
      real SMC(lonr,lsoil,lats_node_r),STC(lonr,lsoil,lats_node_r),
     &     HICE(lonr,lats_node_r),FICE(lonr,lats_node_r),  ! FOR SEA-ICE - NOV04
     &     TISFC(lonr,lats_node_r),

! li added for NST components
     &     xt     (lonr,lats_node_r), xs   (lonr,lats_node_r),
     &     xu     (lonr,lats_node_r), xv   (lonr,lats_node_r),
     &     xz     (lonr,lats_node_r), zm   (lonr,lats_node_r),
     &     xtts   (lonr,lats_node_r), xzts (lonr,lats_node_r),
     &     dt_cool(lonr,lats_node_r), z_c  (lonr,lats_node_r),
     &     c_0    (lonr,lats_node_r), c_d  (lonr,lats_node_r),
     &     w_0    (lonr,lats_node_r), w_d  (lonr,lats_node_r),
     &     d_conv (lonr,lats_node_r), ifd  (lonr,lats_node_r),
     &     Tref   (lonr,lats_node_r), Qrain(lonr,lats_node_r),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     &     TSEA  (lonr,lats_node_r),SHELEG(lonr,lats_node_r),
     &     TG3   (lonr,lats_node_r),
     &     ZORL  (lonr,lats_node_r),CV    (lonr,lats_node_r),
     &     CVB   (lonr,lats_node_r),
     &     CVT   (lonr,lats_node_r),ALVSF (lonr,lats_node_r),
     &     ALVWF (lonr,lats_node_r),
     &     ALNSF (lonr,lats_node_r),ALNWF (lonr,lats_node_r),
     &     SLMSK (lonr,lats_node_r),
     &     VFRAC (lonr,lats_node_r),CANOPY(lonr,lats_node_r),
     &     F10M  (lonr,lats_node_r),
     &     VTYPE (lonr,lats_node_r),STYPE (lonr,lats_node_r),
     &     FACSF (lonr,lats_node_r),
     &     FACWF (lonr,lats_node_r),UUSTAR(lonr,lats_node_r),
     &     FFMM  (lonr,lats_node_r),
     &     FFHH  (lonr,lats_node_r)
!lu [+5L]: add (tprcp,srflag),(slc,snwdph,snoalb,slope,shdmin,shdmax)
     +,    TPRCP (lonr,lats_node_r),SRFLAG(lonr,lats_node_r)
     +,    SLC    (lonr,lsoil,lats_node_r)
     +,    SNWDPH (lonr,lats_node_r)
     +,    SNOALB (lonr,lats_node_r),SLOPE (lonr,lats_node_r)
     +,    SHDMIN (lonr,lats_node_r),SHDMAX(lonr,lats_node_r)
     +,    SNCOVR (lonr,lats_node_r)

      real , allocatable :: SMC1(:,:,:),STC1(:,:,:),
     &  HICE1(:,:),FICE1(:,:),TISFC1(:,:),                   ! FOR SEA-ICE - XW Nov04

! li added for NST components
     &  xt1(:,:),xs1(:,:),xu1(:,:),xv1(:,:),xz1(:,:),zm1(:,:),
     &  xtts1(:,:),xzts1(:,:),dt_cool1(:,:),z_c1(:,:),
     &  c_01(:,:),c_d1(:,:),w_01(:,:),w_d1(:,:),
     &  d_conv1(:,:),ifd1(:,:),Tref1(:,:),Qrain1(:,:),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     &  TSEA1(:,:),SHELEG1(:,:),TG31(:,:),
     &  ZORL1(:,:),CV1(:,:),CVB1(:,:),
     &  CVT1(:,:),ALVSF1(:,:),ALVWF1(:,:),
     &  ALNSF1(:,:),ALNWF1(:,:),SLMSK1(:,:),
     &  VFRAC1(:,:),CANOPY1(:,:),F10M1(:,:),
     &  VTYPE1(:,:),STYPE1(:,:),FACSF1(:,:),
     &  FACWF1(:,:),UUSTAR1(:,:),FFMM1(:,:),
     &  FFHH1(:,:)
!lu [+3L]: add (tprcp1,srflag1),(slc1,snwdph1,slope1,shdmin1,shdmax1,snoalb1)
     +, TPRCP1(:,:),SRFLAG1(:,:)
     +, SLC1(:,:,:),SNWDPH1(:,:),SLOPE1(:,:)
     +, SHDMIN1(:,:),SHDMAX1(:,:),SNOALB1(:,:), SNCOVR1(:,:)

      logical first
      data first/.true./
      save   first,SMC1,STC1,TSEA1,SHELEG1,TG31,ZORL1,CV1,CVB1,CVT1
      save   HICE1,FICE1,TISFC1                        ! FOR SEA-ICE - XW Nov04
!                                                      ! FOR NST   - XL Dec0r97
      save   xt1,xs1,xu1,xv1,xz1,zm1,xtts1,xzts1,
     &       dt_cool1,z_c1,c_01,c_d1,w_01,w_d1,
     &       d_conv1,ifd1,Tref1,Qrain1
      save   ALVSF1,ALVWF1,ALNSF1,ALNWF1,SLMSK1,VFRAC1,CANOPY1,F10M1
      save   VTYPE1,STYPE1,FACSF1,FACWF1,UUSTAR1,FFMM1,FFHH1
!lu [+2L]: save (tprcp1,srflag1),(slc1,snwdph1,slope1,shdmin1,shdmax1,snoalb1)
      save   TPRCP1,SRFLAG1
      save   SLC1,SNWDPH1,SLOPE1,SHDMIN1,SHDMAX1,SNOALB1,SNCOVR1

      integer i,j,k

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      print *,' enter fixwr '                                   ! hmhj
      if (first) then
      print *,' enter fixwr in first allocate'
        allocate (SMC1(lonr,lsoil,lats_node_r))
        allocate (STC1(lonr,lsoil,lats_node_r))
        allocate (HICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (FICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (TISFC1(lonr,lats_node_r))             ! FOR SEA-ICE - XW Nov04
                                                        ! FOR NST     - XL Dec09
        allocate (xt1(lonr,lats_node_r))
        allocate (xs1(lonr,lats_node_r))
        allocate (xu1(lonr,lats_node_r))
        allocate (xv1(lonr,lats_node_r))
        allocate (xz1(lonr,lats_node_r))
        allocate (zm1(lonr,lats_node_r))
        allocate (xtts1(lonr,lats_node_r))
        allocate (xzts1(lonr,lats_node_r))

        allocate (dt_cool1(lonr,lats_node_r))
        allocate (z_c1(lonr,lats_node_r))
        allocate (c_01(lonr,lats_node_r))
        allocate (c_d1(lonr,lats_node_r))
        allocate (w_01(lonr,lats_node_r))
        allocate (w_d1(lonr,lats_node_r))
        allocate (d_conv1(lonr,lats_node_r))
        allocate (ifd1(lonr,lats_node_r))
        allocate (Tref1(lonr,lats_node_r))
        allocate (Qrain1(lonr,lats_node_r))

        allocate (TSEA1(lonr,lats_node_r))
        allocate (SHELEG1(lonr,lats_node_r))
        allocate (TG31(lonr,lats_node_r))
        allocate (ZORL1(lonr,lats_node_r))
        allocate (CV1(lonr,lats_node_r))
        allocate (CVB1(lonr,lats_node_r))
        allocate (CVT1(lonr,lats_node_r))
        allocate (ALVSF1(lonr,lats_node_r))
        allocate (ALVWF1(lonr,lats_node_r))
        allocate (ALNSF1(lonr,lats_node_r))
        allocate (ALNWF1(lonr,lats_node_r))
        allocate (SLMSK1(lonr,lats_node_r))
        allocate (VFRAC1(lonr,lats_node_r))
        allocate (CANOPY1(lonr,lats_node_r))
        allocate (F10M1(lonr,lats_node_r))
        allocate (VTYPE1(lonr,lats_node_r))
        allocate (STYPE1(lonr,lats_node_r))
        allocate (FACSF1(lonr,lats_node_r))
        allocate (FACWF1(lonr,lats_node_r))
        allocate (UUSTAR1(lonr,lats_node_r))
        allocate (FFMM1(lonr,lats_node_r))
        allocate (FFHH1(lonr,lats_node_r))
Clu [+8L]: allocate (tprcp,srflag),(slc,snwdph,slope,shdmin,shdmax,snoalb)
        allocate (TPRCP1(lonr,lats_node_r))
        allocate (SRFLAG1(lonr,lats_node_r))
        allocate (SLC1(lonr,lsoil,lats_node_r))
        allocate (SNWDPH1(lonr,lats_node_r))
        allocate (SLOPE1(lonr,lats_node_r))
        allocate (SHDMIN1(lonr,lats_node_r))
        allocate (SHDMAX1(lonr,lats_node_r))
        allocate (SNOALB1(lonr,lats_node_r))
        allocate (SNCOVR1(lonr,lats_node_r))
        first = .false.
      endif
!
      if(iflag == 1) then
       print *,'fixwr,iflag=1,asve 3hr phys'
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              smc1(i,k,j) = smc(i,k,j)
              stc1(i,k,j) = stc(i,k,j)
              slc1(i,k,j) = slc(i,k,j)        !! Clu [+1L]: slc -> slc1
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            hice1(i,j)    = hice(i,j)                   ! FOR SEA-ICE - XW Nov04
            fice1(i,j)    = fice(i,j)                   ! FOR SEA-ICE - XW Nov04
            tisfc1(i,j)   = tisfc(i,j)                  ! FOR SEA-ICE - XW Nov04
                                                        ! For NST
            xs1(i,j)      = xs(i,j)
            xu1(i,j)      = xu(i,j)
            xv1(i,j)      = xv(i,j)
            xz1(i,j)      = xz(i,j)
            zm1(i,j)      = zm(i,j)
            xtts1(i,j)    = xtts(i,j)
            xzts1(i,j)    = xzts(i,j)

            dt_cool1(i,j) = dt_cool(i,j)
            z_c1(i,j)     = z_c(i,j)
            c_01(i,j)     = c_0(i,j)
            c_d1(i,j)     = c_d(i,j)
            w_01(i,j)     = w_0(i,j)
            w_d1(i,j)     = w_d(i,j)
            d_conv1(i,j)  = d_conv(i,j)
            Tref1(i,j)    = Tref(i,j)
            Qrain1(i,j)   = Qrain(i,j)

            tsea1(i,j)    = tsea(i,j)
            sheleg1(i,j)  = sheleg(i,j)
            tg31(i,j)     = tg3(i,j)
            zorl1(i,j)    = zorl(i,j)
            cv1(i,j)      = cv(i,j)
            cvb1(i,j)     = cvb(i,j)
            cvt1(i,j)     = cvt(i,j)
            alvsf1(i,j)   = alvsf(i,j)
            alvwf1(i,j)   = alvwf(i,j)
            alnsf1(i,j)   = alnsf(i,j)
            alnwf1(i,j)   = alnwf(i,j)
            slmsk1(i,j)   = slmsk(i,j)
            vfrac1(i,j)   = vfrac(i,j)
            canopy1(i,j)  = canopy(i,j)
            f10m1(i,j)    = f10m(i,j)
            vtype1(i,j)   = vtype(i,j)
            stype1(i,j)   = stype(i,j)
            facsf1(i,j)   = facsf(i,j)
            facwf1(i,j)   = facwf(i,j)
            uustar1(i,j)  = uustar(i,j)
            ffmm1(i,j)    = ffmm(i,j)
            ffhh1(i,j)    = ffhh(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            tprcp1(i,j)   = tprcp(i,j)
            srflag1(i,j)  = srflag(i,j)
            snwdph1(i,j)  = snwdph(i,j)
            slope1(i,j)   = slope(i,j)
            shdmin1(i,j)  = shdmin(i,j)
            shdmax1(i,j)  = shdmax(i,j)
            snoalb1(i,j)  = snoalb(i,j)
            sncovr1(i,j)  = sncovr(i,j)
          enddo
        enddo
        print *,'save,tsea(65,58)=',tsea(65,58)
      elseif(iflag == 2) then
       print *,'fixwr,iflag=2,set 3hr phys'
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              smc(i,k,j) = smc1(i,k,j)
              stc(i,k,j) = stc1(i,k,j)
              slc(i,k,j) = slc1(i,k,j)          !! Clu [+1L]: slc1 -> slc
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            hice(i,j)    = hice1(i,j)                 ! FOR SEA-ICE - XW Nov04
            fice(i,j)    = fice1(i,j)                 ! FOR SEA-ICE - XW Nov04
            tisfc(i,j)   = tisfc1(i,j)                ! FOR SEA-ICE - XW Nov04

            xs(i,j)      = xs1(i,j)
            xu(i,j)      = xu1(i,j)
            xv(i,j)      = xv1(i,j)
            xz(i,j)      = xz1(i,j)
            zm(i,j)      = zm1(i,j)
            xtts(i,j)    = xtts1(i,j)
            xzts(i,j)    = xzts1(i,j)

            dt_cool(i,j) = dt_cool1(i,j)
            z_c(i,j)     = z_c1(i,j)
            c_0(i,j)     = c_01(i,j)
            c_d(i,j)     = c_d1(i,j)
            w_0(i,j)     = w_01(i,j)
            w_d(i,j)     = w_d1(i,j)
            d_conv(i,j)  = d_conv1(i,j)
            Tref(i,j)    = Tref1(i,j)
            Qrain(i,j)   = Qrain1(i,j)

            tsea(i,j)    = tsea1(i,j)
            sheleg(i,j)  = sheleg1(i,j)
            tg3(i,j)     = tg31(i,j)
            zorl(i,j)    = zorl1(i,j)
            cv(i,j)      = cv1(i,j)
            cvb(i,j)     = cvb1(i,j)
            cvt(i,j)     = cvt1(i,j)
            alvsf(i,j)   = alvsf1(i,j)
            alvwf(i,j)   = alvwf1(i,j)
            alnsf(i,j)   = alnsf1(i,j)
            alnwf(i,j)   = alnwf1(i,j)
            slmsk(i,j)   = slmsk1(i,j)
            vfrac(i,j)   = vfrac1(i,j)
            canopy(i,j)  = canopy1(i,j)
            f10m(i,j)    = f10m1(i,j)
            vtype(i,j)   = vtype1(i,j)
            stype(i,j)   = stype1(i,j)
            facsf(i,j)   = facsf1(i,j)
            facwf(i,j)   = facwf1(i,j)
            uustar(i,j)  = uustar1(i,j)
            ffmm(i,j)    = ffmm1(i,j)
            ffhh(i,j)    = ffhh1(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            tprcp(i,j)   = tprcp1(i,j)
            srflag(i,j)  = srflag1(i,j)
            snwdph(i,j)  = snwdph1(i,j)
            slope(i,j)   = slope1(i,j)
            shdmin(i,j)  = shdmin1(i,j)
            shdmax(i,j)  = shdmax1(i,j)
            snoalb(i,j)  = snoalb1(i,j)
            sncovr(i,j)  = sncovr1(i,j)
          enddo
        enddo
        print *,'after reset tsea(65,58)=',tsea(65,58)
      endif
      print *,' leave fixwr '
      return
      end subroutine dfi_fixwr
