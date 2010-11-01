      SUBROUTINE read_mtn_hprim_oz(SLMSK,HPRIME,NEEDORO,ORO,
     &           iozondp,ozplin,global_lats_r,lonsperlar)
!
!***********************************************************************
!
      use resol_def, ONLY: latr, lonr, nmtvr
      use layout1,   ONLY: me, nodes, lats_node_r
      use ozne_def,  ONLY: latsozp, levozp, timeoz, pl_coeff
      USE machine,   ONLY: kind_io8, kind_io4
      implicit none

!
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
      real (kind=kind_io8) SLMSK(lonr,lats_node_r),
     &  HPRIME(NMTVR,lonr,lats_node_r),ORO(lonr,lats_node_r)
 
      integer iozondp
      real (kind=kind_io8) ozplin(latsozp,levozp,pl_coeff,timeoz)
 
      real(kind=kind_io4) buff1(lonr,latr),buffm(lonr,latr,nmtvr)
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      real(kind=kind_io8) buff2(lonr,lats_node_r)
      integer kmsk0(lonr,latr)
      integer i,j,k,nmtn
      integer needoro
!
      kmsk0=0
!
!     Read HPRIME from file MTNVAR
!     ****************************
      nmtn=24
!jfe  IF (me.eq.0) THEN
      IF (me.eq.0) THEN   
        READ(nmtn) buffm
!!      do k=1,nmtvr
!!        write(200) buffm(:,:,k)
!!      enddo
      ENDIF
      DO k=1,nmtvr
       call split2d_phys(buffm(1,1,k),buffo,global_lats_r)
       CALL interpred_phys(1,kmsk0,buffo,buff2,global_lats_r,
     &                lonsperlar)
       HPRIME(k,:,:)=buff2(:,:)
      ENDDO
 
!my jordan's mb
!sela  print *, ' (*j*)  nmtvr= ',nmtvr, 'reading hprime'
!my      DO j=1,lats_node_r
!my      DO i=1,lonr
!my      DO k=1,NMTVR
!my        IF(SLMSK(i,j).NE.1.) HPRIME(k,i,j) = 0.
!my      ENDDO
!my      ENDDO
!my      ENDDO
 

 
      IF (iozondp.eq.1) CALL readoz_disprd(ozplin)
!
!     reading the grib orography and scattering the data
!
      if(needoro.eq.1) then

      IF( me==0) then
        CALL ORORD(101,lonr,latr,buff1)
      endif
      call split2d_phys(buff1,buffo,global_lats_r)
      CALL interpred_phys(1,kmsk0,buffo,oro,global_lats_r,lonsperlar)
      endif
      RETURN
      END


      SUBROUTINE read_sfc_nemsio(sfc_fld,NEEDORO,nread,
     &                    cfile,global_lats_r,lonsperlar)
!
!***********************************************************************
!
!      use sfcio_module, ONLY: sfcio_head, sfcio_data, sfcio_realfill,
!     &                        sfcio_srohdc, sfcio_axdata
      use resol_def,    ONLY: latr, latr2, lonr, lsoil
      use layout1,      ONLY: me, nodes, lats_node_r
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use namelist_soilveg ,       only: salp_data, snupx
      use physcons,     only : tgice => con_tice
      USE machine,      ONLY: kind_io4, kind_io8
      use module_nemsio
      implicit none
!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)

      integer jump
      integer needoro

      real(kind=kind_io4) buff1(lonr*latr),buff2(lonr,latr,LSOIL)
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      real(kind=kind_io8) buff3(lonr,lats_node_r)
      integer nread,i,j,k,ij,idate7(7),lonsfc,latsfc,lplsfc(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr),kmskcv(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2,timef,rsnow
      real(4) fhour4
      type(nemsio_gfile) gfile_in
      integer iret, vegtyp,lonb4,latb4,nsoil4,ivs4
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()

      if(me==0) print *,' nread=',nread,' cfile=',cfile
      call nemsio_init()
!
      call nemsio_open(gfile_in,trim(cfile),'read',iret=iret)
!
      IF (me==0) THEN

        call nemsio_getheadvar(gfile_in,'fhour',fhour4,iret=iret)
        call nemsio_getheadvar(gfile_in,'lonb',lonb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'latb',latb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'nsoil',nsoil4,iret=iret)
        call nemsio_getheadvar(gfile_in,'ivs',ivs4,iret=iret)
        call nemsio_getheadvar(gfile_in,'idate',idate7,iret=iret)

!        PRINT 99,nread,head%fhour,head%idate,
!     &           head%lonb,head%latb,head%lsoil,head%ivs,iret
        PRINT 99,nread,fhour4,idate7(1:4),
     &           lonb4,latb4,nsoil4,ivs4,iret
99      FORMAT(1H ,'in fixio nread=',i3,2x,'HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4),4x,'lonsfc,latsfc,lsoil,ivssfc,iret=',5i8)

        if(iret.ne.0) goto 5000
        if(lonb4.ne.lonr) goto 5000
        if(latb4.ne.latr) goto 5000
        if(nsoil4.ne.lsoil) goto 5000

      ENDIF

      kmsk=0
!
      if(me==0) call nemsio_readrecv(gfile_in,'tmp','sfc',1,buff1,
     &    iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TSEA,
     &    global_lats_r,lonsperlar)

      DO K=1, LSOIL

        if(me==0) call nemsio_readrecv(gfile_in,'smc','soil layer',k,
     &      buff1,iret=iret)
        call split2d_phys(buff1, buffo,global_lats_r)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        sfc_fld%SMC(k,:,:)=buff3(:,:)
      ENDDO

      if(me==0) call nemsio_readrecv(gfile_in,'weasd','sfc',1,buff1,
     &    iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHELEG,
     &               global_lats_r,lonsperlar)

      DO K = 1, LSOIL
        if(me==0) call nemsio_readrecv(gfile_in,'stc','soil layer',k,
     &      buff1,iret=iret)
        call split2d_phys(buff1, buffo,global_lats_r)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        sfc_fld%STC(k,:,:)=buff3(:,:)
      ENDDO

      if(me==0) call nemsio_readrecv(gfile_in,'tg3','sfc',1,buff1,
     & iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TG3,
     &    global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'sfcr','sfc',1,buff1,
     &   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ZORL,
     &    global_lats_r,lonsperlar)

      sfc_fld%cv  = 0
      sfc_fld%cvb = 0
      sfc_fld%cvt = 0

      if(me==0) call nemsio_readrecv(gfile_in,'alvsf','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%alvsf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVSF,
     &               global_lats_r,lonsperlar)
      if(me==0) call nemsio_readrecv(gfile_in,'alvwf','sfc',1,buff1,
     &  iret=iret)
!      if(me==0) buff1=data%alvwf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVWF,
     &               global_lats_r,lonsperlar)
      if(me==0) call nemsio_readrecv(gfile_in,'alnsf','sfc',1,buff1,
     &    iret=iret)
!      if(me==0) buff1=data%alnsf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNSF,
     &               global_lats_r,lonsperlar)
      if(me==0) call nemsio_readrecv(gfile_in,'alnwf','sfc',1,buff1,
     &    iret=iret)
!      if(me==0) buff1=data%alnwf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNWF,
     &               global_lats_r,lonsperlar)

!     The mask cannot be interpolated
      if(me==0) call nemsio_readrecv(gfile_in,'land','sfc',1,buff1,
     &    iret=iret)
!      if(me==0) buff1=data%slmsk
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLMSK,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'veg','sfc',1,buff1,
     &    iret=iret)
!      if(me==0) buff1=data%vfrac
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VFRAC,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'cnwat','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%canopy
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%CANOPY,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'f10m','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%f10m
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%F10M,
     &    global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'vtype','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%vtype
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VTYPE,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'sotyp','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%stype
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%STYPE,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'facsf','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%facsf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACSF,
     &               global_lats_r,lonsperlar)
      if(me==0) call nemsio_readrecv(gfile_in,'facwf','sfc',1,buff1,
     &     iret=iret)
!      if(me==0) buff1=data%facwf
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACWF,
     &               global_lats_r,lonsperlar)

!szunyogh 06/16/99
      if(me==0) call nemsio_readrecv(gfile_in,'fricv','sfc',1,buff1,
     &     iret=iret)
!        if(me==0) buff1=data%uustar
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%UUSTAR,
     &               global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'ffmm','sfc',1,buff1,
     &     iret=iret)
!        if(me==0) buff1=data%ffmm
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFMM,
     &                  global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'ffhh','sfc',1,buff1,
     &    iret=iret)
!        if(me==0) buff1=data%ffhh
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFHH,
     &                  global_lats_r,lonsperlar)

!c-- XW: FOR SEA-ICE Nov04
!    Sea-ice (hice/fice) was added to the surface files.

      if(me==0) call nemsio_readrecv(gfile_in,'icetk','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%hice
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%HICE,
     &                  global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'icec','sfc',1,buff1,
     &    iret=iret)
!         if(me==0) buff1=data%fice
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FICE,
     &                  global_lats_r,lonsperlar)

      if(me==0) call nemsio_readrecv(gfile_in,'tisfc','sfc',1,buff1,
     &    iret=iret)
!         if(me==0) buff1=data%tisfc
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TISFC,
     &                  global_lats_r,lonsperlar)
         if (sfc_fld%tisfc(1,1) < 0.0)  then
           DO j=1,lats_node_r
             DO i=1,LONR
                sfc_fld%TISFC(i,j)= sfc_fld%TSEA(i,j)
                IF(sfc_fld%SLMSK(i,j) >=  2. .AND.
     &             sfc_fld%FICE(i,j)  >= 0.5) THEN
                   sfc_fld%TISFC(i,j) = (sfc_fld%TSEA(i,j)
     &            -tgice*(1.-sfc_fld%FICE(i,j))) / sfc_fld%FICE(i,j)
                  sfc_fld%TISFC(i,j)=MIN(sfc_fld%TISFC(i,j),tgice)
                ENDIF
             ENDDO
           ENDDO
         endif

!c-- XW: END SEA-ICE

!lu   11/10/2004
!*     surface files for GFS/Noah contain 8 additional records:
!*     tprcp, srflag, snwdph, slc, shdmin, shdmax, slope, snoalb

      if(me==0) call nemsio_readrecv(gfile_in,'tprcp','sfc',1,buff1,
     &    iret=iret)
!         if(me==0) buff1=data%tprcp
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TPRCP,
     &                  global_lats_r,lonsperlar)

!* srflag
      if(me==0) call nemsio_readrecv(gfile_in,'crain','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%srflag
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SRFLAG,
     &                  global_lats_r,lonsperlar)

!* snwdph
      if(me==0) call nemsio_readrecv(gfile_in,'snod','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%snwdph
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNWDPH,
     &                  global_lats_r,lonsperlar)

!* slc
         DO K=1, LSOIL
!         if(me==0) buff1=data%slc(:,:,k)
      if(me==0) call nemsio_readrecv(gfile_in,'slc','soil layer',k,
     &   buff1,iret=iret)
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,buff3,
     &    global_lats_r,lonsperlar)
         sfc_fld%SLC(k,:,:) = buff3(:,:)
         ENDDO

!* shdmin
      if(me==0) call nemsio_readrecv(gfile_in,'shdmin','sfc',1,buff1,
     &    iret=iret)
!         if(me==0) buff1=data%shdmin
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMIN,
     &                  global_lats_r,lonsperlar)

!* shdmax
      if(me==0) call nemsio_readrecv(gfile_in,'shdmax','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%shdmax
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMAX,
     &                  global_lats_r,lonsperlar)

!* slope
      if(me==0) call nemsio_readrecv(gfile_in,'sltyp','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%slope
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLOPE,
     &                  global_lats_r,lonsperlar)

!* snoalb
      if(me==0) call nemsio_readrecv(gfile_in,'salbd','sfc',1,buff1,
     &     iret=iret)
!         if(me==0) buff1=data%snoalb
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNOALB,
     &                  global_lats_r,lonsperlar)
!     print *,' snoalb=',sfc_fld%snoalb(1,:)
!lu [+67L]: the addition of 8 Noah records ends here .........................

       if(needoro.eq.1) then
         if(me==0) then
           call nemsio_readrecv(gfile_in,'orog','sfc',1,buff1,iret=iret)
!           buff1=data%orog
           needoro=1
           if(all(buff1.ne.-9999.)) needoro=0
           print *,'read sfc orography'
         endif
         call split2d_phys(buff1, buffo,global_lats_r)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%ORO,
     &                  global_lats_r,lonsperlar)
         call skip(needoro)
       endif
!
!Wei initialize snow fraction(sheleg is in mm)
      DO j=1,lats_node_r
        DO i=1,LONR
          sfc_fld%SNCOVR(i,j) = 0.0
          if (sfc_fld%slmsk(i,j) > 0.001 .AND. 
     &        ABS(sfc_fld%VTYPE(i,j)) >= 0.5 ) then
            vegtyp = sfc_fld%VTYPE(i,j)
            RSNOW  = 0.001*sfc_fld%SHELEG(i,j)/SNUPX(vegtyp)
            IF (0.001*sfc_fld%SHELEG(i,j) < SNUPX(vegtyp)) THEN
              sfc_fld%SNCOVR(i,j) = 1.0 - ( EXP(-SALP_DATA*RSNOW)
     &                                    - RSNOW*EXP(-SALP_DATA))
            ELSE
              sfc_fld%SNCOVR(i,j) = 1.0
            ENDIF
!           if (i == 1)
!    &       print*,SNUPX(vegtyp),SALP_DATA,sfc_fld%SNCOVR(i,j),
!    &       '************debug',sfc_fld%SHELEG(i,j),vegtyp,' j=',j
!    &,      ' snoalb1=',sfc_fld%snoalb(i,j)
!
          endif
        ENDDO
       ENDDO
!

       IF (me==0) then
!         call sfcio_axdata(data,iret)
         t2=timef()
         print *,'FIXIO TIME ',t2-t1,t1,t2
       endif
!
      call nemsio_close(gfile_in,iret=iret)
!
      call nemsio_finalize()
!
      RETURN
 5000 PRINT *, ' error in input in routine read_sfc'
      STOP
      END
!
!***********************************************************************
!
      subroutine interpred_phys(iord,kmsk,f,fi,global_lats_r,lonsperlar)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: ipt_lats_node_r, lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!!
      integer              global_lats_r(latr)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(in):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(out):: fi(lonr,lats_node_r)
      integer j,lons,lat
!!
      do j=1,lats_node_r
          lat=global_lats_r(ipt_lats_node_r-1+j)
          lons=lonsperlar(lat)
          if(lons.ne.lonr) then
            call intlon_phys(iord,1,1,lonr,lons,
     &                  kmsk(1,j),f(1,j),fi(1,j))
cjfe        fi(lons+1:lonr,j)=-9999.e9
            fi(lons+1:lonr,j)=0.
          else
            fi(:,j)=f(:,j)
          endif
        enddo
      end subroutine
c
c***********************************************************************
c
      subroutine intlon_phys(iord,imon,imsk,m1,m2,k1,f1,f2)
      use machine, ONLY: kind_io8
      implicit none
      integer,intent(in):: iord,imon,imsk,m1,m2
      integer,intent(in):: k1(m1)
      real (kind=kind_io8),intent(in):: f1(m1)
      real (kind=kind_io8),intent(out):: f2(m2)
      integer i2,in,il,ir
      real (kind=kind_io8) r,x1
      r=real(m1)/real(m2)
      do i2=1,m2
         x1=(i2-1)*r
         il=int(x1)+1
         ir=mod(il,m1)+1
          if(iord.eq.2.and.(imsk.eq.0.or.k1(il).eq.k1(ir))) then
            f2(i2)=f1(il)*(il-x1)+f1(ir)*(x1-il+1)
          else
            in=mod(nint(x1),m1)+1
            f2(i2)=f1(in)
          endif
      enddo
      end subroutine
c
c**********************************************************************
c
      SUBROUTINE readoz_disprd(ozplin)
 
      use ozne_def, ONLY: latsozp, levozp, timeoz, pl_coeff, kozpl
      USE machine,  ONLY: kind_phys, kind_io4
      implicit none
!!
      integer n,k,kk,i
      real (kind=kind_phys) ozplin(latsozp,levozp,pl_coeff,timeoz)
      real(kind=kind_io4) tempin(latsozp)
!
      DO I=1,timeoz
        do n=1,pl_coeff
          DO k=1,levozp
            READ(kozpl) tempin
            ozplin(:,k,n,i) = tempin(:)
          ENDDO
        enddo
      ENDDO
 
      RETURN
      END
c
c***********************************************************************
c
      SUBROUTINE ORORD(LUGB,IORO,JORO,ORO)
!
      use layout1, ONLY: me
      USE machine, ONLY: kind_io4, kind_io8
      implicit none
!!
      integer lugb, ioro, joro, kpdoro, ior, jor, i,k
      CHARACTER*80 FNOROG
!
      real (kind=kind_io4) oro(ioro,joro)
      real (kind=kind_io8) orog(ioro,joro), blnm, bltm
      logical gausm
!
      FNOROG = 'orography'
      kpdoro = 8
      IOR    = IORO
      JOR    = JORO
      CALL FIXRDG(LUGB,IOR,JOR,FNOROG,
     &            KPDORO,OROG,GAUSM,BLNM,BLTM,me)
!
      if (ior .ne. ioro .or. jor .ne. joro) then
         print *,' orography file not o.k. run aborted'
         call abort
      endif
      ORO = OROG
!
      RETURN
      END
c
c***********************************************************************
c
      subroutine split2d_phys(x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,     ONLY: latr, lonr
      use layout1,       ONLY: me, nodes, lats_node_r, ipt_lats_node_r
      use mpi_def,       ONLY: info, mpi_r_io, mpi_comm_all
      USE machine,       ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_io4) tmp(lonr,latr)
      integer global_lats_r(latr)
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
      integer ifld/0/
      save ifld
      real t1,t2,t3,t4,timef,ta,tb
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
      XL=0.
      ifld=ifld+1
      IF (me==0) THEN
!
!         Sending the data
!         ----------------
!-- do not need to send data, all processores read the data
         tmp=0.
         do j=1,latr
            do i=1,lonr
              tmp(i,j)=X(i,j)
            enddo
         enddo
      ENDIF
      call mpi_bcast
     1 (tmp,lonr*latr,MPI_R_IO,0,MPI_COMM_ALL,info)
       call mpi_barrier(mpi_comm_all,info)
!-- get subdomain of data
        do j=1,lats_node_r
           lat=global_lats_r(ipt_lats_node_r-1+j)
           do i=1,lonr
              xl(i,j)=tmp(i,lat)
           enddo
        enddo
      return
      end
c
c***********************************************************************
c
      SUBROUTINE skip(jump)
 
c*************************************************************************
 
      use resol_def
      use layout1
      use mpi_def
      implicit none
 
      integer jump,ipe
 
      ipe=0
 
      CALL MPI_BCAST(jump,1,MPI_INTEGER,ipe,MPI_COMM_ALL,info)
 
      RETURN
      END
!
c
c***********************************************************************
c
      SUBROUTINE EXCHA(lats_nodes_r,global_lats_r,X1,X2,Y1,Y2)
c
c***********************************************************************
c
      use resol_def,  ONLY: latr
      use layout1,    ONLY: nodes, lats_node_r_max, lats_node_r,
     &                      ipt_lats_node_r
      use mpi_def,    ONLY: mc_comp, mpi_r_def
      USE machine,    ONLY: kind_io8
      implicit none
 
      integer n,i,j,ierr,ilat,lat,node,nsend
      integer              global_lats_r(latr)
      integer              lats_nodes_r(nodes)
      real(kind=kind_io8) X1(lats_node_r),X2(lats_node_r)
      real(kind=kind_io8) Y1(latr),Y2(latr)
cjfe  real(kind=kind_mpi) tmps(2,lats_node_r_max,nodes)
cjfe  real(kind=kind_mpi) tmpr(2,lats_node_r_max,nodes)
      real(kind=kind_io8) tmps(2,lats_node_r_max,nodes)
      real(kind=kind_io8) tmpr(2,lats_node_r_max,nodes)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      if (nodes.ne.1) then
        do node=1,nodes
          do i=1,lats_node_r
           lat=global_lats_r(ipt_lats_node_r-1+i)
           tmps(1,i,node)=X1(I)
           tmps(2,i,node)=X2(I)
          enddo
        enddo
!!
        nsend=2*lats_node_r_max
cjfe    call mpi_alltoall(tmps,nsend,MPI_R_MPI,
cjfe x                     tmpr,nsend,MPI_R_MPI,
cjfe x                     MC_COMP,ierr)
        call mpi_alltoall(tmps,nsend,MPI_R_DEF,
     x                     tmpr,nsend,MPI_R_DEF,
     x                     MC_COMP,ierr)
!!
        ilat=1
        do node=1,nodes
          do i=1,lats_nodes_r(node)
             lat=global_lats_r(ilat)
             Y1(lat)=tmpr(1,i,node)
             Y2(lat)=tmpr(2,i,node)
             ilat=ilat+1
          enddo
        enddo
!!
      ELSE
        Y1=X1
        Y2=X2
      ENDIF
!!
      RETURN
      END
c
c***********************************************************************
c
      SUBROUTINE SUMLAT(n,X,nodes)
c
c***********************************************************************
c
      use mpi_def,   ONLY: MC_COMP, MPI_R_DEF, info, mpi_sum
      USE machine,   ONLY: kind_io8, kind_io4
      implicit none
 
      integer n,i,j,np,mr,nodes
      real(kind=kind_io8) X(n),Y(N)
      real(kind=kind_io4) Z(n)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      if (nodes.ne.1) then
        DO i=1,n
          Y(i)=X(i)
        ENDDO
        CALL mpi_allreduce(Y,X,n,MPI_R_DEF,MPI_SUM,
     &                    MC_COMP   ,info)
      endif
        DO i=1,n
          Z(i)=X(i)
        ENDDO
        DO i=1,n
          X(i)=Z(i)
        ENDDO
!!
      RETURN
      END
c
c***********************************************************************
c
      subroutine unsplit2d_phys(ioproc,x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max,
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: info, mpi_comm_all, liope, mpi_r_io,
     &                       stat
      USE machine,     ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_io4) tmp(lonr,latr+2)
      integer global_lats_r(latr),ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,ncc
      data ncc/0/
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
      X=0.
      maxfld=50
      ifldu=ifldu+1
!!
!jw all fcst node need to send data
!jw IF (me.ne.ioproc) THEN
c
c         Sending the data
c         ----------------
!jw         tmp=0.
!jw         tmp(lonr,latr+1)=ipt_lats_node_r
!jw         tmp(lonr,latr+2)=lats_node_r
!jw         do j=1,lats_node_r
!jw            do i=1,lonr
!jw              tmp(i,j)=XL(i,j)
!jw            enddo
!jw         enddo
!jw         if (.NOT.LIOPE) then
!jw           nodesr=nodes
!jw         else
!jw           nodesr=nodes+1
!jw         endif
!jw         msgtag=1000+(me+1)*nodesr*maxfld+ifldu
!jw          call MPI_SEND(tmp(lonr,latr+1),1,MPI_R_IO,ioproc,
!jw     &                  msgtag,MPI_COMM_ALL,info)
!jw          call MPI_SEND(tmp(lonr,latr+2),1,MPI_R_IO,ioproc,
!jw     &                  msgtag,MPI_COMM_ALL,info)
!jw         illen=tmp(lonr,latr+2)
c send the local grid domain
!jw         CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO,ioproc,
!jw     &                  msgtag,MPI_COMM_ALL,info)
!jw      ELSE
!!
!!     for pes ioproc
!jw        if (.NOT.LIOPE) then
!jw          nproct=nodes
!jw          do j=1,lats_node_r
!jw             lat=global_lats_r(ipt_lats_node_r-1+j)
!jw             do i=1,lonr
!jw                x(i,lat)=XL(i,j)
!jw             enddo
!jw          enddo
!jw        else
!jw          nproct=nodes-1
!jw        endif
!jw        DO proc=1,nproct
!jw         if (proc.ne.ioproc+1) then
!jw         msgtag=1000+proc*nodes*maxfld+ifldu
!jw          CALL mpi_recv(tmp(lonr,latr+1),1,MPI_R_IO,proc-1,
!jw     &                msgtag,MPI_COMM_ALL,stat,info)
!jw          CALL mpi_recv(tmp(lonr,latr+2),1,MPI_R_IO,proc-1,
!jw     &                msgtag,MPI_COMM_ALL,stat,info)
!jw         illen=tmp(lonr,latr+2)
!jw          CALL mpi_recv(tmp(1,1),illen*lonr ,MPI_R_IO,proc-1,
!jw     &                msgtag,MPI_COMM_ALL,stat,info)
!jw         if (.NOT.LIOPE) then
!jw           ipt_lats_node_rl=tmp(lonr,latr+1)
!jw           lats_nodes_rl=tmp(lonr,latr+2)
!jw         else
!jw           ipt_lats_node_rl=tmp(lonr,lats_node_r_max+1)
!jw           lats_nodes_rl=tmp(lonr,lats_node_r_max+2)
!jw         endif
!jw         do j=1,lats_nodes_rl
!jw           lat=global_lats_r(ipt_lats_node_rl-1+j)
!jw           do i=1,lonr
!jw              x(i,lat)=tmp(i,j)
!jw           enddo
!jw         enddo
!jw         endif   !(proc.ne.ioproc+1)
!jw        enddo
!!
!jw      ENDIF
!jw         ncc=ncc+1
 
!!
      return
      end
c
c***********************************************************************
c
      subroutine uninterpred(iord,kmsk,f,fi,global_lats_r,lonsperlar)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!!
      integer              global_lats_r(latr)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(out):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(in):: fi(lonr,lats_node_r)
      integer j,lons,lat
!!
      do j=1,lats_node_r
          lat=global_lats_r(ipt_lats_node_r-1+j)
          lons=lonsperlar(lat)
          if(lons.ne.lonr) then
            call intlon_phys(iord,1,1,lons,lonr,
     &                  kmsk(1,j),fi(1,j),f(1,j))
          else
            f(:,j)=fi(:,j)
          endif
        enddo
      end subroutine



      subroutine uninterprez(iord,kmsk,f,fi,global_lats_r,lonsperlar,  
     &    buff_mult_piecea)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io4,kind_io8
      implicit none
!!
      integer,intent(in):: global_lats_r(latr)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(out):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(in):: fi(lonr,lats_node_r)
       real(kind=4) f4(lonr,lats_node_r)
      integer j,lons,lat
      integer i,ubound
!
      real(kind=kind_io4),intent(inout) :: buff_mult_piecea
     &  (1:lonr,1:lats_node_r)
!!
      do j=1,lats_node_r
          lat=global_lats_r(ipt_lats_node_r-1+j)
          lons=lonsperlar(lat)
          if(lons.ne.lonr) then
            call intlon_phys(iord,1,1,lons,lonr,
     &                  kmsk(1,j),fi(1,j),f(1,j))
          f4(:,j)=fi(:,j)
          else
            f(:,j)=fi(:,j)
            f4(:,j)=fi(:,j)
          endif
        enddo
      do j=1,lats_node_r
      do i=1,lonr
      buff_mult_piecea(i,j)=f (i,j)
      end do
      end do
      end subroutine



       subroutine unsplit2z(ioproc,ngridx,ngridt,x,global_lats_r)
c
c***********************************************************************
c
      use resol_def,   ONLY: lonr,latr
      use mod_state,   ONLY: ivar_global_a, buff_mult_pieces
      use layout1,     ONLY: me, nodes_comp
      use mpi_def,     ONLY: liope
      USE machine,     ONLY: kind_io4
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      real(kind=kind_io4) tmp(lonr,latr+2)
      integer global_lats_r(latr),ipt_lats_node_rl,nodesr,ngridx,ngridt
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,nd1,nd2
       character*8 cna
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
      write(cna,985)600+ngridx
 985   format('fort.',i3)
      X=0.
      maxfld=50
      ifldu=ifldu+1
!!
      IF (me.ne.ioproc) THEN
            continue
      ELSE
!!
!!     for pes ioproc
        nproct=nodes_comp
        nd1=0
        DO proc=1,nproct
          ipt_lats_node_rl=ivar_global_a(1,proc)
          lats_nodes_rl=ivar_global_a(2,proc)
          nd2=nd1+lonr*lats_nodes_rl*(ngridx-1)
          do j=1,lats_nodes_rl
           lat=global_lats_r(ipt_lats_node_rl-1+j)
           do i=1,lonr
             x(i,lat)=buff_mult_pieces(nd2+i+(j-1)*lonr)
           enddo
          enddo
          nd1=nd1+lonr*lats_nodes_rl*ngridt
        enddo

!!
      ENDIF
!!
      return
      end
 
c
c***********************************************************************
c
      subroutine unsplit2d_phys_r(ioproc,x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max, 
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: liope, info, stat, mpi_comm_all, 
     &                       mpi_r_io_r
      USE machine,     ONLY: kind_ior, kind_io8
      implicit none
!!
      real(kind=kind_ior) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_ior) tmp(lonr,latr+2)
      integer global_lats_r(latr),ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,ncc
      data ncc/0/
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
!     X=0.               ! commented by moorthi on 20051117
      maxfld=50
      ifldu=ifldu+1
!!
      IF (me.ne.ioproc) THEN
c
c         Sending the data
c         ----------------
         tmp=0.
         tmp(lonr,latr+1)=ipt_lats_node_r
         tmp(lonr,latr+2)=lats_node_r
         do j=1,lats_node_r
            do i=1,lonr
              tmp(i,j)=XL(i,j)
            enddo
         enddo
         if (.NOT.LIOPE) then
           nodesr=nodes
         else
           nodesr=nodes+1
         endif
         msgtag=1000+(me+1)*nodesr*maxfld+ifldu
          call MPI_SEND(tmp(lonr,latr+1),1,MPI_R_IO_R,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
          call MPI_SEND(tmp(lonr,latr+2),1,MPI_R_IO_R,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
         illen=tmp(lonr,latr+2)
c send the local grid domain
         CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO_R,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
!!
!!     for pes ioproc
        x = 0.0               ! added by Moorthi on 2005111700
        if (.NOT.LIOPE) then
          nproct=nodes
          do j=1,lats_node_r
             lat=global_lats_r(ipt_lats_node_r-1+j)
             do i=1,lonr
                x(i,lat)=XL(i,j)
             enddo
          enddo
        else
          nproct=nodes-1
        endif
        DO proc=1,nproct
         if (proc.ne.ioproc+1) then
         msgtag=1000+proc*nodes*maxfld+ifldu
          CALL mpi_recv(tmp(lonr,latr+1),1,MPI_R_IO_R,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
          CALL mpi_recv(tmp(lonr,latr+2),1,MPI_R_IO_R,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
         illen=tmp(lonr,latr+2)
          CALL mpi_recv(tmp(1,1),illen*lonr ,MPI_R_IO_R,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
         if (.NOT.LIOPE) then
           ipt_lats_node_rl=tmp(lonr,latr+1)
           lats_nodes_rl=tmp(lonr,latr+2)
         else
           ipt_lats_node_rl=tmp(lonr,lats_node_r_max+1)
           lats_nodes_rl=tmp(lonr,lats_node_r_max+2)
         endif
         do j=1,lats_nodes_rl
           lat=global_lats_r(ipt_lats_node_rl-1+j)
           do i=1,lonr
              x(i,lat)=tmp(i,j)
           enddo
         enddo
         endif   !(proc.ne.ioproc+1)
        enddo
!!
      ENDIF
         ncc=ncc+1
 
!!
      return
      end
c
c***********************************************************************
c
      subroutine split2d_phys_r(x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,      ONLY: latr, lonr
      use layout1,        ONLY: me, lats_node_r, ipt_lats_node_r, nodes
      use mpi_def,        ONLY: liope, mpi_comm_all, info,mpi_r_io_r
      USE machine,        ONLY: kind_ior, kind_io8
      implicit none
!!
      real(kind=kind_ior) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_ior) tmp(lonr,latr)
      integer global_lats_r(latr)
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifld/0/
      save ifld
      real t1,t2,t3,t4,timef,ta,tb
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
      XL=0.
!     maxfld=50
      ifld=ifld+1
!!
      IF (me.eq.0) THEN
        ta=timef()
        t3=ta
c        DO proc=1,nodes-1
         do proc=1,1
c
c         Sending the data
c         ----------------
         tmp=0.
         do j=1,latr
            do i=1,lonr
              tmp(i,j)=X(i,j)
            enddo
         enddo
!Moor    msgtag=1000+proc*nodes*maxfld+ifld
         t1=timef()
!sela    print *,' GWVX BROADCASTING FROM ',nodes-1
         call mpi_bcast
     1 (tmp,lonr*latr,MPI_R_IO_R,nodes-1,MPI_COMM_ALL,info)
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
c         CALL mpi_send(tmp,lonr*latr,MPI_R_IO_R,proc-1,msgtag,
c     &                  MPI_COMM_ALL,info)
         t2=timef()
!sela    print 102,t2-t1
 
 102    format(' SEND TIME ',f10.5)
        enddo
        t4=timef()
      ELSE
        if (.NOT.LIOPE) then
          nodesr=nodes
        else
          nodesr=nodes+1
        endif
!Moor   msgtag=1000+(me+1)*nodesr*maxfld+ifld
!sela    print *,' GWVX BROADCASTREC  FROM ',nodesr-1
         call mpi_bcast
     1 (tmp,lonr*latr,MPI_R_IO_R,nodesr-1,MPI_COMM_ALL,info)
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
!sela    print *,'GWVX IPT ',ipt
c        CALL mpi_recv(tmp,lonr*latr,MPI_R_IO_R,nodesr-1,
c     &                msgtag,MPI_COMM_ALL,stat,info)
        do j=1,lats_node_r
           lat=global_lats_r(ipt_lats_node_r-1+j)
           do i=1,lonr
              xl(i,j)=tmp(i,lat)
           enddo
        enddo
!!
      ENDIF
!!
!!     for pes nodes-1
      if (.NOT.LIOPE) then
        if (me.eq.nodes-1) then
          do j=1,lats_node_r
             lat=global_lats_r(ipt_lats_node_r-1+j)
             do i=1,lonr
                xl(i,j)=X(i,lat)
             enddo
          enddo
        endif
      endif
!!
      tb=timef()
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
 
!sela  if(icolor.eq.2.and.me.eq.nodes-1)print 103,tb-ta,t4-t3
 103  format(' GLOBAL AND SEND TIMES  split2d_phys',2f10.5)
      return
      end

!
c***********************************************************************
c
      subroutine split2d_rst(x,xl,fieldsize,global_lats_r,lonsperlar)
c
c***********************************************************************
c
      use resol_def,      ONLY: latr, lonr
      use layout1,        ONLY: me, lats_node_r, ipt_lats_node_r, nodes
      use mpi_def,        ONLY: liope, mpi_comm_all, info,mpi_r_io_r
      USE machine,        ONLY: kind_ior, kind_io8
      implicit none
!!
!!
      integer,intent(in) :: fieldsize,global_lats_r(latr),
     &                      lonsperlar(latr)
      real(kind=kind_ior),intent(in) :: x(fieldsize)
      real (kind=kind_io8),intent(inout) :: xl(lonr,lats_node_r)
      integer j,lat,i,lon
!      real t1,t2,t3,t4,timef,ta,tb
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
!--- get subdomain of data
       do j=1,lats_node_r
           lat=global_lats_r(ipt_lats_node_r-1+j)
           if(lat/=1) then
             lon=sum(lonsperlar(1:lat-1))
           else
             lon=0
           endif
!
           do i=1,lonsperlar(lat)
              xl(i,j)=X(lon+i)
           enddo
      enddo
!!

!sela  if(icolor.eq.2.and.me.eq.nodes-1)print 103,tb-ta,t4-t3
 103  format(' GLOBAL AND SEND TIMES  split2d_phys',2f10.5)
      return
      end subroutine split2d_rst


!***********************************************************************
!
      SUBROUTINE read_sfc_r(cfile,sfc_fld,phy_f2d,phy_f3d,num_p3d,
     &           num_p2d,NGPTC,NBLCK,global_lats_r,lonsperlar,NEEDORO)
!
!***********************************************************************
!
      use resol_def,      ONLY: latr, lonr, latr2, lsoil,levs
      use layout1,        ONLY: me, nodes, lats_node_r,ipt_lats_node_r
      USE machine,        ONLY: kind_ior, kind_io8, kind_rad

      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use namelist_soilveg ,       only: salp_data, snupx
      use physcons,                only : tgice => con_tice
      use module_nemsio
!
      implicit none
!
      character(*),intent(in) :: cfile
      TYPE(Sfc_Var_Data),intent(inout) :: sfc_fld
      integer,intent(in)            :: global_lats_r(latr)
      integer,intent(in)            :: lonsperlar(latr)
      integer,intent(in)            :: num_p2d,num_p3d,NGPTC,NBLCK
      real(kind=kind_rad),intent(inout) ::
     &    phy_f2d(lonr,lats_node_r,num_p2d),
     &    phy_f3d(NGPTC,LEVS,NBLCK,lats_node_r,num_p3d)
      integer,intent(inout) :: needoro
!
      integer jump

      real(kind=kind_io8) buff3(lonr,lats_node_r)
!
      real(kind=kind_ior),allocatable :: buff1(:)
!
      integer i,j,k,im,jm,idate(4),lplsfc(latr2)
      real t1,t2,timef,rsnow
!---
      type(nemsio_gfile) :: gfile
      integer iret, vegtyp,fieldsize,iblk,il,lons_lat,njeff,l,lat,lon
      character*2 nump2d,nump3d
      character(255) varname
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()
!
      call nemsio_init()
!
      call nemsio_open(gfile,trim(cfile),'read',iret=iret)
      print *,'after nemsio_open, iret=',iret
      if(iret/=0) then
        PRINT *, ' ERROR in input routine read_sfc_r'
        return
      endif
!
      call nemsio_getfilehead(gfile,dimx=im,dimy=jm,iret=iret)
      fieldsize=im*jm
      allocate(buff1(fieldsize))
!
!-- tsea
      call nemsio_readrecv(gfile,'tmp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%TSEA,fieldsize,global_lats_r,
     &  lonsperlar)
!-- smc
      DO K=1, LSOIL
        call nemsio_readrecv(gfile,'smc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1, sfc_fld%smc(k,:,:),fieldsize,
     &    global_lats_r,lonsperlar)
        print *,'read inrst,smc=',sfc_fld%smc(k,1:5,1:5)
      ENDDO

!-- sheleg
      call nemsio_readrecv(gfile,'weasd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%SHELEG,fieldsize,global_lats_r,
     &  lonsperlar)
!--stc
      DO K = 1, LSOIL
        call nemsio_readrecv(gfile,'stc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1, sfc_fld%stc(k,:,:),fieldsize,
     &    global_lats_r,lonsperlar)
        print *,'read inrst,stc=',sfc_fld%stc(k,1:5,1:5)
      ENDDO

!--tg3
      call nemsio_readrecv(gfile,'tg3','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tg3,fieldsize,global_lats_r,
     &  lonsperlar)
        print *,'read inrst,tg3=',sfc_fld%tg3(1:3,1:3)
!--zorl
      call nemsio_readrecv(gfile,'sfcr','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%zorl,fieldsize,global_lats_r,
     &  lonsperlar)
!
      sfc_fld%cv  = 0
      sfc_fld%cvb = 0
      sfc_fld%cvt = 0
        print *,'read inrst,cwafter cvt'

!-- alvsf
      call nemsio_readrecv(gfile,'alvsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alvsf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- alvwf
      call nemsio_readrecv(gfile,'alvwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alvwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- alnsf
      call nemsio_readrecv(gfile,'alnsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alnsf,fieldsize,global_lats_r,
     &  lonsperlar)
!--alnwf
      call nemsio_readrecv(gfile,'alnwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alnwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slmsk
      call nemsio_readrecv(gfile,'land','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%slmsk,fieldsize,global_lats_r,
     &  lonsperlar)

!-- vfrac
      call nemsio_readrecv(gfile,'veg','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%vfrac,fieldsize,global_lats_r,
     &  lonsperlar)
!-- canopy
      call nemsio_readrecv(gfile,'cnwat','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%canopy,fieldsize,global_lats_r,
     &  lonsperlar)
!-- f10m
      call nemsio_readrecv(gfile,'f10m','10 m above gnd',1,buff1,
     &   iret=iret)
      call split2d_rst(buff1,sfc_fld%f10m,fieldsize,global_lats_r,
     &  lonsperlar)
!--vtype
      call nemsio_readrecv(gfile,'vtype','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%vtype,fieldsize,global_lats_r,
     &  lonsperlar)
!-- stype
      call nemsio_readrecv(gfile,'sotyp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%stype,fieldsize,global_lats_r,
     &  lonsperlar)
!-- facsf
      call nemsio_readrecv(gfile,'facsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%facsf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- facwf
      call nemsio_readrecv(gfile,'facwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%facwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- uustar (fricv)
      call nemsio_readrecv(gfile,'fricv','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%uustar,fieldsize,global_lats_r,
     &  lonsperlar)
!-- ffhh
      call nemsio_readrecv(gfile,'ffhh','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%ffhh,fieldsize,global_lats_r,
     &  lonsperlar)
!-- ffmm
      call nemsio_readrecv(gfile,'ffmm','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%ffmm,fieldsize,global_lats_r,
     &  lonsperlar)
!-- hice
      call nemsio_readrecv(gfile,'icetk','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%hice,fieldsize,global_lats_r,
     &  lonsperlar)
!-- fice
      call nemsio_readrecv(gfile,'icec','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%fice,fieldsize,global_lats_r,
     &  lonsperlar)
!-- tisfc
      call nemsio_readrecv(gfile,'tisfc','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tisfc,fieldsize,global_lats_r,
     &  lonsperlar)
        print *,'read inrst,tisfc=',sfc_fld%tisfc(1:3,1:3)
      if (sfc_fld%tisfc(1,1) < 0.0) then
        DO j=1,lats_node_r
          DO i=1,LONR
             sfc_fld%TISFC(i,j) = sfc_fld%TSEA(i,j)
             IF(sfc_fld%SLMSK(i,j) >=  2. .AND.
     &          sfc_fld%FICE(i,j)  >= 0.5) THEN
                sfc_fld%TISFC(i,j) = (sfc_fld%TSEA(i,j)
     &         -tgice*(1.-sfc_fld%FICE(i,j))) / sfc_fld%FICE(i,j)
                sfc_fld%TISFC(i,j) = MIN(sfc_fld%TISFC(i,j),tgice)
             ENDIF
          ENDDO
        ENDDO
      endif
!-- tprcp
      call nemsio_readrecv(gfile,'tprcp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tprcp,fieldsize,global_lats_r,
     &  lonsperlar)
!-- srflag (crain)
      call nemsio_readrecv(gfile,'crain','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%srflag,fieldsize,global_lats_r,
     &  lonsperlar)
!-- snwdph
      call nemsio_readrecv(gfile,'snod','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%SNWDPH,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slc
      DO K=1, LSOIL
        call nemsio_readrecv(gfile,'slc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1,sfc_fld%slc(k,:,:),fieldsize,
     &    global_latS_r,lonsperlar)
        print *,'read inrst,slc=',sfc_fld%slc(k,1:3,1:3)
      ENDDO
!-- shdmin
      call nemsio_readrecv(gfile,'shdmin','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%shdmin,fieldsize,global_lats_r,
     &  lonsperlar)
!-- shdmax
      call nemsio_readrecv(gfile,'shdmax','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%shdmax,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slope (sltyp)
      call nemsio_readrecv(gfile,'sltyp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%slope,fieldsize,global_lats_r,
     &  lonsperlar)
!-- salbd
      call nemsio_readrecv(gfile,'salbd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%SNOALB,fieldsize,global_lats_r,
     &  lonsperlar)
        print *,'read inrst,snoalb=',sfc_fld%snoalb(1:3,1:3)
!-- orog
      if(needoro.eq.1) then
        call nemsio_readrecv(gfile,'orog','sfc',1,buff1,iret=iret)
        needoro=1
        if(any(buff1.eq.-9999.)) needoro=0
        print *,'read sfc orography'
        call split2d_rst(buff1,sfc_fld%oro,fieldsize,global_lats_r,
     &  lonsperlar)
        call skip(needoro)
      endif
        print *,'read inrst,after orog'
!jw read sncovr from rstart file
!-- read in snow cover from restart file
      sfc_fld%SNCOVR = 0.0
      call nemsio_readrecv(gfile,'sncovr','sfc',1,buff1,iret=iret)
      if(iret==0)
     &call split2d_rst(buff1,sfc_fld%sncovr,fieldsize,global_lats_r,
     &  lonsperlar)
        print *,'read inrst,snoalb=',sfc_fld%sncovr(38,3),
     &    sfc_fld%SHELEG(38,3)
!
!-- num_p2d
      DO K=1, num_p2d
        write(nump2d,'(I2.2)')k
        varname='phyf2d_'//nump2d
        call nemsio_readrecv(gfile,trim(varname),'sfc',1,buff1,
     &    iret=iret)
        print *,'read inrst,',trim(varname),'iret=',iret
        call split2d_rst(buff1,phy_f2d(:,:,k),fieldsize,global_lats_r,
     &    lonsperlar)
      ENDDO
!
!-- num_p3d
      DO K=1, num_p3d
        write(nump3d,'(I2.2)')k
        varname='phyf3d_'//nump3d
        DO L=1, levs
          call nemsio_readrecv(gfile,trim(varname),'mid layer',L,
     &      buff1,iret=iret)
        print *,'read inrst,phy_p3d,',trim(varname),'iret=',iret
          call split2d_rst(buff1,buff3,fieldsize,global_lats_r,
     &    lonsperlar)
!
          do j=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+j)
            lons_lat = lonsperlar(lat)
            iblk=0
            il=1
            do lon=1,lons_lat,NGPTC
              NJEFF=MIN(NGPTC,lons_lat-lon+1)
              iblk=iblk+1
              do i=1,NJEFF
                phy_f3d(i,l,iblk,j,k)=buff3(il,j)
                il=il+1
              enddo
            enddo
          enddo
!
        ENDDO
      ENDDO

      call nemsio_close(gfile)
      call nemsio_finalize()
!
      t2=timef()
      print *,'FIXIO TIME ',t2-t1,t1,t2
!
      RETURN

      STOP
      END
!
!***********************************************************************
!
