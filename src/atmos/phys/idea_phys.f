!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_phys(im,ix,levs,prsi,prsl,                        &
     &adu,adv,adt,adr,nbdsw,nbdlw,ntrac,dtp,lat,                        &
     &solhr,slag,sdec,cdec,sinlat,coslat,                               &
     &xlon,xlat,oro,cozen,swhb,hlwb,dt6dt,                              &
     &thermodyn_id,sfcpress_id,gen_coord_hybrid)
!-----------------------------------------------------------------------
! add temp, wind changes due to viscosity and thermal conductivity
! also solar heating
! Apr 06 2012   Henry Juang, initial implement for NEMS
!-----------------------------------------------------------------------
      use physcons,  amo2=>con_amo2,avgd => con_avgd
      use idea_composition
!
      implicit none
! Argument
      integer, intent(in) :: im  ! number of data points in adt (first dim)
      integer, intent(in) :: ix  ! max data points in adt (first dim)
      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: lat    ! latitude index
      integer, intent(in) :: ntrac  ! number of tracer
      integer, intent(in) :: nbdsw  ! number of sw band
      integer, intent(in) :: nbdlw  ! number of lw band
      real,    intent(in) :: dtp    ! time step in second
      real, intent(in)    :: prsi(ix,levs+1) ! pressure
      real, intent(in)    :: prsl(ix,levs)   ! pressure
      real, intent(in)    :: hlwb(ix,levs,nbdlw)   ! long wave rad (K/s)
      real, intent(in)    :: swhb(ix,levs,nbdsw)   ! short wave rad (K/s)
      real, intent(in)    :: cozen(im)   ! time avg(1 hour) cos zenith angle
      real, intent(in)    :: oro(im)   ! surface height (m)
      real, intent(in) :: solhr,slag,sdec,cdec ! for solar zenith angle
      real, intent(in) :: xlon(im),xlat(im),coslat(im),sinlat(im)
      real, intent(inout) :: adr(ix,levs,ntrac) ! tracer
      real, intent(inout) :: adt(ix,levs)    ! temperature
      real, intent(inout) :: adu(ix,levs)    ! real u
      real, intent(inout) :: adv(ix,levs)    ! real v
      real, intent(inout) :: dt6dt(ix,levs,6)! diagnostic array 
      integer thermodyn_id, sfcpress_id
      logical gen_coord_hybrid
! Local variables
      real cp(ix,levs),cospass(im),dt(ix,levs),rtime1,hold1,n(ix,levs) 
      real  o_n(ix,levs),o2_n(ix,levs),n2_n(ix,levs),o3_n(ix,levs),     &
     & am(ix,levs),dudt(ix,levs),dvdt(ix,levs),dtdt(ix,levs),xmu(im),   &
     & dtco2c(ix,levs),dtco2h(ix,levs)                                  &
     &,dth2oh(ix,levs),dth2oc(ix,levs),dto3(ix,levs),rho(ix,levs)       &
     &,wvl(ix,levs),wvm(ix,levs),wvs(ix,levs),zg(ix,levs)               &
     &,amin,amax,grav(ix,levs)                                          &
     &,prslk(ix,levs),prsik(ix,levs+1),phil(ix,levs),phii(ix,levs+1)
! solar
      real utsec,sda,maglat(im),maglon(im),btot(im),                    &
     &dipang(im),essa(im),dlat,dlon
      integer i,k,dayno,j1,j2
! change to real windl !hmhj already real wind
!hmhj do i=1,im
!hmhj   adu(i,1:levs)=adu(i,1:levs)/coslat(i)
!hmhj   adv(i,1:levs)=adv(i,1:levs)/coslat(i)
!hmhj enddo ! i
! get phil geopotential from temperature
!hmhj call GET_PHI_gc_h(im,ix,levs,ntrac,adt,adr,prsi,phii,phil)
      call get_phi(im,ix,levs,ntrac,adt,adr,                            &
     &             thermodyn_id, sfcpress_id,                           &
     &             gen_coord_hybrid,                                    &
     &             prsi,prsik,prsl,prslk,phii,phil)
! get height
      call phi2z(im,ix,levs,phil,oro,zg,grav)
!     print*,'wwwz',zg(1,1:150)
!     print*,'wwwg',grav(1,1:150)
!     print*,'wwwp',phil(1,1:150),oro(1)
! 
! get composition at layers (/cm3) and rho (kg/m3)
      call idea_tracer(im,ix,levs,ntrac,2,grav,prsi,prsl,adt,adr,       &
     &dtp,o_n,o2_n,n2_n,n,rho,am)
! calculate cp
      call getcp_idea(im,ix,levs,ntrac,adr,cp,                          &
     &                thermodyn_id,gen_coord_hybrid)
! dissipation
      call idea_phys_dissipation(im,ix,levs,grav,prsi,prsl,             &
     &adu,adv,adt,o_n,o2_n,n2_n,dtp,cp,dt6dt)
!
! get cos solar zenith angle (instant)
      call presolar(im,ix,solhr,slag,                                   &
     &              sinlat,coslat,sdec,cdec,xlon,xlat                   &
     &              ,cospass,dayno,utsec,sda                            &
     &              ,maglat,maglon,btot,dipang,essa)
! get solar heating and NO cooling then get temp adjustment
      call idea_sheat(im,ix,levs,adt,dt,cospass,o_n,o2_n,n2_n,rho,      &
     &cp,lat,dayno,prsl,zg,grav,am,maglat,dt6dt)
!     rtime1=3600.*6.
      do i=1,im
        do k=1,levs
          adt(i,k)=adt(i,k)+dt(i,k)*dtp
!         dt3dt(i,k,1)=dt(i,k)*rtime1
        enddo !k
      enddo ! i
!
! ion_drag
! change to /m3
      o_n=o_n*1.e6
      o2_n=o2_n*1.e6
      n2_n=n2_n*1.e6
      n=n*1.e6
      call idea_ion(solhr,cospass,zg,o_n,o2_n,n2_n,cp,                  &
     &adu,adv,adt,dudt,dvdt,dtdt,rho,xlat,xlon,ix,im,levs,              &
     &dayno,utsec,sda,maglon,maglat,btot,dipang,essa) 
      do i=1,im
!      dlat=xlat(i)*180./3.14159
!      dlon=xlon(i)*180./3.14159
!      if(abs(dlat-60.).le.1..and.abs(dlon-270.).le.1.) then
!      print*,'www0',solhr,dudt(i,140)*dtp,dvdt(i,140)*dtp,             &
!    &dtdt(i,140)*dtp,adu(i,140),adv(i,140)
!      endif
       do k=1,levs
         adu(i,k)=adu(i,k)+dtp*dudt(i,k)
         adv(i,k)=adv(i,k)+dtp*dvdt(i,k)
         adt(i,k)=adt(i,k)+dtdt(i,k)*dtp
       enddo
      enddo
! change u,V back !hmhj no need to change back, they are real wind
!hmhj do i=1,im
!hmhj   adu(i,1:levs)=adu(i,1:levs)*coslat(i)
!hmhj   adv(i,1:levs)=adv(i,1:levs)*coslat(i)
!hmhj enddo ! i
! radiation
! co2 cooling, heating
      call idea_co2(im,ix,levs,nlev_co2,ntrac,grav,cp,adr,adt,          &
     &              dtco2c,cospass,dtco2h)
!hmhj&'/mtb/save/wx20fw/fcst07rd',dtco2c,cospass,dtco2h)
! h2o cooling heating 110-41 down ward
      call idea_h2o(im,ix,levs,nlev_h2o,nlevc_h2o,ntrac,grav,cp,        &
     &adr,adt,dth2oh,cospass,dth2oc)
         dt6dt(1:im,1:levs,4)=dtco2c
!        dt6dt(1:im,1:levs,5)=dth2oc
!        dt6dt(1:im,1:levs,6)=dth2oh
!     dth2oc=0.
!     dth2oh=0.
! o2 o3 heating
      call o3pro(im,ix,levs,ntrac,adr,am,n,o3_n)
      call idea_o2_o3(im,ix,levs,cospass,adt,o2_n,o3_n,rho,cp,          &
     &zg,grav,dto3)
! get xmu
      do i=1,im
        if(cospass(i).gt.0.0001.and.cozen(i).gt.0.0001) then
         xmu(i)=cospass(i)/cozen(i)
        else
         xmu(i)=0.
        endif
      enddo
! dt6dt
!     do i=1,im
!     do k=1,levs
!     dt6dt(i,k,2)=dtco2c(i,k)
!     dt6dt(i,k,3)=dtco2h(i,k)
!     dt6dt(i,k,4)=dth2oc(i,k)
!     dt6dt(i,k,5)=dth2oh(i,k)
!     dt6dt(i,k,6)=dto3(i,k)
!     enddo
!     enddo
! merge
      call rad_merge(im,ix,levs,nbdsw,nbdlw,hlwb,swhb,wvl,              &
     &wvm,wvs,xmu,dtco2c,dtco2h,dth2oh,dth2oc,dto3,dt6dt)
      do i=1,im
        do k=1,levs
          adt(i,k)=adt(i,k)+(wvl(i,k)+wvm(i,k)+wvs(i,k))*dtp
        enddo !k
      enddo ! i
      return
      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getcp_idea(im,ix,levs,ntrac,adr,xcp,                   & 
     &                      thermodyn_id,gen_coord_hybrid)
!
      use tracer_const
!hmhj use resol_def , only: thermodyn_id
!hmhj use namelist_def , only: gen_coord_hybrid
!
      implicit none
      integer, intent(in) :: im  ! number of data points in adr (first dim)
      integer, intent(in) :: ix  ! max data points in adr (first dim)
      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: ntrac  ! number of tracer
      real, intent(in) :: adr(ix,levs,ntrac) ! tracer kg/kg
      real, intent(out) :: xcp(ix,levs) !CP (J/kg/k)
      integer thermodyn_id
      logical gen_coord_hybrid
!
! local
      real sumq(ix,levs),work1
      integer i,j,k
        sumq=0.0					
        xcp=0.0					
      if( gen_coord_hybrid .and. thermodyn_id.eq.3 ) then		
        do i=1,ntrac			
         if( cpi(i).ne.0.0 ) then					
         do k=1,levs							
          do j=1,im							
            work1=adr(j,k,i)						
            sumq(j,k)=sumq(j,k)+work1				
            xcp(j,k)=xcp(j,k)+cpi(i)*work1				
          enddo								
         enddo								
         endif								
        enddo							
        do k=1,levs							
         do j=1,im							
           xcp(j,k)=(1.-sumq(j,k))*cpi(0)+xcp(j,k)
         enddo						
        enddo	
      else
        do i=4,ntrac			
         do k=1,levs							
          do j=1,im							
            work1=adr(j,k,i)						
            sumq(j,k)=sumq(j,k)+work1					
            xcp(j,k)=xcp(j,k)+cpi(i)*work1				
          enddo								
         enddo								
        enddo							
        do k=1,levs							
         do j=1,im							
           xcp(j,k)=(1.-sumq(j,k))*cpi(0)+xcp(j,k)				
         enddo								
        enddo								
      endif
      return
      end
      subroutine rad_merge(im,ix,levs,nbdsw,nbdlw,hlwb,swhb,wvl,        &
     &wvm,wvs,xmu,dtco2c,dtco2h,dth2oh,dth2oc,dto3,dt6dt)
!
      use idea_composition
!     use module_radsw_parameters,   only : NBDSW
!     use module_radlw_parameters,   only : NBDLW
      implicit none
      integer, intent(in) :: im  ! number of data points in hlw,dt..(first dim)
      integer, intent(in) :: ix  ! max data points in hlw,... (first dim)
      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: nbdsw   ! number of sw band
      integer, intent(in) :: nbdlw   ! number of lw band
      real, intent(in) :: hlwb(ix,levs,nbdlw) ! GFS lw rad (K/s) 16band
      real, intent(in) :: swhb(ix,levs,nbdsw) ! GFS sw rad (K/s) 11band
      real, intent(in) :: xmu(im) ! mormalized cos zenith angle
      real, intent(in) :: dtco2c(ix,levs)  ! idea co2 cooling(K/s)
      real, intent(in) :: dtco2h(ix,levs)  ! idea co2 heating(K/s)
      real, intent(in) :: dth2oc(ix,levs)  ! idea h2o cooling(K/s)
      real, intent(in) :: dth2oh(ix,levs)  ! idea h2o heating(K/s)
      real, intent(in) :: dto3(ix,levs)    ! idea o3 heating(K/s)
      real, intent(out) :: wvl(ix,levs)  ! GFS idea combined lw rad
      real, intent(out) :: wvm(ix,levs)  ! GFS idea combined mid rad
      real, intent(out) :: wvs(ix,levs)  ! GFS idea combined sw rad
      real, intent(inout) :: dt6dt(ix,levs,6)  
!     local
      real dx
      integer i,k,j,levb(2),levt(2)
      levb(1)=k81
      levb(2)=k47
      levt(1)=k87
      levt(2)=k64
!     data levb/81,47/
!     data levt/87,64/
!
      do i=1,im
        do k=1,levs
           wvl(i,k)=0.
           wvm(i,k)=0.
           wvs(i,k)=0.
! lw 1-11
           do j=1,11
              wvl(i,k)=wvl(i,k)+hlwb(i,k,j)
           enddo
! mid (lw 12-16 sw 9-11)
           do j=12,16
              wvm(i,k)=wvm(i,k)+hlwb(i,k,j)
           enddo
           do j=9,11
              wvm(i,k)=wvm(i,k)+swhb(i,k,j)*xmu(i)
           enddo
! sw 1-8
           do j=1,5
              wvs(i,k)=wvs(i,k)+swhb(i,k,j)*xmu(i)*ef(k)
           enddo
           do j=6,8
              wvs(i,k)=wvs(i,k)+swhb(i,k,j)*xmu(i)
           enddo
        enddo
!       dt6dt(i,1:levs,1)=wvl(i,1:levs)
!       dt6dt(i,1:levs,2)=dth2oc(i,1:levs)+dtco2c(i,1:levs)
!       dt6dt(i,1:levs,3)=wvm(i,1:levs)
!       dt6dt(i,1:levs,4)=dth2oh(i,1:levs)+dtco2h(i,1:levs)
!       dt6dt(i,1:levs,5)=wvs(i,1:levs)
!       dt6dt(i,1:levs,6)=dto3(i,1:levs)
! merge wvl,wvs
        dx=prlog(levt(1))-prlog(levb(1))
        do k=levb(1)+1,levt(1)-1
           wvl(i,k)=(wvl(i,k)*(prlog(levt(1))-prlog(k))+                &
     &     (dth2oc(i,k)+dtco2c(i,k))*(prlog(k)-prlog(levb(1))))/dx
           wvs(i,k)=(wvs(i,k)*(prlog(levt(1))-prlog(k))+                &
     &            dto3(i,k)*(prlog(k)-prlog(levb(1))))/dx
        enddo
        do k=levt(1),levs
           wvl(i,k)=dtco2c(i,k)+dth2oc(i,k)
           wvs(i,k)=dto3(i,k)
        enddo
! merge wvm
        dx=prlog(levt(2))-prlog(levb(2))
        do k=levb(2)+1,levt(2)-1                                        &
           wvm(i,k)=(wvm(i,k)*(prlog(levt(2))-prlog(k))+
     &     (dth2oh(i,k)+dtco2h(i,k))*(prlog(k)-prlog(levb(2))))/dx
        enddo
        do k=levt(2),levs
           wvm(i,k)=dtco2h(i,k)+dth2oh(i,k)
        enddo
!       dt6dt(i,1:levs,3)=wvl(i,1:levs)
!       dt6dt(i,1:levs,4)=wvm(i,1:levs)
!       dt6dt(i,1:levs,5)=wvs(i,1:levs)
      enddo
      return
      end
      subroutine getmax(ain,n1,n,m,amin,j1,amax,j2)
      real ain(n1,m)
      amin=1.e36
      amax=-1.e36
      i1=500
      j1=500
      i2=500
      j2=500
      do i=1,n
      do j=1,m
      if(amin.gt.ain(i,j)) then
      amin=ain(i,j)
      i1=i
      j1=j
      endif
      if(amax.lt.ain(i,j)) then
      amax=ain(i,j)
      i2=i
      j2=j
      endif
      enddo
      enddo
      return
      end
      subroutine getmax2(ain,ain1,n1,n,m,amax,j2)
      real ain(n1,m),ain1(n1,m)
      amax=-1.e36
      i1=500
      j1=500
      i2=500
      j2=500
      do i=1,n
      do j=1,m
      sq=sqrt(ain(i,j)**2+ain1(i,j)**2)
      if(amax.lt.sq) then
      amax=sq
      i2=i
      j2=j
      endif
      enddo
      enddo
      return
      end
      subroutine phi2z(im,ix,levs,phi,soro,z,grav)

! Subroutine to calculate geometric height and gravity from geopotential
! in a hydrostatic atmosphere, assuming a spherically symmetric planet
! and Newton's gravity.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! File history

! Feb 26, 2010: Rashid Akmaev
! Loosely based on Hojun Wang's phi2hgt but generalized to rid of
! recursive calculations, include surface orography, and calculate 
! gravity.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Define constants
! - Earth radius (m) and 
! - gravity at sea level (m/s**2) 

! If used with GFS/WAM codes "use" this module
      use physcons, only: re => con_rerth, g0 => con_g

      implicit none

! If the module is not available, comment out the "use" line above and
! uncomment this line
!      real, parameter:: re = 6.3712e+6, g0 = 9.80665e+0

      real, parameter:: g0re = g0*re, g0re2 = g0*re**2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine parameters
! INPUTS
! - array dimensions (following GFS conventios): first actual, first 
! maximum, number of levels

      integer, intent(in):: im,ix,levs

! - geopotential (m**2/s**2)
! - surface orography (m)

      real, intent(in):: phi(ix,levs),soro(im)

! OUTPUTS
! - height (m)
      
      real, intent(out):: z(ix,levs)

! Optional output
! - gravity (m/s**2)

!     real, intent(out), optional:: grav(ix,levs)
      real, intent(out):: grav(ix,levs)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Local variables

      integer:: i,l
      real:: phis(im)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculate surface geopotential

      do i = 1,im
         phis(i)=g0re*soro(i)/(re+soro(i))
      enddo

! Calculate height

      z(:,:) = 0.
      do l = 1,levs
         do i = 1,im
            z(i,l)=re*(phis(i)+phi(i,l))/(g0re-(phis(i)+phi(i,l)))
         enddo
      enddo

! ***Optionally*** calculate gravity

!     if(present(grav)) then
         grav(:,:) = 0.
         do l = 1,levs
            do i = 1,im
               grav(i,l)=g0re2/(re+z(i,l))**2
            enddo
         enddo
!     endif
      end subroutine phi2z
