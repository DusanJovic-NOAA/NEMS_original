      subroutine setindxh2o(latd,nlats,gaul,
     &                      jindx1_h,jindx2_h,ddy_h)
!
! May 2015 Shrinivas Moorthi - Prepare for H2O interpolation
!
      use machine , only : kind_phys
      use h2o_def , only : jh2o => latsh2o, h2o_lat, h2o_time
!
      implicit none
!
      integer                     latd,nlats
      integer, dimension(latd) :: jindx1_h, jindx2_h
      real(kind=kind_phys)     :: gaul(nlats),ddy_h(latd)
      integer i,j,lat
 
!
         do j=1,nlats
           jindx2_h(j) = jh2o + 1
           do i=1,jh2o
             if (gaul(j) < h2o_lat(i)) then
               jindx2_h(j) = i
               exit
             endif
           enddo
           jindx1_h(j) = max(jindx2_h(j)-1,1)
           jindx2_h(j) = min(jindx2_h(j),jh2o)
           if (jindx2_h(j) /= jindx1_h(j)) then
             ddy_h(j) = (gaul(j)            - h2o_lat(jindx1_h(j)))
     &              / (h2o_lat(jindx2_h(j)) - h2o_lat(jindx1_h(j)))
           else
             ddy_h(j) = 1.0
           endif
!     print *,' j=',j,' gaul=',gaul(j),' jindx12=',jindx1(j),
!    &jindx2(j),' h2o_lat=',h2o_lat(jindx1(j)),h2o_lat(jindx2(j))
!    &,' ddy=',ddy(j)
!csela if(me.eq.0) print*,'1st ddy(j,1) ddy(j,2),j=',ddy(j,1),ddy(j,2),j
 
         enddo
 
csela do j=1,nlats
csela if(me.eq.0) print*,'x1(j,1) jindx1(j,2)',jindx1(j,1),jindx1(j,2),j
csela if(me.eq.0) print*,'x2(j,1) jindx2(j,2)',jindx2(j,1),jindx2(j,2),j
csela enddo
csela do j=1,nlats
csela  if(me.eq.0) print*,'ddy(j,1) ddy(j,2)',ddy(j,1),ddy(j,2)
csela enddo
cyt   if(me.eq.0) print*,'completed setindxh2o for nasa prod. and diss'
 
      return
      end
!
!**********************************************************************
!
      subroutine h2ointerpol(me,latd,nlats,idate,fhour,
     &                       jindx1_h,jindx2_h,h2oplin,h2oplout,ddy_h)
!
! May 2015 Shrinivas Moorthi - Prepare for H2O interpolation
!
      use machine , only : kind_phys
      use h2o_def
      implicit none
      integer             j,j1,j2,l,latd,nc,n1,n2
      real(kind=kind_phys) fhour,tem, tx1, tx2
!
 
      integer  jindx1_h(latd), jindx2_h(latd)
      integer  me,idate(4),nlats
      integer  idat(8),jdat(8)
!
      real(kind=kind_phys) h2oplin(latsh2o,levh2o,h2o_coeff,timeh2o)
      real(kind=kind_phys) ddy_h(latd)
      real(kind=kind_phys) h2oplout(levh2o,latd,h2o_coeff)
      real(kind=kind_phys) rinc(5), rjday
      integer              jdow, jdoy, jday
      real(4)              rinc4(5)
      integer              w3kindreal, w3kindint
!
      idat    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc    = 0.
      rinc(2) = fhour
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4 = rinc
        CALL W3MOVDAT(RINC4,IDAT,JDAT)
      else
        CALL W3MOVDAT(RINC,IDAT,JDAT)
      endif
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      if (rjday < h2o_time(1)) rjday = rjday+365.
!
      n2 = timeh2o + 1
      do j=1,timeh2o
        if (rjday < h2o_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)      n1 = n1 + timeh2o
      if (n2 > timeh2o) n2 = n2 - timeh2o

!
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
!    &,'h2o_time=',h2o_time(n1),h2o_time(n2)
!

      tx1 = (h2o_time(n2) - rjday) / (h2o_time(n2) - h2o_time(n1))
      tx2 = 1.0 - tx1
!
      do nc=1,h2o_coeff
        do l=1,levh2o
          do j=1,nlats
            j1  = jindx1_h(j)
            j2  = jindx2_h(j)
            tem = 1.0 - ddy_h(j)
            h2oplout(l,j,nc) =
     &        tx1*(tem*h2oplin(j1,l,nc,n1)+ddy_h(j)*h2oplin(j2,l,nc,n1))
     &      + tx2*(tem*h2oplin(j1,l,nc,n2)+ddy_h(j)*h2oplin(j2,l,nc,n2))
          enddo
        enddo
      enddo
!
      return
      end
