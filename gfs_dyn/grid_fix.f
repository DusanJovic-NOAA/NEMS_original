      subroutine interpred(iord,kmsk,f,fi,global_lats_a,lonsperlat)
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
!!
      integer              global_lats_a(latg)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonf,lats_node_a)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_io8),intent(in):: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(out):: fi(lonf,lats_node_a)
      integer j,lons,lat
!!
      do j=1,lats_node_a
          lat=global_lats_a(ipt_lats_node_a-1+j)
          lons=lonsperlat(lat)
          if(lons.ne.lonf) then
            call intlon(iord,1,1,lonf,lons,
     &                  kmsk(1,j),f(1,j),fi(1,j))
cjfe        fi(lons+1:lonf,j)=-9999.e9
            fi(lons+1:lonf,j)=0.
          else
            fi(:,j)=f(:,j)
          endif
        enddo
      end subroutine
c
c***********************************************************************
c
      subroutine intlon(iord,imon,imsk,m1,m2,k1,f1,f2)
      use gfs_dyn_machine
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
      subroutine split2d(x,xl,global_lats_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!!
      real(kind=kind_io4) x(lonf,latg)
      real (kind=kind_io8) xl(lonf,lats_node_a)
      real(kind=kind_io4) tmp(lonf,latg)
      integer global_lats_a(latg)
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
!     maxfld=50
      ifld=ifld+1
!!
      IF (icolor.eq.2.and.me.eq.nodes-1) THEN
        ta=timef()
        t3=ta
c        DO proc=1,nodes-1
         do proc=1,1
c
c         Sending the data
c         ----------------
         tmp=0.
         do j=1,latg
            do i=1,lonf
              tmp(i,j)=X(i,j)
            enddo
         enddo
!Moor    msgtag=1000+proc*nodes*maxfld+ifld
         t1=timef()
!sela    print *,' GWVX BROADCASTING FROM ',nodes-1
         call mpi_bcast
     1 (tmp,lonf*latg,MPI_R_IO,nodes-1,MPI_COMM_ALL,info)
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
c         CALL mpi_send(tmp,lonf*latg,MPI_R_IO,proc-1,msgtag,
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
     1 (tmp,lonf*latg,MPI_R_IO,nodesr-1,MPI_COMM_ALL,info)
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
!sela    print *,'GWVX IPT ',ipt
c        CALL mpi_recv(tmp,lonf*latg,MPI_R_IO,nodesr-1,
c     &                msgtag,MPI_COMM_ALL,stat,info)
        do j=1,lats_node_a
           lat=global_lats_a(ipt_lats_node_a-1+j)
           do i=1,lonf
              xl(i,j)=tmp(i,lat)
           enddo
        enddo
!!
      ENDIF
!!
!!     for pes nodes-1
      if (.NOT.LIOPE) then
        if (me.eq.nodes-1) then
          do j=1,lats_node_a
             lat=global_lats_a(ipt_lats_node_a-1+j)
             do i=1,lonf
                xl(i,j)=X(i,lat)
             enddo
          enddo
        endif
      endif
!!
      tb=timef()
         call mpi_comm_rank(MPI_COMM_ALL,i,info)
 
!sela  if(icolor.eq.2.and.me.eq.nodes-1)print 103,tb-ta,t4-t3
 103  format(' GLOBAL AND SEND TIMES  SPLIT2D',2f10.5)
      return
      end
