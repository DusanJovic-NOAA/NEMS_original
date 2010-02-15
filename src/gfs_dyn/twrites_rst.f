      subroutine twrites_rst(fname,IOPROC,fhour,idate,
     x           si,ls_nodes,max_ls_nodes,
     x           gze_ls,qe_ls,tee_ls,die_ls,zee_ls,rqe_ls,
     x           gzo_ls,qo_ls,teo_ls,dio_ls,zeo_ls,rqo_ls)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  write spectral variables for restart
!-------------------------------------------------------------------
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def	
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use module_nemsio
!
      implicit none
!
      character(*),intent(in) :: fname
      integer,intent(in) :: ioproc
      real(kind=kind_evod),intent(in) :: fhour
      integer,intent(in) :: idate(4)
!
      real(kind=kind_evod),intent(in):: si(levp1)
      integer,intent(in)             :: ls_nodes(ls_dim,nodes)
      integer,intent(in)             :: max_ls_nodes(nodes)
!
      real(kind=kind_evod),intent(in) :: gze_ls(len_trie_ls,2)
     &,                     qe_ls(len_trie_ls,2)
     &,                    tee_ls(len_trie_ls,2,levs)
     &,                    die_ls(len_trie_ls,2,levs)
     &,                    zee_ls(len_trie_ls,2,levs)
     &,                    rqe_ls(len_trie_ls,2,levh)
!
     &,                    gzo_ls(len_trio_ls,2)
     &,                     qo_ls(len_trio_ls,2)
     &,                    teo_ls(len_trio_ls,2,levs)
     &,                    dio_ls(len_trio_ls,2,levs)
     &,                    zeo_ls(len_trio_ls,2,levs)
     &,                    rqo_ls(len_trio_ls,2,levh)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer              ierr,j,k,l,lenrec,locl,n,node,idate7(7)
!
      integer              indjoff
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      integer              ixgr
!
      real(kind=kind_ior), target ::   buf(lnt2)
      real(kind=kind_ior),allocatable ::   Z_R(:)
      real(kind=kind_ior)   tmps(4+nodes+jcap1*nodes)
      real(kind=kind_ior)   tmpr(3+nodes+jcap1*(nodes-1))
!
      type(nemsio_gfile) gfile
      integer nmetavarr
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:)
      character(16),allocatable :: varrname(:)
      real,allocatable :: varrval(:)
      integer iret, idvt
      integer il,ilen,i,msgtag,ls_diml,nrec
      integer nfhour,nfminute,nfsecondn,nfsecondd,nframe,jrec,nmeta
      logical first
      save first,nmetavarr,recname,reclevtyp, 
     &     reclev,varrname,varrval,nmeta,nrec
      save Z_R
      data first /.true./
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
      integer              joff, jj
!
      include 'function2'
!
      joff(n,l)=(jcap1)*(jcap2)-(jcap1-l)*(jcap2-l)+2*(n-l)
!
      real(kind=kind_mpi_r),allocatable :: trieo_ls_node (:,:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_ls_nodes(:,:,:,:)
!
      real(kind=kind_mpi_r),allocatable :: trieo_gz_node (:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_gz_nodes(:,:,:)
!
      integer      lan,lat,iblk,lons_lat,lon,NJEFF,nn,lv
!
!---------------------------------------------------------------------
!
!      print *,' enter twrites_rst ' 

      call mpi_comm_size(MPI_COMM_ALL,i,ierr)
!
!-- allocate
      allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     x                            2, 3*levs+1*levh+1 ) )
!
!-- compute gze_ls only once
      if(first ) then
        allocate(trieo_gz_node(len_trie_ls_max+len_trio_ls_max,2))
        do j=1,len_trie_ls
          trieo_gz_node(j,1) = gze_ls(j,1)
          trieo_gz_node(j,2) = gze_ls(j,2)
        enddo
        do j=1,len_trio_ls
          trieo_gz_node(j+len_trie_ls_max,1) = gzo_ls(j,1)
          trieo_gz_node(j+len_trie_ls_max,2) = gzo_ls(j,2)
        enddo
      endif
!
!-- collect data
      do j=1,len_trie_ls
        trieo_ls_node(j,1,kwq) = qe_ls(j,1)
        trieo_ls_node(j,2,kwq) = qe_ls(j,2)
      enddo
!
      do j=1,len_trio_ls
        trieo_ls_node(j+len_trie_ls_max,1,kwq) = qo_ls(j,1)
        trieo_ls_node(j+len_trie_ls_max,2,kwq) = qo_ls(j,2)
      enddo
!
      do k=1,levs
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwte+  k-1) = tee_ls(j,1,k)
          trieo_ls_node(j,2,kwte+  k-1) = tee_ls(j,2,k)
          trieo_ls_node(j,1,kwdz+  k-1) = die_ls(j,1,k)
          trieo_ls_node(j,2,kwdz+  k-1) = die_ls(j,2,k)
          trieo_ls_node(j,1,kwdz+levs+k-1) = zee_ls(j,1,k)
          trieo_ls_node(j,2,kwdz+levs+k-1) = zee_ls(j,2,k)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwte+  k-1) = teo_ls(j,1,k)
          trieo_ls_node(jj,2,kwte+  k-1) = teo_ls(j,2,k)
          trieo_ls_node(jj,1,kwdz+  k-1) = dio_ls(j,1,k)
          trieo_ls_node(jj,2,kwdz+  k-1) = dio_ls(j,2,k)
          trieo_ls_node(jj,1,kwdz+levs+k-1) = zeo_ls(j,1,k)
          trieo_ls_node(jj,2,kwdz+levs+k-1) = zeo_ls(j,2,k)
        enddo
      enddo
!
      do k=1,levh
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwrq+  k-1) = rqe_ls(j,1,k)
          trieo_ls_node(j,2,kwrq+  k-1) = rqe_ls(j,2,k)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwrq+  k-1) = rqo_ls(j,1,k)
          trieo_ls_node(jj,2,kwrq+  k-1) = rqo_ls(j,2,k)
        enddo
      enddo
!
!-- collect data to ioproc
      if ( me .eq. ioproc ) then
!       write(0,*)'ALLOC PARMS TWRITE ',len_trie_ls_max+len_trio_ls_max,
!     &      2,3*levs+1*levh+1, nodes,1
 
         allocate ( trieo_ls_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, 3*levs+1*levh+1, nodes ),
     &              stat=ierr )
         if(first) then
           allocate ( trieo_gz_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, nodes ),stat=ierr )
         endif
      endif
      if (ierr .ne. 0) then
        write (0,*) ' GWX trieo_ls_nodes allocate failed'
        call mpi_abort(mpi_comm_all,ierr,i)
       endif
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
      if(nodes >1 )then

        lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * 
     &   (3*levs+1*levh+1)
!
        t1=rtc()
       call mpi_gather( trieo_ls_node , lenrec, MPI_R_MPI_R,
     x                 trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     x                 ioproc, MPI_COMM_ALL, ierr)
       call mpi_barrier(MPI_COMM_ALL,ierr)
       if(first) then
         lenrec=(len_trie_ls_max+len_trio_ls_max) * 2 
         call mpi_gather( trieo_gz_node , lenrec, MPI_R_MPI_R,
     x                 trieo_gz_nodes, lenrec, MPI_R_MPI_R,
     x                 ioproc, MPI_COMM_ALL, ierr)
       endif
      else
       trieo_ls_nodes(:,:,:,1)=trieo_ls_node(:,:,:)
       if(first) then
        trieo_gz_nodes(:,:,1)=trieo_gz_node(:,:)
       endif
      endif
      deallocate ( trieo_ls_node  )
      if(first) then
        deallocate ( trieo_gz_node  )
        if ( me .eq. ioproc ) then
!
         allocate(Z_R(lnt2) )
!
         do node=1,nodes
          jbasev=0
          do locl=1,max_ls_nodes(node)
            l=ls_nodes(locl,node)
            indev1 = indlsev(L,L)
            if (mod(L,2).eq.mod(jcap+1,2)) then
              indev2 = indlsev(jcap-1,L)
            else
              indev2 = indlsev(jcap  ,L)
            endif
            indjoff=joff(l,l)
            do indev = indev1 , indev2
              Z_R(indjoff+1) = trieo_gz_nodes(indev,1,node)
              Z_R(indjoff+2) = trieo_gz_nodes(indev,2,node)
              indjoff=indjoff+4
            end do
            jbasev=jbasev+(jcap+3-l)/2
           end do
!
           jbasod=len_trie_ls_max
           do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
!
              if ( l .ne. jcap ) then  ! fix for out of bound error
                indod1 = indlsod(L+1,L)
                if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
                else
                  indod2 = indlsod(jcap-1,L)
                endif
                indjoff=joff(l+1,l)
                do indod = indod1 , indod2
                  Z_R(indjoff+1) = trieo_gz_nodes(indod,1,node)
                  Z_R(indjoff+2) = trieo_gz_nodes(indod,2,node)
                  indjoff=indjoff+4
                end do
                jbasod=jbasod+(jcap+2-l)/2
              endif  ! fix for out of bound error
            end do
         end do
         deallocate(trieo_gz_nodes)
        endif    !end ioproc
      endif
!
      t2=rtc()
      call mpi_barrier(MPI_COMM_ALL,ierr)
!
!-- write out data
      IF (me.eq.ioproc) THEN
 
        print *,' in TWRITES fhour=',fhour
!
        if (first) then
!
          nmeta=5
          nrec=3*levs+1*levh+2
          ntrac=levh/levs
!          print *,'write spec rst, nrec=',nrec,'ntrac=',ntrac
!
!-- field infomation
          allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
          recname(1)='gz'
          recname(2)='pres'
          recname(3:levs+2)='tmp'
          recname(levs+3:2*levs+2)='di'
          recname(2*levs+3:3*levs+2)='ze'
          recname(3*levs+3:3*levs+levh+2)='rq'
          reclevtyp(1)='sfc'
          reclevtyp(2)='sfc'
          reclevtyp(3:3*levs+2)='mid layer'
          reclevtyp(3*levs+3:3*levs+levh+2)='tracer layer'
          reclev(1)=1
          reclev(2)=1
          do i=1,levs
            reclev(i+2)=i
            reclev(i+2+levs)=i
            reclev(i+2+2*levs)=i
          enddo
          do i=1,levh
            reclev(i+2+3*levs)=i
          enddo
!
          nmetavarr=1
          allocate(varrname(nmetavarr),varrval(nmetavarr))
          varrname(1:nmetavarr)=(/'fhour  '/)
!
!endof first
        endif
!
        idate7=0.;idate7(7)=1
        idate7(1)=idate(4)
        idate7(2:3)=idate(2:3)
        idate7(4)=idate(1)
        nfhour=int(fhour)
        nfminute=int((fhour-nfhour)*60)
        nfsecondn=int(((fhour-nfhour)*60-nfminute)*60)
        nfsecondd=1
!
        varrval(1)=fhour
!
        call nemsio_init()
!
        call nemsio_open(gfile,fname,'write',iret,modelname='GFS',   
     &    gdatatype='bin8',idate=idate7,nfhour=nfhour,nfminute=nfminute,
     &    nfsecondn=nfsecondn,nfsecondd=nfsecondd,dimx=lnt2,dimy=1,
     &    dimz=levs,nmeta=nmeta,nrec=nrec,jcap=jcap,ntrac=ntrac,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev,
     &    extrameta=.true.,
     &    nmetavarr=nmetavarr,varrname=varrname,varrval=varrval)
!        print *,'after nemsio_open,',trim(fname),'iret=',iret
!
!--- write out data
!
       jrec=1
       call nemsio_writerec(gfile,1,Z_R,iret=iret)
!       print *,'after snemsio_writerec gz,',maxval(Z_r),minval(Z_R),
!     &   'iret=',iret
!
       do k=1,3*levs+1*levh+1
         jrec=k+1
         do node=1,nodes
!
           jbasev=0
           do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
              indev1 = indlsev(L,L)
              if (mod(L,2).eq.mod(jcap+1,2)) then
                indev2 = indlsev(jcap-1,L)
              else
                indev2 = indlsev(jcap  ,L)
              endif
              indjoff=joff(l,l)
              do indev = indev1 , indev2
                buf(indjoff+1) = trieo_ls_nodes(indev,1,k,node)
                buf(indjoff+2) = trieo_ls_nodes(indev,2,k,node)
                indjoff=indjoff+4
              end do
              jbasev=jbasev+(jcap+3-l)/2
            end do
!
            jbasod=len_trie_ls_max
            do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
!
              if ( l .ne. jcap ) then  ! fix for out of bound error
                indod1 = indlsod(L+1,L)
                if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
                else
                  indod2 = indlsod(jcap-1,L)
                endif
                indjoff=joff(l+1,l)
                do indod = indod1 , indod2
                  buf(indjoff+1) = trieo_ls_nodes(indod,1,k,node)
                  buf(indjoff+2) = trieo_ls_nodes(indod,2,k,node)
                  indjoff=indjoff+4
                end do

                jbasod=jbasod+(jcap+2-l)/2
              endif  ! fix for out of bound error
            end do
          end do
!
          call nemsio_writerec(gfile,jrec,buf,iret=iret)
!        print *,'after snemsio_writerec,jrec=',jrec,'buf=',maxval(buf),
!     &     minval(buf),'iret=',iret
        end do
!
        t4=rtc  ()
!sela print *, ' DISK TIME FOR SIG TWRITEO WRT ',t4-t3
!
        deallocate ( trieo_ls_nodes )
!
        call nemsio_close(gfile)
        call nemsio_finalize()
!
!
      endif   !me.eq.ioproc
!!
      if(first) then
          first = .false.
      endif
      call mpi_barrier(MPI_COMM_ALL,ierr)

      return
      end
