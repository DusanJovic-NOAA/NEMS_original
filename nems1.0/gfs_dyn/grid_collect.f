      SUBROUTINE grid_collect
     &        (zsg,psg,uug,vvg,teg,rqg,dpg,
     &        global_lats_a,lonsperlat)
!
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!!
      INTEGER              global_lats_a(latg)
      INTEGER              lonsperlat(latg)
!!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) teg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
!
      real(kind=kind_io8) buffo(lonf,lats_node_a)
      real(kind=kind_io8) buffi(lonf,lats_node_a)
      integer kmsk(lonf,lats_node_a),kmskcv(lonf,lats_node_a)
      integer k,il
      integer ubound
      integer icount
       integer  ierr
!
      data  icount/0/
      integer maxlats_comp
!
      print *,' enter grid collect '

      ngridg=1
      if(allocated(buff_mult_pieceg)) then
         continue
      else
         allocate(buff_mult_pieceg(lonf,ngrids_gg,lats_node_a))
      endif
!
      kmsk = 0
!
      buffi(:,:) = zsg(:,:)
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
!
      buffi(:,:) = psg(:,:)
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
!
      do k=1,levs
        buffi(:,:) = dpg(:,:,k)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
      enddo
!
      do k=1,levs
        buffi(:,:) = uug(:,:,k)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
      enddo
!
      do k=1,levs
        buffi(:,:) = vvg(:,:,k)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
      enddo
!
      do k=1,levs
        buffi(:,:) = teg(:,:,k)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
      enddo
!
      if (levh .gt. 0) then
        do k=1,levh
          buffi(:,:) = rqg(:,:,k)
          CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
        enddo
      endif
!
      print *,' finished gridc_collect for  ngridg=',ngridg
  999 continue
      ngrid=1
      return
      end

       subroutine atmgg_move(ioproc)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!
      integer ipt_lats_node_al,nodesr
      integer lats_nodes_al
c     integer lats_nodes_a(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,ubound
      integer icount
      data icount/0/
      integer maxlats_comp
c  allocate the data structures
c
      if(icount .eq. 0) then
         allocate(ivarg_global(10))
         allocate(ivarg_global_a(10,nodes))
         ivarg_global(1)=ipt_lats_node_a
         ivarg_global(2)= lats_node_a
         ivarg_global(3)=lats_node_a_max
         call mpi_gather(ivarg_global,10,MPI_INTEGER,
     &       ivarg_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
         icount=icount+1
      endif
!     print *,' icount=',icount
!!
      if(allocated(buff_mult_piecesg)) then
          continue
      else
          maxlats_comp=lats_node_a_max
          if(.not. liope .or. me .ne. ioproc) then
            continue
          else
c           maxlats_comp=ivarg_global_a(3,ioproc)
            maxlats_comp=ivarg_global_a(3,1)
          endif
!         print *,' INDEX FOR MAXLAT SET ',ioproc
cgwv watch this!!
          print *,' allocating ', lonf,maxlats_comp,nodes
          allocate
     1    (buff_mult_piecesg(lonf,ngrids_gg,maxlats_comp,nodes))
          print *,' allocated', lonf,ngrids_gg,maxlats_comp,nodes
      endif


c
c   SENDLOOP of grids from comp processors to I/O task.  The
c   I/O task may or may not be a comp task also.  The
c   send logic on that task is different for these two cases
c
c  big send
c     if(me .gt. -1) return
!
!
      IF (ME .ne. ioproc) THEN    !   Sending the data
         msgtag=me
         illen=lats_node _a
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_pieceg,illen*lonf*ngrids_gg,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP .ne. MPI_COMM_NULL) then
!
c iotask is also a compute task.  send is replaced with direct
c  array copy
           buff_mult_piecesg(:,:,1:lats_node_a,ioproc+1)=
     1     buff_mult_pieceg(:,:,1:lats_node_a)
!                              END COMPUTE TASKS PORTION OF LOGIC
        endif
!
c  END COMPUTE TASKS PORTION OF LOGIC
c  receiving part of I/O task
!
!!
!!      for pes ioproc
        DO proc=1,nodes_comp
          if (proc.ne.ioproc+1) then
            msgtag=proc-1
            illen=ivarg_global_a(2,proc)
!           print *,' pux target ',ubound(buff_mult_piecesg)
            CALL mpi_recv(buff_mult_piecesg(1,1,1,proc),
     1        illen*lonf*ngrids_gg
     1        ,MPI_R_IO,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
          endif
        enddo
      ENDIF
      call mpi_barrier(mpi_comm_all,ierr)
!!
      return
      end
      SUBROUTINE atmgg_wrt(IOPROC,cfile,xhour,idate
     &,                  global_lats_a,lonsperlat,pdryini)
cc
      use gfsio_module
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfsio_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_io_header
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rk => con_rocp
      implicit none
!!
      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   p0=100000.0, p0i=1.0/p0
      real (kind=kind_io4), parameter  :: zero4=0.0
      integer IOPROC
      character*40 cfile, tracer
      real(kind=kind_grid) xhour, pdryini
      integer idate(4),k,il, ngridgg, nt
!
      INTEGER              global_lats_a(latg),   lonsperlat(latg)
!!
!     real(kind=kind_evod) gencode,ppid,realform
      real(kind=kind_io4) yhour, pdryini4
      real(kind=kind_io4), allocatable :: vcoord4(:,:)
      real(kind=kind_io4)              :: pup(lonf*latg)
     &,                                   pdn(lonf*latg)
     &,                                   plyr(lonf*latg)
     &,                                   pupk(lonf*latg)
     &,                                   pdnk(lonf*latg)
!     real(kind=kind_io8), allocatable :: buff(:)
!     real(kind=kind_io4) yhour, pdryini4, vcoord4(levp1,3)
      integer iret, ks,iorder_l,i
!     integer iret, ks,irealf,iorder, idusr
      character * 8, allocatable :: recname(:), reclevtyp(:)
      integer,       allocatable :: reclev(:)
      logical first
      save first, recname, reclevtyp, reclev, vcoord4
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build upper air fields in to buff_mult
!
      ngridg=1
      do ngridgg=1,ngrids_gg
!       print *,' inside atmgg_wrt calling unsp ngridgg=',ngridgg
        call unsplit2g(ioproc,buff_multg(1,ngridgg),global_lats_a)
      enddo
!     print *,' finished ngrid loop ngrids_gg=',ngrids_gg
!    Building upper air  field is done
!
      if (me.eq.ioproc) then
!
!     print *,' me=',me,' In the atmgg_wrt'
        if (first) then
          first = .false.
          allocate (recname(levs*(4+ntrac)+2),
     &              reclevtyp(levs*(4+ntrac)+2),
     &              reclev(levs*(4+ntrac)+2))
!    &              buff(lonf*lats_nodes_a))
          recname(1)   = 'hgt'
          recname(2)   = 'pres'
          reclevtyp(1) = 'sfc'
          reclevtyp(2) = 'sfc'
          reclev(1)    = 1
          reclev(2)    = 1
          do k=1,levs
            recname(k+2)        = 'dpres'
            recname(k+2+levs)   = 'ugrd'
            recname(k+2+levs*2) = 'vgrd'
            recname(k+2+levs*3) = 'tmp'
            recname(k+2+levs*4) = 'spfh'
          enddo
          do nt=2,ntrac
              write(tracer,'("tracer",i2)') nt
!     print *,' TRACER=',tracer,' ntrac=',ntrac,' ntoz=',ntoz,
!    &' ntcw=',ntcw,' nt=',nt
            if (nt == ntoz) then
              do k=1,levs
                recname(k+2+levs*(3+ntoz)) = 'o3mr'
              enddo
            elseif (nt == ntcw) then
              do k=1,levs
                recname(k+2+levs*(3+ntcw)) = 'clwmr'
              enddo
            else
              write(tracer,'("tracer",i2)') nt
!     print *,' TRACER=',tracer
              do k=1,levs
                recname(k+2+levs*(3+nt)) = trim(tracer)
              enddo
            endif
          enddo
          do nt=1,ntrac+4
            do k=1,levs
              reclevtyp(k+2+(nt-1)*levs) = 'layer'
              reclev(k+2+(nt-1)*levs)    = k
            enddo
          enddo
!
          idpp  = 0
          idusr = 0
          idrun = 0
!         gfile_out = gfile_in
          if (gen_coord_hybrid) then                                      ! hmhj
            idvc    = vertcoord_id
            idvm    = thermodyn_id*10 + sfcpress_id    ! 1: ln(ps) 2:ps   ! hmhj
            idsl    = 2    ! idsl=2 for middle of layer                   ! hmhj
            nvcoord = 3
            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1                                                  ! hmhj
              vcoord4(k,1) = ak5(k)*1000.                                 ! hmhj
              vcoord4(k,2) = bk5(k)                                       ! hmhj
              vcoord4(k,3) = ck5(k)*1000.                                 ! hmhj
            enddo                                                         ! hmhj
          else if (hybrid) then                                           ! hmhj
            idvc    = 2                        ! for hybrid vertical coord.
            nvcoord = 2
            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1
              vcoord4(k,1) = ak5(levp1+1-k)*1000.
              vcoord4(k,2) = bk5(levp1+1-k)
              print 190,k,vcoord4(k,1),vcoord4(k,2)
190           format('in gfsio k=',i2,'  ak5r4=',f13.6,'  bk5r4=',e13.5)
            enddo
          else
            idvc    = 1    ! for sigma vertical coord. (default)
            nvcoord = 1
            allocate (vcoord4(levp1,nvcoord))
            vcoord4(:,1) = si (:)
          endif
        endif
!
        pdryini4 = pdryini
        iorder_l = 2
        irealf   = 2
        yhour    = xhour
        idvt    = (ntoz-1) + 10 * (ntcw-1)
!       idusr    = usrid

      print *,' calling gfsio_open lonf=',lonf,' latg=',latg
     &,' idate=',idate,' yhour=',yhour
!
        call gfsio_open(gfile_out,trim(cfile),'write',iret,
     &    version=ivsupa,fhour=yhour,idate=idate,nrec=2+levs*(4+ntrac),
     &    latb=latg,lonb=lonf,levs=levs,jcap=jcap,itrun=itrun,
     &    iorder=iorder_l,irealf=irealf,igen=igen,latf=latg,lonf=lonf,
     &    latr=latg,lonr=lonf,ntrac=ntrac,icen2=icen2,iens=iens,
     &    idpp=idpp,idsl=idsl,idvc=idvc,idvm=idvm,idvt=idvt,idrun=idrun,
     &    idusr=idusr,pdryini=pdryini4,ncldt=ncldt,nvcoord=nvcoord,
     &    vcoord=vcoord4,recname=recname,reclevtyp=reclevtyp,
     &    reclev=reclev)
!
!     print *,' after calling gfsio_open iret=',iret
!     if (yhour .gt. 5.99) call mpi_quit(3333)
!
!     print *,' buff_multg=',buff_multg(lonb*latb/2,:)
!        buff(:) = buff_multg(:,1)
         call gfsio_writerecv(gfile_out,'hgt','sfc',1,buff_multg(:,1),
     &                                                iret)
!     print *,' hgt ',' iret=',iret,' hgt=',buff_multg(1000,1)
!        buff(:) = buff_multg(:,2)
         call gfsio_writerecv(gfile_out,'pres','sfc',1,buff_multg(:,2),
     &                                                iret)
!     print *,' pres ',' iret=',iret,' pres=',buff_multg(1000,2)
!
         ks = 2
         do k=1,levs
!          call gfsio_getrechead(gfile_out,ngrid,vname,vlevtyp,vlev,iret)
!          print *,jrec,vname,vlevtyp,vlev
!          buff(:) = buff_multg(:,k+ks)
           call gfsio_writerecv(gfile_out,'dpres','layer',k,
     &                          buff_multg(:,k+ks), iret)
!     print *,' dp k=',k,' iret=',iret,' dp=',buff_multg(1000,k+ks)
         enddo
!
!     write out the layer mean pressure
!
         pdn(:) = buff_multg(:,2)
         pdnk   = (pdn*p0i) ** rk
         do k=1,levs
           pup(:) = max(pdn(:)-buff_multg(:,2+k), zero4)
!     print *,' pup=',pup(1),' pdn=',pdn(1),' k=',k,' pdnk=',pdnk(1)
           if (idvc == 3 .and. mod(idvm,10) == 2) then
             plyr = 0.5 * (pup + pdn)
           else
             do i=1,lonf*latg
               pupk(i) = (pup(i)*p0i) ** rk
               plyr(i) = p0*((pdnk(i)*pdn(i)-pupk(i)*pup(i)) /
     &                   (rk1*(pdn(i)-pup(i)))) ** rkr
               pdn(i)  = pup(i)
               pdnk(i) = pupk(i)
             enddo
           endif
!     print *,' pupk=',pupk(1),' plyr=',plyr(1),' k=',k,' rkr=',rkr
           call gfsio_writerecv(gfile_out,'pres','layer',k,plyr, iret)
         enddo
!
         ks = ks + levs
         do k=1,levs
!          buff(:) = buff_multg(:,k+ks)
           call gfsio_writerecv(gfile_out,'ugrd','layer',k,
     &                          buff_multg(:,k+ks), iret)
!     print *,' u k=',k,' iret=',iret,' u=',buff_multg(1000,k+ks)
         enddo
         ks = ks + levs
         do k=1,levs
!          buff(:) = buff_multg(:,k+ks)
           call gfsio_writerecv(gfile_out,'vgrd','layer',k,
     &                          buff_multg(:,k+ks), iret)
!     print *,' v k=',k,' iret=',iret,' v=',buff_multg(1000,k+ks)
         enddo
         ks = ks + levs
         do k=1,levs
!          buff(:) = buff_multg(:,k+ks)
           call gfsio_writerecv(gfile_out,'tmp','layer',k,
     &                          buff_multg(:,k+ks), iret)
!     print *,' T k=',k,' iret=',iret,' T=',buff_multg(1000,k+ks)
         enddo
         ks = ks + levs
         do k=1,levs
!          buff(:) = buff_multg(:,k+ks)
           call gfsio_writerecv(gfile_out,'spfh','layer',k,
     &                          buff_multg(:,k+ks), iret)
!     print *,' Q k=',k,' iret=',iret,' Q=',buff_multg(1000,k+ks)
         enddo
         if (ntoz .gt. 0) then
           ks = 2 + levs*(ntoz+3)
           do k=1,levs
!            buff(:) = buff_multg(:,k+ks)
             call gfsio_writerecv(gfile_out,'o3mr','layer',k,
     &                          buff_multg(:,k+ks), iret)
           enddo
         endif
         if (ntcw .gt. 0) then
           ks = 2 + levs*(ntcw+3)
           do k=1,levs
!            buff(:) = buff_multg(:,k+ks)
             call gfsio_writerecv(gfile_out,'clwmr','layer',k,
     &                          buff_multg(:,k+ks), iret)
           enddo
         endif
         if (ntrac .gt. ntcw) then
           do nt=ntcw+1,ntrac
             ks = 2 + levs*(nt+3)
             write(tracer,'("tracer",i2)') nt
!            print *,' tracer=',tracer
             do k=1,levs
!              buff(:) = buff_multg(:,k+ks)
               call gfsio_writerecv(gfile_out,trim(tracer),'layer',k,
     &                          buff_multg(:,k+ks), iret)
             enddo
           enddo
         endif
      endif
!
!     print *,' return code before closing iret=',iret
      call gfsio_close(gfile_out,iret)
!     print *,' return code after closing iret=',iret
!     if (allocated(vcoord4)) deallocate(vcoord4)
!     print *,' after all atmgg writes iret=',iret
      return
      end
!
!

      subroutine uninterpreg(iord,kmsk,f,fi,global_lats_a,lonsperlat)
!!
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      implicit none
!!
      integer              global_lats_a(latg)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonf,lats_node_a)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_io8),intent(out):: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(in):: fi(lonf,lats_node_a)
!      real(kind=4) f4(lonf,lats_node_a)
      integer j,lons,lat
      integer i,ubound
!!
      do j=1,lats_node_a
         lat=global_lats_a(ipt_lats_node_a-1+j)
         lons=lonsperlat(lat)
         if(lons.ne.lonf) then
           call intlon(iord,1,1,lons,lonf,
     &                 kmsk(1,j),fi(1,j),f(1,j))
!          f4(:,j)=fi(:,j)
         else
            f(:,j)=fi(:,j)
!           f4(:,j)=fi(:,j)
         endif
      enddo
!     print *,' ngridg=',ngridg
      do j=1,lats_node_a
        do i=1,lonf
          buff_mult_pieceg(i,ngridg,j) = f (i,j)
        end do
      end do
      ngridg=ngridg+1
      end subroutine
       subroutine unsplit2g(ioproc,x,global_lats_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!!
      real(kind=kind_io4) x(lonf,latg)
      real(kind=kind_io4) tmp(lonf,latg+2)
      integer global_lats_a(latg),ipt_lats_node_al,nodesr
      integer lats_nodes_al
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen
       character*8 cna
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
!     write(cna,985)600+ngridg
!985   format('fort.',i3)
      X=0.
!     maxfld=50
      ifldu=ifldu+1
!!
      IF (me.ne.ioproc) THEN
            continue
      ELSE
!!
!!     for pes ioproc
c        if (.NOT.LIOPE) then
c            continue
c        else
c          nproct=nodes-1
c        endif
           nproct=nodes_comp
!     print *,' NGRIDG=',ngridg,' ifldu=',ifldu,' nproct=',nproct
        DO proc=1,nproct
c         if (proc.ne.ioproc+1) then
c         if (.NOT.LIOPE) then
c             continue
c         else
            ipt_lats_node_al=ivarg_global_a(1,proc)
            lats_nodes_al=ivarg_global_a(2,proc)
c         endif
         do j=1,lats_nodes_al
           lat=global_lats_a(ipt_lats_node_al-1+j)
           do i=1,lonf
c              x(i,lat)=tmp(i,j)
              x(i,lat)=buff_mult_piecesg(i,ngridg,j,proc)
           enddo
         enddo
c         endif   !(proc.ne.ioproc+1)
        enddo
!!
c        call baclose(563,i)
c         print *,cna,' UNSPLITFCLOSE  ',i
c        call baopenw(563,cna,i)
c         print *,cna,' UNSPLITF OPEN  ',i
      ENDIF
        ngridg=ngridg+1
!!
      return
      end
