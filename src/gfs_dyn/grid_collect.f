      SUBROUTINE grid_collect
     &        (zsg,psg,uug,vvg,teg,rqg,dpg,
     &        global_lats_a,lonsperlat)
!
      use gfs_dyn_resol_def
!jw
      use gfs_dyn_mod_state, only: buff_mult_pieceg,ngrid,ngridg
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
!jw
      use gfs_dyn_coordinate_def, only: idvc,idvm
      use gfs_dyn_physcons, rk => con_rocp
      implicit none

!jws
      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   p0=100000.0, p0i=1.0/p0
      real (kind=kind_io4), parameter  :: zero4=0.0
      integer i,j
!jwe
!
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
!jws
      real(kind=kind_io4) :: pup(lonf,lats_node_a)
     &,                      pdn(lonf,lats_node_a)
     &,                      pupk(lonf,lats_node_a)
     &,                      pdnk(lonf,lats_node_a)
!jwe
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
!jw         allocate(buff_mult_pieceg(lonf,ngrids_gg,lats_node_a))
         allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
      endif
!
      kmsk = 0
!
      buffi(:,:) = zsg(:,:)
!jw      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
!jw
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &  buff_mult_pieceg(1,1,1) )
      write(0,*)'in grid collect, buff_zsg=',
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,1)),
     & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,1))
!
      buffi(:,:) = psg(:,:)
!jw      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2) )
      write(0,*)'in grid collect, buff_psg=',
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2)),
     & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2))

!
      do k=1,levs
        buffi(:,:) = dpg(:,:,k)
!jw        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2+k) )
       write(0,*)'buff_dp layer=',k,
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+k)),
     &  minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+k))
      enddo
!!
!***  write out the layer mean pressure
         pdn(1:lonf,1:lats_node_a) =
     &      buff_mult_pieceg(1:lonf,1:lats_node_a,2)
         pdnk   = (pdn*p0i) ** rk
!         write(0,*)'before pup,size=', size(buff_mult_pieceg,1),
!     &      size(buff_mult_pieceg,2),'pdn=',size(pdn,1),size(pdn,2),
!     &  'pup=', minval(buff_mult_pieceg(1:lonf,1:lats_node_a,3)),
!     &      maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,3))
         do k=1,levs
           pup(1:lonf,1:lats_node_a) = max(pdn(1:lonf,1:lats_node_a)
     &       -buff_mult_pieceg(1:lonf,1:lats_node_a,2+k), zero4)
!      write(0,*)'pup= pdn=',pdn(1,1),' k=',k,' pdnk=',
!     &    pdnk(1,1),' dp=',buff_mult_pieceg(1,1,2+k),'idvc=',idvc,
!     &    'idvm=',idvm,'pup=',maxval(pup(1:lonf,1:lats_node_a)),
!     &    minval(pup(1:lonf,1:lats_node_a)),'pdn=',
!     &     minval(pdn(1:lonf,1:lats_node_a)),
!     &    maxval(pdn(1:lonf,1:lats_node_a))
           if (idvc == 3 .and. mod(idvm,10) == 2) then
             buff_mult_pieceg(1:lonf,1:lats_node_a,2+levs+k) = 
     &    0.5 * (pup(1:lonf,1:lats_node_a) + pdn(1:lonf,1:lats_node_a))
             pdn(1:lonf,1:lats_node_a)  = pup(1:lonf,1:lats_node_a)
           else
             do j=1,lats_node_a
             do i=1,lonf
               pupk(i,j) = (pup(i,j)*p0i) ** rk
               buff_mult_pieceg(i,j,2+levs+k) = p0*((pdnk(i,j)*pdn(i,j)-
     &            pupk(i,j)*pup(i,j)) /(rk1*(pdn(i,j)-pup(i,j)))) ** rkr
             pdn(i,j)  = pup(i,j)
             pdnk(i,j) = pupk(i,j)
             enddo
             enddo
           endif
       write(0,*)'buff_p layer=',k,
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+levs+k)),
     &  minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+levs+k))
!
         enddo
!
      do k=1,levs
        buffi(:,:) = uug(:,:,k)
!jw     CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2+2*levs+k) )
       write(0,*)'buff_u layer=',k,
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+2*levs+k)),
     &  minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+2*levs+k))
      enddo
!
      do k=1,levs
        buffi(:,:) = vvg(:,:,k)
!jw        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2+3*levs+k) )
      enddo
!
      do k=1,levs
        buffi(:,:) = teg(:,:,k)
!jw        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2+4*levs+k) )
       write(0,*)'buff_t layer=',k,'teg=',maxval(teg(:,:,k)),
     &  minval(teg(:,:,k)), 'buff=',
     &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+4*levs+k)),
     &  minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+4*levs+k))
      enddo
!
      if (levh .gt. 0) then
        do k=1,levh
          buffi(:,:) = rqg(:,:,k)
!jw          CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat)
          CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     & buff_mult_pieceg(1,1,2+5*levs+k) )
        enddo
      endif
!
      write(0,*)' finished gridc_collect for  ngridg=',ngridg
  999 continue
!jw      ngrid=1
!jw      ngridg=1
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
      integer illen,ubound,nd1
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
          if(me .ne. ioproc) then
            continue
          else
            maxlats_comp=ivarg_global_a(3,1)
          endif
      endif
cgwv watch this!!
          write(0,*)' allocating ', lonf,maxlats_comp,nodes,ngridg
      if(me .eq. ioproc) then
          allocate(buff_mult_piecesg(lonf*latg*ngrids_gg))
!jw     1    (buff_mult_piecesg(lonf,ngrids_gg,maxlats_comp,nodes))
          write(0,*)' allocated', lonf,ngrids_gg,maxlats_comp,nodes
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
        if(nodes_comp==1) then
        buff_mult_piecesg(1:lonf*lats_node_a*ngrids_gg)=
     &   reshape(buff_mult_pieceg(1:lonf,1:lats_node_a,1:ngrids_gg),
     &     (/lonf*lats_node_a*ngrids_gg/) )
!                              END COMPUTE TASKS PORTION OF LOGIC
        else
!
c  END COMPUTE TASKS PORTION OF LOGIC
c  receiving part of I/O task
!
!!      for pes ioproc
        nd1=0
        DO proc=1,nodes_comp
          illen=ivarg_global_a(2,proc)
          if (proc.ne.ioproc+1) then
            msgtag=proc-1
!           print *,' pux target ',ubound(buff_mult_piecesg)
            CALL mpi_recv(buff_mult_piecesg(nd1+1),
     1        illen*lonf*ngrids_gg
     1        ,MPI_R_IO,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
          else
            buff_mult_piecesg(nd1+1:nd1+illen*lonf*ngrids_gg)=
     &   reshape(buff_mult_pieceg(1:lonf,1:illen,1:ngrids_gg),
     &     (/lonf*illen*ngrids_gg/) )
          endif
          nd1=nd1+illen*lonf*ngrids_gg
        enddo
       endif 
       endif
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
      use gfs_dyn_tracer_const, only : cpi,ri
      use gfs_dynamics_output, only : DYN_INT_STATE_ISCALAR,
     & DYN_INT_STATE_RSCALAR,DYN_INT_STATE_LSCALAR, DYN_INT_STATE_1D_I,
     & DYN_INT_STATE_2D_I,DYN_INT_STATE_1D_R,DYN_INT_STATE_2D_R,
     & DYN_INT_STATE_3D_R_DIAB,DYN_INT_STATE_3D_R_ADIAB,
     & DYN_INT_STATE_4D_R
      use nemsio_module
!
      implicit none
!!
      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   p0=100000.0, p0i=1.0/p0
      real (kind=kind_io4), parameter  :: zero4=0.0
      integer IOPROC
      character*40 cfile, tracer
      real(kind=kind_grid) xhour, pdryini
      integer idate(4),k,il, ngridgg, nt,idate7(7)
!
      INTEGER              global_lats_a(latg),   lonsperlat(latg)
!!
!     real(kind=kind_evod) gencode,ppid,realform
      real(kind=kind_io4) yhour, pdryini4
      real(kind=kind_io4), allocatable :: vcoord4(:,:,:)
      real(kind=kind_io4)              :: pup(lonf*latg)
     &,                                   pdn(lonf*latg)
     &,                                   plyr(lonf*latg)
     &,                                   pupk(lonf*latg)
     &,                                   pdnk(lonf*latg)
!     real(kind=kind_io8), allocatable :: buff(:)
!     real(kind=kind_io4) yhour, pdryini4, vcoord4(levp1,3)
      integer iret, ks,iorder_l,i
!     integer iret, ks,irealf,iorder, idusr
      character * 16, allocatable :: recname(:), reclevtyp(:)
      integer,       allocatable :: reclev(:)
!jws
      integer j,ndim3,N2DR,kount,nrec
      integer ::nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr
      character(16),allocatable :: variname(:),varrname(:),
     &    varlname(:),aryiname(:),aryrname(:)
      integer,allocatable :: varival(:),aryilen(:),aryrlen(:),
     &    aryival(:,:)
      real(kind=kind_io4),allocatable    :: varrval(:),aryrval(:,:)
      logical,allocatable    :: varlval(:)
      REAL(KIND=KIND_io4) ,allocatable :: buff_multg(:,:)
      character * 16 :: recname1, reclevtyp1
      integer reclev1
      real(kind=kind_io4),allocatable :: tmp(:)
      type(nemsio_gfile) :: gfileout
!jwe
      logical first
      save first, recname, reclevtyp, reclev, vcoord4
      save nrec,nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr,
     &     variname,varrname,varlname,aryiname,aryrname,
     &     varival,aryilen,aryrlen,aryival,aryrval,varrval,varlval
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build upper air fields in to buff_mult
!
      allocate(buff_multg(lonf*latg,ngrids_gg))
      do ngridgg=1,ngrids_gg
!       print *,' inside atmgg_wrt calling unsp ngridgg=',ngridgg
        call unsplit2g(ioproc,ngridgg,ngrids_gg,buff_multg(1,ngridgg),
     &        global_lats_a)
      enddo
!     print *,' finished ngrid loop ngrids_gg=',ngrids_gg
!    Building upper air  field is done
!
      if (me.eq.ioproc) then
!
!     print *,' me=',me,' In the atmgg_wrt'
        if (first) then
          first = .false.
!***jw
!--- get names from module output 
          nrec=4+ntrac
          kount=size(DYN_INT_STATE_ISCALAR,2)
!for integer var::
          nmetavari=0
          do i=1,kount
           if(trim(DYN_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SIG') 
     &     nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          do i=1,kount
           if(trim(DYN_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SIG')then
           variname(i)=trim(DYN_INT_STATE_ISCALAR(1,i))
!jw ???
           if(i==1) varival(i)=latg
           if(i==2) varival(i)=lonf
           if(i==3) varival(i)=levs
           if(i==4) varival(i)=jcap
           if(i==5) varival(i)=ntoz
           if(i==6) varival(i)=ntcw
           if(i==7) varival(i)=ncldt
           if(i==8) varival(i)=vertcoord_id
           if(i==9) varival(i)=thermodyn_id
           if(i==10) varival(i)=sfcpress_id
           if(i==11) varival(i)=ienst
           if(i==12) varival(i)=iensi
           if(i==13) varival(i)=itrun
           if(i==14) varival(i)=icen2
           endif
          enddo
!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(DYN_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SIG')
     &     nmetavarr=nmetavarr+1
          enddo
          allocate(varrname(nmetavarr),varrval(nmetavarr))
          do i=1,kount
           if(trim(DYN_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SIG')then
            varrname(i)=trim(DYN_INT_STATE_RSCALAR(1,i))
            if(i==1) varrval(i)=pdryini
            if(i==2) varrval(i)=xhour
           endif
          enddo
!for logical var::
          nmetavarl=0
          do i=1,kount
           if(trim(DYN_INT_STATE_LSCALAR(2,i)).eq.'OGFS_SIG')
     &     nmetavarl=nmetavarl+1
          enddo
          allocate(varlname(nmetavarl),varlval(nmetavarl))
          do i=1,kount
           if(trim(DYN_INT_STATE_LSCALAR(2,i)).eq.'OGFS_SIG')then
            varlname(i)=trim(DYN_INT_STATE_LSCALAR(1,i))
            if(i==1) varlval(i)=hybrid
            if(i==2) varlval(i)=gen_coord_hybrid
            if(i==3) varlval(i)=adiabatic
           endif
          enddo
!for 1D integer array
          nmetaaryi=0
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_I(2,i)).eq.'OGFS_SIG')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_I(2,i)).eq.'OGFS_SIG') then
            aryiname(i)=trim(DYN_INT_STATE_1D_I(1,i))
            if(i==1) aryilen(i)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi))
          aryival(:,1)=idate(:)
!for 1D real array
          nmetaaryr=0
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_R(2,i)).eq.'OGFS_SIG')
     &     nmetaaryr=nmetaaryr+1
          enddo
          allocate(aryrname(nmetaaryr),aryrlen(nmetaaryr))
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_R(2,i)).eq.'OGFS_SIG') then
            aryrname(i)=trim(DYN_INT_STATE_1D_R(1,i))
            if(i==1) aryrlen(i)=size(ak5)
            if(i==2) aryrlen(i)=size(bk5)
            if(i==3) aryrlen(i)=size(ck5)
            if(i==4) aryrlen(i)=size(si)
            if(i==5) aryrlen(i)=size(cpi)
            if(i==6) aryrlen(i)=size(ri)
           endif
          enddo
          allocate(aryrval(maxval(aryrlen),nmetaaryr))
          aryrval(1:aryrlen(1),1)=ak5(:)
          aryrval(1:aryrlen(2),2)=bk5(:)
          aryrval(1:aryrlen(3),3)=ck5(:)
          aryrval(1:aryrlen(4),4)=si(:)
          aryrval(1:aryrlen(5),5)=cpi(:)
          aryrval(1:aryrlen(6),6)=ri(:)
!
!for record name, levtyp and lev          
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          N2DR=0
          do i=1,kount
           if(trim(DYN_INT_STATE_2D_R(2,i)).eq.'OGFS_SIG') then
            recname(i)=trim(DYN_INT_STATE_2D_R(1,i))
            reclevtyp(i)='sfc'
            reclev(i)=1
            N2DR=N2DR+1
           endif
          enddo
          do j=1,Kount
           if(trim(DYN_INT_STATE_3D_R_DIAB(2,j)).eq.'OGFS_SIG') then
            if(trim(DYN_INT_STATE_3D_R_DIAB(3,j)).eq.'levs') then
             do i=1,levs
               recname(N2DR+1)=trim(DYN_INT_STATE_3D_R_DIAB(1,j))
               reclevtyp(N2DR+1)='mid layer'
               reclev(N2DR+1)=i
               N2DR=N2DR+1
             enddo
            endif
           endif
          enddo
          if(ntrac>ntcw) then
           do j=1,Kount
            if(trim(DYN_INT_STATE_4D_R(2,j)).eq.'OGFS_SIG') then
             if(trim(DYN_INT_STATE_4D_R(3,j))=='levs')then
              NDIM3=levs
             elseif(trim(DYN_INT_STATE_4D_R(3,j))=='levsp1')then
              NDIM3=levs+1
             endif
             write(tracer,'("_",i2)') j

             do i=1,NDIM3
                recname(N2DR+1)=trim(DYN_INT_STATE_4D_R(1,i))//TRACER
                reclevtyp(N2DR+1)='mid layer'
                reclev(N2DR+1)=i
                N2DR=N2DR+1
             enddo
            Endif
           enddo
          endif
          write(0,*)'gfs sig file, total records=',nrec, 'N2DR=',N2DR
!
          idpp  = 0
          idusr = 0
          idrun = 0
          ALLOCATE(VCOORD4(levs+1,3,2))
          vcoord4=0.
          if (gen_coord_hybrid) then                                      ! hmhj
            idvc    = vertcoord_id
            idvm    = thermodyn_id*10 + sfcpress_id    ! 1: ln(ps) 2:ps   ! hmhj
            idsl    = 2    ! idsl=2 for middle of layer                   ! hmhj
            nvcoord = 3
!            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1                                                  ! hmhj
              vcoord4(k,1,1) = ak5(k)*1000.                                 ! hmhj
              vcoord4(k,2,1) = bk5(k)                                       ! hmhj
              vcoord4(k,3,1) = ck5(k)*1000.                                 ! hmhj
            enddo                                                         ! hmhj
          else if (hybrid) then                                           ! hmhj
            idvc    = 2                        ! for hybrid vertical coord.
            nvcoord = 2
!            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1
              vcoord4(k,1,1) = ak5(levp1+1-k)*1000.
              vcoord4(k,2,1) = bk5(levp1+1-k)
!              print 190,k,vcoord4(k,1),vcoord4(k,2)
190           format('in gfsio k=',i2,'  ak5r4=',f13.6,'  bk5r4=',e13.5)
            enddo
          else
            idvc    = 1    ! for sigma vertical coord. (default)
            nvcoord = 1
!jw            allocate (vcoord4(levp1,nvcoord))
            vcoord4(:,1,1) = si (:)
          endif
!end first
        endif
!
        pdryini4 = pdryini
        iorder_l = 2
        irealf   = 2
        yhour    = xhour
        idvt    = (ntoz-1) + 10 * (ntcw-1)
        idate7=0
        idate7(1:4)=idate(1:4)
!       idusr    = usrid
!
      print *,' calling gfsio_open lonf=',lonf,' latg=',latg
     &,' idate=',idate,' yhour=',yhour
!
        call nemsio_open(gfileout,trim(cfile),'write',iret=iret,
     &    modelname='gfs',gdatatype='grib',version=ivsupa,
     &    nfhour=int(yhour),idate=idate7,nrec=nrec,
     &    dimx=latg,dimy=lonf,dimz=levs,ncldt=ncldt,
     &    vcoord=vcoord4,extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,nmetavarl=nmetavarl,
     &    nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,varlname=varlname,varlval=varlval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    aryrname=aryrname,aryrlen=aryrlen,aryrval=aryrval,
     &    idsl=idsl,idvc=idvc,idvm=idvm,
     &    ntrac=ntrac,
     &    recname=recname,reclevtyp=reclevtyp,
     &    reclev=reclev)
!
      print *,' after calling gfsio_open iret=',iret
!     if (yhour .gt. 5.99) call mpi_quit(3333)
!
!     print *,' buff_multg=',buff_multg(lonb*latb/2,:)
      allocate(tmp(lonf*latg) )
       do k=1,nrec
         call nemsio_getrechead(gfileout,k,recname1,reclevtyp1,
     &     reclev1,iret)
         tmp(1:lonf*latg)=buff_multg(1:lonf*latg,k)
         call nemsio_writerec(gfileout,k,tmp,iret=ireT)
         write(0,*)'ig file,k=',k,'recname=',recname1,'reclevtyp=',
     &     reclevtyp1,'reclev=',reclev,'data=',maxval(buff_multg(:,k)),
     &     minval(buff_multg(:,k))
       enddo
       deallocate(tmp)
       deallocate(buff_multg)
!     print *,' return code before closing iret=',iret

      call  nemsio_close(gfileout,iret)
!
!endif ioproc
      endif
!     print *,' return code after closing iret=',iret
!     if (allocated(vcoord4)) deallocate(vcoord4)
!     print *,' after all atmgg writes iret=',iret
      return
      end
!
!

      subroutine uninterpreg(iord,kmsk,f,fi,global_lats_a,lonsperlat, 
     &    buff_mult)
!!
      use gfs_dyn_resol_def
!jw      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      implicit none
!!
      integer              global_lats_a(latg)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonf,lats_node_a)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_io8),intent(out):: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(in):: fi(lonf,lats_node_a)
      real(kind=kind_io4),intent(inout)::buff_mult(lonf,lats_node_a_max)
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
!jw          buff_mult_pieceg(i,ngridg,j) = f (i,j)
          buff_mult(i,j) = f (i,j)
        end do
      end do
!jw      ngridg=ngridg+1
      end subroutine
       subroutine unsplit2g(ioproc,ngridx,ngridt,x,global_lats_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_mod_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!!
      real(kind=kind_io4) x(lonf*latg)
      real(kind=kind_io4) tmp(lonf,latg+2)
      integer global_lats_a(latg),ipt_lats_node_al,nodesr
      integer lats_nodes_al
      integer maxfld,ioproc,nproct,ngridx,ngridt
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,nd1,nd2
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
       nproct=nodes_comp
!     print *,' NGRIDG=',ngridg,' ifldu=',ifldu,' nproct=',nproct
!for all the other pes
       nd1=0
       DO proc=1,nproct
         ipt_lats_node_al=ivarg_global_a(1,proc)
         lats_nodes_al=ivarg_global_a(2,proc)
         nd2=lats_nodes_al*lonf*(ngridx-1)
         do j=1,lats_nodes_al
           lat=global_lats_a(ipt_lats_node_al-1+j)
           do i=1,lonf
             x(i+(lat-1)*lonf)=buff_mult_piecesg(nd1+nd2+i+(j-1)*lonf)
           enddo
         enddo
         nd1=nd1+lats_nodes_al*lonf*ngridt
       enddo
!
      ENDIF
!!
      return
      end
