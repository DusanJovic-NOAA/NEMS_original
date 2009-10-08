      module gfs_dyn_mod_state
!
c new module to supply domain information
c to the GFS output routines called by
c wrtout.
!
! May 2009 Jun Wang, modified to use write grid component
!
      use gfs_dyn_machine
      use gfs_dyn_resol_def
      use gfsio_module
      use gfsio_def
      implicit none
!
      real(kind=kind_io4), allocatable,target :: buff_mult_pieceg(:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecesg(:)
!
      real(kind=kind_io4), allocatable :: buff_mult_piece(:,:,:),
     1                                    buff_mult_pieces(:,:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecef(:,:,:),
     1                                    buff_mult_piecesf(:,:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecea(:,:,:),
     1                                    buff_mult_piecesa(:,:,:,:)
      integer , allocatable :: ivar_global(:),ivar_global_a(:,:)
     &,                        ivarg_global(:),ivarg_global_a(:,:)
!
      integer ngrid ,ngrida,ngridg
      save ngrid,ngrida,buff_mult_piece,buff_mult_pieces,ivar_global
     &,    ngridg,buff_mult_pieceg,buff_mult_piecesg,ivarg_global
      end module gfs_dyn_mod_state

      subroutine wrtout_dynamics(phour,fhour,zhour,idate,
     &                  TRIE_LS,TRIO_LS,grid_gr,
     &                  sl,si,
     &    ls_node,ls_nodes,max_ls_nodes,
     &    lats_nodes_a,global_lats_a,lonsperlat,nblck,
     &    colat1,cfhour1,
     &    epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,
     &    pdryini)
!!
!! write out only grid values for gfsio
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_gg_def
!     use sigio_module
!     use sigio_r_module
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, cp => con_cp 
     &                    , rd => con_rd, fv => con_fvirt
     &                    , rkappa => con_rocp
      implicit none
cc
      CHARACTER(16) :: CFHOUR1         ! for the ESMF Export State Creation
      integer ixgr
      real(kind=kind_evod) phour,fhour,zhour
cc
      integer              idate(4),nblck,km,iostat,no3d,ks
      logical lfnhr
      real colat1, lat, lan
      real(kind=8) t1,t2,t3,t4,t5,ta,tb,tc,td,te,tf,rtc,tx,ty
      real timesum
cc
      real(kind=kind_evod) sl(levs), si(levp1)
cc
      integer              ls_node(ls_dim,3)
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)

      real(kind=kind_evod)   tfac(lonf,levs), sumq(lonf,levs)
      real(kind=kind_evod)   tki(lonf,levs+1)
      real(kind=kind_evod)   tkrt0, tx2(levs), tem
      real(kind=kind_evod), parameter :: one=1.0, cb2pa=1000.0
      real(kind=kind_evod), parameter :: qmin=1.e-10
      real(kind=kind_evod)  tx1
      integer               lons_lat,nn,kk,nnl
cc
      integer               ierr,i,j,k,l,lenrec,locl,n,node
      integer               nosig,nfill,jlonf
      character*16 cosfc
      data timesum/0./
cc
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
     &,                    TRIO_LS(LEN_TRIO_LS,2,lotls)
      REAL(KIND=KIND_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!!
      character CFHOUR*40,CFORM*40
      integer jdate(4),nzsig,ndigyr,ndig,kh,ioproc
!!
      REAL (KIND=KIND_grid) pdryini
      INTEGER              GLOBAL_lats_a(latg),   lonsperlat(latg)
!
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!!
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
!!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
      real(kind=kind_grid) ttg(lonf,lats_node_a,levs)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!!
      real(kind=kind_mpi),allocatable :: trieo_ls_nodes_buf(:,:,:,:,:)
      real(kind=kind_mpi),allocatable :: trieo_ls_node(:,:,:)
      save trieo_ls_nodes_buf,trieo_ls_node
      real(kind=8) tba,tbb,tbc,tbd
      integer iret
!
      t3=rtc()
!jw      call mpi_barrier(mpi_comm_all,ierr)
      call mpi_barrier(mc_comp,ierr)
      t4=rtc()
      tba=t4-t3
!jw      if(nodes_comp .lt. 1 .or. nodes_comp .gt. nodes) then
!jw        print *, '  NODES_COMP UNDEFINED, CANNOT DO I.O '
!jw        call mpi_finalize()
!jw         stop 333
!jw      endif
!
      ioproc=nodes_comp-1
      if(allocated ( trieo_ls_node)) then
        continue
      else
        allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     x                            2, 3*levs+1*levh+1 ) )
      endif
      t3=rtc()
!jw      call shapeset (ls_nodes,max_ls_nodes,pdryini)
!jw      call MPI_BARRIER(mpi_comm_all,ierr)
!jw      call MPI_BARRIER(mpi_comp,ierr)
      t4=rtc()
      tbb=t4-t3
       
      if ( allocated (trieo_ls_nodes_buf) )then
        continue
      else
        allocate( trieo_ls_nodes_buf ( len_trie_ls_max+len_trio_ls_max,
     x                               2, 3*levs+1*levh+1, nodes,1 ) )
      endif
      t1=rtc()

cc

!!
      JDATE=IDATE
      ndigyr=4
      IF(NDIGYR.EQ.2) THEN
        JDATE(4)=MOD(IDATE(4)-1,100)+1
      ENDIF

csela set lfnhr to false for writing one step output etc.
      lfnhr=.true.    ! no output
!      lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      lfnhr=3600*abs(fhour-nint(fhour)).le.1
      IF(LFNHR) THEN
        KH=NINT(FHOUR)
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS=NINT(FHOUR*3600)
        KH=KS/3600
        KM=(KS-KH*3600)/60
        KS=KS-KH*3600-KM*60
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      if( nfill(ens_nam) == 0 ) then
      CFHOUR = CFHOUR(1:nfill(CFHOUR))
      else
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      endif
      print *,' in wrtout_dynamics cfhour=',cfhour,' ens_nam=',ens_nam
cjfe
      nosig=61
!!
      t3=rtc()
      call MPI_BARRIER(mpi_comm_all,ierr)
      t4=rtc()
!
C*** BUILD STATE ON EACH NODE ********
c build state on each node.   COMP tasks only
c assemble upair state first then sfc state,
c then (only if liope)  flux state.
!
      t3=rtc()
      if(mc_comp .ne. MPI_COMM_NULL) then

          do lan=1,lats_node_a
            jlonf = (lan-1)*lonf
            zsg(1:lonf,lan) = grid_gr(jlonf+1:jlonf+lonf,g_gz)
          enddo
          do k=1,levh
            do lan=1,lats_node_a
              jlonf = (lan-1)*lonf
              rqg(1:lonf,lan,k)=
     &        grid_gr(jlonf+1:jlonf+lonf,g_rq-1+k)
            enddo
          enddo

          do lan=1,lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            tx1      = one / coslat_a(lat)
            jlonf = (lan-1)*lonf

            if (gen_coord_hybrid) then
              psg(1:lons_lat,lan) = grid_gr(jlonf+1:jlonf+lons_lat,g_q)
            else
              psg(1:lons_lat,lan) = 
     &        exp(grid_gr(jlonf+1:jlonf+lons_lat,g_q))
            endif

            if (gen_coord_hybrid) then        ! for general sigma-thera-p hybrid
              tki(:,1)       = 0.0
              tki(:,levs+1)  = 0.0
              do k=2,levs
                do i=1,lons_lat
                  tkrt0 = ( grid_gr(i+jlonf,g_tt-1+k-1)
     &                     +grid_gr(i+jlonf,g_tt-1+k) )
     &                      /(thref(k-1)+thref(k))
                  tki (i,k) = ck5(k)*tkrt0**rkappa
                enddo
              enddo
              do k=1,levs
                do i=1,lons_lat
                  dpg(i,lan,k) = ak5(k)-ak5(k+1)+(bk5(k)-bk5(k+1))
     &                     * psg(i,lan) + tki(i,k) - tki(i,k+1)
                enddo
              enddo
            elseif (hybrid) then              ! for sigma-p hybrid (ECWMF)
              do k=1,levs
                kk = levs - k + 1
                do i=1,lons_lat
                  dpg(i,lan,k) = ak5(kk+1)-ak5(kk)
     &                     + (bk5(kk+1)-bk5(kk)) * psg(i,lan)
                enddo
              enddo
            else                 ! For sigma coordinate
              do k=1,levs
                do i=1,lons_lat
                  dpg(i,lan,k) = (si(k) - si(k+1)) * psg(i,lan)
                enddo
              enddo
            endif
            if (thermodyn_id == 3) then
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = 0.0
                  sumq(i,k) = 0.0
                enddo
              enddo
              do nn=1,ntrac
                nnl = (nn-1)*levs
                if (cpi(nn) .ne. 0.0) then
                  do k=1,levs
                    do i=1,lons_lat
                      sumq(i,k) = sumq(i,k) + rqg(i,lan,nnl+k)
                      tfac(i,k) = tfac(i,k) + cpi(nn)*rqg(i,lan,nnl+k)
                    enddo
                  enddo
                endif
              enddo
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = (one-sumq(i,k))*cpi(0) + tfac(i,k)
                enddo
              enddo
            else
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = one + fv*max(rqg(i,lan,k),qmin)
                enddo
              enddo
            endif
            do k=1,levs
              do i=1,lons_lat
                uug(i,lan,k) = grid_gr(i+jlonf,g_uu-1+k) * tx1
                vvg(i,lan,k) = grid_gr(i+jlonf,g_vv-1+k) * tx1
                ttg(i,lan,k) = grid_gr(i+jlonf,g_tt-1+k) / tfac(i,k)
              enddo
            enddo
            do k=1,levs
              do i=1,lons_lat
                dpg(i,lan,k) = cb2pa*dpg(i,lan,k)
              enddo
            enddo
            do i=1,lons_lat
              psg(i,lan) = cb2pa*psg(i,lan)
            enddo

          enddo


      endif                 ! comp node
!
c  done with state build
c  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
c  NODES (currently done with a send to the I/O task_
c  send state to I/O task.  All tasks
!
        call grid_collect (zsg,psg,uug,vvg,ttg,rqg,dpg,
     &                         global_lats_a,lonsperlat)
!jw      if (.not.quilting ) then
!jw         call atmgg_move(ioproc)
!
c ioproc only
!jw         CFHOUR1 = CFHOUR          !for the ESMF Export State Creation
!jw         ta=rtc()
!jw         if(me .eq. ioproc) then
!jw           CFORM = 'SIG.F'//CFHOUR
!jw           print *,' calling atmgg_wrt fhour=',fhour
!jw     &,                     ' cform=',cform,' idate=',idate
!jw           call atmgg_wrt(IOPROC,CFORM,fhour,idate
!jw     &,                global_lats_a,lonsperlat,pdryini)
!jw           print *,' returning fromatmgg_wrt=',fhour
!jw         endif
!jw      endif
!
!jw      tc=rtc()
!jw      if(me .eq. 0) t2=rtc()
cgwv  t2=rtc()
!jw      t3=rtc()
!jw      if(MC_COMP   .ne. MPI_COMM_NULL) then
!jw        call mpi_barrier(mc_comp,info)
!jw      endif
!
!      write(0,*)'me=',me,'ioproc=',ioproc,'fhour=',fhour
      if(me .eq. ioproc)  call wrtlog_dynamics(phour,fhour,idate)
!jw      tb=rtc()
!jw      tf=tb-ta
!jw      t2=rtc()
!jw 1011 format(' WRTOUT_DYNAMICS TIME ',f10.4)
!jw      timesum=timesum+(t2-t1)
!jw 1012 format(
!jw     1 ' WRTOUT_DYNAMICS TIME ALL TASKS  ',f10.4,f10.4,
!jw     1 ' state, send, io  iobarr, (beginbarr),
!jw     1 spectbarr,open, openbarr )  ' ,
!jw     1  8f9.4)
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SUBROUTINE wrt_restart_dynamics(TRIE_LS,TRIO_LS,grid_gr,
!    &        SI,SL,fhour,idate,
!    &        igen,pdryini,
!    x        ls_node,ls_nodes,max_ls_nodes,
!    &        global_lats_a,lonsperlat,SNNP1EV,SNNP1OD,
!    &        ngptc, nblck, ens_nam)
!c
!     use gfs_dyn_resol_def
!     use gfs_dyn_layout1
!     use gfs_dyn_mpi_def
!     use sigio_module
!     use sigio_r_module
!     implicit none
!c
!     real(kind=kind_evod) fhour
!     character (len=*)  :: ens_nam
cc
!     integer              idate(4), ixgr
!     INTEGER              LS_NODE (LS_DIM*3)
!     integer              ls_nodes(ls_dim,nodes)
!     integer              max_ls_nodes(nodes)
!
!     REAL(KIND=KIND_EVOD) SNNP1EV(LEN_TRIE_LS)
!     REAL(KIND=KIND_EVOD) SNNP1OD(LEN_TRIO_LS)
!!
!     integer              ngptc, nblck
!!
!     real(kind=kind_evod) sl(levs)
!     real(kind=kind_evod) si(levp1)
!c
!     REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
!     REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,lotls)
!     REAL(KIND=KIND_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!!
!     integer igen
!!!
!     INTEGER              GLOBAL_lats_a(latg)
!     INTEGER              lonsperlat(latg)
!     integer IOPROC, IPRINT
!     integer needoro, iret, nfill
!
!c
!!
!     real runid,usrid
!     integer n3,n4,nflop
!     character*20 cflop,sigr51, sigr52
!     real pdryini
!     integer nn
!!
!     IPRINT = 0
!
!     sigr51 = 'SIGR1' // ens_nam(1:nfill(ens_nam))
!     sigr52 = 'SIGR2' // ens_nam(1:nfill(ens_nam))
!     print *,' sigr51=',sigr51,' sigr52=',sigr52
!    &,'ens_nam=',ens_nam(1:nfill(ens_nam))
!
!     n3=51
!     call sigio_rwopen(n3,sigr51,iret)
!     rewind(n3)
!     IF (icolor.eq.2) then
!        IOPROC=nodes-1
!     else
!        IOPROC=nodes
!     endif
!
!
!       CALL TWRITEEO(n3,ioproc,FHOUR,idate,
!    X                TRIE_LS(1,1,P_ZQ), TRIE_LS(1,1,P_QM ),
!    X                TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_DIM),
!    X                TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_RM),
!    X                TRIE_LS(1,1,P_GZ),
!    X                TRIO_LS(1,1,P_ZQ), TRIO_LS(1,1,P_QM ),
!    X                TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_DIM),
!    X                TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_RM),
!    X                TRIO_LS(1,1,P_GZ),
!    X                SL,SI,pdryini,
!    X                LS_NODES,MAX_LS_NODES,ixgr,
!    &                global_lats_a,lonsperlat,nblck)
!
!     IF (icolor.eq.2.and.me.eq.ioproc) print *,' closed ',n3
!
!     n4=52
!     call sigio_rwopen(n4,sigr52,iret)
!     rewind(n4)
!     IF (icolor.eq.2) then
!        IOPROC=nodes-1
!     else
!        IOPROC=nodes
!     endif
!     ixgr = 0
!       CALL TWRITEEO(n4,ioproc,FHOUR,idate,
!    X                TRIE_LS(1,1,P_ZQ), TRIE_LS(1,1,P_Q ),
!    X                TRIE_LS(1,1,P_TE), TRIE_LS(1,1,P_DI),
!    X                TRIE_LS(1,1,P_ZE), TRIE_LS(1,1,P_RQ),
!    X                TRIE_LS(1,1,P_GZ),
!    X                TRIO_LS(1,1,P_ZQ), TRIO_LS(1,1,P_Q ),
!    X                TRIO_LS(1,1,P_TE), TRIO_LS(1,1,P_DI),
!    X                TRIO_LS(1,1,P_ZE), TRIO_LS(1,1,P_RQ),
!    X                TRIO_LS(1,1,P_GZ),
!    X                SL,SI,pdryini,
!    X                LS_NODES,MAX_LS_NODES,ixgr,
!    &                global_lats_a,lonsperlat,nblck)
!jfe
!
!     nflop=53
!     IF (icolor.eq.2) then
!        IOPROC=nodes-1
!     else
!        IOPROC=nodes
!     endif
!
!     return
!     end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtlog_dynamics(phour,fhour,idate)
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      implicit none

      integer idate(4),ndigyr,nolog
      integer ks,kh,km,ndig,nfill
      character CFHOUR*40,CFORM*40
      logical lfnhr
      real phour,fhour
c
c     CREATE CFHOUR

csela set lfnhr to false for writing one step output etc.
      lfnhr=.true.    ! no output
ccmr  lfnhr=.false.   !    output
!      lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      lfnhr=3600*abs(fhour-nint(fhour)).le.1
      IF(LFNHR) THEN
        KH=NINT(FHOUR)
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS=NINT(FHOUR*3600)
        KH=KS/3600
        KM=(KS-KH*3600)/60
        KS=KS-KH*3600-KM*60
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      if( nfill(ens_nam) == 0 ) then
      CFHOUR = CFHOUR(1:nfill(CFHOUR))
      else
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      endif
!      print *,' in wrtlog_dynamics cfhour=',cfhour,' ens_nam=',ens_nam

      nolog=99
      OPEN(NOlog,FILE='LOG.F'//CFHOUR,FORM='FORMATTED')
      write(nolog,100)fhour,idate
100   format(' completed mrf fhour=',f10.3,2x,4(i4,2x))
      CLOSE(NOlog)

      RETURN
      END


      subroutine  shapeset (ls_nodes,max_ls_nodes,pdryini)
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      implicit none
!
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
cc
      integer              ierr,j,k,l,lenrec,locl,n,node
cc
      integer              indjoff
      integer              indev
      integer              indod
cc
      real(kind=kind_evod) gencode,order,ppid,realform
      real(kind=kind_evod) subcen,tracers,trun,vcid,vmid,vtid
cc
      real(kind=kind_evod) dummy(201-levp1-levs)
      real(kind=kind_evod) ensemble(2),dummy2(18)
cc
      real(kind=kind_io4)   tmps(4+nodes+jcap1*nodes)
      real(kind=kind_io4)   tmpr(3+nodes+jcap1*(nodes-1))
      REAL (KIND=KIND_grid) pdryini
cc
      INTEGER              GLOBAL_lats_a(latg)
      INTEGER                 lonsperlat(latg)
cc
      integer  il,ilen,i,msgtag,ls_diml,nodesl,ioproc, itmpr
                                                                                                        
c  Now define shape of the coefficients array
c  as a function of node. This will define how
c  to assemble the few wavenumbers on each node
c  into a full coefficient array.
c
       IOPROC=nodes
       IF (LIOPE) then
 199    format(' GWVX MAX_LS_NODES ',i20)
        if (me.eq.0.or. me .eq. ioproc) then
        tmps=0.
        tmps(1)=PDRYINI
        tmps(2:nodes_comp+1)=max_ls_nodes(1:nodes_comp)
        tmps(nodes_comp+2)=ls_dim
        tmps(nodes_comp+3)=len_trie_ls_max
        tmps(nodes_comp+4)=len_trio_ls_max
        il=nodes_comp+4
        do i=1,nodes_comp
        do j=1,ls_dim
           il=il+1
           tmps(il)=ls_nodes(j,i)
        enddo
        enddo
        ilen=4+nodes_comp+jcap1*nodes_comp
        msgtag=2345
        if(me .eq. 0) then
            CALL mpi_send(tmps,ilen,MPI_R_IO,ioproc,
     &                msgtag,MPI_COMM_ALL,info)
           endif
        endif
!
        if (me.eq.ioproc) then
         ilen=4+nodes_comp+jcap1*(nodes_comp)
         msgtag=2345
             CALL mpi_recv(tmpr,ilen,MPI_R_IO,0,
     &                msgtag,MPI_COMM_ALL,stat,info)

          itmpr=3+nodes+jcap1*(nodes-1)
          tmps(1:itmpr) = tmpr(1:itmpr)
          ls_nodes=0
          pdryini=tmps(1)
          max_ls_nodes(1:nodes_comp)=int(tmps(2:nodes_comp+1))
          ls_diml= int(tmps(nodes_comp+2))
          len_trie_ls_max=int(tmps(nodes_comp+3))
          len_trio_ls_max=int(tmps(nodes_comp+4))
           il=nodes_comp+3+1
                                                                                                        
          do i=1,nodes_comp
          do j=1,ls_diml
             il=il+1
             ls_nodes(j,i)=int(tmps(il))
          enddo
          enddo
        endif
      ENDIF

      return
      end
 

      INTEGER FUNCTION nfill(C)
      implicit none
      integer j
      CHARACTER*(*) C
      NFILL=LEN(C)
      DO J=1,NFILL
        IF(C(J:J).EQ.' ') THEN
          NFILL=J-1
          RETURN
        ENDIF
      ENDDO
      RETURN
      END
 
 
