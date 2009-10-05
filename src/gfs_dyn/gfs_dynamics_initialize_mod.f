
! !module: gfs_dynamics_initialize_mod 
!          --- initialize module of the
!              gridded component of the gfs dynamics system.
!              gfs dynamics related initilization 
!
! !description: gfs dynamics gridded component initialize module.
!
! !revision history:
!
!  november 2004  weiyu yang     initial code.
!  january 2006  s. moorthi      update to the new gfs version
!  august 2006   h. juang        add option to run generalized coordinates
!  december 2006 s. moorthi      gfsio included
!  january 2007 h. juang         change for dynamics only
!  May     2008 j. wang          change for gfs wrt grid component
!  Oct     2009 sarah lu         init xlon, xlat, lats_nodes_a_fix
!
!
! !interface:
!
      module gfs_dynamics_initialize_mod
!
!!uses:
!
      use gfs_dynamics_getcf_mod
      use gfs_dyn_machine, only : kind_io4
      use gfsio_module , only : gfsio_init
!jw
      use gfs_dyn_mod_state, only : buff_mult_pieceg
      use gfs_dyn_layout1, only : ipt_lats_node_a
      use gfs_dyn_resol_def, only : adiabatic

      implicit none

      contains

      subroutine gfs_dynamics_initialize(gis_dyn, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      integer,                                    intent(out)   :: rc

      integer 		:: ierr

      integer 		:: j, l, n, ilat, locl, ikey, nrank_all
!!
      integer 		nf0, nf1
      real 		fd2
      integer           indlsev,jbasev,indev
      integer           indlsod,jbasod,indod
!!
      include 'function2'
!!

! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------
      me     = gis_dyn%me
!
      nodes  = gis_dyn%nodes
      nlunit = gis_dyn%nam_gfs_dyn%nlunit

      call compns_dynamics(gis_dyn%deltim, gis_dyn%iret, gis_dyn%ntrac,	&
                           gis_dyn%nxpt,   gis_dyn%nypt, gis_dyn%jintmx,&
                           gis_dyn%jcap,   gis_dyn%levs, gis_dyn%levr, 	&
                           gis_dyn%lonf,   gis_dyn%latg,          	&
                           gis_dyn%ntoz,   gis_dyn%ntcw, gis_dyn%ncld, 	&
                           gis_dyn%spectral_loop,               	&
                     me,   gis_dyn%nam_gfs_dyn%nlunit, 			&
                           gis_dyn%nam_gfs_dyn%gfs_dyn_namelist)
      write(0,*)'after compns_dyn,me=',me,'ntoz=',gis_dyn%ntoz
!
      call get_tracer_const(gis_dyn%ntrac,me,gis_dyn%nam_gfs_dyn%nlunit)
      write(0,*)'after gettracer_con=',gis_dyn%ntrac
!
      ntrac   = gis_dyn%ntrac
      nxpt    = gis_dyn%nxpt
      nypt    = gis_dyn%nypt
      jintmx  = gis_dyn%jintmx
      jcap    = gis_dyn%jcap
      levs    = gis_dyn%levs
      levr    = gis_dyn%levr
      lonf    = gis_dyn%lonf
      latg    = gis_dyn%latg
      ntoz    = gis_dyn%ntoz
      ntcw    = gis_dyn%ntcw
      ncld    = gis_dyn%ncld
      if (gis_dyn%nam_gfs_dyn%total_member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",i2.2)') gis_dyn%nam_gfs_dyn%member_id
      endif
!
      levh   = ntrac*levs
      latgd  = latg+ 2*jintmx 
      jcap1  = jcap+1 
      jcap2  = jcap+2 
      latg2  = latg/2 
      levm1  = levs-1 
      levp1  = levs+1 
      lonfx  = lonf + 1 + 2*nxpt+1 
      lnt    = jcap2*jcap1/2 
      lnuv   = jcap2*jcap1 
      lnt2   = 2*lnt 
      lnt22  = 2*lnt+1 
      lnte   = (jcap2/2)*((jcap2/2)+1)-1 
      lnto   = (jcap2/2)*((jcap2/2)+1)-(jcap2/2) 
      lnted  = lnte 
      lntod  = lnto 

!jw      ngrids_gg       = 2+levs*(4+ntrac)
      ngrids_gg       = 2+levs*(5+ntrac)
      gis_dyn%lnt2    = lnt2

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latgd))
      allocate(lon_dims_ext(latgd))

      allocate(colrad_a(latg2))
      allocate(wgt_a(latg2))
      allocate(wgtcs_a(latg2))
      allocate(rcs2_a(latg2))
      allocate(sinlat_a(latg))
      allocate(coslat_a(latg))

      allocate(am(levs,levs))
      allocate(bm(levs,levs))
      allocate(cm(levs,levs))
      allocate(dm(levs,levs,jcap1))
      allocate(tor(levs))
      allocate(si(levp1))
      allocate(sik(levp1))
      allocate(sl(levs))
      allocate(slk(levs))
      allocate(del(levs))
      allocate(rdel2(levs))
      allocate(ci(levp1))
      allocate(cl(levs))
      allocate(tov(levs))
      allocate(sv(levs))

      allocate(ak5(levp1))
      allocate(bk5(levp1))
      allocate(ck5(levp1)) 
      allocate(thref(levp1))
      allocate(ck(levs))
      allocate(dbk(levs))
      allocate(bkl(levs))
      allocate(amhyb(levs,levs))
      allocate(bmhyb(levs,levs))
      allocate(svhyb(levs))
      allocate(tor_hyb(levs))
      allocate(d_hyb_m(levs,levs,jcap1))
      allocate(dm205_hyb(jcap1,levs,levs))

      allocate(spdmax(levs))

!jw      if (gfsio_out .or. gfsio_in) then
      if (gfsio_in) then
        call gfsio_init(ierr)
      endif
!jw      if (gfsio_out) then
!jw        allocate(buff_multg(lonf*latg,ngrids_gg))
!jw      endif

      allocate(lbasdz(4,2,levs),lbasiz(4,2,levs),detai(levp1), &
               detam(levs),etamid(levs),etaint(levp1),         &
               sinlamg(lonf,latg2),coslamg(lonf,latg2))
!

      allocate(tor_sig(levs), d_m(levs,levs,jcap1),            &
         dm205(jcap1,levs,levs))
         dm205=555555555.
         d_m  =444444444.
!

!hmhj allocate(z(lnt2))	! not used
!hmhj allocate(z_r(lnt2))	! not used
!
      allocate(gis_dyn%lonsperlat(latg))

      if( reduced_grid ) then
        print *,' run with reduced gaussian grid '
        call set_lonsgg(gis_dyn%lonsperlat)
      else
        print *,' run with full gaussian grid '
        do j=1,latg
          gis_dyn%lonsperlat(j) = lonf
        enddo
      endif
!
      p_gz   =         1     !      gze/o(lnte/od,2),
      p_zem  = p_gz   +1     !     zeme/o(lnte/od,2,levs),
      p_dim  = p_zem  +levs  !     dime/o(lnte/od,2,levs),
      p_tem  = p_dim  +levs  !     teme/o(lnte/od,2,levs),
      p_rm   = p_tem  +levs  !      rme/o(lnte/od,2,levh),
      p_qm   = p_rm   +levh  !      qme/o(lnte/od,2),

      p_ze   = p_qm   +1     !      zee/o(lnte/od,2,levs),
      p_di   = p_ze   +levs  !      die/o(lnte/od,2,levs),
      p_te   = p_di   +levs  !      tee/o(lnte/od,2,levs),
      p_rq   = p_te   +levs  !      rqe/o(lnte/od,2,levh),
      p_q    = p_rq   +levh  !       qe/o(lnte/od,2),
      p_dlam = p_q    +1     !  dpdlame/o(lnte/od,2),
      p_dphi = p_dlam +1     !  dpdphie/o(lnte/od,2),
      p_uln  = p_dphi +1     !     ulne/o(lnte/od,2,levs),
      p_vln  = p_uln  +levs  !     vlne/o(lnte/od,2,levs),
      p_zslam= p_vln  +levs  !    zslam/o(lnte/od,2),
      p_zsphi= p_zslam+1     !    zsphi/o(lnte/od,2),
                                                                                
      p_w    = p_zsphi+1     !       we/o(lnte/od,2,levs),
      p_x    = p_w    +levs  !       xe/o(lnte/od,2,levs),
      p_y    = p_x    +levs  !       ye/o(lnte/od,2,levs),
      p_rt   = p_y    +levs  !      rte/o(lnte/od,2,levh),
      p_zq   = p_rt   +levh  !      zqe/o(lnte/od,2),

      lotls  = p_zq 

      g_gz   = 1
      g_uum  = g_gz  + 1	!  for grid point 
      g_vvm  = g_uum + levs	!  for grid point 
      g_ttm  = g_vvm + levs	!  for grid point 
      g_rm   = g_ttm + levs	!  for grid point 
      g_qm   = g_rm  + levh	!  for grid point 

      g_uu   = g_qm  + 1	!  for grid point 
      g_vv   = g_uu  + levs	!  for grid point 
      g_tt   = g_vv  + levs	!  for grid point 
      g_rq   = g_tt  + levs	!  for grid point 
      g_q    = g_rq  + levh	!  for grid point 

      g_u    = g_q   + 1	!  for grid point 
      g_v    = g_u   + levs	!  for grid point 
      g_t    = g_v   + levs	!  for grid point 
      g_rt   = g_t   + levs	!  for grid point 
      g_zq   = g_rt  + levh	!  for grid point 

      g_p    = g_zq  + 1   	!  for grid point 
      g_dp   = g_p   + levs   	!  for grid point 
      g_dpdt = g_dp  + levs   	!  for grid point 

      lotgr  = g_dpdt+ levs - 1
!c
      lots = 5*levs+1*levh+5 
      lotd = 6*levs+2*levh+0 
      lota = 3*levs+1*levh+1 
      lotm = 3*levs+1*levh+1 
!
      kwq  = 0*levs+0*levh+1   !   qe/o_ls
      kwte = 0*levs+0*levh+2   !  tee/o_ls
      kwdz = 1*levs+0*levh+2   !  die/o_ls  zee/o_ls
      kwrq = 3*levs+0*levh+2   !  rqe/o_ls
!
      ksz     =1
      ksd     =ksz+levs
      kst     =ksd+levs
      ksr     =kst+levs
      ksq     =ksr+levh
      ksplam  =ksq+1
      kspphi  =ksplam+1
      ksu     =kspphi+1
      ksv     =ksu+levs
      kzslam  =ksv+levs
      kzsphi  =kzslam+1
!
      ksum    =1
      ksvm    =ksum+levs
      kstm    =ksvm+levs
      ksrm    =kstm+levs
      kspsm   =ksrm+levh
!
      kau     =1
      kav     =kau+levs
      kat     =kav+levs
      kar     =kat+levs
      kaps    =kar+levh
      kazs    =kaps+1
!
      kdtphi  =1
      kdrphi  =kdtphi+levs
      kdtlam  =kdrphi+levh
      kdrlam  =kdtlam+levs
      kdulam  =kdrlam+levh
      kdvlam  =kdulam+levs
      kduphi  =kdvlam+levs
      kdvphi  =kduphi+levs
!
      gis_dyn%p_gz    = p_gz             !      gze/o(lnte/od,2),
      gis_dyn%p_zem   = p_zem            !     zeme/o(lnte/od,2,levs),
      gis_dyn%p_dim   = p_dim            !     dime/o(lnte/od,2,levs),
      gis_dyn%p_tem   = p_tem            !     teme/o(lnte/od,2,levs),
      gis_dyn%p_rm    = p_rm             !      rme/o(lnte/od,2,levh),
      gis_dyn%p_qm    = p_qm             !      qme/o(lnte/od,2),
      gis_dyn%p_zslam = p_zslam          ! hmhj
      gis_dyn%p_zsphi = p_zsphi          ! hmhj
      gis_dyn%p_ze    = p_ze             !      zee/o(lnte/od,2,levs),
      gis_dyn%p_di    = p_di             !      die/o(lnte/od,2,levs),
      gis_dyn%p_te    = p_te             !      tee/o(lnte/od,2,levs),
      gis_dyn%p_rq    = p_rq             !      rqe/o(lnte/od,2,levh),
      gis_dyn%p_q     = p_q              !       qe/o(lnte/od,2),
      gis_dyn%p_dlam  = p_dlam           !  dpdlame/o(lnte/od,2),
      gis_dyn%p_dphi  = p_dphi           !  dpdphie/o(lnte/od,2),
      gis_dyn%p_uln   = p_uln            !     ulne/o(lnte/od,2,levs),
      gis_dyn%p_vln   = p_vln            !     vlne/o(lnte/od,2,levs),
      gis_dyn%p_w     = p_w              !       we/o(lnte/od,2,levs),
      gis_dyn%p_x     = p_x              !       xe/o(lnte/od,2,levs),
      gis_dyn%p_y     = p_y              !       ye/o(lnte/od,2,levs),
      gis_dyn%p_rt    = p_rt             !      rte/o(lnte/od,2,levh),
      gis_dyn%p_zq    = p_zq             !      zqe/o(lnte/od,2)
!
      gis_dyn%g_gz    = g_gz             !      gze/o(lnte/od,2),
      gis_dyn%g_uum   = g_uum            !     uume/o(lnte/od,2,levs),
      gis_dyn%g_vvm   = g_vvm            !     vvme/o(lnte/od,2,levs),
      gis_dyn%g_ttm   = g_ttm            !     teme/o(lnte/od,2,levs),
      gis_dyn%g_rm    = g_rm             !      rme/o(lnte/od,2,levh),
      gis_dyn%g_qm    = g_qm             !      qme/o(lnte/od,2),
      gis_dyn%g_uu    = g_uu             !      uue/o(lnte/od,2,levs),
      gis_dyn%g_vv    = g_vv             !      vve/o(lnte/od,2,levs),
      gis_dyn%g_tt    = g_tt             !      tee/o(lnte/od,2,levs),
      gis_dyn%g_rq    = g_rq             !      rqe/o(lnte/od,2,levh),
      gis_dyn%g_q     = g_q              !       qe/o(lnte/od,2),
      gis_dyn%g_u     = g_u              !       ue/o(lnte/od,2,levs),
      gis_dyn%g_v     = g_v              !       ve/o(lnte/od,2,levs),
      gis_dyn%g_t     = g_t              !       ye/o(lnte/od,2,levs),
      gis_dyn%g_rt    = g_rt             !      rte/o(lnte/od,2,levh),
      gis_dyn%g_zq    = g_zq             !      zqe/o(lnte/od,2)
      gis_dyn%g_p     = g_p              !        p/o(lnte/od,2)
      gis_dyn%g_dp    = g_dp             !       dp/o(lnte/od,2)
      gis_dyn%g_dpdt  = g_dpdt           !     dpdt/o(lnte/od,2)
!
      gis_dyn%lotls = lotls
      gis_dyn%lotgr = lotgr
      gis_dyn%lots = lots 
      gis_dyn%lotd = lotd
      gis_dyn%lota = lota
      gis_dyn%lotm = lotm
!
      allocate(gis_dyn%tee1(levs))

      gis_dyn%lslag=.false.  ! if false eulerian scheme =.true. for semilag

!     print *,' finish dimension in gfs_dynamics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
      if (me == 0) write(0,*) 'io option ,liope :',liope
!
!##jw      call mpi_comm_dup(mpi_comm_all, mpi_comm_all_dup, ierr)
!##jw      call mpi_barrier (mpi_comm_all_dup,               ierr)

!##jw      if (nodes == 1) liope=.false.
!##jw      if (liope) then
!##jw        call mpi_comm_rank(mpi_comm_all_dup,nrank_all,ierr)
!##jw        icolor=1
!##jw        ikey=1
!##jw        nodes_comp=nodes-1
!##jw        if (nrank_all.eq.nodes-1) then
!!  io server
!##jw          write(*,*) 'io server task'
!##jw         icolor=2
!##jw        gis_dyn%kcolor=mpi_undefined
!##jw          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
!##jw          call mpi_comm_split(mpi_comm_all_dup,gis_dyn%kcolor,ikey,mc_comp,ierr)
!##jw        else
!sela     write(*,*) 'compute server task '
!##jw          icolor=mpi_undefined
!##jw          gis_dyn%kcolor=1
!##jw          call mpi_comm_split(mpi_comm_all_dup,gis_dyn%kcolor,ikey,mc_comp,ierr)
!##jw          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
!##jw          call mpi_comm_size(mc_comp,nodes,ierr)
!##jw        endif
!##jw      else
!##jw        icolor=2
!##jw        mc_comp=mpi_comm_all_dup
!##jw        nodes_comp=nodes
!##jw      endif
!!jw
      nodes_comp=nodes
!c
      call f_hpminit(me,"evod")  !jjt hpm stuff
!c
      call f_hpmstart(25,"get_ls_gftlons")
!c
!!
      call synchro
      call init_countperf(latg)
!$$$      time0=timer()
!jfe  call countperf(0,15,0.)
!
      if (me.eq.0) then
      print 100, jcap,levs
100   format (' smf ',i3,i3,' created august 2000 ev od ri ')
      print*,'number of threads is ',num_parthds()
!jw        if (liope) then
!jw          print*,'number of mpi procs is ',nodes
!jw          print*,'number of mpi io procs is 1 (nodes)'
!jw        else
          print*,'number of mpi procs is ',nodes
!jw        endif
      endif
!c
      gis_dyn%cons0    =    0.0d0     !constant
      gis_dyn%cons0p5  =    0.5d0     !constant
      gis_dyn%cons1200 = 1200.d0      !constant
      gis_dyn%cons3600 = 3600.d0      !constant
!c
!jw      if (liope) then
!jw         if (icolor.eq.2) then
!jw           ls_dim = jcap1
!jw         else
!jw           ls_dim = (jcap1-1)/nodes+1
!jw         endif
!jw      else
         ls_dim = (jcap1-1)/nodes+1
!jw      endif
!!
!c
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c
!c
      allocate (      gis_dyn%ls_node (ls_dim*3) )
      allocate (      gis_dyn%ls_nodes(ls_dim,nodes) )
      allocate (  gis_dyn%max_ls_nodes(nodes) )
!c
      allocate (  gis_dyn%lats_nodes_a_fix(nodes))     ! added for mGrid 
!c
      allocate (  gis_dyn%lats_nodes_a(nodes) )
      allocate ( gis_dyn%global_lats_a(latg) )
!c
      allocate (   gis_dyn%lats_nodes_ext(nodes) )
      allocate ( gis_dyn%global_lats_ext(latg+2*jintmx+2*nypt*(nodes-1)) )
!c
!c
      gis_dyn%iprint = 0
      call get_ls_node( me, gis_dyn%ls_node, ls_max_node, gis_dyn%iprint )
!c
!c
      len_trie_ls=0
      len_trio_ls=0
      do locl=1,ls_max_node
           gis_dyn%ls_node(locl+  ls_dim)=len_trie_ls
          gis_dyn%ls_node(locl+2*ls_dim)=len_trio_ls
         l=gis_dyn%ls_node(locl)
         len_trie_ls=len_trie_ls+(jcap+3-l)/2
         len_trio_ls=len_trio_ls+(jcap+2-l)/2
      enddo
!c
!c
      allocate (       gis_dyn%epse  (len_trie_ls) )
      allocate (       gis_dyn%epso  (len_trio_ls) )
      allocate (       gis_dyn%epsedn(len_trie_ls) )
      allocate (       gis_dyn%epsodn(len_trio_ls) )
!c
      allocate (      gis_dyn%snnp1ev(len_trie_ls) )
      allocate (      gis_dyn%snnp1od(len_trio_ls) )
!c
      allocate (       gis_dyn%ndexev(len_trie_ls) )
      allocate (       gis_dyn%ndexod(len_trio_ls) )
!c
      allocate (      gis_dyn%plnev_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%plnod_a(len_trio_ls,latg2) )
      allocate (      gis_dyn%pddev_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%pddod_a(len_trio_ls,latg2) )
      allocate (      gis_dyn%plnew_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%plnow_a(len_trio_ls,latg2) )
!c
      gis_dyn%maxstp=36

 
      if(me.eq.0) 							&
        print*,'from compns_dynamics : iret=',gis_dyn%iret		&
       ,' nsout=',nsout,' nsres=',nsres
      if(gis_dyn%iret.ne.0) then
        if(me.eq.0) print *,' incompatible namelist - aborted in main'
        call mpi_quit(13)
      endif
!!
      gis_dyn%lats_nodes_ext = 0
      call getcon_dynamics(gis_dyn%n3,gis_dyn%n4,			&
           gis_dyn%ls_node,gis_dyn%ls_nodes,gis_dyn%max_ls_nodes,       &
           gis_dyn%lats_nodes_a,gis_dyn%global_lats_a,                  &
           gis_dyn%lonsperlat,gis_dyn%lats_node_a_max,                  &
           gis_dyn%lats_nodes_ext,gis_dyn%global_lats_ext,              &
           gis_dyn%epse,gis_dyn%epso,gis_dyn%epsedn,gis_dyn%epsodn,     &
           gis_dyn%snnp1ev,gis_dyn%snnp1od,				&
           gis_dyn%ndexev,gis_dyn%ndexod,  				&
           gis_dyn%plnev_a,gis_dyn%plnod_a,				&
           gis_dyn%pddev_a,gis_dyn%pddod_a, 				&
           gis_dyn%plnew_a,gis_dyn%plnow_a,gis_dyn%colat1)
!jw
      gis_dyn%lats_node_a=gis_dyn%lats_nodes_a(me+1)
      gis_dyn%ipt_lats_node_a=ipt_lats_node_a
      write(0,*)'after getcon_dynamics,lats_node_a=',gis_dyn%lats_node_a &
       ,'ipt_lats_node_a=',gis_dyn%ipt_lats_node_a,'ngptc=',ngptc
!
!!
      gis_dyn%nblck=lonf/ngptc+1

!
! initialize coord def (xlon,xlat) and lats_nodes_a_fix          
!
      gis_dyn%lats_nodes_a_fix(:) = gis_dyn%lats_node_a_max     
      allocate ( gis_dyn%XLON(lonf,gis_dyn%lats_node_a) )     
      allocate ( gis_dyn%XLAT(lonf,gis_dyn%lats_node_a) )     

      call gfs_dyn_lonlat_para(gis_dyn%global_lats_a,        &   
              gis_dyn%xlon, gis_dyn%xlat, gis_dyn%lonsperlat)  

      call countperf(0,18,0.)
!!
!jw      if (.not.liope.or.icolor.ne.2) then
!!
      call countperf(0,15,0.)
      allocate (      gis_dyn%trie_ls(len_trie_ls,2,lotls) )
      allocate (      gis_dyn%trio_ls(len_trio_ls,2,lotls) )
      allocate (      gis_dyn%grid_gr(lonf*lats_node_a_max,lotgr) )
      allocate (      gis_dyn%pwat   (lonf,lats_node_a) )
      allocate (      gis_dyn%ptot   (lonf,lats_node_a) )
!c
      allocate (     gis_dyn%syn_ls_a(4*ls_dim,gis_dyn%lots,latg2) )
      allocate (     gis_dyn%dyn_ls_a(4*ls_dim,gis_dyn%lotd,latg2) )
!c
      allocate (   gis_dyn%syn_gr_a_1(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%syn_gr_a_2(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%sym_gr_a_2(lonfx*gis_dyn%lotm,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_1(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_2(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_1(lonfx*gis_dyn%lota,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_2(lonfx*gis_dyn%lota,lats_dim_ext) )
      write(0,*)'after allocate lonf=',lonf,'lats_node_a_max=',lats_node_a_max, &
       'ngrids_gg=',ngrids_gg
!!
!jw      endif !(.not.liope.or.icolor.ne.2)
!
!##jws allocate output vars
      allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
      buff_mult_pieceg=0.
      adiabatic=gis_dyn%adiabatic
!##jwe

!!
      allocate (      gis_dyn%fhour_idate(1,5) )
      write(0,*)'after allocate fhour_idate'
!
     print*, ' lats_dim_a=', lats_dim_a, ' lats_node_a=', lats_node_a
     print*, ' lats_dim_ext=', lats_dim_ext,              &
              ' lats_node_ext=', lats_node_ext
!c
      if (gis_dyn%lslag) then
        ilat=lats_node_ext
      else
        ilat=lats_node_a
      endif
      call countperf(1,15,0.)
!!
!c......................................................................
!c
      call countperf(0,15,0.)
      write(0,*)'before f_hpmstop'
      call f_hpmstop(25)
!c
      write(0,*) 'number of latitudes ext. :',lats_node_ext,              &
                  lats_dim_ext,lats_node_a
!!
      call countperf(1,15,0.)
!!
!      print *,' sig_ini=',gis_dyn%nam_gfs_dyn%sig_ini,			  &
!              ' sig_ini2=',gis_dyn%nam_gfs_dyn%sig_ini2 
      call countperf(0,18,0.)
      gis_dyn%pdryini = 0.0

      write(0,*)'before input_fields,ntoz=',ntoz
      call input_fields(gis_dyn%n1, gis_dyn%n2, gis_dyn%pdryini,           &
        gis_dyn%trie_ls, gis_dyn%trio_ls,  gis_dyn%grid_gr ,               &
        gis_dyn%ls_node, gis_dyn%ls_nodes, gis_dyn%max_ls_nodes,           &
        gis_dyn%snnp1ev, gis_dyn%snnp1od,                                  &
        gis_dyn%global_lats_a,  gis_dyn%nblck,   gis_dyn%lonsperlat,       &
        gis_dyn%epse, gis_dyn%epso, gis_dyn%plnev_a, gis_dyn%plnod_a,      &
        gis_dyn%plnew_a, gis_dyn%plnow_a, gis_dyn%lats_nodes_a,	           &
        gis_dyn%nam_gfs_dyn%sig_ini, gis_dyn%nam_gfs_dyn%sig_ini2,         &
        gis_dyn%pwat, gis_dyn%ptot)
      write(0,*)'after input_fields,ntoz=',ntoz
!!
        call countperf(1,18,0.)
!!
 
      tov = 0.0
      if (.not. (hybrid.or.gen_coord_hybrid) ) then                   ! hmhj
       call setsig(si,ci,del,sl,cl,rdel2,tov,me)
       am=-8888888.
       bm=-7777777.
       call amhmtm(del,sv,am)
       call bmdi_sig(ci,bm)
      endif
      call deldifs(gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,		&
                   gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,      	&
                   gis_dyn%epso,gis_dyn%epso,gis_dyn%epso,		&
                   gis_dyn%epso,gis_dyn%epso,gis_dyn%epso,      	& 
                   gis_dyn%cons0,sl,gis_dyn%ls_node,gis_dyn%epse,	&
                   0,hybrid,gen_coord_hybrid,nislfv)  
      write(0,*)'after deldifs,ntoz=',ntoz
 
!c
      call f_hpmstart(26,"step1")
!c
!!
      call countperf(1,18,0.)
!!
      gis_dyn%zhour        = fhour
      gis_dyn% start_step  = .true.
      gis_dyn% reset_step  = .false.
!
!
      end subroutine gfs_dynamics_initialize
!
! =========================================================================
!
      subroutine set_lonsgg(lonsperlat)
      use gfs_dyn_resol_def
      use gfs_dyn_reduce_lons_grid_module, only : gfs_dyn_reduce_grid   ! hmhj
      integer numreduce                                                 ! hmhj
      integer lonsperlat(latg)

      integer lonsperlat_62(94)
      integer lonsperlat_126(190)
      integer lonsperlat_170(256)
      integer lonsperlat_190(288)
      integer lonsperlat_254(384)
      integer lonsperlat_382(576)
      integer lonsperlat_510(766)

      data lonsperlat_62/                                                &
        30,  30,  30,  40,  48,  56,  60,  72,  72,  80,  90,  90,       &
        96, 110, 110, 120, 120, 128, 144, 144, 144, 144, 154, 160,       &
       160, 168, 168, 180, 180, 180, 180, 180, 180, 192, 192, 192,       &
       192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 47*0 /

      data lonsperlat_126      /                                         &
          30,   30,   36,   48,   56,   60,   72,   72,   80,   90,      &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,      &
         160,  180,  180,  180,  192,  192,  210,  210,  220,  220,      &
         240,  240,  240,  240,  240,  252,  256,  280,  280,  280,      &
         280,  288,  288,  288,  288,  308,  308,  308,  320,  320,      &
         320,  320,  330,  330,  360,  360,  360,  360,  360,  360,      &
         360,  360,  360,  360,  360,  360,  384,  384,  384,  384,      &
         384,  384,  384,  384,  384,  384,  384,  384,  384,  384,      &
         384,  384,  384,  384,  384,  384,  384,  384,  384,  384,      &
         384,  384,  384,  384,  384, 95*0 /
 
      data lonsperlat_170 /                                              &
         48,  48,  48,  48,  48,  56,  60,  72,  72,  80,  90,  96,      &
        110, 110, 120, 120, 128, 144, 144, 144, 154, 160, 168, 180,      &
        180, 180, 192, 210, 210, 220, 220, 240, 240, 240, 240, 240,      &
        252, 256, 280, 280, 280, 288, 288, 288, 308, 308, 320, 320,      &
        320, 320, 330, 360, 360, 360, 360, 360, 360, 360, 384, 384,      &
        384, 384, 384, 384, 420, 420, 420, 440, 440, 440, 440, 440,      &
        440, 440, 440, 440, 462, 462, 462, 462, 462, 480, 480, 480,      &
        480, 480, 480, 480, 480, 480, 480, 480, 504, 504, 504, 504,      &
        504, 504, 504, 504, 504, 512, 512, 512, 512, 512, 512, 512,      &
        512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,      &
        512, 512, 512, 512, 512, 512, 512, 512, 128*0 /
 
      data lonsperlat_190 /                                              &
        64,  64,  64,  64,  64,  64,  64,  70,  80,  84,                 &
        88, 110, 110, 110, 120, 126, 132, 140, 144, 154,                 &
       160, 168, 176, 176, 192, 192, 198, 210, 210, 220,                 &
       220, 240, 240, 240, 252, 252, 256, 264, 280, 280,                 &
       280, 288, 308, 308, 308, 320, 320, 320, 330, 336,                 &
       352, 352, 352, 352, 360, 384, 384, 384, 384, 384,                 &
       396, 396, 420, 420, 420, 420, 420, 440, 440, 440,                 &
       440, 440, 448, 448, 462, 462, 462, 480, 480, 480,                 &
       480, 480, 504, 504, 504, 504, 504, 504, 504, 512,                 &
       512, 528, 528, 528, 528, 528, 528, 560, 560, 560,                 &
       560, 560, 560, 560, 560, 560, 560, 560, 560, 560,                 &
       560, 576, 576, 576, 576, 576, 576, 576, 576, 576,                 &
       576, 576, 576, 576, 576, 576, 576, 576, 576, 576,                 &
       576, 576, 576, 576, 576, 576, 576, 576, 576, 576,                 &
       576, 576, 576, 576, 144*   0/
!
      data lonsperlat_254      /                                         &
          64,   64,   64,   64,   64,   64,   72,   72,   80,   90,      &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,      &
         168,  180,  180,  180,  192,  192,  210,  220,  220,  240,      &
         240,  240,  240,  252,  256,  280,  280,  280,  288,  288,      &
         288,  308,  308,  320,  320,  320,  330,  360,  360,  360,      &
         360,  360,  360,  384,  384,  384,  384,  420,  420,  420,      &
         440,  440,  440,  440,  440,  440,  462,  462,  462,  480,      &
         480,  480,  480,  480,  480,  504,  504,  504,  504,  512,      &
         512,  560,  560,  560,  560,  560,  560,  576,  576,  576,      &
         576,  576,  576,  576,  576,  616,  616,  616,  616,  616,      &
         616,  640,  640,  640,  640,  640,  640,  640,  640,  640,      &
         640,  660,  660,  660,  720,  720,  720,  720,  720,  720,      &
         720,  720,  720,  720,  720,  720,  720,  720,  720,  720,      &
         720,  720,  720,  720,  720,  720,  720,  720,  768,  768,      &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,      &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,      &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,      &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,      &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,      &
         768,  768,  192*0/
 
      data lonsperlat_382      /                                         &
         64,  64,  64,  64,  64,  64,  64,  70,  80,  84,                &
         88,  96, 110, 110, 120, 126, 132, 140, 144, 154,                &
        160, 168, 176, 180, 192, 192, 198, 210, 220, 220,                &
        224, 240, 240, 252, 252, 256, 264, 280, 280, 280,                &
        288, 308, 308, 308, 320, 320, 330, 336, 352, 352,                &
        352, 360, 384, 384, 384, 384, 396, 396, 420, 420,                &
        420, 420, 440, 440, 440, 448, 448, 462, 462, 480,                &
        480, 480, 504, 504, 504, 504, 512, 528, 528, 528,                &
        560, 560, 560, 560, 560, 560, 576, 576, 616, 616,                &
        616, 616, 616, 616, 616, 616, 630, 630, 640, 640,                &
        660, 660, 660, 660, 672, 672, 704, 704, 704, 704,                &
        704, 704, 720, 720, 720, 768, 768, 768, 768, 768,                &
        768, 768, 768, 768, 768, 792, 792, 792, 792, 792,                &
        840, 840, 840, 840, 840, 840, 840, 840, 840, 840,                &
        880, 880, 880, 880, 880, 880, 880, 880, 880, 880,                &
        896, 896, 896, 896, 924, 924, 924, 924, 924, 924,                &
        960, 960, 960, 960, 960, 960, 960, 960, 960, 960,                &
        990, 990, 990, 990, 990, 990, 990, 990,1008,1008,                &
       1008,1008,1008,1008,1024,1024,1024,1024,1024,1024,                &
       1056,1056,1056,1056,1056,1056,1056,1056,1056,1056,                &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,                &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,                &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,                &
       1120,1152,1152,1152,1152,1152,1152,1152,1152,1152,                &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,                &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,                &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,                &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,                &
       1152,1152,1152,1152,1152,1152,1152,1152, 288*   0/
 
      data lonsperlat_510      /                                         &
          64,   64,   64,   64,   64,   64,   72,   72,   80,   90,      &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,      &
         168,  180,  180,  180,  192,  210,  210,  220,  220,  240,      &
         240,  240,  240,  252,  256,  280,  280,  288,  288,  288,      &
         308,  308,  320,  320,  320,  330,  360,  360,  360,  360,      &
         360,  384,  384,  384,  384,  420,  420,  440,  440,  440,      &
         440,  440,  440,  462,  462,  462,  480,  480,  480,  480,      &
         504,  504,  504,  504,  512,  512,  560,  560,  560,  560,      &
         576,  576,  576,  576,  576,  576,  616,  616,  616,  616,      &
         640,  640,  640,  640,  640,  640,  640,  660,  720,  720,      &
         720,  720,  720,  720,  720,  720,  720,  720,  720,  720,      &
         720,  768,  768,  768,  768,  768,  768,  768,  768,  840,      &
         840,  840,  840,  840,  840,  840,  840,  880,  880,  880,      &
         880,  880,  880,  880,  880,  880,  880,  924,  924,  924,      &
         924,  924,  924,  924,  960,  960,  960,  960,  960,  960,      &
         960,  960,  960,  960,  960,  990,  990,  990, 1008, 1008,      &
        1008, 1008, 1008, 1024, 1024, 1024, 1024, 1024, 1120, 1120,      &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,      &
        1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,      &
        1152, 1152, 1152, 1152, 1152, 1152, 1232, 1232, 1232, 1232,      &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260, 1260,      &
        1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260,      &
        1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1320,      &
        1320, 1320, 1320, 1386, 1386, 1386, 1386, 1386, 1386, 1386,      &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386,      &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,      &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,      &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,      &
        1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,      &
        1536, 1536, 1536,  383*0/
 
      integer i
      if (jcap .eq. 62) then
         lonsperlat=lonsperlat_62
      endif
      if (jcap .eq. 126) then
         lonsperlat=lonsperlat_126
      endif
      if (jcap .eq. 170) then
         lonsperlat=lonsperlat_170
      endif
      if (jcap .eq. 190) then
         lonsperlat=lonsperlat_190
      endif
      if (jcap .eq. 254) then
         lonsperlat=lonsperlat_254
      endif
      if (jcap .eq. 382) then
         lonsperlat=lonsperlat_382
      endif
      if (jcap .eq. 510) then
         lonsperlat=lonsperlat_510
      endif

      if (jcap .ne. 62 .and. jcap .ne. 126 .and. jcap .ne. 170 .and.     &
          jcap .ne. 190 .and.                                            &
          jcap .ne. 254 .and. jcap .ne. 382 .and. jcap .ne. 510) then
! compute reduced grid using juang 2003
         print*,' non standard resolution  - lonsperlat',     &
                ' computed locally'
         numreduce=4                                                    ! hmhj
         call gfs_dyn_reduce_grid (numreduce,jcap,latg,lonsperlat)      ! hmhj
         print*,' reduced grid is computed - lonsperlat '    		! hmhj
      endif

!     print*,' jcap = ',jcap
!     print*,'min,max of lonsperlat = ',minval(lonsperlat),              &
!             maxval(lonsperlat)

      end subroutine set_lonsgg

      end module gfs_dynamics_initialize_mod
