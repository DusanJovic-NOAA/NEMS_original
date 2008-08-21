
! !module: gfs_physics_initialize_mod 
!          --- initialization module of the gridded component of gfs physics.
!
! !description: gfs physics gridded component initialize module.
!
! !revision history:
!
!  november 2004  weiyu yang     initial code.
!  january 2006  s. moorthi      update to the new gfs version
!  august 2006   h. juang        add option to run generalized coordinates
!  january 2007 h. juang        change for dynamics only
!  july    2007 s. moorthi      change for physics only
!  november 2007 h. juang        continue for physics
!
!
! !interface:
!
      module gfs_physics_initialize_mod

!
!!uses:
!
      use gfs_physics_getcf_mod
      use machine, only : kind_io4

      implicit none

      contains

      subroutine gfs_physics_initialize(gis_phy, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------

      type(gfs_physics_internal_state), pointer, intent(inout) :: gis_phy
      integer,                                   intent(out)   :: rc

      integer 		:: ierr

      integer 		:: l, n, ilat, locl, ikey, nrank_all
!!
      real (kind=kind_io4) blatc4
      real (kind=kind_io4), allocatable :: pl_lat4(:)
      real (kind=kind_io4), allocatable :: pl_pres4(:)
      real (kind=kind_io4), allocatable :: pl_time4(:)
!!
!     include 'function2'
!!

! set up gfs internal state (gis_phy) dimension and values for physics etc
!-------------------------------------------------------------------
      me     = gis_phy%me
      nodes  = gis_phy%nodes

      call compns_physics(gis_phy%deltim, gis_phy%iret,  gis_phy%ntrac,	&
                          gis_phy%nxpt,   gis_phy%nypt,  gis_phy%jintmx,&
                          gis_phy%jcap,   gis_phy%levs,  gis_phy%levr,  &
                          gis_phy%lonr,   gis_phy%latr,                 &
                          gis_phy%ntoz,   gis_phy%ntcw,  gis_phy%ncld,  &
                          gis_phy%lsoil,  gis_phy%nmtvr,                &
                          gis_phy%num_p3d, gis_phy%num_p2d,             &
                          gis_phy%thermodyn_id, gis_phy%sfcpress_id,    &
                          gis_phy%nam_gfs_phy%nlunit, me,               &
                          gis_phy%nam_gfs_phy%gfs_phy_namelist)
!
      CALL set_soilveg(me,gis_phy%nam_gfs_phy%nlunit)
      call set_tracer_const(gis_phy%ntrac,me,gis_phy%nam_gfs_phy%nlunit)
!
      nlunit  = gis_phy%nam_gfs_phy%nlunit
      ntrac   = gis_phy%ntrac
      nxpt    = gis_phy%nxpt
      nypt    = gis_phy%nypt
      jintmx  = gis_phy%jintmx
      jcap    = gis_phy%jcap
      levs    = gis_phy%levs
      levr    = gis_phy%levr
      lonr    = gis_phy%lonr
      latr    = gis_phy%latr
      ntoz    = gis_phy%ntoz
      ntcw    = gis_phy%ntcw
      ncld    = gis_phy%ncld
      lsoil   = gis_phy%lsoil
      nmtvr   = gis_phy%nmtvr
      num_p3d = gis_phy%num_p3d
      num_p2d = gis_phy%num_p2d
      thermodyn_id = gis_phy%thermodyn_id
      sfcpress_id  = gis_phy%sfcpress_id
      if (gis_phy%nam_gfs_phy%Total_Member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",I2.2)') gis_phy%nam_gfs_phy%Member_Id
      endif
!
!     ivssfc  = 200501
      ivssfc  = 200509
      ivssfc_restart  = 200509
      if (ivssfc .gt. ivssfc_restart) ivssfc_restart = ivssfc
      ivsupa  = 0
      if (levs .gt. 99) ivsupa  = 200509
!
      levh   = ntrac*levs
      latrd  = latr + 2*jintmx
      latr2  = latr/2
      levm1  = levs-1 
      levp1  = levs+1 
      lonrx  = lonr + 1 + 2*nxpt + 1
!
      ngrids_sfcc = 32+LSOIL*3   ! No CV, CVB, CVT! includes T2M, Q2M, TISFC
      ngrids_flx  = 66+30        ! additional fields (most related to land surface)
      nfxr        = 27
      ngrids_gg   = 2+LEVS*(4+ntrac)
!
      allocate ( lon_dims_r(latr), stat = ierr )
      allocate ( lon_dims_ext(latrd), stat = ierr )
!
      allocate(colrad_r(latr2), stat = ierr)
      allocate(wgt_r(latr2), stat = ierr)
      allocate(wgtcs_r(latr2), stat = ierr)
      allocate(rcs2_r(latr2), stat = ierr)
      allocate(sinlat_r(latr), stat = ierr)
      allocate(coslat_r(latr), stat = ierr)
!
      allocate(buff_mult(lonr,latr,ngrids_sfcc), stat = ierr)
!
      allocate ( gis_phy% lats_nodes_r(nodes), stat = ierr )
      allocate ( gis_phy% lats_nodes_ext(nodes), stat = ierr )
      allocate ( gis_phy% global_lats_r(latr), stat = ierr )
      allocate ( gis_phy% global_lats_ext(latr), stat = ierr )
      allocate ( gis_phy% lonsperlar(latr), stat = ierr)

      call set_lonsgg(gis_phy%lonsperlar)
!
      g_gz   = 1
      g_ps   = g_gz  + 1
      g_t    = g_ps  + 1     
      g_u    = g_t   + levs
      g_v    = g_u   + levs
      g_q    = g_v   + levs
      g_p    = g_q   + levh
      g_dp   = g_p   + levs
      g_dpdt = g_dp  + levs
       
      lotgr  = g_dpdt+ levs - 1

      gis_phy%g_gz    = g_gz  
      gis_phy%g_ps    = g_ps  
      gis_phy%g_t     = g_t  
      gis_phy%g_u     = g_u  
      gis_phy%g_v     = g_v  
      gis_phy%g_q     = g_q  
      gis_phy%g_p     = g_p  
      gis_phy%g_dp    = g_dp 
      gis_phy%g_dpdt  = g_dpdt 
!
      gis_phy%lotgr = lotgr


      print *,' g_gz ',g_gz
      print *,' g_ps ',g_ps
      print *,' g_t  ',g_t 
      print *,' g_u  ',g_u 
      print *,' g_v  ',g_v 
      print *,' g_q  ',g_q 
      print *,' g_p  ',g_p 
      print *,' g_dp ',g_dp
      print *,' g_dpdt ',g_dpdt
      print *,' lotgr ',lotgr

      if (ras) then
        nrcm = max(nrcmax, nint((nrcmax*gis_phy%deltim)/600.0))
      else
        nrcm = 1
      endif
!
      if (ntoz .le. 0) then      ! Diagnostic ozone
        rewind (kozc)
        read (kozc,end=101) latsozc, levozc, timeozc, blatc4
  101   if (levozc .lt. 10 .or. levozc .gt. 100) then
          rewind (kozc)
          levozc  = 17
          latsozc = 18
          blatc   = -85.0
        else
          blatc   = blatc4
        endif
        latsozp   = 2
        levozp    = 1
        timeoz    = 1
        pl_coeff  = 0
      else                       ! Prognostic Ozone
        rewind (kozpl)
        read (kozpl) pl_coeff, latsozp, levozp, timeoz
        allocate (pl_lat(latsozp), pl_pres(levozp),pl_time(timeoz+1), stat = ierr)
        allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1), stat = ierr)
        rewind (kozpl)
        read (kozpl) pl_coeff, latsozp, levozp, timeoz, pl_lat4, pl_pres4,  &
                     pl_time4
        pl_pres(:) = pl_pres4(:)
        pl_lat(:)  = pl_lat4(:)
        pl_time(:) = pl_time4(:)
        latsozc = 2
        blatc   = 0.0
      endif
      dphiozc = -(blatc+blatc)/(latsozc-1)
!
      if (me .eq. 0) then
        print *,' latsozp=',latsozp,' levozp=',levozp,' timeoz=',timeoz
        print *,' latsozc=',latsozc,' levozc=',levozc,' timeozc=',        &
                  timeozc, 'dphiozc=',dphiozc
        print *,' pl_lat=',pl_lat
        print *,' pl_pres=',pl_pres
        print *,' pl_time=',pl_time
      endif
      pl_pres(:) = log(0.1*pl_pres(:))       ! Natural log of pres in cbars
!
      allocate(gis_phy%OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz), stat = ierr) !OZONE P-L coeffcients
!

      print *,' finished array allocation in gfs_physics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
      if (me == 0) write(*,*) 'io option ,liope :',liope
!
      call mpi_comm_dup(mpi_comm_all, mpi_comm_all_dup, ierr)
      call mpi_barrier (mpi_comm_all_dup,               ierr)

      if (nodes == 1) liope=.false.
      if (liope) then
        call mpi_comm_rank(mpi_comm_all_dup,nrank_all,ierr)
        icolor=1
        ikey=1
        nodes_comp=nodes-1
        if (nrank_all.eq.nodes-1) then
!!  io server
          write(*,*) 'io server task'
          icolor=2
          gis_phy%kcolor=mpi_undefined
          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
          call mpi_comm_split(mpi_comm_all_dup,gis_phy%kcolor,ikey,mc_comp,ierr)
        else
!sela     write(*,*) 'compute server task '
          icolor=mpi_undefined
          gis_phy%kcolor=1
          call mpi_comm_split(mpi_comm_all_dup,gis_phy%kcolor,ikey,mc_comp,ierr)
          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
          call mpi_comm_size(mc_comp,nodes,ierr)
        endif
      else
        icolor=2
        mc_comp=mpi_comm_all_dup
        nodes_comp=nodes
      endif
!!
!c
!     call f_hpminit(me,"evod")  !jjt hpm stuff
!c
!     call f_hpmstart(25,"get_ls_gftlons")
!c
!!
      call synchro
      call init_countperf(latr)
!$$$      time0=timer()
!jfe  call countperf(0,15,0.)
!
      if (me.eq.0) then
      print 100, jcap,levs
100   format (' smf ',i3,i3,' in gfs physics initialize ')
      print*,'number of threads is ',num_parthds()
        if (liope) then
          print*,'number of mpi procs is ',nodes
          print*,'number of mpi io procs is 1 (nodes)'
        else
          print*,'number of mpi procs is ',nodes
        endif
      endif
!c
      gis_phy%cons0    =    0.0d0     !constant
      gis_phy%cons0p5  =    0.5d0     !constant
      gis_phy%cons1200 = 1200.d0      !constant
      gis_phy%cons3600 = 3600.d0      !constant
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!c
      if(gis_phy%iret.ne.0) then
        if(me.eq.0) print *,' incompatible physics namelist -',         &
                            ' aborted in gfs_phy_initilize'
        call mpi_quit(13)
      endif
!!
!     if predicted ozon is desired set jo3=2
      jo3=2          !using predicted ozone in radiation.
!
      gis_phy%lats_nodes_ext = 0

      call getcon_physics(gis_phy%n3,gis_phy%n4,                        &
           gis_phy%lats_nodes_r,gis_phy%global_lats_r,                  &
           gis_phy%lonsperlar,                                          &
           gis_phy%lats_nodes_ext,gis_phy%global_lats_ext,              &
           gis_phy%colat1)

      gis_phy%lats_node_r_max = lats_node_r_max

      call sfcvar_aldata(lonr,lats_node_r,lsoil,gis_phy%sfc_fld,ierr)
      call flxvar_aldata(lonr,lats_node_r,gis_phy%flx_fld,ierr)

      print *,' check after sfc flx var_aldata ' 

      allocate (   gis_phy%grid_gr(lonr*lats_node_r_max,lotgr), stat = ierr )
      ALLOCATE (   gis_phy%XLON(LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%XLAT(LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%COSZDG(LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%CLDCOV(LEVS,LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%SFALB(LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%HPRIME(NMTVR,LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%FLUXR(nfxr,LONR,LATS_NODE_R), stat = ierr)

      gis_phy%NBLCK=LONR/NGPTC+1

      ALLOCATE (   gis_phy%SWH(NGPTC,LEVS,gis_phy%NBLCK,LATS_NODE_R), stat = ierr)
      ALLOCATE (   gis_phy%HLW(NGPTC,LEVS,gis_phy%NBLCK,LATS_NODE_R), stat = ierr)
!
      ALLOCATE (gis_phy%JINDX1(LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%JINDX2(LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%DDY(LATS_NODE_R), stat = ierr)
!
      allocate (gis_phy%phy_f3d(NGPTC,LEVS,gis_phy%NBLCK,lats_node_r,num_p3d), stat = ierr)
      allocate (gis_phy%phy_f2d(lonr,lats_node_r,num_p2d), stat = ierr)
!
      allocate (   gis_phy%fhour_idate(1,5), stat = ierr )
!
!     print *,' check after lots allocates ' 

      if (ldiag3d) then
        call d3d_init(ngptc,gis_phy%nblck,lonr,lats_node_r,levs,pl_coeff)
      endif

!     if (icolor /= 2 .or. .not. liope) then
        if (num_p3d .gt. 0) gis_phy%phy_f3d = 0.0
        if (num_p2d .gt. 0) gis_phy%phy_f2d = 0.0
!     endif
!!
      call countperf(0,18,0.)
!!
      call fix_fields(gis_phy%LONSPERLAR,gis_phy%GLOBAL_LATS_R,           &
        gis_phy%XLON,gis_phy%XLAT,gis_phy%sfc_fld,                        &
        gis_phy%HPRIME,gis_phy%JINDX1,gis_phy%JINDX2,gis_phy%DDY,         &
        gis_phy%OZPLIN,gis_phy%nam_gfs_phy%sfc_ini)

      print *,' finish fix_fields '
!!
      call countperf(1,18,0.)
!!
      gis_phy%zhour=fhour
      gis_phy%FLUXR=0.
!
      call flx_init(gis_phy%flx_fld,ierr)
!
      print *,' finish fix_init '

      if (ldiag3d) then
        call d3d_zero
      endif
!
!
      end subroutine gfs_physics_initialize
!
! =========================================================================
!
      subroutine set_lonsgg(lonsperlat)
      use resol_def
      use reduce_lons_grid_module, only : reduce_grid           ! hmhj
      integer numreduce                                         ! hmhj
      integer lonsperlat(latr)
       
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
         call reduce_grid (numreduce,jcap,latr,lonsperlat)              ! hmhj
         print*,' reduced grid is computed - lonsperlat '               ! hmhj
      endif
                                                                                
      print*,' jcap = ',jcap
      print*,'min,max of lonsperlat = ',minval(lonsperlat),              &
              maxval(lonsperlat)
      end subroutine set_lonsgg
!
      end module gfs_physics_initialize_mod
