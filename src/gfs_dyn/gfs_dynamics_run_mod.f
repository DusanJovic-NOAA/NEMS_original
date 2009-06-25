!
! !module: gfs_dynamics_run_mod --- run module of the grided
!                              component of the gfs dynamics system.
!
! !description: gfs run module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang, updated to the new version of gfs.
!  janusry  2007      hann-ming henry juang for gfs dynamics only
!
!
! !interface:
!
      module gfs_dynamics_run_mod
!
!!uses:
!
      use gfs_dynamics_internal_state_mod

      implicit none

      contains

      subroutine gfs_dynamics_run(gis_dyn, rc)

      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      integer, optional,                          intent(out)   :: rc

!     real , save	:: timestep=0.0
      integer, save     :: kdt_save=0
      integer		rc1, k 
!***********************************************************************
!
!     lsfwd      logical true during a first forward step
!     phour      real forecast hour at the end of the previous timestep
!     fhrot      if =0 read one time level sigmas and te=tem zem=ze ...
!
!     lsout controls freq. of output
!     nsout controls freq. of output in time steps
!     fhout controls freq. of output in hours
!
!     nsres time steps between writing restart files
!
!***********************************************************************

       rc1 = 0.0
! ---------------------------------------------------------------------
! change temperature and pressure back to model grid value at n+1  slot.
!
       write(0,*)' before of common_to_model_vars,kdt=',gis_dyn%kdt,      &
         'nsout=',nsout,'lsout=',gis_dyn%LSOUT,'t=',  &
         maxval(gis_dyn%grid_gr(:,gis_dyn%g_t)),minval(gis_dyn%grid_gr(:,gis_dyn%g_t))
      if( .not. gis_dyn% start_step ) then
!       print *,' change common variables to model variables '
        gis_dyn%lsout = mod(gis_dyn%kdt,nsout).eq.0
        call common_to_model_vars (gis_dyn%grid_gr(1,gis_dyn%g_zq),      &
                                   gis_dyn%grid_gr(1,gis_dyn%g_t ),      &
                                   gis_dyn%grid_gr(1,gis_dyn%g_rt),      &
                                   gis_dyn%grid_gr(1,gis_dyn%g_u ),      &
                                   gis_dyn%grid_gr(1,gis_dyn%g_v ),      &
                                   gis_dyn%global_lats_a,	           &
                                   gis_dyn%lonsperlat,     		   &
                                   gis_dyn%pwat, gis_dyn%ptot )
!       lsfwd=.false.
      else
        print *,' the first time step, model running from internal state.'
!       lsfwd=.true.
      endif
!
! ---------------------------------------------------------------------
! check whether time step is changes to update matrix
!
!     if( gis_dyn%deltim .ne. timestep ) then
!       timestep=gis_dyn%deltim
!       if( lsfwd ) timestep = 0.5 * timestep
!       print *,' call get_cd for deltim=',timestep
!       if(hybrid)then
!         call get_cd_hyb(timestep)
!       else if( gen_coord_hybrid ) then           
!         call get_cd_hyb_gc(timestep)    
!       else
!         call get_cd_sig(am,bm,timestep,tov,sv)
!       endif
!     endif
!
! ---------------------------------------------------------------------
! check whether reset step due to dfi
!
      if( gis_dyn%kdt .lt. kdt_save ) then
        gis_dyn%reset_step = .true.
      endif
      kdt_save = gis_dyn%kdt

! ======================================================================
!                     do one time step with one-loop
! ---------------------------------------------------------------------
      if( gis_dyn% spectral_loop == 1 ) then
!
        call  do_dynamics_one_loop(    gis_dyn% deltim         ,	&
               gis_dyn% kdt           ,gis_dyn% phour          ,	&
               gis_dyn% trie_ls       ,gis_dyn% trio_ls        ,	&
               gis_dyn% grid_gr       ,					&
               gis_dyn% ls_node       ,gis_dyn% ls_nodes       ,	&
               gis_dyn% max_ls_nodes  ,					&
               gis_dyn% lats_nodes_a  ,gis_dyn% global_lats_a  ,	&
               gis_dyn% lonsperlat    ,					&
               gis_dyn% lats_nodes_ext,gis_dyn% global_lats_ext,	&
               gis_dyn% epse          ,gis_dyn% epso           ,	&
               gis_dyn% epsedn        ,gis_dyn% epsodn         ,	&
               gis_dyn% snnp1ev       ,gis_dyn% snnp1od        ,	&
               gis_dyn% ndexev        ,gis_dyn% ndexod         ,	&
               gis_dyn% plnev_a       ,gis_dyn% plnod_a        ,	&
               gis_dyn% pddev_a       ,gis_dyn% pddod_a        ,	&
               gis_dyn% plnew_a       ,gis_dyn% plnow_a        ,	&
               gis_dyn% syn_ls_a      ,gis_dyn% dyn_ls_a       ,	&
               gis_dyn% syn_gr_a_1    ,gis_dyn% dyn_gr_a_1     ,	&
               gis_dyn% anl_gr_a_1    ,gis_dyn% sym_gr_a_2     ,        &
               gis_dyn% syn_gr_a_2    ,gis_dyn% dyn_gr_a_2     ,	&
               gis_dyn% anl_gr_a_2    ,gis_dyn% lslag          ,        &
               gis_dyn% pwat          ,gis_dyn% ptot           ,	&
               gis_dyn% pdryini       ,gis_dyn% nblck          ,	&
               gis_dyn% zhour         ,					&
               gis_dyn% n1            ,gis_dyn% n4             ,	&
               gis_dyn% lsout	      ,					&
               gis_dyn% colat1        ,gis_dyn% cfhour1	       ,	&
               gis_dyn% start_step    ,                                 &
               gis_dyn% reset_step    ,gis_dyn% end_step       )
        write(0,*)'after gfs_dyn_oneloop_run, t=',  &
         maxval(gis_dyn%grid_gr(:,gis_dyn%g_t)),minval(gis_dyn%grid_gr(:,gis_dyn%g_t))

!
! ======================================================================
!                     do one time step with two-loop
! ---------------------------------------------------------------------
      else if( gis_dyn% spectral_loop == 2 ) then
!
        call  do_dynamics_two_loop(    gis_dyn% deltim         ,	&
               gis_dyn% kdt           ,gis_dyn% phour 	       ,	&
               gis_dyn% trie_ls       ,gis_dyn% trio_ls        ,	&
               gis_dyn% grid_gr       ,	 				&
               gis_dyn% ls_node       ,gis_dyn% ls_nodes       ,	&
               gis_dyn% max_ls_nodes  ,					&
               gis_dyn% lats_nodes_a  ,gis_dyn% global_lats_a  ,	&
               gis_dyn% lonsperlat    ,					&
               gis_dyn% lats_nodes_ext,gis_dyn% global_lats_ext,	&
               gis_dyn% epse          ,gis_dyn% epso           ,	&
               gis_dyn% epsedn	      ,gis_dyn% epsodn         ,	&
               gis_dyn% snnp1ev       ,gis_dyn% snnp1od        ,	&
               gis_dyn% ndexev        ,gis_dyn% ndexod         ,	&
               gis_dyn% plnev_a       ,gis_dyn% plnod_a        ,	&
               gis_dyn% pddev_a       ,gis_dyn% pddod_a        ,	&
               gis_dyn% plnew_a       ,gis_dyn% plnow_a        ,	&
               gis_dyn% syn_ls_a      ,gis_dyn% dyn_ls_a       ,	&
               gis_dyn% syn_gr_a_1    ,gis_dyn% dyn_gr_a_1     ,	&
               gis_dyn% anl_gr_a_1    ,gis_dyn% sym_gr_a_2     ,        &
               gis_dyn% syn_gr_a_2    ,gis_dyn% dyn_gr_a_2     ,	&
               gis_dyn% anl_gr_a_2    ,gis_dyn% lslag          ,        &
               gis_dyn% pwat          ,gis_dyn% ptot           ,	&
               gis_dyn% pdryini       ,gis_dyn% nblck          ,	&
               gis_dyn% zhour         ,					&
               gis_dyn% n1            ,gis_dyn% n4             ,	&
               gis_dyn% lsout	      ,					&
               gis_dyn% colat1        ,gis_dyn% cfhour1	,		&	
               gis_dyn%start_step     ,                                 &
               gis_dyn%reset_step     ,gis_dyn% end_step       )

      else

        print *,' number of spectral loop is wrong. it is ',		&
                gis_dyn%spectral_loop
        call mpi_quit(99)
        stop

      endif
!
      write(0,*)'in gfs_dynamics,after one loop,grid_gr(zq)=',             &
       maxval(gis_dyn%grid_gr(:,gis_dyn%g_zq)),minval(gis_dyn%grid_gr(:,gis_dyn%g_zq)), &
       't=',maxval(gis_dyn%grid_gr(:,gis_dyn%g_t)),minval(gis_dyn%grid_gr(:,gis_dyn%g_t))
! =======================================================================
!                   end of one- or two-loop computation
! =======================================================================
!
! check whether wind speed is out of bound.
!
!jw      if (.not.liope.or.icolor.ne.2) then
        do k=1,levs
          if(spdmax(k).gt.0. .and. spdmax(k).lt.1000.) then
            continue
          else
            print *,'unphysical maximum speed',spdmax(k),' me=',me
            call mpi_quit(7)
            stop
          endif
        enddo
!jw      endif


! =======================================================================

! update hour
!
      gis_dyn% phour      = fhour

! ---------------------------------------------------------------------
! prepare n+1 grid value back to common temperature and pressure
!
      call model_to_common_vars (gis_dyn% grid_gr(1,gis_dyn%g_zq),       &
                                 gis_dyn% grid_gr(1,gis_dyn%g_t ),       &
                                 gis_dyn% grid_gr(1,gis_dyn%g_rt),       &
                                 gis_dyn% grid_gr(1,gis_dyn%g_u ),       &
                                 gis_dyn% grid_gr(1,gis_dyn%g_v ),       & 
                                 gis_dyn% grid_gr(1,gis_dyn%g_p ),       & 
                                 gis_dyn% grid_gr(1,gis_dyn%g_dp ),      & 
                                 gis_dyn% grid_gr(1,gis_dyn%g_dpdt ),    & 
                                 gis_dyn%global_lats_a,	                   &
                                 gis_dyn%lonsperlat)
      write(0,*)'after n+1grid value,t=',  &
         maxval(gis_dyn%grid_gr(:,gis_dyn%g_t)),minval(gis_dyn%grid_gr(:,gis_dyn%g_t))

!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c
      if(present(rc)) then
          rc = rc1
      end if

      end subroutine gfs_dynamics_run

      end module gfs_dynamics_run_mod
