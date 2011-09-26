!
! !module: gfs_physics_run_mod --- run module of the grided
!                              component of the gfs physics.
!
! !description: gfs run module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang, updated to the new version of gfs.
!  janusry  2007      hann-ming henry juang for gfs dynamics only
!  july     2007      shrinivas moorthi for gfs physics only
!  november 2007      hann-ming henry juang continue for gfs physics
!  october  2009      jun wang add nsout option
!  oct 11   2009      sarah lu, grid_gr replaced by grid_fld
!  oct 17   2009      sarah lu, q is replaced by tracers(1)
!  dec 08   2009      sarah lu, add g3d_fld to do_physics_one_step
!                     calling argument
!  July     2010      Shrinivas Moorthi - Updated for new physics and added nst
!                     eliminated calls to common_vars
!  jul 21  2010       sarah lu, add g2d_fld to do_physics_one_step
!                     calling argument
!  Aug 03  2010       jun wang, fix lsout for dfi
!  Aug 25  2010       Jun Wang, add zhour_dfi for filtered dfi fields output
!  Oct 18  2010       s. moorthi added fscav to do tstep
!  Dec 23  2010       Sarah Lu, setup fscav from gfs_phy_tracer 
!
! !interface:
!
      module gfs_physics_run_mod
!
!!uses:
!
      use gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state
      USE date_def,                       ONLY: fhour
      USE namelist_physics_def,           ONLY: nsout,ldfi,ndfi

      implicit none

      contains

      subroutine gfs_physics_run(gis_phy, rc)

      type(gfs_physics_internal_state), pointer, intent(inout) :: gis_phy
      integer, optional,                          intent(out)   :: rc

      real , save	:: timestep=0.0
      integer		rc1, k , i1, i2

!***********************************************************************
!
!     lsfwd      logical true during a first forward step
!     lssav      logical true during a step for which
!                diagnostics are accumulated
!     lscca      logical true during a step for which convective clouds
!                are calculated from convective precipitation rates
!     phour      real forecast hour at the end of the previous timestep
!
!     lsout controls freq. of output
!     nsout controls freq. of output in time steps
!     fhout controls freq. of output in hours
!     nszer time steps between zeroing fluxes
!
!***********************************************************************

!       print *,' enter gfs_physics_run '
!
       rc1 = 0
!
!      print *,' uug=',gis_phy%grid_gr(1,gis_phy%g_u:gis_phy%g_u+gis_phy%levs-1)
!      print *,' pg=',gis_phy%grid_gr(1,gis_phy%g_p:gis_phy%g_p+gis_phy%levs-1)
!      print *,' dpg=',gis_phy%grid_gr(1,gis_phy%g_dp:gis_phy%g_dp+gis_phy%levs-1)
!      call common_to_physics_vars(gis_phy%grid_fld%ps,      &
!                                  gis_phy%grid_fld%t ,      &
!*                                 gis_phy%grid_fld%tracers(1)%flds ,  &
!                                  gis_phy%grid_fld%u ,      &
!                                  gis_phy%grid_fld%v ,      &
!                                  gis_phy%grid_fld%p ,      &
!                                  gis_phy%grid_fld%dp ,     &
!                                  gis_phy%grid_fld%dpdt ,   &
!                                  gis_phy%global_lats_r,                &
!                                  gis_phy%lonsperlar)

!      print *,' end of common_to_physics_vars '

!
! ---------------------------------------------------------------------
! ======================================================================
!                     do one physics time step
! ---------------------------------------------------------------------
       if (.not.ldfi) then
        gis_phy%LSOUT = MOD(gis_phy%kdt ,NSOUT).EQ.0  .or. gis_phy%kdt==1
       else
        gis_phy%LSOUT = (MOD(gis_phy%kdt ,NSOUT).EQ.0 .and.              &
           (gis_phy%kdt<=ndfi/2.or.gis_phy%kdt>ndfi)).or. gis_phy%kdt==1
       endif
!
!       write(0,*)' end of common_to_physics_vars,kdt=',gis_phy%kdt,      &
!         'nsout=',nsout,'lsout=',gis_phy%LSOUT,'zhour=',gis_phy%ZHOUR,   &
!         'ldfi=',ldfi,'ndfi=',ndfi,gis_phy%kdt<=ndfi/2,gis_phy%kdt>ndfi, &
!             gis_phy%kdt<=ndfi/2.or.gis_phy%kdt>ndfi
!       if(gis_phy%kdt>=96.and.gis_phy%kdt<=98.or.gis_phy%kdt>=4.and.gis_phy%kdt<=6) then
!       print *,'be phys one,kdt=',gis_phy%kdt,'ps=',maxval(gis_phy%grid_fld%ps), &
!        minval(gis_phy%grid_fld%ps),'t=',maxval(gis_phy%grid_fld%t), &
!        minval(gis_phy%grid_fld%t),'spfh=',maxval(gis_phy%grid_fld%tracers(1)%flds),  &
!        minval(gis_phy%grid_fld%tracers(1)%flds)
!       endif
        
!!  gfs_phy_tracer%fscav(1:ntrac) is for all tracers
!!  fscav(1:ntrac-ncld-1) is for ozone + aerosl species

        if ( gis_phy%lgocart ) then
           i1 = gis_phy%gfs_phy_tracer%ntrac_met+1   ! 1st chemical tracer (excluding o3)
           i2 = gis_phy%gfs_phy_tracer%ntrac         ! last chemical tracer (excluding o3)
           gis_phy%fscav(2:gis_phy%ntrac-gis_phy%ncld-1) =      &
                  gis_phy%gfs_phy_tracer%fscav(i1:i2)
       
           if ( gis_phy%kdt .eq. 2 ) then
	      print *, 'LU_TRC:',i1,i2,                        &
                  gis_phy%gfs_phy_tracer%fscav(i1:i2)
           endif
        endif
 
!
! ======================================================================
        call do_physics_one_step(                                         &
                 gis_phy%deltim,   gis_phy%kdt,     gis_phy%phour,        &
!*               gis_phy%grid_gr,  gis_phy%sfc_fld, gis_phy%flx_fld,      &
                 gis_phy%grid_fld, gis_phy%sfc_fld, gis_phy%flx_fld,      &
                 gis_phy%nst_fld,  gis_phy%g3d_fld, gis_phy%g2d_fld,      &
                 gis_phy%lats_nodes_r,   gis_phy%global_lats_r,           &
                 gis_phy%lonsperlar,                                      &
                 gis_phy%XLON,    gis_phy%XLAT,    gis_phy%COSZDG,        &
                 gis_phy%HPRIME,  gis_phy%SWH,     gis_phy%HLW,           &
                 gis_phy%FLUXR,   gis_phy%SFALB,                          &
                 gis_phy%SLAG,    gis_phy%SDEC,    gis_phy%CDEC,          &
                 gis_phy%OZPLIN,  gis_phy%JINDX1,  gis_phy%JINDX2,        &
                 gis_phy%DDY,                                             &
                 gis_phy%phy_f3d, gis_phy%phy_f2d, gis_phy%NBLCK,         &
                 gis_phy%ZHOUR,   gis_phy%ZHOUR_DFI,                      &
                 gis_phy%N3,      gis_phy%N4,                             &
                 gis_phy%LSOUT,   gis_phy%COLAT1,  gis_phy%CFHOUR1,       &
                 gis_phy%fscav )

!       if(gis_phy%kdt>=96.and.gis_phy%kdt<=98.or.gis_phy%kdt>=4.and.gis_phy%kdt<=6) then
!       if(gis_phy%kdt<=1) then
!        print *,'af phys one,kdt=',gis_phy%kdt,'ps=',maxval(gis_phy%grid_fld%ps), &
!        minval(gis_phy%grid_fld%ps),'t=',maxval(gis_phy%grid_fld%t), &
!        minval(gis_phy%grid_fld%t),'spfh=',maxval(gis_phy%grid_fld%tracers(1)%flds),  &
!        minval(gis_phy%grid_fld%tracers(1)%flds)
!       endif
!                        
! =======================================================================
!
! update hour
!
      gis_phy%phour      = fhour

! --------------------------------------------------------------------------
!      call physics_to_common_vars(gis_phy%grid_fld%ps,      &
!                                  gis_phy%grid_fld%t ,      &
!*                                 gis_phy%grid_fld%q ,      &
!                                  gis_phy%grid_fld%tracers(1)%flds, &
!                                  gis_phy%grid_fld%u ,      &
!                                  gis_phy%grid_fld%v ,      &
!                                  gis_phy%grid_fld%p ,      &
!                                  gis_phy%grid_fld%dp ,     &
!                                  gis_phy%grid_fld%dpdt ,   &
!                                  gis_phy%global_lats_r,                &
!                                  gis_phy%lonsperlar)
! --------------------------------------------------------------------------
!
      if(present(rc)) then
          rc = rc1
      end if

      end subroutine gfs_physics_run

      end module gfs_physics_run_mod
