      SUBROUTINE input_fields_rst(gread,gread2,cread, cread2, 
     &                 PDRYINI,TRIE_LS,TRIO_LS,
     &                 grid_gr,LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 SNNP1EV,SNNP1OD,global_lats_a,lonsperlat,
     &                 epse,epso,plnev_a,plnod_a,plnew_a,plnow_a,
     &                 lats_nodes_a)
!!
!program log
!  20100205  J. WANG     Read in input restart files without computing 
!                        pwat nad ptot
!  20100908  J. WANG     remove gfsio module
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, cp => con_cp , rd => con_rd
      use gfs_dyn_coordinate_def
      IMPLICIT NONE
!!
      REAL(KIND=kind_grid) PDRYINI
      CHARACTER (len=*)   :: CREAD, CREAD2,gread,gread2
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTLS)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTLS)
      REAL(KIND=KIND_GRID) GRID_GR(lonf,lats_node_a_max,lotgr)
      REAL(KIND=KIND_EVOD) SNNP1EV(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) SNNP1OD(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) EPSE   (LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) EPSO   (LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) PLNEV_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOD_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNEW_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOW_a(LEN_TRIE_LS,latg2)
!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
      real(kind=kind_grid) ttg(lonf,lats_node_a,levs)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!
      integer nblck
      integer global_lats_a(latg), lonsperlat(latg)
!
      INTEGER              LS_NODE (LS_DIM*3)
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES(NODES)
      integer            lats_nodes_a(nodes)
!
      INTEGER              IERR,IPRINT,J,JDT,K,L,LOCL,N,i
      REAL(KIND=KIND_EVOD) TEE1(LEVS)
      REAL(KIND=KIND_EVOD)  YE1(LEVS)
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS2=2.0,
     &                                   CONS600=600.0
      LOGICAL LSLAG
      integer		lan, lat, lons_lat, jlonf,nnl,nn,kk,lon
!
!------------------------------------------------------------------
!
      if(me.eq.0) PRINT  9876,FHOUR,idate
 9876 FORMAT(1H ,'FHOUR IN input_fields ',F6.2,
     & ' idate no yet read in',4(1x,i4))
      IPRINT = 0
!$$$  IF ( ME .EQ. 0 ) IPRINT = 1
!
!--- n time step grid 
! 
      print *,'in inputfile_rst,gread=',gread,'gread2=',gread2,
     & 'cread=',cread,'cread2=',cread2
      if (me .eq. 0) write(0,*)' gread=',gread,'ntoz=',ntoz
        CALL TREADG_nemsio(gread,FHOUR,IDATE,
     &               Zsg, psg, ttg, uug, vvg, rqg,
     X               PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat) 
       print *,'after treadg_nemsio,mypsg n-1=',psg(5,3)

      do j=1,lats_node_a
        grid_gr(1:lonf,j,g_gz) = zsg(1:lonf,j)
        grid_gr(1:lonf,j,g_qm) = psg(1:lonf,j)
      enddo
      do k=1,levs
        do j=1,lats_node_a
          grid_gr(1:lonf,j,g_ttm+k-1) = ttg(1:lonf,j,k)
          grid_gr(1:lonf,j,g_uum+k-1) = uug(1:lonf,j,k)
          grid_gr(1:lonf,j,g_vvm+k-1) = vvg(1:lonf,j,k)
        enddo
      enddo
      do k=1,levh
        do j=1,lats_node_a
          grid_gr(1:lonf,j,g_rm +k-1) = rqg(1:lonf,j,k)
        enddo
      enddo
!
!---------------------------------------------------------------
!--------------------------------------------------------
!
!--- n time step grid 
! 
      if (me .eq. 0) write(0,*)' gread2=',gread2
          CALL TREADG_nemsio(gread2,fhour,idate,
     &                 zsg, psg, ttg, uug, vvg, rqg,
     X                 PDRYINI,IPRINT,
     &                 global_lats_a,lats_nodes_a,lonsperlat) 
       print *,'after treadg_nemsio,mypsg n=',psg(1:5,3)
!
      do j=1,lats_node_a
!       grid_gr(1:lonf,j,g_gz) = zsg(1:lonf,j)
        grid_gr(1:lonf,j,g_q ) = psg(1:lonf,j)
      enddo
      do k=1,levs
        do j=1,lats_node_a
          grid_gr(1:lonf,j,g_tt +k-1) = ttg(1:lonf,j,k)
          grid_gr(1:lonf,j,g_uu +k-1) = uug(1:lonf,j,k)
          grid_gr(1:lonf,j,g_vv +k-1) = vvg(1:lonf,j,k)
        enddo
      enddo
      do k=1,levh
        do j=1,lats_node_a
          grid_gr(1:lonf,j,g_rq +k-1) = rqg(1:lonf,j,k)
        enddo
      enddo
!
! =======================================================================
!
!-------------------------------------------------------------------------
!---read in n-1 time step spectral file
      if (me .eq. 0) write(0,*)'in input, sread1, cread=',trim(cread)
          CALL TREADS_nemsio(cread,FHOUR,IDATE,
     X                 TRIE_LS(1,1,P_GZ), TRIE_LS(1,1,P_QM ),
     X                 TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_DIM),
     X                 TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_RM),
     X                 TRIO_LS(1,1,P_GZ), TRIO_LS(1,1,P_QM ),
     X                 TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_DIM),
     X                 TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_RM),
     X                 LS_NODE,LS_NODES,MAX_LS_NODES,
     X                 SNNP1EV,SNNP1OD,
     &                 epse, epso, plnew_a, plnow_a, 
     &                 plnev_a, plnod_a)
        print *,'after treads n-1, P-GZ=',P_GZ,'p_QM=',P_QM,'P_TEM=',
     &  P_TEM, 'P_DIM=',P_DIM,'P_ZEM=',P_ZEM,'P_RM=',P_RM
!---
      fhini=fhour
      if(me.eq.0) PRINT 9877, FHOUR
 9877 FORMAT(1H ,'FHOUR AFTER TREAD',F6.2)

      if (me .eq. 0) write(0,*)' fhini=',fhini,'last_fcst_pe=',
     &     last_fcst_pe,'fhrot=',fhrot
      if (me<=last_fcst_pe) then
        CALL RMS_spect(TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_DIM),
     X             TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_ZEM),
     X             TRIE_LS(1,1,P_RM ),
     X             TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_DIM),
     X             TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_ZEM),
     X             TRIO_LS(1,1,P_RM ),
     X             LS_NODES,MAX_LS_NODES)
      endif

!
!--- N time step spectral 
      if (me .eq. 0) write(0,*)'in input, sread2, cread=',trim(cread2)
          CALL TREADS_nemsio(cread2,FHOUR,IDATE,
     X                 TRIE_LS(1,1,P_GZ), TRIE_LS(1,1,P_Q ),
     X                 TRIE_LS(1,1,P_TE), TRIE_LS(1,1,P_DI),
     X                 TRIE_LS(1,1,P_ZE), TRIE_LS(1,1,P_RQ),
     X                 TRIO_LS(1,1,P_GZ), TRIO_LS(1,1,P_Q ),
     X                 TRIO_LS(1,1,P_TE), TRIO_LS(1,1,P_DI),
     X                 TRIO_LS(1,1,P_ZE), TRIO_LS(1,1,P_RQ),
     X                 LS_NODE,LS_NODES,MAX_LS_NODES,
     X                 SNNP1EV,SNNP1OD,
     &                 epse, epso, plnew_a, plnow_a,
     &                 plnev_a, plnod_a)
        print *,'after treads n, P-GZ=',P_GZ,'p_Q=',P_Q,'P_TE=',P_TE,
     &  'P_DI=',P_DI,'P_ZE=',P_ZE,'P_RQ=',P_RQ
        print *,'after treads n, p_zQ=',P_ZQ,'P_y=',P_y,
     &  'P_x=',P_x,'P_w=',P_w,'P_Rt=',P_Rt
        print *,'trie_ls(3:5,1,2)=',trie_ls(3:5,1,2)
!
        if(me.eq.0) PRINT 9878, FHOUR
 9878   FORMAT(1H ,'FHOUR AFTER TREAD',F6.2)
!
!--------------------------------------------------------------
! fill up n+1 grid_gr in case of internal2export used.
!
        trie_ls(:,:,p_zq)=trie_ls(:,:,p_q )
        trie_ls(:,:,p_y :p_y +levs-1)=trie_ls(:,:,p_te:p_te+levs-1)
        trie_ls(:,:,p_x :p_x +levs-1)=trie_ls(:,:,p_di:p_di+levs-1)
        trie_ls(:,:,p_w :p_w +levs-1)=trie_ls(:,:,p_ze:p_ze+levs-1)
        trie_ls(:,:,p_rt:p_rt+levh-1)=trie_ls(:,:,p_rq:p_rq+levh-1)
        trio_ls(:,:,p_zq)=trio_ls(:,:,p_q )
        trio_ls(:,:,p_y :p_y +levs-1)=trio_ls(:,:,p_te:p_te+levs-1)
        trio_ls(:,:,p_x :p_x +levs-1)=trio_ls(:,:,p_di:p_di+levs-1)
        trio_ls(:,:,p_w :p_w +levs-1)=trio_ls(:,:,p_ze:p_ze+levs-1)
        trio_ls(:,:,p_rt:p_rt+levh-1)=trio_ls(:,:,p_rq:p_rq+levh-1)
!
!--------------------------------------------------------------
! fill up n+1 grid_gr in case of internal2export used.
!
        grid_gr(:,:,g_zq)=grid_gr(:,:,g_q )
        grid_gr(:,:,g_t :g_t +levs-1)=grid_gr(:,:,g_tt:g_tt+levs-1)
        grid_gr(:,:,g_u :g_u +levs-1)=grid_gr(:,:,g_uu:g_uu+levs-1)
        grid_gr(:,:,g_v :g_v +levs-1)=grid_gr(:,:,g_vv:g_vv+levs-1)
        grid_gr(:,:,g_rt:g_rt+levh-1)=grid_gr(:,:,g_rq:g_rq+levh-1)
     
!--------------------------------------------------------
!!
      RETURN
      END
