      SUBROUTINE input_fields(n1,n2, PDRYINI,TRIE_LS,TRIO_LS,grid_gr,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,SNNP1EV,SNNP1OD,
     &                 global_lats_a,nblck,lonsperlat,
     &                 epse,epso,plnev_a,plnod_a,plnew_a,plnow_a,
     &                 lats_nodes_a, cread, cread2,pwat,ptot)
!!
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfsio_module
      use gfsio_def
      use gfs_dyn_mpi_def
      IMPLICIT NONE
!!
 
cmy fix pdryini type
cmy      REAL(KIND=KIND_EVOD) PDRYINI
      REAL(KIND=kind_grid) PDRYINI
      INTEGER              N1,N2
      CHARACTER (len=*)   :: CREAD, CREAD2
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTLS)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTLS)
      REAL(KIND=KIND_GRID) GRID_GR(lonf*lats_node_a_max,lotgr)
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

      REAL(KIND=KIND_GRID) pwat   (lonf,lats_node_a)
      REAL(KIND=KIND_GRID) ptot   (lonf,lats_node_a)
!
      integer nblck
      integer global_lats_a(latg), lonsperlat(latg)
 
cmy bug fix on dimension of ls_node
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
      integer		lan, lat, lons_lat, jlonf
!
      if(me.eq.0) PRINT  9876,N1,N2,FHOUR,idate
 9876 FORMAT(1H ,'N1,N2,FHOUR IN input_fields ',2(I4,1X),F6.2,
     & ' idate no yet read in',4(1x,i4))
      IPRINT = 0
c$$$  IF ( ME .EQ. 0 ) IPRINT = 1
!
      if (me .eq. 0) write(0,*)' cread=',cread,'ntoz=',ntoz
        CALL TREADEO_gfsio(FHOUR,IDATE,
     X               TRIE_LS(1,1,P_GZ ), TRIE_LS(1,1,P_QM ),
     X               TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_DIM),
     X               TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_RM ),
     X               TRIO_LS(1,1,P_GZ ), TRIO_LS(1,1,P_QM ),
     X               TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_DIM),
     X               TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_RM ),
     &               zsg, psg, ttg, uug, vvg, rqg,
     X               LS_NODE,LS_NODES,MAX_LS_NODES,
     X               SNNP1EV,SNNP1OD,PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat, cread,
     &               epse, epso, plnew_a, plnow_a, 
     &               plnev_a, plnod_a, pwat, ptot)

      do j=1,lats_node_a
        jlonf=(j-1)*lonf
        grid_gr(jlonf+1:jlonf+lonf,g_gz) = zsg(1:lonf,j)
        grid_gr(jlonf+1:jlonf+lonf,g_qm) = psg(1:lonf,j)
      enddo
      do k=1,levs
        do j=1,lats_node_a
          jlonf=(j-1)*lonf
          grid_gr(jlonf+1:jlonf+lonf,g_ttm+k-1) = ttg(1:lonf,j,k)
          grid_gr(jlonf+1:jlonf+lonf,g_uum+k-1) = uug(1:lonf,j,k)
          grid_gr(jlonf+1:jlonf+lonf,g_vvm+k-1) = vvg(1:lonf,j,k)
        enddo
      enddo
      do k=1,levh
        do j=1,lats_node_a
          jlonf=(j-1)*lonf
          grid_gr(jlonf+1:jlonf+lonf,g_rm +k-1) = rqg(1:lonf,j,k)
        enddo
      enddo

 
      fhini=fhour
      if(me.eq.0) PRINT 9877, N1,FHOUR
 9877 FORMAT(1H ,'N1,FHOUR AFTER TREAD',1(I4,1X),F6.2)
 
      if (me .eq. 0) write(0,*)' fhini=',fhini,'last_fcst_pe=',
     &     last_fcst_pe,'fhrot=',fhrot
!jw      if (.NOT.LIOPE.or.icolor.ne.2) then
      if (me<=last_fcst_pe) then 
!sela   print*,'liope=',liope,' icolor=',icolor
        CALL RMS_spect(TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_DIM),
     X             TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_ZEM),
     X             TRIE_LS(1,1,P_RM ),
     X             TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_DIM),
     X             TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_ZEM),
     X             TRIO_LS(1,1,P_RM ),
     X             LS_NODES,MAX_LS_NODES)
      endif
!---------------------------------------------------------------
      if(fhini.eq.fhrot) THEN
!set n time level values to n-1 time
! spectral
!       print *,' set time level n to time level n-1 '
        do i=1,len_trie_ls
           trie_ls(i,1,P_q )=trie_ls(i,1,P_qm )
           trie_ls(i,2,P_q )=trie_ls(i,2,P_qm )
        enddo
        do i=1,len_trio_ls
           trio_ls(i,1,P_q )=trio_ls(i,1,P_qm )
           trio_ls(i,2,P_q )=trio_ls(i,2,P_qm )
        enddo
 
        do k=1,levs
          do i=1,len_trie_ls
            trie_ls(i,1,P_te +k-1)=trie_ls(i,1,P_tem +k-1)
            trie_ls(i,2,P_te +k-1)=trie_ls(i,2,P_tem +k-1)
 
            trie_ls(i,1,P_di +k-1)=trie_ls(i,1,P_dim +k-1)
            trie_ls(i,2,P_di +k-1)=trie_ls(i,2,P_dim +k-1)
 
            trie_ls(i,1,P_ze +k-1)=trie_ls(i,1,P_zem +k-1)
            trie_ls(i,2,P_ze +k-1)=trie_ls(i,2,P_zem +k-1)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,P_te +k-1)=trio_ls(i,1,P_tem+k-1)
            trio_ls(i,2,P_te +k-1)=trio_ls(i,2,P_tem+k-1)
 
            trio_ls(i,1,P_di +k-1)=trio_ls(i,1,P_dim+k-1)
            trio_ls(i,2,P_di +k-1)=trio_ls(i,2,P_dim+k-1)
 
            trio_ls(i,1,P_ze +k-1)=trio_ls(i,1,P_zem+k-1)
            trio_ls(i,2,P_ze +k-1)=trio_ls(i,2,P_zem+k-1)
          enddo
        enddo
 
        do k=1,levh
          do i=1,len_trie_ls
            trie_ls(i,1,P_rq +k-1)=trie_ls(i,1,P_rm +k-1)
            trie_ls(i,2,P_rq +k-1)=trie_ls(i,2,P_rm +k-1)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,P_rq +k-1)=trio_ls(i,1,P_rm+k-1)
            trio_ls(i,2,P_rq +k-1)=trio_ls(i,2,P_rm+k-1)
          enddo
        enddo
! grid
        grid_gr(:,g_q )=grid_gr(:,g_qm )
        grid_gr(:,g_tt:g_tt+levs-1)=grid_gr(:,g_ttm:g_ttm+levs-1)
        grid_gr(:,g_uu:g_uu+levs-1)=grid_gr(:,g_uum:g_uum+levs-1)
        grid_gr(:,g_vv:g_vv+levs-1)=grid_gr(:,g_vvm:g_vvm+levs-1)
        grid_gr(:,g_rq:g_rq+levh-1)=grid_gr(:,g_rm :g_rm +levh-1)

!--------------------------------------------------------
      else
!--------------------------------------------------------
        IPRINT = 0
c$$$      IF ( ME .EQ. 0 ) IPRINT = 1
      if (me .eq. 0) write(0,*)' cread2=',cread2
          CALL TREADEO_gfsio(FHOUR,IDATE,
     X                 TRIE_LS(1,1,P_GZ), TRIE_LS(1,1,P_Q ),
     X                 TRIE_LS(1,1,P_TE), TRIE_LS(1,1,P_DI),
     X                 TRIE_LS(1,1,P_ZE), TRIE_LS(1,1,P_RQ),
     X                 TRIO_LS(1,1,P_GZ), TRIO_LS(1,1,P_Q ),
     X                 TRIO_LS(1,1,P_TE), TRIO_LS(1,1,P_DI),
     X                 TRIO_LS(1,1,P_ZE), TRIO_LS(1,1,P_RQ),
     &                 zsg, psg, ttg, uug, vvg, rqg,
     X                 LS_NODE,LS_NODES,MAX_LS_NODES,
     X                 SNNP1EV,SNNP1OD,PDRYINI,IPRINT,
     &                 global_lats_a,lats_nodes_a,lonsperlat, cread2,
     &                 epse, epso, plnew_a, plnow_a, 
     &                 plnev_a, plnod_a, pwat, ptot)

      do j=1,lats_node_a
        jlonf=(j-1)*lonf
!       grid_gr(jlonf+1:jlonf+lonf,g_gz) = zsg(1:lonf,j)
        grid_gr(jlonf+1:jlonf+lonf,g_q ) = psg(1:lonf,j)
      enddo
      do k=1,levs
        do j=1,lats_node_a
          jlonf=(j-1)*lonf
          grid_gr(jlonf+1:jlonf+lonf,g_tt +k-1) = ttg(1:lonf,j,k)
          grid_gr(jlonf+1:jlonf+lonf,g_uu +k-1) = uug(1:lonf,j,k)
          grid_gr(jlonf+1:jlonf+lonf,g_vv +k-1) = vvg(1:lonf,j,k)
        enddo
      enddo
      do k=1,levh
        do j=1,lats_node_a
          jlonf=(j-1)*lonf
          grid_gr(jlonf+1:jlonf+lonf,g_rq +k-1) = rqg(1:lonf,j,k)
        enddo
      enddo

        if(me.eq.0) PRINT 9878, N2,FHOUR
 9878   FORMAT(1H ,'N2,FHOUR AFTER TREAD',1(I4,1X),F6.2)
      endif
!
! fill up n+1 grid_gr in case of internal2export used.
!
        grid_gr(:,g_zq)=grid_gr(:,g_q )
        grid_gr(:,g_t :g_t +levs-1)=grid_gr(:,g_tt:g_tt+levs-1)
        grid_gr(:,g_u :g_u +levs-1)=grid_gr(:,g_uu:g_uu+levs-1)
        grid_gr(:,g_v :g_v +levs-1)=grid_gr(:,g_vv:g_vv+levs-1)
        grid_gr(:,g_rt:g_rt+levh-1)=grid_gr(:,g_rq:g_rq+levh-1)
     
!--------------------------------------------------------
!!
      RETURN
      END
