      module do_dynamics_mod
cc
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
      private
!
      public do_dynamics_gridc2syn
      public do_dynamics_gridt2anl
      public do_dynamics_gridn2anl
      public do_dynamics_gridm2sym
      public do_dynamics_spectupdatewrt
      public do_dynamics_spectupdatexyzq
      public do_dynamics_spectn2c
      public do_dynamics_spectn2m
      public do_dynamics_spectc2n
      public do_dynamics_syn2gridn
      public do_dynamics_gridomega
      public do_dynamics_gridfilter
      public do_dynamics_gridn2c
      public do_dynamics_gridn2m
      public do_dynamics_gridupdate
      public do_dynamics_gridpdp
      public do_dynamics_griddpm
      public do_dynamics_gridcheck

      contains


! --------------------------------------------------------------
      subroutine  do_dynamics_gridc2syn(grid_gr,syn_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
          jlonf = (lan-1)*lonf
          do k=1,levs
            do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksu-2+k)*lon_dim,lan)=grid_gr(ilan,g_uu+k-1)
             syn_gr_a_2(i+(ksv-2+k)*lon_dim,lan)=grid_gr(ilan,g_vv+k-1)
             syn_gr_a_2(i+(kst-2+k)*lon_dim,lan)=grid_gr(ilan,g_tt+k-1)
            enddo
          enddo
          do k=1,levh
            do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksr-2+k)*lon_dim,lan)=grid_gr(ilan,g_rq+k-1)
            enddo
          enddo
          do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksq-1)*lon_dim,lan)=grid_gr(ilan,g_q)
          enddo
        enddo

      return
      end subroutine do_dynamics_gridc2syn
!
! --------------------------------------------------------------
      subroutine  do_dynamics_gridn2anl(grid_gr,anl_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kau-2+k)*lon_dim,lan)=grid_gr(ilan,g_u  +k-1)
            anl_gr_a_2(i+(kav-2+k)*lon_dim,lan)=grid_gr(ilan,g_v  +k-1)
            anl_gr_a_2(i+(kat-2+k)*lon_dim,lan)=grid_gr(ilan,g_t  +k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)=grid_gr(ilan,g_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kaps-1)*lon_dim,lan)=grid_gr(ilan,g_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2anl
!
! --------------------------------------------------------------
      subroutine  do_dynamics_gridm2sym(grid_gr,sym_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) sym_gr_a_2(lonfx*lotm,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            sym_gr_a_2(i+(ksum-2+k)*lon_dim,lan)=grid_gr(ilan,g_uum+k-1)
            sym_gr_a_2(i+(ksvm-2+k)*lon_dim,lan)=grid_gr(ilan,g_vvm+k-1)
            sym_gr_a_2(i+(kstm-2+k)*lon_dim,lan)=grid_gr(ilan,g_ttm+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            sym_gr_a_2(i+(ksrm-2+k)*lon_dim,lan)=grid_gr(ilan,g_rm +k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            sym_gr_a_2(i+(kspsm-1)*lon_dim,lan)=grid_gr(ilan,g_qm)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridm2sym
!
! --------------------------------------------------------------
      subroutine do_dynamics_gridt2anl(grid_gr,anl_gr_a_2,rdt2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real	rdt2
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kau-2+k)*lon_dim,lan)=
     &     (grid_gr(ilan,G_u  +k-1)-grid_gr(ilan,G_uum+k-1))*rdt2
            anl_gr_a_2(i+(kav-2+k)*lon_dim,lan)=
     &     (grid_gr(ilan,G_v  +k-1)-grid_gr(ilan,G_vvm+k-1))*rdt2
            anl_gr_a_2(i+(kat-2+k)*lon_dim,lan)=
     &     (grid_gr(ilan,G_t  +k-1)-grid_gr(ilan,G_ttm+k-1))*rdt2
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)=
     &     (grid_gr(ilan,G_rt +k-1)-grid_gr(ilan,G_rm +k-1))*rdt2
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kaps-1)*lon_dim,lan)=
     &     (grid_gr(ilan,G_zq)-grid_gr(ilan,G_qm))*rdt2
        enddo
      enddo

      return
      end subroutine do_dynamics_gridt2anl
!
!----------------------------------------------------------
      subroutine do_dynamics_spectupdatewrt(trie_ls,trio_ls,dt2)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      real,   intent(in):: dt2

      integer	k,i

      do k=1,levs
         do i=1,len_trie_ls
            trie_ls(i,1,P_w  +k-1)=
     &      trie_ls(i,1,P_zem+k-1)+dt2*trie_ls(i,1,P_w+k-1)
            trie_ls(i,2,P_w  +k-1)=
     &      trie_ls(i,2,P_zem+k-1)+dt2*trie_ls(i,2,P_w+k-1)
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_w  +k-1)=
     &      trio_ls(i,1,P_zem+k-1)+dt2*trio_ls(i,1,P_w+k-1)
            trio_ls(i,2,P_w  +k-1)=
     &      trio_ls(i,2,P_zem+k-1)+dt2*trio_ls(i,2,P_w+k-1)
         enddo
      enddo
      do k=1,levh
         do i=1,len_trie_ls
            trie_ls(i,1,P_rt+k-1)=
     &      trie_ls(i,1,P_rm+k-1)+dt2* trie_ls(i,1,P_rt+k-1)  
            trie_ls(i,2,P_rt+k-1)=
     &      trie_ls(i,2,P_rm+k-1)+dt2* trie_ls(i,2,P_rt+k-1) 
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_rt+k-1)=
     &      trio_ls(i,1,P_rm+k-1)+dt2* trio_ls(i,1,P_rt+k-1)
            trio_ls(i,2,P_rt+k-1)=
     &      trio_ls(i,2,P_rm+k-1)+dt2* trio_ls(i,2,P_rt+k-1)
         enddo
      enddo

      return
      end subroutine do_dynamics_spectupdatewrt
!
!----------------------------------------------------------
      subroutine do_dynamics_spectupdatexyzq(trie_ls,trio_ls,dt2)
 
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      real,   intent(in):: dt2

      integer	k,i

      do k=1,levs                                                       
         do i=1,len_trie_ls                                             
            trie_ls(i,1,P_x  +k-1)=                                     
     &      trie_ls(i,1,P_dim+k-1)+dt2*trie_ls(i,1,P_x+k-1)    
            trie_ls(i,2,P_x  +k-1)=                                     
     &      trie_ls(i,2,P_dim+k-1)+dt2*trie_ls(i,2,P_x+k-1)    
            trie_ls(i,1,P_y  +k-1)=                                     
     &      trie_ls(i,1,P_tem+k-1)+dt2*trie_ls(i,1,P_y+k-1)    
            trie_ls(i,2,P_y  +k-1)=                                     
     &      trie_ls(i,2,P_tem+k-1)+dt2*trie_ls(i,2,P_y+k-1)    
         enddo                                                          
         do i=1,len_trio_ls                                             
            trio_ls(i,1,P_x  +k-1)=                                     
     &      trio_ls(i,1,P_dim+k-1)+dt2*trio_ls(i,1,P_x+k-1)    
            trio_ls(i,2,P_x  +k-1)=                                     
     &      trio_ls(i,2,P_dim+k-1)+dt2*trio_ls(i,2,P_x+k-1)    
            trio_ls(i,1,P_y  +k-1)=                                     
     &      trio_ls(i,1,P_tem+k-1)+dt2*trio_ls(i,1,P_y+k-1)    
            trio_ls(i,2,P_y  +k-1)=                                     
     &      trio_ls(i,2,P_tem+k-1)+dt2*trio_ls(i,2,P_y+k-1)    
         enddo                                                          
      enddo                                                             
!
         do i=1,len_trie_ls                                             
            trie_ls(i,1,P_zq)=                                          
     &      trie_ls(i,1,P_qm)+dt2*trie_ls(i,1,P_zq)            
            trie_ls(i,2,P_zq)=                                          
     &      trie_ls(i,2,P_qm)+dt2*trie_ls(i,2,P_zq)            
         enddo                                                          
         do i=1,len_trio_ls                                             
            trio_ls(i,1,P_zq)=                                          
     &      trio_ls(i,1,P_qm)+dt2*trio_ls(i,1,P_zq)            
            trio_ls(i,2,P_zq)=                                          
     &      trio_ls(i,2,P_qm)+dt2*trio_ls(i,2,P_zq)            
         enddo                                                          

      return 
      end subroutine do_dynamics_spectupdatexyzq

!--------------------------------------------
      subroutine do_dynamics_spectn2c(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

      DO K=1,LEVS
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_DI+K-1)=TRIE_LS(J,1,P_X+K-1)
         TRIE_LS(J,2,P_DI+K-1)=TRIE_LS(J,2,P_X+K-1)
         TRIE_LS(J,1,P_ZE+K-1)=TRIE_LS(J,1,P_W+K-1)
         TRIE_LS(J,2,P_ZE+K-1)=TRIE_LS(J,2,P_W+K-1)
         TRIE_LS(J,1,P_TE+K-1)=TRIE_LS(J,1,P_Y+K-1)
         TRIE_LS(J,2,P_TE+K-1)=TRIE_LS(J,2,P_Y+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_DI+K-1)=TRIO_LS(J,1,P_X+K-1)
         TRIO_LS(J,2,P_DI+K-1)=TRIO_LS(J,2,P_X+K-1)
         TRIO_LS(J,1,P_ZE+K-1)=TRIO_LS(J,1,P_W+K-1)
         TRIO_LS(J,2,P_ZE+K-1)=TRIO_LS(J,2,P_W+K-1)
         TRIO_LS(J,1,P_TE+K-1)=TRIO_LS(J,1,P_Y+K-1)
         TRIO_LS(J,2,P_TE+K-1)=TRIO_LS(J,2,P_Y+K-1)
       ENDDO
      ENDDO
      DO K=1,LEVH
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_RQ+K-1)=TRIE_LS(J,1,P_RT+K-1)
         TRIE_LS(J,2,P_RQ+K-1)=TRIE_LS(J,2,P_RT+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_RQ+K-1)=TRIO_LS(J,1,P_RT+K-1)
         TRIO_LS(J,2,P_RQ+K-1)=TRIO_LS(J,2,P_RT+K-1)
       ENDDO
      ENDDO
      DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_Q)=TRIE_LS(J,1,P_ZQ)
         TRIE_LS(J,2,P_Q)=TRIE_LS(J,2,P_ZQ)
      ENDDO
      DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_Q)=TRIO_LS(J,1,P_ZQ)
         TRIO_LS(J,2,P_Q)=TRIO_LS(J,2,P_ZQ)
      ENDDO

      return
      end subroutine do_dynamics_spectn2c

!--------------------------------------------
      subroutine do_dynamics_spectn2m(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

      DO K=1,LEVS
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_DIM+K-1)=TRIE_LS(J,1,P_X+K-1)
         TRIE_LS(J,2,P_DIM+K-1)=TRIE_LS(J,2,P_X+K-1)
         TRIE_LS(J,1,P_ZEM+K-1)=TRIE_LS(J,1,P_W+K-1)
         TRIE_LS(J,2,P_ZEM+K-1)=TRIE_LS(J,2,P_W+K-1)
         TRIE_LS(J,1,P_TEM+K-1)=TRIE_LS(J,1,P_Y+K-1)
         TRIE_LS(J,2,P_TEM+K-1)=TRIE_LS(J,2,P_Y+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_DIM+K-1)=TRIO_LS(J,1,P_X+K-1)
         TRIO_LS(J,2,P_DIM+K-1)=TRIO_LS(J,2,P_X+K-1)
         TRIO_LS(J,1,P_ZEM+K-1)=TRIO_LS(J,1,P_W+K-1)
         TRIO_LS(J,2,P_ZEM+K-1)=TRIO_LS(J,2,P_W+K-1)
         TRIO_LS(J,1,P_TEM+K-1)=TRIO_LS(J,1,P_Y+K-1)
         TRIO_LS(J,2,P_TEM+K-1)=TRIO_LS(J,2,P_Y+K-1)
       ENDDO
      ENDDO
      DO K=1,LEVH
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_RM+K-1)=TRIE_LS(J,1,P_RT+K-1)
         TRIE_LS(J,2,P_RM+K-1)=TRIE_LS(J,2,P_RT+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_RM+K-1)=TRIO_LS(J,1,P_RT+K-1)
         TRIO_LS(J,2,P_RM+K-1)=TRIO_LS(J,2,P_RT+K-1)
       ENDDO
      ENDDO
      DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_QM)=TRIE_LS(J,1,P_ZQ)
         TRIE_LS(J,2,P_QM)=TRIE_LS(J,2,P_ZQ)
      ENDDO
      DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_QM)=TRIO_LS(J,1,P_ZQ)
         TRIO_LS(J,2,P_QM)=TRIO_LS(J,2,P_ZQ)
      ENDDO

      return
      end subroutine do_dynamics_spectn2m

!--------------------------------------------
      subroutine do_dynamics_spectc2n(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

      DO K=1,LEVS
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_X+K-1)=TRIE_LS(J,1,P_DI+K-1)
         TRIE_LS(J,2,P_X+K-1)=TRIE_LS(J,2,P_DI+K-1)
         TRIE_LS(J,1,P_W+K-1)=TRIE_LS(J,1,P_ZE+K-1)
         TRIE_LS(J,2,P_W+K-1)=TRIE_LS(J,2,P_ZE+K-1)
         TRIE_LS(J,1,P_Y+K-1)=TRIE_LS(J,1,P_TE+K-1)
         TRIE_LS(J,2,P_Y+K-1)=TRIE_LS(J,2,P_TE+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_X+K-1)=TRIO_LS(J,1,P_DI+K-1)
         TRIO_LS(J,2,P_X+K-1)=TRIO_LS(J,2,P_DI+K-1)
         TRIO_LS(J,1,P_W+K-1)=TRIO_LS(J,1,P_ZE+K-1)
         TRIO_LS(J,2,P_W+K-1)=TRIO_LS(J,2,P_ZE+K-1)
         TRIO_LS(J,1,P_Y+K-1)=TRIO_LS(J,1,P_TE+K-1)
         TRIO_LS(J,2,P_Y+K-1)=TRIO_LS(J,2,P_TE+K-1)
       ENDDO
      ENDDO
      DO K=1,LEVH
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_RT+K-1)=TRIE_LS(J,1,P_RQ+K-1)
         TRIE_LS(J,2,P_RT+K-1)=TRIE_LS(J,2,P_RQ+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_RT+K-1)=TRIO_LS(J,1,P_RQ+K-1)
         TRIO_LS(J,2,P_RT+K-1)=TRIO_LS(J,2,P_RQ+K-1)
       ENDDO
      ENDDO
      DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_ZQ)=TRIE_LS(J,1,P_Q)
         TRIE_LS(J,2,P_ZQ)=TRIE_LS(J,2,P_Q)
      ENDDO
      DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_ZQ)=TRIO_LS(J,1,P_Q)
         TRIO_LS(J,2,P_ZQ)=TRIO_LS(J,2,P_Q)
      ENDDO

      return
      end subroutine do_dynamics_spectc2n

!--------------------------------------------
      subroutine do_dynamics_syn2gridn(syn_gr_a_2,grid_gr,
     &                                 global_lats_a,lonsperlat,nislfv)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      integer,intent(in):: nislfv

      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_u+k-1)= syn_gr_a_2(i+(ksu-2+k)*lon_dim,lan)
            grid_gr(ilan,G_v+k-1)= syn_gr_a_2(i+(ksv-2+k)*lon_dim,lan)
            grid_gr(ilan,G_t+k-1)= syn_gr_a_2(i+(kst-2+k)*lon_dim,lan)
          enddo
        enddo
!hmhj test
        if( nislfv.le.1 ) then
! ---------------------
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rt+k-1)= syn_gr_a_2(i+(ksr-2+k)*lon_dim,lan)
          enddo
        enddo
        endif
! --------------------
        do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_zq)= syn_gr_a_2(i+(ksq-1)*lon_dim,lan)
        enddo
      enddo

      return
      end subroutine do_dynamics_syn2gridn

!--------------------------------------------
      subroutine do_dynamics_gridomega(syn_gr_a_2,dyn_gr_a_2,
     &                                 grid_gr,rcs2,
     &                                 global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      real(kind=kind_evod) dyn_gr_a_2(lonfx*lotd,lats_dim_ext)
      real(kind=kind_grid) rcs2(latg2)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real(kind=kind_grid)  ugr (lonf,levs), vgr (lonf,levs)
      real(kind=kind_grid)  gtv (lonf,levs), gd  (lonf,levs)
      real(kind=kind_grid)  gtvx(lonf,levs), gtvy(lonf,levs)
      real(kind=kind_grid)  gphi(lonf)     , glam(lonf)     , gq  (lonf)
      real(kind=kind_grid)  vvel(lonf,levs)

      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf

        do k=1,levs
          do i=1,lons_lat
            ugr (i,k)= syn_gr_a_2(i+(ksu-2+k)*lon_dim,lan)
            vgr (i,k)= syn_gr_a_2(i+(ksv-2+k)*lon_dim,lan)
            gd  (i,k)= syn_gr_a_2(i+(ksd-2+k)*lon_dim,lan)
            gtv (i,k)= syn_gr_a_2(i+(kst-2+k)*lon_dim,lan)
            gtvx(i,k)= dyn_gr_a_2(i+(kdtlam-2+k)*lon_dim,lan)
            gtvy(i,k)= dyn_gr_a_2(i+(kdtphi-2+k)*lon_dim,lan)
          enddo
        enddo
        do i=1,lons_lat
          gq  (i)=syn_gr_a_2(i+(ksq   -1)*lon_dim,lan)
          gphi(i)=syn_gr_a_2(i+(kspphi-1)*lon_dim,lan)
          glam(i)=syn_gr_a_2(i+(ksplam-1)*lon_dim,lan)
        enddo

        if( gen_coord_hybrid ) then 
          call omega_gch(lons_lat,lonf,rcs2(min(lat,latg-lat+1)),
     &                   gq,gphi,glam,gtv,gtvx,gtvy,gd,ugr,vgr,vvel)
        else if( hybrid )then 
          call omega_hyb(lons_lat,lonf,rcs2(min(lat,latg-lat+1)),
     &                  gq,gphi,glam,gd,ugr,vgr,vvel)
        else
          call omega_sig(lons_lat,lonf,rcs2(min(lat,latg-lat+1)),
     &                  gq,gphi,glam,gd,ugr,vgr,vvel)
        endif

        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,g_dpdt+k-1)=vvel(i,k)
          enddo
        enddo

      enddo

      return
      end subroutine do_dynamics_gridomega

!-------------------------------------------------------------------
      subroutine do_dynamics_gridfilter(grid_gr,filta,filtb,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: filta, filtb

      integer 	lan,lat,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uum+k-1)=grid_gr(ilan,G_uu +k-1) *filta+
     &     (grid_gr(ilan,G_uum+k-1)+grid_gr(ilan,G_u  +k-1))*filtb
            grid_gr(ilan,G_vvm+k-1)=grid_gr(ilan,G_vv +k-1) *filta+
     &     (grid_gr(ilan,G_vvm+k-1)+grid_gr(ilan,G_v  +k-1))*filtb
            grid_gr(ilan,G_ttm+k-1)=grid_gr(ilan,G_tt +k-1) *filta+
     &     (grid_gr(ilan,G_ttm+k-1)+grid_gr(ilan,G_t  +k-1))*filtb
            grid_gr(ilan,G_uu +k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vv +k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_tt +k-1)=grid_gr(ilan,G_t  +k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rm +k-1)=grid_gr(ilan,G_rq +k-1) *filta+
     &     (grid_gr(ilan,G_rm +k-1)+grid_gr(ilan,G_rt +k-1))*filtb
            grid_gr(ilan,G_rq +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_qm)=grid_gr(ilan,G_q )
          grid_gr(ilan,G_q )=grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridfilter

!-------------------------------------------------------------------
      subroutine do_dynamics_gridn2c(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer	lan,lat,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uu +k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vv +k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_tt +k-1)=grid_gr(ilan,G_t  +k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rq +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_q )=grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2c

!-------------------------------------------------------------------
      subroutine do_dynamics_gridn2m(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer	lan,lat,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uum+k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vvm+k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_ttm+k-1)=grid_gr(ilan,G_t  +k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rm +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_qm)=grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2m

! -------------------------------------------------------------------
      subroutine do_dynamics_gridupdate(grid_gr,anl_gr_a_2,dt2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: dt2

      integer 	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_u  +k-1)=grid_gr(ilan,G_uum+k-1)+
     &      anl_gr_a_2(i+(kau-2+k)*lon_dim,lan)*dt2
            grid_gr(ilan,G_v  +k-1)=grid_gr(ilan,G_vvm+k-1)+
     &      anl_gr_a_2(i+(kav-2+k)*lon_dim,lan)*dt2
            grid_gr(ilan,G_t  +k-1)=grid_gr(ilan,G_ttm+k-1)+
     &      anl_gr_a_2(i+(kat-2+k)*lon_dim,lan)*dt2
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rt +k-1)=grid_gr(ilan,G_rm +k-1)+
     &      anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)*dt2
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_zq)=grid_gr(ilan,G_qm)+
     &      anl_gr_a_2(i+(kaps-1)*lon_dim,lan)*dt2
        enddo

      enddo
 
      return
      end subroutine do_dynamics_gridupdate
!
! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridpdp(grid_gr,
     &                                  global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real(kind=kind_grid)  gtv (lonf,levs)
      real(kind=kind_grid)  gq  (lonf)
      real(kind=kind_grid)  prsl(lonf,levs), dprs(lonf,levs)

      integer 	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            gtv(i,k) = grid_gr(ilan,G_t  +k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            gq(i) = grid_gr(ilan,G_zq)
        enddo

        if( gen_coord_hybrid ) then 
          call gch2press(lons_lat,lonf,gq, gtv, prsl, dprs)
        else if( hybrid )then 
          call hyb2press(lons_lat,lonf,gq, prsl, dprs)
        else
          call sig2press(lons_lat,lonf,gq, prsl, dprs)
        endif

        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,g_p   +k-1)=prsl(i,k)
            grid_gr(ilan,g_dp  +k-1)=dprs(i,k)
          enddo
        enddo

      enddo
 
      return
      end subroutine do_dynamics_gridpdp
!
! -------------------------------------------------------------------
!
      subroutine do_dynamics_griddpm(grid_gr,
     &                                  global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real(kind=kind_grid)  gtv (lonf,levs)
      real(kind=kind_grid)  gq  (lonf)
      real(kind=kind_grid)  prsl(lonf,levs), dprs(lonf,levs)

      integer 	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            gtv(i,k) = grid_gr(ilan,G_ttm+k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            gq(i) = grid_gr(ilan,G_qm)
        enddo

        if( gen_coord_hybrid ) then 
          call gch2press(lons_lat,lonf,gq, gtv, prsl, dprs)
        else if( hybrid )then 
          call hyb2press(lons_lat,lonf,gq, prsl, dprs)
        else
          call sig2press(lons_lat,lonf,gq, prsl, dprs)
        endif

        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,g_dp  +k-1)=dprs(i,k)
          enddo
        enddo

      enddo
 
      return
      end subroutine do_dynamics_griddpm
!
! -------------------------------------------------------------------
      subroutine do_dynamics_gridcheck(grid_gr,
     &                                 global_lats_a,lonsperlat,chr)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      character*(*) chr

      integer 	lan,lat,lons_lat,k

      print *,' check: g_ttm g_tt g_t ',g_ttm,g_tt,g_t
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        print *,' gridcheck: lan lat lons_lat ',lan,lat,lons_lat
        do k=1,levs
          print *,' check grid of ttm tt t at k=',k
          call mymaxmin(grid_gr(1,g_ttm+k-1),lons_lat,lonf,1,chr)
          call mymaxmin(grid_gr(1,g_tt +k-1),lons_lat,lonf,1,chr)
          call mymaxmin(grid_gr(1,g_t  +k-1),lons_lat,lonf,1,chr)
        enddo
      enddo
 
      return
      end subroutine do_dynamics_gridcheck
!
      end module do_dynamics_mod
