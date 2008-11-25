      subroutine gridrk_to_spect
     &    (trie_ls,trio_ls,grid_gr,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epse,epso,plnew_a,plnow_a)
!!
!! hmhj - this routine do spectral to grid transform in model partial reduced grid
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      implicit none
!!
!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!!
      real(kind=kind_evod) anl_gr_a_1(lonfx*levs,lats_dim_a)
      real(kind=kind_evod) anl_gr_a_2(lonfx*levs,lats_dim_a)
!
      integer              ls_node(ls_dim,3)
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
!
      real(kind=kind_evod)  epse(len_trie_ls)
      real(kind=kind_evod)  epso(len_trio_ls)
!
      real(kind=kind_evod)  plnew_a(len_trie_ls,latg2)
      real(kind=kind_evod)  plnow_a(len_trio_ls,latg2)
!
      integer              i,j,k,n
      integer              l,lan,lat,ilan,jlonf
      integer              lon_dim,lons_lat
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!
      logical 	lslag
!
      include 'function2'
!
!--------------------------------------------------------------------
!
      lslag   = .false.
      
      do k=1,levs
        trie_ls(:,:,p_rk+k-1) = 0.0
        trio_ls(:,:,p_rk+k-1) = 0.0
      enddo
!
!--------------------------------------------------------------------
      do lan=1,lats_node_a
        lon_dim = lon_dims_a(lan)
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
             ilan=i+jlonf
             anl_gr_a_2(i+(k-1)*lon_dim,lan)=grid_gr(ilan,g_rk+k-1)
          enddo
        enddo
      enddo
!
      do lan=1,lats_node_a
!
         lon_dim = lon_dims_a(lan)
!
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         call grid2four_thread(anl_gr_a_2(1,lan),anl_gr_a_1(1,lan),
     &                  lon_dim,lons_lat,lonfx,levs)
!
      enddo
!
      dimg=0
      call four2fln(lslag,lats_dim_a,levs,levs,anl_gr_a_1,
     x              ls_nodes,max_ls_nodes,
     x              lats_nodes_a,global_lats_a,lon_dims_a,
     x              lats_node_a,ipt_lats_node_a,dimg,
     x              lat1s_a,lonfx,latg,latg2,
     x              trie_ls(1,1,p_rk), trio_ls(1,1,p_rk),
     x              plnew_a, plnow_a,
     x              ls_node,0)
!
!
!     print *,' exit gridrk_to_spect '
!!
      return
      end
