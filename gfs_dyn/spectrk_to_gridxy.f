      subroutine spectrk_to_gridxy
     x    (trie_ls,trio_ls,
     &     pyn_gr_a_1,pyn_gr_a_2,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,lonsperlat,
     x     epse,epso,epsedn,epsodn,
     x     snnp1ev,snnp1od,plnev_a,plnod_a)
!
! H.-M. H. Juang:  compute dpx dpy phix and phiy
c
#include "f_hpm.h"
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      implicit none
!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
!
      integer              ls_node(ls_dim,3)
!
      integer              ls_nodes(ls_dim,nodes)
!
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
!
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
!
      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
!
      real(kind=kind_evod) pyn_gr_a_1(lonfx*lotk,lats_dim_a)
      real(kind=kind_evod) pyn_gr_a_2(lonfx*lotk,lats_dim_a)
!
      integer              i,j,k,kx,ky,l
      integer              lan,lat,lotx
      integer              lon_dim,lons_lat,n,node
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      real(kind=kind_evod) cons0,cons2     !constant
!
      logical lslag
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
      lslag=.false.
      lotx = levs * 2
! ................................................................
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epse,epso,ls_node)
!$omp+private(k)
      do k=1,levs

      call delnpe(trie_ls(1,1,P_rk   +k-1),
     x            trio_ls(1,1,P_rkphi+k-1),
     x            trie_ls(1,1,P_rklam+k-1),
     x            epse,epso,ls_node)
!
      call delnpo(trio_ls(1,1,P_rk   +k-1),
     x            trie_ls(1,1,P_rkphi+k-1),
     x            trio_ls(1,1,P_rklam+k-1),
     x            epse,epso,ls_node)

      enddo
!!
      dimg=0
      CALL countperf(0,1,0.)
      call f_hpmstart(8,"ga sumflna")
!
      call sumflna(trie_ls(1,1,P_rklam),
     x            trio_ls(1,1,P_rklam),
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            lotx,ls_node,latg2,
     x            lslag,lats_dim_a,lotk,
     x            pyn_gr_a_1,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_a,global_lats_a,
     x            lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
     x            lonsperlat,lonfx,latg)
!
!11111111111111111111111111111111111111111111111111111111111111111111
      do lan=1,lats_node_a  
 
         lon_dim = lon_dims_a(lan)
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         CALL FOUR2GRID_thread(pyn_gr_a_1(1,lan),pyn_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx,lan,me)

      enddo 
!
      return
      end
