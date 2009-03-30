      module gfs_dyn_layout1
      implicit none
      
cc
      integer           nodes, nodes_comp,nodes_io,
     x                  me,
     x                  ls_dim,
     x                  ls_max_node,
     x                  lats_dim_a,
     x                  lats_dim_ext,
     x                  lats_node_a,
     x                  lats_node_a_max,
     x                  lats_node_ext,
     x                  ipt_lats_node_a,
     x                  ipt_lats_node_ext,
     x                  len_trie_ls,
     x                  len_trio_ls,
     x                  len_trie_ls_max,
     x                  len_trio_ls_max,
     x                  me_l_0
cc
      INTEGER ,ALLOCATABLE :: lat1s_a(:),
     .  lon_dims_a(:),lon_dims_ext(:)

!hmhj nislfv
      integer   lonfull,lonhalf,lonpart,mylonlen
      integer   latfull,lathalf,latpart,mylatlen

      integer, allocatable :: lonstr(:),lonlen(:)
      integer, allocatable :: latstr(:),latlen(:)
      real, allocatable :: cosglat(:)
      real, allocatable :: gglat(:),gglati(:)
      real, allocatable :: gglon(:),ggloni(:)

      end module gfs_dyn_layout1
