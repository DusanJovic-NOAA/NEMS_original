      module layout1
      implicit none
      save
cc
      integer           nodes, nodes_comp,nodes_io,
     x                  me,
     x                  ls_dim,
     x                  ls_max_node,
     x                  lats_dim_r,
     x                  lats_dim_ext,
     x                  lats_node_r,
     x                  lats_node_r_max,
     x                  lats_node_ext,
     x                  ipt_lats_node_r,
     x                  ipt_lats_node_ext,
     x                  me_l_0
cc
      integer ,allocatable :: lon_dims_r(:),lon_dims_ext(:)
      end module layout1
