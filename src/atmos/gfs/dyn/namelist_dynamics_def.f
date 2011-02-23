      module namelist_dynamics_def

!
! program lot
! 20 Feb 2011:    Henry Juang add some options for NDSL
!
      use gfs_dyn_machine
      implicit none
      
      integer nsres,nsout,igen,ngptc,num_reduce
      real(kind=kind_evod) fhrot,fhmax,fhout,fhres,fhini,fhdfi
      real(kind=kind_evod) filta,ref_temp
      logical lsfwd
      logical shuff_lats_a,reshuff_lats_a
!jw      logical hybrid,gen_coord_hybrid,zflxtvd,explicit
      logical,target :: hybrid,gen_coord_hybrid
      logical zflxtvd,explicit


      logical nemsio_in, nemsio_out
      logical reduced_grid, semi_implicit_temp_profile
      logical mass_dp, process_split, linear_grid
      logical ndslfv
      character*20 ens_nam
!
      end module namelist_dynamics_def
