      module namelist_dynamics_def

      use gfs_dyn_machine
      implicit none
      save
      integer nsres,nsout,igen,ngptc
      real(kind=kind_evod) fhrot,fhmax,fhout,fhres,fhini
      real(kind=kind_evod) filta,ref_temp
      logical lsfwd
      logical shuff_lats_a,reshuff_lats_a
      logical hybrid,gen_coord_hybrid,zflxtvd,explicit
      logical gfsio_in, gfsio_out
      character*20 ens_nam
!
      end module namelist_dynamics_def
