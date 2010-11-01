      module namelist_dynamics_def
!
!---- revision history
!--- Sep 8 2010    J. Wang   change gfsio to nemsio
!
!
      use gfs_dyn_machine
      implicit none
      
      integer nsres,nsout,igen,ngptc
      real(kind=kind_evod) fhrot,fhmax,fhout,fhres,fhini,fhdfi
      real(kind=kind_evod) filta,ref_temp
      logical lsfwd
      logical shuff_lats_a,reshuff_lats_a
!jw      logical hybrid,gen_coord_hybrid,zflxtvd,explicit
      logical,target :: hybrid,gen_coord_hybrid
      logical zflxtvd,explicit


      logical nemsio_in, nemsio_out
      logical reduced_grid
      integer nislfv
      character*20 ens_nam
!
      end module namelist_dynamics_def
