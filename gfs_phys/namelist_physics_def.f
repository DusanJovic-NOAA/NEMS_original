      module namelist_physics_def
      use machine
      implicit none
      save
      integer nszer,nsres,nslwr,nsout,nsswr,nscyc,igen,jo3,ngptc
     &,       lsm,ens_mem,ncw(2)
      real(kind=kind_evod) fhswr,fhlwr,fhrot,fhseg,fhmax,fhout,fhres,
     & fhzer,fhini,fhcyc,crtrh(3),flgmin(2),ccwf
      logical ldiag3d,ras,zhao_mic,sashal,newsas
      logical lsfwd,lssav,lscca,lsswr,lslwr
      logical shuff_lats_r,reshuff_lats_r
      logical hybrid,gen_coord_hybrid,zflxtvd
      logical pre_rad,random_xkt2,old_monin,cnvgwd 
      logical restart, gfsio_in, gfsio_out
      character*20 ens_nam
!
!     Radiation control parameters
!
      integer isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw
      end module namelist_physics_def
