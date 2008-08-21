      module gfs_dyn_deldifs_def
      use gfs_dyn_machine
      implicit none
      save
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNE(:),DNO(:),
     . SF(:),RTRD(:),RTHK(:),BKLY(:),CKLY(:)			! hmhj
      end module gfs_dyn_deldifs_def
