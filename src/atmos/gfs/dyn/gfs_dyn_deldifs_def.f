      module gfs_dyn_deldifs_def
      use gfs_dyn_machine
      implicit none
      
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNE(:),DNO(:),
     . SF(:),RTRD(:),RTHK(:),BKLY(:),CKLY(:)			! hmhj
!hmhj idea Apr 06 2012
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNcE(:),DNcO(:)
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNvE(:),DNvO(:)
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNdE(:),DNdO(:)
      end module gfs_dyn_deldifs_def
