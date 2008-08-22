      module gfs_dyn_vert_def
      use gfs_dyn_machine
      implicit none

      REAL(KIND=KIND_EVOD) ,ALLOCATABLE :: am(:,:),bm(:,:),cm(:,:),
     . dm(:,:,:),tor(:), si(:),sl(:),del(:),rdel2(:),ci(:),
     . cl(:),tov(:),sv(:)
      real(kind=kind_evod), allocatable :: slk(:), sik(:)
      end module gfs_dyn_vert_def
