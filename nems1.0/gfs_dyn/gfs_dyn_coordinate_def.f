      module gfs_dyn_coordinate_def
      use gfs_dyn_machine
      implicit none
      
       real(kind=kind_evod) , allocatable ::
     . AK5(:),BK5(:),CK5(:),CK(:),DBK(:),bkl(:),   		! hmhj
     . AMHYB(:,:),BMHYB(:,:),SVHYB(:),tor_hyb(:),
     . D_HYB_m(:,:,:),THREF(:),dm205_hyb(:,:,:)			! hmhj
       real(kind=kind_evod) vertcoord_id,eps_si			! hmhj

!
      real(kind=kind_evod) , allocatable :: vcoord(:,:)
      integer nvcoord, idsl, idvc, idvm
      end module gfs_dyn_coordinate_def
