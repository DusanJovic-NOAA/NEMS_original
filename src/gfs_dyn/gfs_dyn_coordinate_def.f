      module gfs_dyn_coordinate_def
      use gfs_dyn_machine
      implicit none
      
       real(kind=kind_evod) , allocatable ::
!jw     . AK5(:),BK5(:),CK5(:),CK(:),DBK(:),bkl(:),   		! hmhj
     .                      CK(:),DBK(:),bkl(:),   		! hmhj
     . AMHYB(:,:),BMHYB(:,:),SVHYB(:),tor_hyb(:),
     . D_HYB_m(:,:,:),THREF(:),dm205_hyb(:,:,:)			! hmhj
!jw       real(kind=kind_evod) vertcoord_id,eps_si			! hmhj
!jws
       real(kind=kind_evod) eps_si                              ! hmhj
       integer(kind=kind_io4),target :: vertcoord_id            ! hmhj
       real(kind=kind_evod),allocatable,target :: AK5(:),BK5(:),CK5(:)
       integer,target :: idsl, idvc, idvm
!jwe

!
      real(kind=kind_evod) , allocatable :: vcoord(:,:)
      integer nvcoord
      end module gfs_dyn_coordinate_def
