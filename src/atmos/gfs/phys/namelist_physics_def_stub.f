      module namelist_physics_def
      implicit none
      
      integer lsm
      real(kind=8),public:: fhswr,fhcyc
      logical,public:: ldiag3d,ras,sashal
      logical,public:: lssav,lscca,lsswr,lslwr,ldfi
      logical,public:: pre_rad 
!
!     Radiation control parameters
!
      integer ,public::isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw
      end module namelist_physics_def
