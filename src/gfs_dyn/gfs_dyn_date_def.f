      module gfs_dyn_date_def
      use gfs_dyn_machine
      implicit none
      
!jw      integer idate(4)
      integer,target :: idate(4)                                           !jwang
!jw      real(kind=kind_evod) fhour,shour,thour,z00
      real(kind=kind_evod)       shour,thour,z00                           !jwang
      real(kind=kind_evod),target :: fhour                                 !jwang

      REAL(KIND=KIND_EVOD) ,ALLOCATABLE :: spdmax(:)

      end module gfs_dyn_date_def
