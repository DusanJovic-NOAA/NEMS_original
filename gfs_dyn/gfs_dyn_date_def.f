      module gfs_dyn_date_def
      use gfs_dyn_machine
      implicit none
      save
      integer idate(4)
      real(kind=kind_evod) fhour,shour,thour,z00
      REAL(KIND=KIND_EVOD) ,ALLOCATABLE :: spdmax(:)

      end module gfs_dyn_date_def
