    module module_timers
      implicit none
      public timef
      contains
      function timef()
       real(kind=8):: timef,mpi_wtime
       timef=MPI_Wtime()
      end function timef
    end module module_timers
