      SUBROUTINE mpi_quit(iret)
      use gfs_dyn_mpi_def
      implicit none
!
      integer iret
 
      write(*,*) 'CALL stop mpi_quit ',iret
      CALL MPI_ABORT(MC_COMP,iret,info)
 
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
