      module gfs_dyn_io_header
      use gfs_dyn_machine
      implicit none
      save
      integer              ifin
      integer              icen
      integer              icen2
      integer              ienst
      integer              iensi
      integer              itrun
!
      integer lonb, latb, iens(5), idpp, idvt, idrun
     &,       idusr, ncldt, irealf, iorder
!!
      character*8   lab(4)
      end module gfs_dyn_io_header
