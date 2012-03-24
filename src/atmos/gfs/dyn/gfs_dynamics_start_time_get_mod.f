#include "../../../ESMFVersionDefine.h"

      subroutine gfs_dynamics_start_time_get(				&
                 yy, mm, dd, hh, mns, sec, kfhour, n1, n2,    		&
                 grib_inp, cfile, cfile2, rc)

! this subroutine gets and calculates the start time from reading the
! sigma file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!  Feb   2010      jun wang   read data from nemsio file
!  Sep   2010      jun wang   remove gfsio option
!  Dec   2010      jun wang   change to nemsio library
!  Feb   2011      sarah lu   change to read nfhour
!
!uses:
!
      USE esmf_mod, ONLY: esmf_success

      use gfs_dyn_machine,      only: kind_io4, kind_evod
      use gfs_dyn_date_def,     only: idate,idate7
      use nemsio_module

      implicit none

!
! arguments:
!-----------

      integer,                intent(in)  :: grib_inp
      integer,                intent(out) :: yy, mm, dd, hh, mns, sec
      integer,                intent(out) :: n1, n2
      integer,                intent(out) :: kfhour
      integer,                intent(out) :: rc     ! return code

      integer                        :: rc1 = esmf_success
      integer                        :: nfhour 
      type(nemsio_gfile) :: nfile
      character (len=*)              :: cfile, cfile2
      integer iret, khour

      n1    = 11
      n2    = 12
 
      print *,' grib_inp=',grib_inp,' cfile=',cfile
!        write(0,*)'in dyn_start_time'
        call nemsio_init()
        call nemsio_open(nfile,trim(cfile),'read',iret=iret)
!        print *,'nemsio open instart iret=',iret
        call nemsio_getfilehead(nfile,idate=idate7,nfhour=kfhour,   & 
     &     iret=iret)                                            
        call nemsio_close(nfile)
        call nemsio_finalize()
!
        idate(1)=idate7(4)
        idate(2:3)=idate7(2:3)
        idate(4)=idate7(1)
        yy     = idate7(1)
        mm     = idate7(2)
        dd     = idate7(3)
        hh     = idate7(4)
        mns    = idate7(5)
        if(idate7(7)/=0) then
          sec    = idate7(6)*1./idate7(7)
        else
          sec    = 0
        endif
!      endif

!      print *,' fhour=',fhour,' idate=',idate7,' iret=',iret
      if (iret .ne. 0) call mpi_quit(5555)
      print *,' idate=',idate,' kfhour=',kfhour

      rc = rc1

      end subroutine gfs_dynamics_start_time_get
