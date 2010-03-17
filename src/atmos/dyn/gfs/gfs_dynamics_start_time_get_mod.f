      subroutine gfs_dynamics_start_time_get(				&
                 yy, mm, dd, hh, mns, sec, kfhour, n1, n2,    		&
                 grib_inp, cfile, cfile2, rc)

! this subroutine gets and calculates the start time from reading the
! sigma file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!  Feb   2010      jun wang   read data from nemsio file
!
!uses:
!
      use esmf_mod,     only: esmf_success
      use gfs_dyn_machine,      only: kind_io4, kind_evod
      use gfs_dyn_date_def,     only: idate,idate7
!     use sigio_module
!     use sigio_r_module
      use gfsio_module
      use gfsio_def
      use module_nemsio

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
      real(kind = kind_evod)         :: fhour
      real(kind = kind_io4)          :: fhour4
      type(nemsio_gfile) :: nfile
      character (len=*)              :: cfile, cfile2
      integer iret, khour

      n1    = 11
      n2    = 12
 
      print *,' grib_inp=',grib_inp,' cfile=',cfile
      call gfsio_open(gfile_in,trim(cfile),'read',iret)
      if(iret==0) then
        print *,'after gfsio open, iert=',iret
        call gfsio_getfilehead(gfile_in,iret=iret,idate=idate,fhour=fhour4)
        print *,'after gfsio getfile head, fhour=',iret
        call gfsio_close(gfile_in)
        yy     = idate(4)
        mm     = idate(2)
        dd     = idate(3)
        hh     = idate(1)
        mns    = 0
        sec    = 0
        fhour = fhour4
      else
!        
        call gfsio_close(gfile_in)
!
        call nemsio_init()
        call nemsio_open(nfile,trim(cfile),'read',iret=iret)
        print *,'after nemsio_open, iret=',iret
        call nemsio_getheadvar(nfile,'idate',idate7,iret=iret)
        call nemsio_getheadvar(nfile,'fhour',fhour,iret=iret)
        print *,'after nemsio,idate=',idate7,'fhour=',fhour
        call nemsio_close(nfile)
        call nemsio_finalize()
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
      endif

      print *,' fhour=',fhour,' idate=',idate7,' iret=',iret
      if (iret .ne. 0) call mpi_quit(5555)
      kfhour = nint(fhour)
      print *,' idate=',idate,' fhour=',fhour,' kfhour=',kfhour

      rc = rc1

      end subroutine gfs_dynamics_start_time_get
