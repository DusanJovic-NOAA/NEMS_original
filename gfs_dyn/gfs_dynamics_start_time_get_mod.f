      subroutine gfs_dynamics_start_time_get(				&
                 yy, mm, dd, hh, mns, sec, kfhour, n1, n2,    		&
                 grib_inp, fhrot, cfile, cfile2, rc)

! this subroutine gets and calculates the start time from reading the
! sigma file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!
!uses:
!
      use esmf_mod,     only: esmf_success
      use gfs_dyn_machine,      only: kind_io4, kind_evod
      use gfs_dyn_date_def,     only: idate
!     use sigio_module
!     use sigio_r_module
      use gfsio_module
      use gfsio_def

      implicit none

!
! arguments:
!-----------

      integer,                intent(in)  :: grib_inp
      integer,                intent(out) :: yy, mm, dd, hh, mns, sec
      integer,                intent(out) :: n1, n2
      integer,                intent(out) :: kfhour
 !logical,                intent(in)  :: gfsio_in
      real(kind = kind_evod), intent(in)  :: fhrot
      integer,                intent(out) :: rc     ! return code

      integer                        :: rc1 = esmf_success
      real(kind = kind_evod)         :: fhour
      real(kind = kind_io4)          :: fhour4
!     type(sigio_head) head
      character (len=*)              :: cfile, cfile2
      integer iret, khour

      n1    = 11
      n2    = 12
 
      print *,' grib_inp=',grib_inp,' n1=',n1
!     if (grib_inp .le. 0) then
!       call sigio_rropen(n1,cfile,iret)
!       call sigio_rrhead(n1,head,iret)

!       if(head%fhour /= fhrot) then
!          call sigio_rropen(n2,cfile2,iret)
!          call sigio_rrhead(n2,head,iret)
!       end if
!
!       fhour  = head%fhour
!       idate  = head%idate
!     else
        print *,' grib_inp=',grib_inp,' cfile=',cfile
        call gfsio_open(gfile_in,trim(cfile),'read',iret)
        call gfsio_getfilehead(gfile_in,iret=iret,idate=idate,fhour=fhour4)
        fhour = fhour4

        print *,' fhour4=',fhour4,' idate=',idate,' iret=',iret
        if (iret .ne. 0) call mpi_quit(5555)
!     endif
      yy     = idate(4)
      mm     = idate(2)
      dd     = idate(3)
      hh     = idate(1)
      mns    = 0
      sec    = 0
      kfhour = nint(fhour)
      print *,' idate=',idate,' fhour=',fhour,' kfhour=',kfhour

      rc = rc1

      end subroutine gfs_dynamics_start_time_get
