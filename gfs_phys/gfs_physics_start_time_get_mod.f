      subroutine gfs_physics_start_time_get(				&
                 yy, mm, dd, hh, mns, sec, kfhour, n1, cfile, rc)

! this subroutine gets and calculates the start time from reading the
! surface file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!        2006      modified version 
!  november 2007   henry juang
!
!uses:
!
      use esmf_mod,     only: esmf_success
      use machine,      only: kind_io4, kind_evod
      use date_def,     only: idate
      use sfcio_module

      implicit none

!
! arguments:
!-----------

      integer,                intent(out) :: yy, mm, dd, hh, mns, sec
      integer,                intent(out) :: n1
      integer,                intent(out) :: kfhour
      integer,                intent(out) :: rc     ! return code

      integer                        :: rc1 = esmf_success
      real(kind = kind_evod)         :: fhour
      real(kind = kind_io4)          :: fhour4
      type(sfcio_head) head
      character (len=*)              :: cfile
      integer iret, khour

      n1    = 13
 
      call sfcio_sropen(n1,cfile,iret)
      call sfcio_srhead(n1,head,iret)
      call sfcio_sclose(n1,iret)

      fhour  = head%fhour
      idate  = head%idate

      yy     = idate(4)
      mm     = idate(2)
      dd     = idate(3)
      hh     = idate(1)
      mns    = 0
      sec    = 0
      kfhour = nint(fhour)
      print *,' idate=',idate,' fhour=',fhour,' kfhour=',kfhour

      rc = rc1

      end subroutine gfs_physics_start_time_get
