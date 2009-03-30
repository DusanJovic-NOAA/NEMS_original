!
! !description: atmos gridded component error messages
!
! !revision history:
!
!  january 2007 	hann-ming henry juang
!
!
! !interface:
!
      module atmos_err_msg_mod

!
!!uses:
!
      use esmf_mod 

      implicit none
      logical, parameter :: lprint = .false.

      contains

      subroutine atmos_err_msg_int(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      integer,           intent(in) :: val
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',val,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass ',msg,' ',val
      end if
      return
      end subroutine atmos_err_msg_int

      subroutine atmos_err_msg_val(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      real,              intent(in) :: val
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',val,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass ',msg,' ',val
      end if
      return
      end subroutine atmos_err_msg_val

      subroutine atmos_err_msg_var(rc1,msg,chr,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      character (len=*), intent(in) :: chr
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',chr,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass ',msg,' ',chr
      end if
      return
      end subroutine atmos_err_msg_var

      subroutine atmos_err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      if(esmf_logmsgfounderror(rc1, msg)) then
          rc  = esmf_failure
          print*, 'error happened for ',msg, ' rc = ', rc1
          rc1 = esmf_success
      else
          if(lprint) print*, 'pass ',msg
      end if
      return
      end subroutine atmos_err_msg

      subroutine atmos_err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal
      integer, intent(inout)        :: rc
      character (len=*), intent(in) :: msg
      if(rcfinal == esmf_success) then
          if(lprint) print*, "final pass: ",msg
      else
          print*, "final fail: ",msg
      end if
!     if(present(rc)) then
          rc = rcfinal
!     end if
      return
      end subroutine atmos_err_msg_final

      end module atmos_err_msg_mod
