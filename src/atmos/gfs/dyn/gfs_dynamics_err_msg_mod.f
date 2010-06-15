!
! !description: gfs dynamics gridded component error messages
!
! !revision history:
!
!  january 2007 	hann-ming henry juang
!
!
! !interface:
!
      module gfs_dynamics_err_msg_mod

!
!!uses:
!
      use esmf_mod                 ! the esmf library.

      implicit none
      logical,parameter::lprint=.false.
      contains

      subroutine gfs_dynamics_err_msg_var(rc1,msg,var,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      character (len=*), intent(in) :: var
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened in dynamics for ',msg,' ',var,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass in dynamics for ',msg,' ',var
      end if
      end subroutine gfs_dynamics_err_msg_var

      subroutine gfs_dynamics_err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      if(esmf_logmsgfounderror(rc1, msg)) then
          rc  = esmf_failure
          print*, 'error happened in dynamics for ',msg, ' rc = ', rc1
          rc1 = esmf_success
      else
          if(lprint) print*, 'pass in dynamics for ',msg
      end if
      end subroutine gfs_dynamics_err_msg

      subroutine gfs_dynamics_err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal
      integer, intent(inout)        :: rc
      character (len=*), intent(in) :: msg
      if(rcfinal == esmf_success) then
         if(lprint) print*, "final pass in dynamics for ",msg
      else
         print*, "final fail in dynamics for ",msg
      end if
      rc = rcfinal
      end subroutine gfs_dynamics_err_msg_final

      end module gfs_dynamics_err_msg_mod
