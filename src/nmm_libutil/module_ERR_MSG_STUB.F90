!
! !description: error messages
!
! !revision history:
!
!  january 2007 	hann-ming henry juang
!
!
! !interface:
!
      module module_err_msg

!
!!uses:
!
      use esmf_mod 

      implicit none

      private
      public :: err_msg,message_check


      contains

      subroutine err_msg_int(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      integer,           intent(in) :: val
      end subroutine err_msg_int

      subroutine err_msg_val(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      real,              intent(in) :: val
      end subroutine err_msg_val

      subroutine err_msg_var(rc1,msg,chr,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      character (len=*), intent(in) :: chr
      end subroutine err_msg_var

      subroutine err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      end subroutine err_msg

      subroutine err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal
      integer, intent(inout)        :: rc
      character (len=*), intent(in) :: msg
      end subroutine err_msg_final

      end module module_err_msg
