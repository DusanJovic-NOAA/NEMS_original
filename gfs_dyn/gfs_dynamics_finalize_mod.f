!
! !module: gfs_dynamics_finalize_mod 
!          --- finalize module of the grided
!              component of the gfs dynamics system.
!
! !description: gfs finalize module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  february 2006      shrinivas moorthi - removed some comments
!  january  2007      hann-ming henry juang -- modify to be dynamics only
!
!
! !interface:
!
      module gfs_dynamics_finalize_mod
!
!!uses:
!
      use gfs_dynamics_internal_state_mod
      use gfsio_module , only : gfsio_finalize

      implicit none

      contains

      subroutine gfs_dynamics_finalize(gis_dyn, rc)

      type(gfs_dynamics_internal_state)                :: gis_dyn
      integer, optional,                 intent(out)   :: rc

!
!***********************************************************************
!
      call countperf(0,15,0.)
      call synchro
      call countperf(1,15,0.)
!!
      call end_countperf()
      call write_countperf(nodes,me,fhour)
!!
!     if(me.eq.0) then
!       call w3tage('gsm     ')
!     endif
!!
!c    jjt terminate hpm ... this is needed
      call f_hpmterminate(me)
!
      if (gfsio_out .or. gfsio_in) then
        call gfsio_finalize()
      endif

      if(present(rc)) then
      rc = 0
      end if

      end subroutine gfs_dynamics_finalize

      end module gfs_dynamics_finalize_mod
