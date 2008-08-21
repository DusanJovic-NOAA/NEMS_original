      module gfs_dyn_tracer_const
      use gfs_dyn_machine , only : kind_grid
      implicit none
      SAVE

      real(kind=kind_grid) ri(0:20),cpi(0:20)
      integer, parameter :: num_tracer=3

      contains
! -------------------------------------------------------------------   
      subroutine get_tracer_const (ntrac,me,nlunit)
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_physcons , only : rd => con_rd , cpd => con_cp
      implicit none
      integer ntrac,me,nlunit
      namelist /tracer_constant/ ri,cpi

!     print *,' enter get_tracer_const',ntrac,num_tracer
c
      if( ntrac.ne.num_tracer ) then
        if( me.eq.0 ) then
          write(*,*) ' Error ; inconsistent number of tracer '
          write(*,*) ' ntrac=',ntrac,' num_tracer=',num_tracer
        endif
        call abort
      endif

      ri=0.0
      cpi=0.0
      ri(0)=rd
      cpi(0)=cpd

      rewind(nlunit)
      read(nlunit, tracer_constant)
      write(*, tracer_constant)

!     print *,' done get_tracer_const'

      return
      end subroutine get_tracer_const

      end module gfs_dyn_tracer_const
