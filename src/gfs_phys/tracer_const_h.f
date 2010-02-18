      module tracer_const
      use machine , only : kind_phys
      implicit none

! !revision history:
!
!  09Feb2010   Sarah Lu, ri/cpi changed to allocatable


!     real(kind=kind_phys) ri(0:20),cpi(0:20)
      real(kind=kind_phys), allocatable ::  ri(:),cpi(:)
      integer, parameter :: num_tracer=3

      contains
! -------------------------------------------------------------------   
      subroutine set_tracer_const (ntrac,me,nlunit)
      use machine , only : kind_phys
      use physcons , only : rd => con_rd , cpd => con_cp
      implicit none
      integer ntrac,me,nlunit
      namelist /tracer_constant/ ri,cpi

c
      if( ntrac.ne.num_tracer ) then
        if( me.eq.0 ) then
          write(0,*) ' Error ; inconsistent number of tracer '
          write(0,*) ' ntrac=',ntrac,' num_tracer=',num_tracer
        endif
        call abort
      endif

!
!! This routine is now called by NMMB only                   (Sarah Lu)
!! For GFS core, CPI/RI is passed in from DYN export state
!! The allocation below is to support NMMB+GFS_physics package
      if (.not. allocated(ri)) then
        allocate( ri(0:num_tracer))
        allocate(cpi(0:num_tracer))
      endif
!
      ri=0.0
      cpi=0.0
      ri(0)=rd
      cpi(0)=cpd

      rewind(nlunit)
      read(nlunit, tracer_constant)
      write(0, tracer_constant)

      return
      end subroutine set_tracer_const

      end module tracer_const
