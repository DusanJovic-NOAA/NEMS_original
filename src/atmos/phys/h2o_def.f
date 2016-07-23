      module h2o_def
      use machine , only : kind_phys
      implicit none
      save
      integer, parameter :: kh2opltc=29
      integer latsh2o, levh2o, timeh2o,  h2o_coeff
      real (kind=kind_phys), allocatable :: h2o_lat(:), h2o_pres(:)
     &,                                     h2o_time(:)
!     integer  nlevpl, nlevt, nlevc
!     real (kind=kind_phys), allocatable :: prefpl(:), preft(:), prefc(:)
!    &, h2opl_lat(:), h2ot_lat(:), h2oc_lat(:)
      end module h2o_def
