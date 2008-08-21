      module gg_def
      use machine

      implicit none
      save
      REAL(KIND=KIND_EVOD) ,ALLOCATABLE ::  colrad_r(:),wgt_r(:),
     & wgtcs_r(:),rcs2_r(:),sinlat_r(:),coslat_r(:)
      end module gg_def
