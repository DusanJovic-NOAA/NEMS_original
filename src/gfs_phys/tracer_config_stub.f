!
! module: tracer_config_stub
!
! ! Revision history:
!   Oct 17 2009   Sarah Lu, First version
! -------------------------------------------------------------------------
!
      module gfs_phy_tracer_config
      implicit none

      type    gfs_phy_tracer_type
        character*20,    pointer     :: vname(:)    ! variable name
        integer                  :: ntrac
        integer                  :: ntrac_met
        integer                  :: ntrac_chem
      endtype gfs_phy_tracer_type
!
      contains

! -------------------------------------------------------------------   
      subroutine tracer_config_init (gfs_phy_tracer,ntrac,      &
     &                               ntoz,ntcw,ncld,me)

      implicit none
      integer         ::  me, ntoz,ntcw,ncld, ntrac
      type (gfs_phy_tracer_type)    ::  gfs_phy_tracer

      return

      end subroutine tracer_config_init


      end module gfs_phy_tracer_config
