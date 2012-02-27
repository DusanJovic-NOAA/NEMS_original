#include "../../../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#else
#define ESMF_520r
#endif

!
      module atmos_phy_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!
! stub for atmos_phy_chem_cpl_comp_mod 
!
!! Code Revision:
!! 10Aug 2011     Jun Wang,   create stub for atmos_phy_chem_cpl_comp_mod
!------------------------------------------------------------------------------

#ifdef ESMF_520r
      USE esmf
#else
      USE esmf_mod
#endif
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      end module atmos_phy_chem_cpl_comp_mod
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
