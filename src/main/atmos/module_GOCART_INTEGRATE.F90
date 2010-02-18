!-----------------------------------------------------------------------
!
      MODULE MODULE_GOCART_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PRIMARY INTEGRATION RUNSTREAM OF GOCART
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-02-04  Lu    First crack.
!-----------------------------------------------------------------------

!
      USE ESMF_MOD
      USE MODULE_ERR_MSG

!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!
      PUBLIC :: GOCART_INTEGRATE
!
!-----------------------------------------------------------------------

      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GOCART_INTEGRATE(                                     &
                                  GC_GFS_CHEM,                         &
                                  GC_PHY2CHEM_CPL,                     &
                                  EXP_GFS_PHY,                         &
                                  IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                  CLOCK_ATM, RC_LOOP                     )

!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-02-04  Lu    First crack.
!
!-----------------------------------------------------------------------

      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_CHEM                 !<-- The GOCART grid component
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_PHY2CHEM_CPL              !<-- The Phy-to-Chem coupler component
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_GFS_PHY,                 &  !<-- The export states for Physics component
                                        IMP_GFS_CHEM,EXP_GFS_CHEM       !<-- The imp/exp states for Chemistry component
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM Component's ESMF Clock

      INTEGER,INTENT(OUT) :: RC_LOOP                                   !<-- Return code

! Locals
      INTEGER             :: RC=ESMF_SUCCESS  

       print *, 'LU_TST: enter GOCART_INTEGRATE'

!-----------------------------------------------------------------------
!***  Couple Physics export state to Chemistry import State
!-----------------------------------------------------------------------
       print *, 'LU_TST: GOCART_INTEGRATE: run GC_PHY2CHEM_CP'
       CALL ESMF_CplCompRun(cplcomp     = GC_PHY2CHEM_CPL            &
                            ,importstate= EXP_GFS_PHY                &
                            ,exportstate= IMP_GFS_CHEM               &
                            ,clock      = CLOCK_ATM                  &
                            ,rc         = RC)
!
       CALL ERR_MSG(RC,'couple phy_exp-to-chem_imp',RC_LOOP)
       print *, 'LU_TST: GOCART_INTEGRATE: exit GC_PHY2CHEM_CP', RC

!-----------------------------------------------------------------------
!***  Execute the Run step of the Chemistry Component
!-----------------------------------------------------------------------
       print *, 'LU_TST: GOCART_INTEGRATE: skip GC_GFS_CHEM'
!	print *, 'LU_TST: GOCART_INTEGRATE: run GC_GFS_CHEM'
!          CALL ESMF_GridCompRun(gridcomp   =GC_GFS_CHEM                 &
!                               ,importstate=IMP_GFS_CHEM                &
!                               ,exportstate=EXP_GFS_CHEM                &
!                               ,clock      =CLOCK_ATM                   &
!                               ,rc         =RC)
!
!          CALL ERR_MSG(RC,'execute chemistry',RC_LOOP)
!	print *, 'LU_TST: GOCART_INTEGRATE: exit GC_GFS_CHEM', RC

!-----------------------------------------------------------------------
!***  The chem-to-dyn coupler is not needed, exit now
!-----------------------------------------------------------------------
       print *, 'LU_TST: GOCART_INTEGRATE: GC_CHEM2PHY_CPL not needed'

       print *, 'LU_TST: exit GOCART_INTEGRATE'
!
      END SUBROUTINE GOCART_INTEGRATE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GOCART_INTEGRATE
!
!-----------------------------------------------------------------------


