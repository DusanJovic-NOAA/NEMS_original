!-----------------------------------------------------------------------
!
      MODULE MODULE_GOCART_ROUTINES
!
!-----------------------------------------------------------------------
!
!*** The stub version needed for NMMB core 
!***
!*** THIS MODULE CONTAINS THE ROUTINES TO SETUP, INITIALIZE, AND RUN
!*** THE AEROSOL MODULE (GOCART)
!***
!*** THE SETUP AND INIT ROUTINES ARE CALLED FROM GFS_ATM_INIT 
!*** THE INTEGRATE ROUTINE IS CALLED FROM GFS_INTEGRATE
!*** 
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-03-05  Lu    - Create the module
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: GOCART_SETUP, GOCART_INIT
!
      CONTAINS


!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      SUBROUTINE GOCART_SETUP ( GC_GFS_CHEM                            &
                               ,IMP_GFS_CHEM                           &
                               ,EXP_GFS_CHEM                           &
                               ,GC_PHY2CHEM_CPL                        &
                               ,GC_CHEM2PHY_CPL                        &
                               ,CHEMISTRY_ON                           &
                               ,RC_SETUP                               &
                                )

!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_GFS_CHEM                    !<-- The gocart import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_CHEM                    !<-- The gocart export state
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_PHY2CHEM_CPL                 !<-- Phy to Chem coupler component
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_CHEM2PHY_CPL                 !<-- Chem to Phy coupler component
      TYPE(ESMF_Logical), INTENT(OUT)   :: CHEMISTRY_ON                    !<-- The option to activate gocart
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_SETUP                        !<-- Return code for the SETUP step
!
      RC_SETUP =ESMF_SUCCESS  
!
      END SUBROUTINE GOCART_SETUP


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GOCART_INIT ( GC_GFS_CHEM                              &
                              ,EXP_GFS_PHY                              &
                              ,IMP_GFS_CHEM                             &
                              ,EXP_GFS_CHEM                             &
                              ,GC_PHY2CHEM_CPL                          &
                              ,GC_CHEM2PHY_CPL                          &
                              ,CLOCK_ATM                                &
                              ,RC_INIT                                  &
                              )

!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_PHY                     !<-- The ATM Run step's export
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_GFS_CHEM                    !<-- The ATM Run step's import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_CHEM                    !<-- The ATM Run step's export
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_PHY2CHEM_CPL                 !<-- The Phy to Chem coupler component
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_CHEM2PHY_CPL                 !<-- The Phy to Chem coupler component
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_INIT                         !<-- Return code for the INIT step
!
      RC_INIT =ESMF_SUCCESS  
!
      END SUBROUTINE GOCART_INIT

!
      END MODULE MODULE_GOCART_ROUTINES

!-----------------------------------------------------------------------
