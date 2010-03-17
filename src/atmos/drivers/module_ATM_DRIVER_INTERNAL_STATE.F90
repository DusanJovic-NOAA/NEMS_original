!---------------------------------------------------------------------------
!
! !MODULE: MODULE_ATM_DRIVER_INTERNAL_STATE --- Internal state definition 
!                                               of the ATM_DRIVER component.
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
      MODULE MODULE_ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
      USE ESMF_Mod
!
!---------------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: ATM_DRIVER_INTERNAL_STATE                                   &
               ,WRAP_ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      TYPE ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
        TYPE(ESMF_GridComp),DIMENSION(:),POINTER :: ATM_GRID_COMP              !<-- NMM ATM gridded components of all domains
!
        TYPE(ESMF_State),DIMENSION(:),POINTER :: IMP_STATE_ATM                 !<-- The import state of the NMM ATM components
        TYPE(ESMF_State),DIMENSION(:),POINTER :: EXP_STATE_ATM                 !<-- The export state of the NMM ATM components
!
!---------------------------------------------------------------------------
!
      END TYPE ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!***  This state is supported by C pointers but not by F90 pointers
!***  therefore we use this "WRAP".
!---------------------------------------------------------------------------
!
      TYPE WRAP_ATM_DRIVER_INTERNAL_STATE
        TYPE(ATM_DRIVER_INTERNAL_STATE),POINTER :: ATM_DRV_INT_STATE
      END TYPE WRAP_ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      END MODULE MODULE_ATM_DRIVER_INTERNAL_STATE
!
!---------------------------------------------------------------------------
