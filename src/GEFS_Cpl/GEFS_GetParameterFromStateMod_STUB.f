!----------------------------------------------------------------------
! !MODULE: GEFS_GetParameterFromStateMod
!        --- Get required parameters from the GEFS Coupler ESMF import state
!            for the ensemble coupler to do the spectral transform
!            for the stochastic perturbation scheme, the second step.
!
! !DESCRIPTION: Get all required parameters from the GEFS Cpl ESMF import state.
!
! !REVISION HISTORY:
!
!  May      2007     Weiyu Yang Initial code.
!  March    2009     Weiyu Yang Modified for the NEMS model.
!
!
! !INTERFACE:
!

 MODULE GEFS_GetParameterFromStateMod

 USE ESMF_Mod
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 CONTAINS

 SUBROUTINE GEFS_GetParameterFromState(State, Int_State, rc)

 TYPE(ESMF_State),                      INTENT(inout) :: State
 TYPE(GEFS_Cpl_InternalState), POINTER, INTENT(inout) :: Int_State
 INTEGER, OPTIONAL,                     INTENT(out)   :: rc

 INTEGER                                              :: rc1, rcfinal

 END SUBROUTINE GEFS_GetParameterFromState

 END MODULE GEFS_GetParameterFromStateMod

