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

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

! One by one get the parameters from the GFS ESMF export state.
!--------------------------------------------------------------
 CALL ESMF_AttributeGet(State, 'NTRAC', Int_State%ntrac, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get ntrac from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting ntrac from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'MPI_R_MPI_R', Int_State%MPI_R_MPI_R, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get MPI_R_MPI_R from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting MPI_R_MPI_R from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'JCAP', Int_State%jcap, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get JCAP from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting JCAP from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'NODES_COMP', Int_State%nodes_comp, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get NODES_COMP from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting NODES_COMP from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'ME_COMP', Int_State%me_comp, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get ME_COMP from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting ME_COMP from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'MC_COMP', Int_State%MC_COMP, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get MC_COMP from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting MC_COMP from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'LATS_NODE_A', Int_State%lats_node_a, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get LATS_NODE_A from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting LATS_NODE_A from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'IPT_LATS_NODE_A', Int_State%ipt_lats_node_a, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get IPT_LATS_NODE_A from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting IPT_LATS_NODE_A from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'LONF', Int_State%lonf, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get LONF from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting LONF from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'LATG', Int_State%latg, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get LATG from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting LATG from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 NULLIFY (Int_State%global_lats_a)
 NULLIFY (Int_State%lonsperlat   )
 ALLOCATE(Int_State%global_lats_a(Int_State%latg))
 ALLOCATE(Int_State%lonsperlat   (Int_State%latg))

 CALL ESMF_AttributeGet(State, 'GLOBAL_LATS_A', Int_State%latg, Int_State%global_lats_a, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get GLOBAL_LATS_A from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting GLOBAL_LATS_A from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeGet(State, 'LONSPERLAT', Int_State%latg, Int_State%lonsperlat, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get LONSPERLAT from the GEFS Cpl import state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting LONSPERLAT from the GEFS Cpl import state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: GEFS_GetParameterFromStateMod.f"
 ELSE
     PRINT*, "FAIL: GEFS_GetParameterFromStateMod.f"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE GEFS_GetParameterFromState

 END MODULE GEFS_GetParameterFromStateMod

