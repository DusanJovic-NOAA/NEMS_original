!
! !MODULE: ENS_CplState_ESMFMod --- Run module of the ESMF grided
!                                   component of the EARTH ensemble coupler.
!
! !DESCRIPTION: ENS coupler run module.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!  May      2008     Weiyu Yang updated to use the ESMF 3.1.0r library.
!
!
! !INTERFACE:
!
 MODULE ENS_CplState_ESMFMod
!
!!USES:
!
 USE ESMF_Mod
 USE ENS_Cpl_InternalState_ESMFMod
 USE Lib_ESMFStateAddGetMod

 TYPE(ESMF_Grid), SAVE                       :: grid1
 TYPE(ESMF_Grid), SAVE                       :: grid2

 IMPLICIT none

 CONTAINS

 SUBROUTINE ENS_Cpl_ESMFImportState2InternalState(impENS, Cpl_Int_State, rc)

! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_State),            INTENT(in)    :: impENS
 TYPE(ENS_Cpl_InternalState), INTENT(inout) :: Cpl_Int_State

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER, OPTIONAL,           INTENT(out)   :: rc

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 CHARACTER(ESMF_MAXSTR)                     :: name
 INTEGER                                    :: rc1
 INTEGER                                    :: rcfinal

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

 CALL GetF90ArrayFromState(impENS, 'pps', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- PS.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- PS, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'tt', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- T.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- T, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T"
 CALL DistributeForStep1(Cpl_Int_State%t, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'uu', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- U.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- U, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U"
 CALL DistributeForStep1(Cpl_Int_State%u, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'vv', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- V.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- V, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V"
 CALL DistributeForStep1(Cpl_Int_State%v, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'sshum', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- q.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- q, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "q"
 CALL DistributeForStep1(Cpl_Int_State%q, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- q.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- q, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ooz', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- OZ.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- OZ, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "OZ"
 CALL DistributeForStep1(Cpl_Int_State%oz, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- OZ.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- OZ, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ccld', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- CLW.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- CLW, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "CLW"
 CALL DistributeForStep1(Cpl_Int_State%clw, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- CLW.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- CLW, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'psm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- PS(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- PS(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%psm, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'tm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- T(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- T(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%tm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'um', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- U(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- U(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%um, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'vm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- V(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- V(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%vm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'shumm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- q(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- q(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "q(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%qm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- q(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- q(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ozm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- OZ(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- OZ(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "OZ(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%ozm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- OZ(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- OZ(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'cldm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- CLW(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- CLW(-DT), rc = ', rc1
     END IF

 name = "CLW(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%clwm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- CLW(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- CLW(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'pps6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- PS(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- PS(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps6, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'tt6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- T(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- T(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%t6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'uu6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- U(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- U(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%u6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'vv6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- V(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- V(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%v6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'sshum6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- q(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- q(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "q(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%q6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- q(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- q(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ooz6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- OZ(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- OZ(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "OZ(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%oz6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- OZ(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- OZ(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ccld6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- CLW(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- CLW(-6H), rc = ', rc1
     END IF

 name = "CLW(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%clw6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- CLW(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- CLW(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'psm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- PS(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- PS(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps6m, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'tm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- T(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- T(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%t6m,  name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'um6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- U(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- U(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%u6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'vm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- V(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- V(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%v6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'shumm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- q(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- q(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "q(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%q6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- q(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- q(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'ozm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- OZ(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- OZ(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "OZ(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%oz6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- OZ(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- OZ(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impENS, 'cldm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impENS to the Cpl Internal State -- CLW(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impENS to the Cpl Internal State -- CLW(-6H-DT), rc = ', rc1
     END IF

 name = "CLW(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%clw6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- CLW(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- CLW(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: ENS_Cpl_ESMFImportState2InternalState"
 ELSE
     PRINT*, "FAIL: ENS_Cpl_ESMFImportState2InternalState"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE ENS_Cpl_ESMFImportState2InternalState





 SUBROUTINE ENS_Cpl_InternalState2ESMFExportState(impENS, expENS, Cpl_Int_State, rc)
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_State),            INTENT(inout) :: impENS
 TYPE(ESMF_State),            INTENT(inout) :: expENS
 TYPE(ENS_Cpl_InternalState), INTENT(inout) :: Cpl_Int_State

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER, OPTIONAL,           INTENT(out)   :: rc

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_VM)                              :: vm
 TYPE(ESMF_Array)                           :: ESMFArray1
 TYPE(ESMF_Array)                           :: ESMFArray2
 TYPE(ESMF_DistGrid)                        :: DistGrid1 
 TYPE(ESMF_DistGrid)                        :: DistGrid2   
 INTEGER                                    :: arraysize_1
 INTEGER                                    :: arraysize_2
 INTEGER                                    :: arraysize_max
 INTEGER                                    :: arraysize_1_max
 INTEGER                                    :: i, j
 INTEGER                                    :: rc1
 INTEGER                                    :: rcfinal
 INTEGER, DIMENSION(:), POINTER             :: arraysize_1_pointer
 INTEGER, DIMENSION(:), POINTER             :: arraysize_1_gather
 LOGICAL                                    :: first

 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: t_wk

 SAVE first

 DATA first/.true./

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

 IF(first) THEN
     ALLOCATE(arraysize_1_pointer(1))
     ALLOCATE(arraysize_1_gather (Cpl_Int_State%nodes))

     CALL ESMF_VMGetGlobal(vm, rc = rc1)

! Get the ESMF grid of the 2-D array.
!------------------------------------
     CALL ESMF_StateGet(impENS, 'pps', ESMFArray1, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF Array1")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF Array1, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL ESMF_ArrayGet(ESMFArray1, distgrid = DistGrid1, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF DistGrid1")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF DistGrid1, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     grid1 = ESMF_GridCreate(name = "Cpl Grid1", distgrid = DistGrid1, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF Grid1")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF Grid1, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Get the ESMF grid of the 3-D array.
!------------------------------------
     CALL ESMF_StateGet(impENS, 'tt', ESMFArray2, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF Array2")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF Array2, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL ESMF_ArrayGet(ESMFArray2, distgrid = DistGrid2, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF DistGrid2")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF DistGrid2, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     grid2 = ESMF_GridCreate(name = "Cpl Grid2", distgrid = DistGrid2, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF Grid2")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF Grid2, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Get the array sizes.
!---------------------
     CALL ESMF_ArrayGet(ESMFArray2, 0, t_wk, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get t for its Size")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting t for its Size, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     arraysize_1 = SIZE(t_wk, 1)
     arraysize_2 = SIZE(t_wk, 2)

     arraysize_1_pointer(1)    = arraysize_1 
     Cpl_Int_State%arraysize_1 = arraysize_1 
     Cpl_Int_State%arraysize_2 = arraysize_2 
     NULLIFY(t_wk)

! Allocate the writing out arrays at the current time step.
!----------------------------------------------------------
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tw  )) &
         ALLOCATE(Cpl_Int_State%tw  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%qw  )) &
         ALLOCATE(Cpl_Int_State%qw  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ozw )) &
         ALLOCATE(Cpl_Int_State%ozw (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clww)) &
         ALLOCATE(Cpl_Int_State%clww(arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%uw  )) &
         ALLOCATE(Cpl_Int_State%uw  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vw  )) &
         ALLOCATE(Cpl_Int_State%vw  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psw )) &
         ALLOCATE(Cpl_Int_State%psw (arraysize_1, 1))

! Allocate the writing out arrays at the last time step.
!-------------------------------------------------------
     IF(.NOT. ASSOCIATED(Cpl_Int_State%twm  )) &
         ALLOCATE(Cpl_Int_State%twm  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%qwm  )) &
         ALLOCATE(Cpl_Int_State%qwm  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ozwm )) &
         ALLOCATE(Cpl_Int_State%ozwm (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clwwm)) &
         ALLOCATE(Cpl_Int_State%clwwm(arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%uwm  )) &
         ALLOCATE(Cpl_Int_State%uwm  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vwm  )) &
         ALLOCATE(Cpl_Int_State%vwm  (arraysize_1, arraysize_2))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%pswm )) &
         ALLOCATE(Cpl_Int_State%pswm (arraysize_1, 1))

     DO j = 1, arraysize_2
         DO i = 1, arraysize_1
             Cpl_Int_State%tw   (i, j) = 0.0
             Cpl_Int_State%qw   (i, j) = 0.0
             Cpl_Int_State%clww (i, j) = 0.0
             Cpl_Int_State%uw   (i, j) = 0.0
             Cpl_Int_State%vw   (i, j) = 0.0
             Cpl_Int_State%twm  (i, j) = 0.0
             Cpl_Int_State%qwm  (i, j) = 0.0
             Cpl_Int_State%clwwm(i, j) = 0.0
             Cpl_Int_State%uwm  (i, j) = 0.0
             Cpl_Int_State%vwm  (i, j) = 0.0
         END DO
     END DO

     DO i = 1, arraysize_1
         Cpl_Int_State%psw (i, 1) = 0.0
         Cpl_Int_State%pswm(i, 1) = 0.0
     END DO

! Get the size of the one piece for the alltoall -- ARRAY_TOT_SIZ2.
!------------------------------------------------------------------
     CALL ESMF_VMAllGather(vm,                  &
                           arraysize_1_pointer, &
                           arraysize_1_gather,  &
                           1,                   &
                           rc     = rc1)
     arraysize_1_max = 0
     DO i = 1, Cpl_Int_State%nodes
         IF(arraysize_1_max < arraysize_1_gather(i)) THEN
             arraysize_1_max = arraysize_1_gather(i)
         END IF
     END DO
     arraysize_max = arraysize_1_max * arraysize_2

     Cpl_Int_State%ARRAY_ONE_SIZ2 = arraysize_1_max / Cpl_Int_State%nodes
     IF(MOD(arraysize_1_max, Cpl_Int_State%nodes) /= 0) THEN
         Cpl_Int_State%ARRAY_ONE_SIZ2 = Cpl_Int_State%ARRAY_ONE_SIZ2 + 1
     END IF

     Cpl_Int_State%ARRAY_TOT_SIZ2 = arraysize_max / Cpl_Int_State%nodes
     IF(MOD(arraysize_max, Cpl_Int_State%nodes) /= 0) THEN
         Cpl_Int_State%ARRAY_TOT_SIZ2 = Cpl_Int_State%ARRAY_TOT_SIZ2 + 1
     END IF

! Get the whole size of the array -- ARRAY_TOT_SIZ3.
!---------------------------------------------------
     Cpl_Int_State%ARRAY_ONE_SIZ3 = Cpl_Int_State%ARRAY_ONE_SIZ2 * Cpl_Int_State%nodes
     Cpl_Int_State%ARRAY_TOT_SIZ3 = Cpl_Int_State%ARRAY_TOT_SIZ2 * Cpl_Int_State%nodes

     IF(.NOT. ASSOCIATED(Cpl_Int_State%work3)) &
         ALLOCATE(Cpl_Int_State%work3(Cpl_Int_State%ARRAY_ONE_SIZ3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%work5)) &
         ALLOCATE(Cpl_Int_State%work5(Cpl_Int_State%ARRAY_TOT_SIZ3))

! Get the size of the array for one ensemble member -- ARRAY_TOT_SIZ4.
!---------------------------------------------------------------------
     Cpl_Int_State%ARRAY_ONE_SIZ4 = Cpl_Int_State%ARRAY_ONE_SIZ3 / Cpl_Int_State%Total_member
     Cpl_Int_State%ARRAY_TOT_SIZ4 = Cpl_Int_State%ARRAY_TOT_SIZ3 / Cpl_Int_State%Total_member

! Allocate all input related arrays.
!-----------------------------------
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps  )) &
         ALLOCATE(Cpl_Int_State%ps  (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psm )) &
         ALLOCATE(Cpl_Int_State%psm (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps6 )) &
         ALLOCATE(Cpl_Int_State%ps6 (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps6m)) &
         ALLOCATE(Cpl_Int_State%ps6m(Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps_mean))  &
         ALLOCATE(Cpl_Int_State%ps_mean(Cpl_Int_State%ARRAY_ONE_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psm_mean)) &
         ALLOCATE(Cpl_Int_State%psm_mean(Cpl_Int_State%ARRAY_ONE_SIZ4))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%t   )) &
         ALLOCATE(Cpl_Int_State%t   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tm  )) &
         ALLOCATE(Cpl_Int_State%tm  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%t6  )) &
         ALLOCATE(Cpl_Int_State%t6  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%t6m )) &
         ALLOCATE(Cpl_Int_State%t6m (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%u   )) &
         ALLOCATE(Cpl_Int_State%u   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%um  )) &
         ALLOCATE(Cpl_Int_State%um  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%u6  )) &
         ALLOCATE(Cpl_Int_State%u6  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%u6m )) &
         ALLOCATE(Cpl_Int_State%u6m (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%v   )) &
         ALLOCATE(Cpl_Int_State%v   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vm  )) &
         ALLOCATE(Cpl_Int_State%vm  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%v6  )) &
         ALLOCATE(Cpl_Int_State%v6  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%v6m )) &
         ALLOCATE(Cpl_Int_State%v6m (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%q   )) &
         ALLOCATE(Cpl_Int_State%q   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%qm  )) &
         ALLOCATE(Cpl_Int_State%qm  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%q6  )) &
         ALLOCATE(Cpl_Int_State%q6  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%q6m )) &
         ALLOCATE(Cpl_Int_State%q6m (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%oz  )) &
         ALLOCATE(Cpl_Int_State%oz  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ozm )) &
         ALLOCATE(Cpl_Int_State%ozm (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%oz6 )) &
         ALLOCATE(Cpl_Int_State%oz6 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%oz6m)) &
         ALLOCATE(Cpl_Int_State%oz6m(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%clw )) &
         ALLOCATE(Cpl_Int_State%clw (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clwm)) &
         ALLOCATE(Cpl_Int_State%clwm(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clw6)) &
         ALLOCATE(Cpl_Int_State%clw6(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clw6m)) &
         ALLOCATE(Cpl_Int_State%clw6m(Cpl_Int_State%ARRAY_TOT_SIZ4,Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps_step1  )) &
         ALLOCATE(Cpl_Int_State%ps_step1  (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psm_step1 )) &
         ALLOCATE(Cpl_Int_State%psm_step1 (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%t_step1   )) &
         ALLOCATE(Cpl_Int_State%t_step1   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tm_step1  )) &
         ALLOCATE(Cpl_Int_State%tm_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%u_step1   )) &
         ALLOCATE(Cpl_Int_State%u_step1   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%um_step1  )) &
         ALLOCATE(Cpl_Int_State%um_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%v_step1   )) &
         ALLOCATE(Cpl_Int_State%v_step1   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vm_step1  )) &
         ALLOCATE(Cpl_Int_State%vm_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%q_step1   )) &
         ALLOCATE(Cpl_Int_State%q_step1   (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%qm_step1  )) &
         ALLOCATE(Cpl_Int_State%qm_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%oz_step1  )) &
         ALLOCATE(Cpl_Int_State%oz_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ozm_step1 )) &
         ALLOCATE(Cpl_Int_State%ozm_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%clw_step1 )) &
         ALLOCATE(Cpl_Int_State%clw_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clwm_step1)) &
         ALLOCATE(Cpl_Int_State%clwm_step1(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%t_mean   )) ALLOCATE(Cpl_Int_State%t_mean   (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%u_mean   )) ALLOCATE(Cpl_Int_State%u_mean   (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%v_mean   )) ALLOCATE(Cpl_Int_State%v_mean   (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%q_mean   )) ALLOCATE(Cpl_Int_State%q_mean   (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%oz_mean  )) ALLOCATE(Cpl_Int_State%oz_mean  (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clw_mean )) ALLOCATE(Cpl_Int_State%clw_mean (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tm_mean  )) ALLOCATE(Cpl_Int_State%tm_mean  (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%um_mean  )) ALLOCATE(Cpl_Int_State%um_mean  (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vm_mean  )) ALLOCATE(Cpl_Int_State%vm_mean  (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%qm_mean  )) ALLOCATE(Cpl_Int_State%qm_mean  (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%ozm_mean )) ALLOCATE(Cpl_Int_State%ozm_mean (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%clwm_mean)) ALLOCATE(Cpl_Int_State%clwm_mean(Cpl_Int_State%ARRAY_TOT_SIZ4))
! Start to link the current time step arrays to the ESMF export state.
!---------------------------------------------------------------------
     CALL AddF90ArrayToState(expENS, grid2, "tt",   Cpl_Int_State%tw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add t to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding t to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "sshum",   Cpl_Int_State%qw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add q to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding q to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "ooz",  Cpl_Int_State%ozw,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add oz to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding oz to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "ccld", Cpl_Int_State%clww, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add clw to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding clw to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "uu",   Cpl_Int_State%uw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add u to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding u to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "vv",   Cpl_Int_State%vw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add v to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding v to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid1, "pps",  Cpl_Int_State%psw,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add ps to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding ps to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Start to link the last time step arrays to the ESMF export state.
!------------------------------------------------------------------
     CALL AddF90ArrayToState(expENS, grid2, "tm",   Cpl_Int_State%twm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add tm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding tm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "shumm",   Cpl_Int_State%qwm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add shumm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding shumm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "ozm",  Cpl_Int_State%ozwm,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add ozm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding ozm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "cldm", Cpl_Int_State%clwwm, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add cldm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding cldm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "um",   Cpl_Int_State%uwm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add um to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding um to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid2, "vm",   Cpl_Int_State%vwm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add vm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding vm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expENS, grid1, "psm",  Cpl_Int_State%pswm,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add psm to expENS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding psm to expENS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Finish the first = .true. run.
!-------------------------------
     first = .false.
 END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: ENS_Cpl_InternalState2ESMFExportState"
 ELSE
     PRINT*, "FAIL: ENS_Cpl_InternalState2ESMFExportState"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE ENS_Cpl_InternalState2ESMFExportState

 END MODULE ENS_CplState_ESMFMod
