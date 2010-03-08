!
! !MODULE: GEFS_CplState_ESMFMod --- Run module of the ESMF grided
!                                    component of the GFS ensemble coupler.
!
! !DESCRIPTION: GEFS coupler run module.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!  May      2008     Weiyu Yang updated to use the ESMF 3.1.0r library.
!  November 2009     Weiyu Yang, NEMS GEFS model.
!
!
! !INTERFACE:
!
 MODULE GEFS_CplState_ESMFMod
!
!!USES:
!
 USE ESMF_Mod
 USE GEFS_Cpl_InternalState_ESMFMod
 USE Lib_ESMFStateAddGetMod

 TYPE(ESMF_Grid), SAVE                       :: mgrid

 IMPLICIT none

 CONTAINS

 SUBROUTINE GEFS_Cpl_ESMFImportState2InternalState(impGEFS, Cpl_Int_State, rc)

! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_State),             INTENT(in)    :: impGEFS
 TYPE(GEFS_Cpl_InternalState), INTENT(inout) :: Cpl_Int_State

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER, OPTIONAL,            INTENT(out)   :: rc

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_Field)                            :: Field
 TYPE(ESMF_FieldBundle)                      :: Bundle
 CHARACTER(ESMF_MAXSTR)                      :: name
 INTEGER                                     :: i, j
 INTEGER                                     :: rc1
 INTEGER                                     :: rcfinal

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

 CALL GetF90ArrayFromState(impGEFS, 'pps', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- PS.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- PS, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'tt', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- T.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- T, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T"
 CALL DistributeForStep1(Cpl_Int_State%t, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'uu', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- U.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- U, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U"
 CALL DistributeForStep1(Cpl_Int_State%u, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'vv', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- V.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- V, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V"
 CALL DistributeForStep1(Cpl_Int_State%v, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_StateGet(impGEFS, 'tracers', Bundle, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Retrieve Bundle from impGEFS in GEFS Cpl.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Retrieving Bundle from impGEFS in GEFS Cpl, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 DO j = 2, 5
     DO i = 1, Cpl_Int_State%ntrac
         CALL ESMF_FieldBundleGet(Bundle,                    &
                                  trim(Cpl_Int_State%vname(i, j)), &
                                  Field,                     &
                                  rc = rc1)

         IF(ESMF_LogMsgFoundAllocError(rc1, "Retrieve Field from Bundle in GEFS Cpl.")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Retrieving Field from Bundle in GEFS Cpl, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

         NULLIFY(Cpl_Int_State%work2)
         CALL ESMF_FieldGet(Field, FArray = Cpl_Int_State%work2, localDE = 0, rc = rc1)

         IF(ESMF_LogMsgFoundAllocError(rc1, "Retrieve FArray from Field in GEFS Cpl.")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Retrieving FArray from Field in GEFS Cpl, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

         SELECT CASE(j)
             CASE(2) 
                 SELECT CASE(i)
                     CASE(1)
                         name = 'q'
                     CASE(2)
                         name = 'OZ'
                     CASE(3)
                         name = 'CLW'
                 END SELECT
                 CALL DistributeForStep1(Cpl_Int_State%tracer  (:, :, i), name, Cpl_Int_State, rc1)
             CASE(3) 
                 SELECT CASE(i)
                     CASE(1)
                         name = 'q(-DT)'
                     CASE(2)
                         name = 'OZ(-DT)'
                     CASE(3)
                         name = 'CLW(-DT)'
                 END SELECT
                 CALL DistributeForStep1(Cpl_Int_State%tracerm (:, :, i), name, Cpl_Int_State, rc1)
             CASE(4) 
                 SELECT CASE(i)
                     CASE(1)
                         name = 'q(-6H)'
                     CASE(2)
                         name = 'OZ(-6H)'
                     CASE(3)
                         name = 'CLW(-6H)'
                 END SELECT
                 CALL DistributeForStep1(Cpl_Int_State%tracer6 (:, :, i), name, Cpl_Int_State, rc1)
             CASE(5) 
                 SELECT CASE(i)
                     CASE(1)
                         name = 'q(-6H-DT)'
                     CASE(2)
                         name = 'OZ(-6H-DT)'
                     CASE(3)
                         name = 'CLW(-6H-DT)'
                 END SELECT
                 CALL DistributeForStep1(Cpl_Int_State%tracer6m(:, :, i), name, Cpl_Int_State, rc1)
         END SELECT

         IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- tracer.")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- tracer, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF
     END DO
 END DO

 CALL GetF90ArrayFromState(impGEFS, 'psm', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- PS(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- PS(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%psm, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'tm', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- T(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- T(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%tm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'um', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- U(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- U(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%um, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'vm', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- V(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- V(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-DT)"
 CALL DistributeForStep1(Cpl_Int_State%vm, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'pps6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- PS(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- PS(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps6, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'tt6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- T(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- T(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%t6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'uu6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- U(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- U(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%u6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'vv6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- V(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- V(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-6H)"
 CALL DistributeForStep1(Cpl_Int_State%v6, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-6H).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-6H), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'psm6', Cpl_Int_State%work1, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- PS(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- PS(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL DistributeForStep1_1(Cpl_Int_State%ps6m, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- PS(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- PS(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'tm6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- T(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- T(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "T(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%t6m,  name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- T(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- T(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'um6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- U(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- U(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "U(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%u6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- U(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- U(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL GetF90ArrayFromState(impGEFS, 'vm6', Cpl_Int_State%work2, 0, rc = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Get From impGEFS to the Cpl Internal State -- V(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting From impGEFS to the Cpl Internal State -- V(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 name = "V(-6H-DT)"
 CALL DistributeForStep1(Cpl_Int_State%v6m, name, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, "Distribute for the Step1 Scheme -- V(-6H-DT).")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Distributing for the Step1 Scheme -- V(-6H-DT), rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: GEFS_Cpl_ESMFImportState2InternalState"
 ELSE
     PRINT*, "FAIL: GEFS_Cpl_ESMFImportState2InternalState"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE GEFS_Cpl_ESMFImportState2InternalState





 SUBROUTINE GEFS_Cpl_InternalState2ESMFExportState(impGEFS, expGEFS, Cpl_Int_State, rc)
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_State),             INTENT(inout) :: impGEFS
 TYPE(ESMF_State),             INTENT(inout) :: expGEFS
 TYPE(GEFS_Cpl_InternalState), INTENT(inout) :: Cpl_Int_State

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER, OPTIONAL,            INTENT(out)   :: rc

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_VM)                               :: vm
 TYPE(ESMF_FieldBundle)                      :: Bundle
 TYPE(ESMF_Field)                            :: ESMFField
 INTEGER                                     :: arraysize_1
 INTEGER                                     :: arraysize_2
 INTEGER                                     :: arraysize_3
 INTEGER                                     :: arraysize_4
 INTEGER                                     :: arraysize_max
 INTEGER                                     :: arraysize_1_max
 INTEGER                                     :: i, j, k, l
 INTEGER                                     :: rc1
 INTEGER                                     :: rcfinal
 INTEGER, DIMENSION(:), POINTER              :: arraysize_1_pointer
 INTEGER, DIMENSION(:), POINTER              :: arraysize_1_gather
 LOGICAL                                     :: first

 REAL(KIND=KIND_EVOD), DIMENSION(:, :, :), POINTER :: t_wk
 REAL(KIND=KIND_EVOD), DIMENSION(:, :, :), POINTER :: FArr3D

 SAVE first

 DATA first/.true./

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

 IF(first) THEN
     ALLOCATE(arraysize_1_pointer(1))
     ALLOCATE(arraysize_1_gather (Cpl_Int_State%nodes))

     CALL ESMF_VMGetGlobal(vm, rc = rc1)

! Get the ESMF grid.
!-------------------
     CALL ESMF_StateGet(impGEFS, 'tt', ESMFField, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF Field")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF Field, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL ESMF_FieldGet(ESMFField, grid = mgrid, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get ESMF mgrid")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting the ESMF mgrid, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Get the array sizes.
!---------------------
     CALL ESMF_FieldGet(ESMFField, FArray = t_wk, localDE = 0, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Get t for its Size")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Getting t for its Size, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     arraysize_1 = SIZE(t_wk, 1)
     arraysize_2 = SIZE(t_wk, 2)
     arraysize_3 = SIZE(t_wk, 3)
     arraysize_4 = arraysize_1 * arraysize_2

     Cpl_Int_State%lonsize = Cpl_Int_State%lonf
     Cpl_Int_State%latsize = arraysize_1 / Cpl_Int_State%lonsize

     arraysize_1_pointer(1)    = arraysize_4 
     Cpl_Int_State%arraysize_1 = arraysize_1 
     Cpl_Int_State%arraysize_2 = arraysize_2 
     Cpl_Int_State%arraysize_3 = arraysize_3 
     Cpl_Int_State%arraysize_4 = arraysize_4 
     NULLIFY(t_wk)

! Allocate the writing out arrays at the current time step.
!----------------------------------------------------------
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tw )) &
         ALLOCATE(Cpl_Int_State%tw(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%uw )) &
         ALLOCATE(Cpl_Int_State%uw(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vw )) &
         ALLOCATE(Cpl_Int_State%vw(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psw)) &
         ALLOCATE(Cpl_Int_State%psw(arraysize_1, arraysize_2))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracerw)) &
         ALLOCATE(Cpl_Int_State%tracerw(arraysize_1, arraysize_2, arraysize_3, Cpl_Int_State%ntrac))

! Allocate the writing out arrays at the last time step.
!-------------------------------------------------------
     IF(.NOT. ASSOCIATED(Cpl_Int_State%twm )) &
         ALLOCATE(Cpl_Int_State%twm(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%uwm )) &
         ALLOCATE(Cpl_Int_State%uwm(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vwm )) &
         ALLOCATE(Cpl_Int_State%vwm(arraysize_1, arraysize_2, arraysize_3))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%pswm)) &
         ALLOCATE(Cpl_Int_State%pswm(arraysize_1, arraysize_2))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracerwm)) &
         ALLOCATE(Cpl_Int_State%tracerwm(arraysize_1, arraysize_2, arraysize_3, Cpl_Int_State%ntrac))

     DO k = 1, arraysize_3
         DO j = 1, arraysize_2
             DO i = 1, arraysize_1
                 Cpl_Int_State%tw (i, j, k) = 0.0
                 Cpl_Int_State%uw (i, j, k) = 0.0
                 Cpl_Int_State%vw (i, j, k) = 0.0
                 Cpl_Int_State%twm(i, j, k) = 0.0
                 Cpl_Int_State%uwm(i, j, k) = 0.0
                 Cpl_Int_State%vwm(i, j, k) = 0.0
             END DO
         END DO
     END DO

     DO l = 1, Cpl_Int_State%ntrac
         DO k = 1, arraysize_3
             DO j = 1, arraysize_2
                 DO i = 1, arraysize_1
                     Cpl_Int_State%tracerw (i, j, k, l) = 0.0
                     Cpl_Int_State%tracerwm(i, j, k, l) = 0.0
                 END DO
             END DO
         END DO
     END DO

     DO j = 1, arraysize_2
         DO i = 1, arraysize_1
             Cpl_Int_State%psw (i, j) = 0.0
             Cpl_Int_State%pswm(i, j) = 0.0
         END DO
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
     arraysize_max = arraysize_1_max * arraysize_3

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

     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracer  )) &
         ALLOCATE(Cpl_Int_State%tracer  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracerm )) &
         ALLOCATE(Cpl_Int_State%tracerm (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracer6 )) &
         ALLOCATE(Cpl_Int_State%tracer6 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracer6m)) &
         ALLOCATE(Cpl_Int_State%tracer6m(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%ps_step1))  &
         ALLOCATE(Cpl_Int_State%ps_step1 (Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%psm_step1)) &
         ALLOCATE(Cpl_Int_State%psm_step1(Cpl_Int_State%ARRAY_ONE_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%t_step1  )) &
         ALLOCATE(Cpl_Int_State%t_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tm_step1 )) &
         ALLOCATE(Cpl_Int_State%tm_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%u_step1  )) &
         ALLOCATE(Cpl_Int_State%u_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%um_step1 )) &
         ALLOCATE(Cpl_Int_State%um_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%v_step1  )) &
         ALLOCATE(Cpl_Int_State%v_step1  (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vm_step1 )) &
         ALLOCATE(Cpl_Int_State%vm_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracer_step1 )) &
         ALLOCATE(Cpl_Int_State%tracer_step1 (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracerm_step1)) &
         ALLOCATE(Cpl_Int_State%tracerm_step1(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%Total_member, Cpl_Int_State%ntrac))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%t_mean )) ALLOCATE(Cpl_Int_State%t_mean (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%u_mean )) ALLOCATE(Cpl_Int_State%u_mean (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%v_mean )) ALLOCATE(Cpl_Int_State%v_mean (Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tm_mean)) ALLOCATE(Cpl_Int_State%tm_mean(Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%um_mean)) ALLOCATE(Cpl_Int_State%um_mean(Cpl_Int_State%ARRAY_TOT_SIZ4))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%vm_mean)) ALLOCATE(Cpl_Int_State%vm_mean(Cpl_Int_State%ARRAY_TOT_SIZ4))

     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracer_mean )) &
         ALLOCATE(Cpl_Int_State%tracer_mean (Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%ntrac))
     IF(.NOT. ASSOCIATED(Cpl_Int_State%tracerm_mean)) &
         ALLOCATE(Cpl_Int_State%tracerm_mean(Cpl_Int_State%ARRAY_TOT_SIZ4, Cpl_Int_State%ntrac))

! Start to link the current time step arrays to the ESMF export state.
!---------------------------------------------------------------------
     CALL AddF90ArrayToState(expGEFS, mgrid, "tt",   Cpl_Int_State%tw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add t to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding t to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "uu",   Cpl_Int_State%uw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add u to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding u to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "vv",   Cpl_Int_State%vw,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add v to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding v to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "pps",  Cpl_Int_State%psw,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add ps to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding ps to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     Bundle = ESMF_FieldBundleCreate(name = 'tracers_GEFS', grid = mgrid, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Create Empty Fieldbundle")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Hapened When Creating Empty Fieldbundle, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     DO k = 1, Cpl_Int_State%ntrac
         NULLIFY(FArr3D)
         FArr3D => Cpl_Int_State%tracerw (:, :, :, k)
         ESMFField = ESMF_FieldCreate(mgrid, FArr3D, &
             name = trim(Cpl_Int_State%vname(k, 2)), rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Create ESMF Field")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Creating ESMF Field, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

         CALL ESMF_FieldBundleAdd(Bundle, ESMFField, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add the ESMF Field into Bundle")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding the ESMF Field into Bundle, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

         NULLIFY(FArr3D)
         FArr3D => Cpl_Int_State%tracerwm(:, :, :, k)
         ESMFField = ESMF_FieldCreate(mgrid, FArr3D, &
             name = trim(Cpl_Int_State%vname(k, 3)), rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Create ESMF Field")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Creating ESMF Field, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

         CALL ESMF_FieldBundleAdd(Bundle, ESMFField, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add the ESMF Field into Bundle")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding the ESMF Field into Bundle, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF
     END DO

     CALL ESMF_StateAdd(expGEFS, Bundle, rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "ESMF State Adds the Bundle")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When ESMF State Adds the Bundle, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Start to link the last time step arrays to the ESMF export state.
!------------------------------------------------------------------
     CALL AddF90ArrayToState(expGEFS, mgrid, "tm",   Cpl_Int_State%twm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add tm to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding tm to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "um",   Cpl_Int_State%uwm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add um to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding um to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "vm",   Cpl_Int_State%vwm,   rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add vm to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding vm to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

     CALL AddF90ArrayToState(expGEFS, mgrid, "psm",  Cpl_Int_State%pswm,  rc = rc1)

         IF(ESMF_LogMsgFoundError(rc1, "Add psm to expGEFS")) THEN
             rcfinal = ESMF_FAILURE
             PRINT*, 'Error Happened When Adding psm to expGEFS, rc = ', rc1
             rc1 = ESMF_SUCCESS
         END IF

! Finish the first = .true. run.
!-------------------------------
     first = .false.
 END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: GEFS_Cpl_InternalState2ESMFExportState"
 ELSE
     PRINT*, "FAIL: GEFS_Cpl_InternalState2ESMFExportState"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE GEFS_Cpl_InternalState2ESMFExportState

 END MODULE GEFS_CplState_ESMFMod
