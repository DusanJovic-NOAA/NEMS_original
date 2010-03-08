 SUBROUTINE DistributeForStep1_1(Array, Cpl_Int_State, rc)

 USE machine, ONLY: KIND_EVOD
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 INCLUDE 'mpif.h'

 TYPE(GEFS_Cpl_InternalState),                                  INTENT(inout) :: Cpl_Int_State
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%ARRAY_ONE_SIZ3), INTENT(inout) :: Array
 INTEGER,                                                       INTENT(out)   :: rc

 INTEGER                                                                      :: i

 END SUBROUTINE DistributeForStep1_1





 SUBROUTINE DistributeForStep1(Array, Name, Cpl_Int_State, rc)

 USE machine, ONLY: KIND_EVOD
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 INCLUDE 'mpif.h'

 CHARACTER(ESMF_MAXSTR),                                        INTENT(in)    :: Name
 TYPE(GEFS_Cpl_InternalState),                                  INTENT(inout) :: Cpl_Int_State
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%ARRAY_TOT_SIZ3), INTENT(inout) :: Array
 INTEGER,                                                       INTENT(out)   :: rc

 INTEGER                                                                      :: i, j, k

 END SUBROUTINE DistributeForStep1





 SUBROUTINE DistributeBackFromStep1_1(Input, Output, Cpl_Int_State, rc)

 USE machine, ONLY: KIND_EVOD
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 INCLUDE 'mpif.h'

 TYPE(GEFS_Cpl_InternalState),                                  INTENT(inout) :: Cpl_Int_State
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%ARRAY_ONE_SIZ3), INTENT(in)    :: Input
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%arraysize_1),    INTENT(inout) :: Output
 INTEGER,                                                       INTENT(inout) :: rc

 INTEGER                                                                      :: rcfinal
 INTEGER                                                                      :: i

 END SUBROUTINE DistributeBackFromStep1_1





 SUBROUTINE DistributeBackFromStep1(Input, Output, Name, Cpl_Int_State, rc)

 USE machine, ONLY: KIND_EVOD
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 INCLUDE 'mpif.h'

 CHARACTER(ESMF_MAXSTR),                                                                INTENT(in)    :: Name
 TYPE(GEFS_Cpl_InternalState),                                                          INTENT(inout) :: Cpl_Int_State
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%ARRAY_ONE_SIZ3),                         INTENT(in)    :: Input
 REAL(KIND=KIND_EVOD), DIMENSION(Cpl_Int_State%arraysize_1, Cpl_Int_State%arraysize_2), INTENT(inout) :: Output
 INTEGER,                                                                               INTENT(out)   :: rc

 INTEGER                                                                                              :: rcfinal
 INTEGER                                                                                              :: i, j, k

 END SUBROUTINE DistributeBackFromStep1
