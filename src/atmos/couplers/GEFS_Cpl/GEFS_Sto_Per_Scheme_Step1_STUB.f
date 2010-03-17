 SUBROUTINE GEFS_Sto_Per_Scheme_Step1(cst)

! March 2009 Weiyu Yang  Modified for the NEMS model.
!----------------------------------------------------

! This subroutine is used to compute the first step of the 
! stochastic perturbation scheme, in which X_i_dot = T_i + S_i 
! and S_i ~ SUM(W_i,j P_j.
!-------------------------------------------------------------


! In the current code, assume each ensemble member uses the same 
! number of the processors.
!---------------------------------------------------------------
 USE GEFS_Cpl_InternalState_ESMFMod

 IMPLICIT none

 TYPE(GEFS_Cpl_InternalState),                         INTENT(inout) :: cst

! Working arrays.
!----------------
 INTEGER                                             :: tb
 INTEGER                                             :: i, j, k
 REAL(KIND=KIND_EVOD), DIMENSION(cst%ARRAY_ONE_SIZ4) :: pstm, psmtm
 REAL(KIND=KIND_EVOD), DIMENSION(cst%ARRAY_TOT_SIZ4) :: ttm, tmtm, utm, umtm, vtm, vmtm,      &
                                                        qtm, qmtm, oztm, ozmtm, clwtm, clwmtm
 END SUBROUTINE GEFS_Sto_Per_Scheme_Step1





 SUBROUTINE GEFS_Sto_Per_Scheme_Step1_2(impGEFS, cst, rc)

! This subroutine is used to compute the first step of the 
! stochastic perturbation scheme, in which X_i_dot = T_i + S_i 
! and S_i ~ SUM(W_i,j P_j.
!-------------------------------------------------------------

! In the current code, assume each ensemble member uses the same 
! number of the processors.
!---------------------------------------------------------------

 USE ESMF_Mod
 USE GEFS_Cpl_InternalState_ESMFMod
 USE Lib_ESMFStateAddGetMod

 IMPLICIT none

 TYPE(ESMF_State),             INTENT(inout) :: impGEFS
 TYPE(GEFS_Cpl_InternalState), INTENT(inout) :: cst
 INTEGER,                      INTENT(out)   :: rc

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 INTEGER                                     :: i, j
 INTEGER                                     :: rc1

 END SUBROUTINE GEFS_Sto_Per_Scheme_Step1_2
