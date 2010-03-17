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

 TYPE(ESMF_Grid), SAVE                       :: grid1
 TYPE(ESMF_Grid), SAVE                       :: grid2

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
 CHARACTER(ESMF_MAXSTR)                      :: name
 INTEGER                                     :: rc1
 INTEGER                                     :: rcfinal

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
 TYPE(ESMF_Array)                            :: ESMFArray1
 TYPE(ESMF_Array)                            :: ESMFArray2
 TYPE(ESMF_DistGrid)                         :: DistGrid1 
 TYPE(ESMF_DistGrid)                         :: DistGrid2   
 INTEGER                                     :: arraysize_1
 INTEGER                                     :: arraysize_2
 INTEGER                                     :: arraysize_max
 INTEGER                                     :: arraysize_1_max
 INTEGER                                     :: i, j
 INTEGER                                     :: rc1
 INTEGER                                     :: rcfinal
 INTEGER, DIMENSION(:), POINTER              :: arraysize_1_pointer
 INTEGER, DIMENSION(:), POINTER              :: arraysize_1_gather
 LOGICAL                                     :: first

 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: t_wk

 SAVE first

 DATA first/.true./

 END SUBROUTINE GEFS_Cpl_InternalState2ESMFExportState

 END MODULE GEFS_CplState_ESMFMod
