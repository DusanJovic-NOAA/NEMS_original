 SUBROUTINE GEFS_bcst_global(var, peid, rc)

!----------------------------------------------------------------------
! SUBROUTINE bcst_global
!
! This subroutine broadcasts the inputted variable var to all PEs and all
! ensemble members.  And output contains all related index parameters.
!
! DESCRIPTION: VM Broadcast tool software.
!
! REVISION HISTORY:
!
!  Setpember 2007     Weiyu Yang Initial code.
!
!
! INTERFACE:
!   
!   var    -- inputted single variable.
!   peid   -- PE ID of var.
!   vm     -- the global ESMF VM.
!
 USE ESMF_Mod
 USE machine

 REAL(KIND = kind_evod)                    :: var 
 REAL(ESMF_KIND_R8), DIMENSION(:), POINTER :: var_work 
 INTEGER                                   :: peid
 TYPE(ESMF_VM)                             :: vm
 INTEGER                                   :: rc

 END SUBROUTINE GEFS_bcst_global





 SUBROUTINE GEFS_bcst_global_i4(var, peid, rc)

!----------------------------------------------------------------------
! SUBROUTINE bcst_global_i4
!
! This subroutine broadcasts the inputted variable var to all PEs and all
! ensemble members.  And output contains all related index parameters.
!
! DESCRIPTION: VM Broadcast tool software.
!
! REVISION HISTORY:
!
!  Setpember 2007     Weiyu Yang Initial code.
!
!
! INTERFACE:
!   
!   var    -- inputted single variable.
!   peid   -- PE ID of var.
!   vm     -- the global ESMF VM.
!
 USE ESMF_Mod
 USE machine

 INTEGER                                   :: var 
 REAL(ESMF_KIND_I4), DIMENSION(:), POINTER :: var_work 
 INTEGER                                   :: peid
 TYPE(ESMF_VM)                             :: vm
 INTEGER                                   :: rc

 END SUBROUTINE GEFS_bcst_global_i4
