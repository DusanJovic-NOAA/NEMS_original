 SUBROUTINE GEFS_Sto_Per_Scheme_Step2(Int_State, USES1, GRIDUSE, &
     Jul_Day, slat1, slat2, rc)

!DHOU, 10/17/2007  Added arguments 7 and 8, i.e. KEMAX(nreg) and nreg for regional rescaling
!DHOU, 10/17/2007  Added argument 6 i.e. RS_GLOBAL for global rescaling factor 
!DHOU  09/11/2007  Added Arguments
!USES1, 0=Zero-out, 1=no-change, 2=replaced with S2, 3=replace S1 and s2 with 0.5(s1+S2).
!           For arrays related to state 1,zem to qm.
!GRIDUSE, 1=S1, 2=S2; Convert the state S1/S2 from spectral to grid and  and back to spec.

! This subroutine is used to compute the second step of the
! stochastic perturbation scheme, in which it carries out the spectral 
! transform computation into the Gaussian grid space arrays, then
! computes the second step of the stochastic perturbation scheme that
! considering the local weighting influences.
!---------------------------------------------------------------------

! !REVISION HISTORY:
!
!  May 2007       Weiyu Yang Initial code for wave-grid conversion for model state .
!  Nov 2007       Dingchen Hou ddopted the code for global/regional rescaling as well 
!                 as conversion, for model state or its perturbation.    
!  Mar 2009       Weiyu Yang Modified for the NEMS model.
!--------------------------------------------------------

 USE ESMF_Mod
 USE GEFS_Cpl_InternalState_ESMFMod
 USE machine,  ONLY: kind_evod, kind_phys, kind_rad

 INCLUDE 'mpif.h'

 IMPLICIT none

 TYPE(GEFS_Cpl_InternalState), INTENT(inout) :: Int_State
 INTEGER,                      INTENT(out)   :: rc
 INTEGER                                     :: USES1, GRIDUSE
 INTEGER,                      INTENT(in)    :: Jul_Day
 INTEGER                                     :: ireg, imem, k500
 REAL(KIND = kind_evod)                      :: slat1, slat2

 INTEGER                                     :: i, j, k, i1
 INTEGER                                     :: lat
 INTEGER                                     :: lon_lat
 INTEGER                                     :: jlonsize

 REAL(KIND = kind_evod)                      :: keavg
 REAL(KIND = kind_evod), PARAMETER           :: pi = 3.1415926535897931

 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: t
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: q
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: oz
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: clw
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: u
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: v

 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: ps
 REAL(KIND=KIND_EVOD), DIMENSION(:, :), POINTER :: KER

 REAL(KIND = kind_evod)                      :: tm
 REAL(KIND = kind_evod)                      :: qm,   ozm,  clwm
 REAL(KIND = kind_evod)                      :: um,   vm
 REAL(KIND = kind_evod)                      :: psm

 REAL(KIND = kind_evod)                      :: ts
 REAL(KIND = kind_evod)                      :: qs,   ozs,  clws 
 REAL(KIND = kind_evod)                      :: us,   vs
 REAL(KIND = kind_evod)                      :: pss


 INTEGER                                     :: rc1
 INTEGER                                     :: rcfinal

 END SUBROUTINE GEFS_Sto_Per_Scheme_Step2





 SUBROUTINE GET_SCALING_FACTORS(KE,nlat,nmember,slat,KEMAX,KER,factor1,nregion,slat1,slat2)
 USE machine,  ONLY: kind_evod
 INTEGER nlat,nmember,nregion
 REAL(KIND = kind_evod)                      :: KE(nlat,nmember)
 REAL(KIND = kind_evod)                      :: factor1(nregion,nmember) 
 REAL(KIND = kind_evod)                      :: slat(nlat),slat1,slat2 
 REAL(KIND = kind_evod)                      :: KEMAX(nregion) 
 REAL(KIND = kind_evod)                      :: KER(nregion,nmember) 
 REAL(KIND = kind_evod)                      :: WEIGHT(3) 
 REAL(KIND = kind_evod)                      :: coslat
 INTEGER i,j,k

 END
