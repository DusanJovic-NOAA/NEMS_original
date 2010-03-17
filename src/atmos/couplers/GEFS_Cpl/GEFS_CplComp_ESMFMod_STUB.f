!----------------------------------------------------------------------
! !MODULE: GEFS_CplComp_ESMFMod
!        --- ESMF coupler gridded component of the GFS Ensemble 
!            Forecast Operational system. 
!
! !DESCRIPTION: GFS coupler gridded component main module.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!                           
!
!  March 2007         Dingchen Hou added Stochatic Perturbation Combination Coefficient array.
!  January to November 2007      Dingchen and Weiyu Yang
!                     Added Broadcasting procedure and Global variables/arrays
!  November 2007      Dingchen, added minimum documentation, mainly for the arrays added during 2007
!  March    2009      Weiyu yang, modified for the NEMS model.

! !INTERFACE:
!
 MODULE GEFS_CplComp_ESMFMod
 
!!USES:
!------

! Define the ESMF internal state and all routines related to run 
! the GFS ensemble coupler grid component.
!---------------------------------------------------------------
 USE GEFS_Cpl_ESMFMod

 IMPLICIT none

#include "ESMF_LogMacros.inc"

 PRIVATE   ! By default data is private to this module
!
! !PUBLIC TYPES:
!---------------

 PUBLIC GEFS_CplCompSetServices

!EOP
!-------------------------------------------------------------------------


 CONTAINS


!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GEFS_CplCompSetServices --- Set services for GFS Ensemble 
!                                       Coupler Gridded Component.
! 
! !INTERFACE:
!
 SUBROUTINE GEFS_CplCompSetServices(CplGEFS, rc)

! !ARGUMENTS:
!------------

 TYPE(ESMF_CplComp),  INTENT(inout) :: CplGEFS ! gridded component
 INTEGER,             INTENT(out)   :: rc      ! return code
     
! !DESCRIPTION: Set services (register) for the GFS Ensemble Coupler
!               Grid Component.
!         
!EOP         
!----------------------------------------------------------------------
  
 INTEGER                            :: rc1     = ESMF_SUCCESS
 rc = ESMF_SUCCESS

! REGISTER SERVICES FOR THIS COMPONENT
! ------------------------------------

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETINIT,  Cpl_Initialize, &
                                 ESMF_SINGLEPHASE, rc1)

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETRUN,   Cpl_Run,        &
                                 ESMF_SINGLEPHASE, rc1)

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETFINAL, Cpl_Finalize,   &
                                 ESMF_SINGLEPHASE, rc1)

 END SUBROUTINE GEFS_CplCompSetServices





!----------------------------------------------------------------------
!BOP
! !ROUTINE:  Cpl_Initialize --- initialize routine to initialize 
!                               and set up the GFS ensemble coupler.
!
! !DESCRIPTION: This subroutine initializes the GFS ensemble coupler
!               before the main running routine.
!
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:
!

 SUBROUTINE Cpl_Initialize(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

!
! !INPUT/OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------------

 TYPE(ESMF_CplComp), INTENT(inout)     :: CplGEFS 
 TYPE(ESMF_State),   INTENT(inout)     :: impGEFS
 TYPE(ESMF_State),   INTENT(inout)     :: expGEFS
 TYPE(ESMF_Clock),   INTENT(inout)     :: clock

!
! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------

 INTEGER,             INTENT(out)       :: rcfinal

! !EOP
!------------------------------------------------------------------------- 
 
 END SUBROUTINE Cpl_Initialize





!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cpl_Run --- Main grid component routine to run the GFS 
!                       ensemble coupler.
!
! !DESCRIPTION: This subroutine will run the most part computations 
!               of the GFS ensemble coupler.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:
!

 SUBROUTINE Cpl_Run(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

!
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_CplComp), INTENT(inout)     :: CplGEFS   
 TYPE(ESMF_State),   INTENT(in)        :: impGEFS 
 
! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 TYPE(ESMF_Clock),   INTENT(inout)     :: clock
 TYPE(ESMF_State),   INTENT(inout)     :: expGEFS
 INTEGER,            INTENT(out)       :: rcfinal 
!
!EOP
!-------------------------------------------------------------------------

!
 END SUBROUTINE Cpl_Run





!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cpl_Finalize --- finalizing routine to finish the 
!                            GFS ensemble coupler.
!
! !DESCRIPTION: This subroutine will finish the GFS ensemble coupler
! !             and will release the memory space.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:

 SUBROUTINE Cpl_Finalize(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

!
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_CplComp), INTENT(inout)  :: CplGEFS
 TYPE(ESMF_State),   INTENT(inout)  :: impGEFS
 TYPE(ESMF_State),   INTENT(inout)  :: expGEFS
 TYPE(ESMF_Clock),   INTENT(inout)  :: clock

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER,            INTENT(out)    :: rcfinal


 END SUBROUTINE Cpl_Finalize

 END MODULE GEFS_CplComp_ESMFMod
