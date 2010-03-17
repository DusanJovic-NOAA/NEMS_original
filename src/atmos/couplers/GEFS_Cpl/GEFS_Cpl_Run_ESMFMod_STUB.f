!
! !MODULE: GEFS_Cpl_Run_ESMFMod --- Run module of the ESMF grided
!                                   component of the GFS ensemble coupler.
!
! !DESCRIPTION: GEFS coupler run module.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!  March    2009     Weiyu Yang, modified for the NEMS model.
!
!
! !INTERFACE:
!
 MODULE GEFS_Cpl_Run_ESMFMod
!
!!USES:
!
 USE ESMF_Mod
 USE GEFS_Cpl_InternalState_ESMFMod

 INCLUDE 'mpif.h'

 IMPLICIT none

 CONTAINS

 SUBROUTINE GEFS_Cpl_Run(impState, clock, cst, rc)

! This subroutine is used to create the new initial conditions for the 
! next GFS ensemble forecast run from the last GFS run inputs and outputs.
! The new created initial condition fields need to be put back. 
!-------------------------------------------------------------------------

 TYPE(ESMF_State),                      INTENT(inout) :: impState
 TYPE(ESMF_Clock),                      INTENT(inout) :: clock
 TYPE(GEFS_Cpl_InternalState), POINTER, INTENT(inout) :: cst
 INTEGER, OPTIONAL,                     INTENT(out)   :: rc



! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_Time)                                      :: currTime
 TYPE(ESMF_TimeInterval)                              :: runDuration
 INTEGER                                              :: hh
 INTEGER                                              :: rc1
 INTEGER                                              :: rcfinal
     
 INTEGER                                              :: i, j, k, l
 INTEGER                                              :: year, month, day, hour
 INTEGER                                              :: Jul_Day
 CHARACTER(ESMF_MAXSTR)                               :: name

 END SUBROUTINE GEFS_Cpl_Run





 SUBROUTINE Cal_Jul_Day(Jul_Day,year,month,day,hour,hh)
 INTEGER Jul_Day,year,month,day,hour,hh,m 
!year,month,day and hour: for the time at which the forecast is valid. 
!hh is the forecast lead time in hours, at NEXT stop of integration.
 INTEGER dd(12)
 DATA dd/31,28,31,30,31,30,31,31,30,31,30,31/

 END SUBROUTINE Cal_Jul_Day

 END MODULE GEFS_Cpl_Run_ESMFMod
