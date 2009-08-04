!-----------------------------------------------------------------------
!
      MODULE MODULE_NMM_FWD_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.)
!***  IN MODULE_MAIN_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!
      PUBLIC :: NMM_FWD_INTEGRATE       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!-----------------------------------------------------------------------
      INCLUDE '../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_FWD_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                            &
				  ,HALFDFIINTVAL                        &
                                  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
				  ,MYPE                                 &
				  ,NUM_TRACERS_MET                      &
				  ,NUM_TRACERS_CHEM                     &
				  ,NTIMESTEP                            &
				  ,NPE_PRINT)                            
!
!-----------------------------------------------------------------------
!
      USE MODULE_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP             !<-- The ATM gridded component
!
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
!
      TYPE(ESMF_Clock),INTENT(INOUT) 	     :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
      INTEGER(KIND=KINT),INTENT(IN)          :: NPE_PRINT
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE                    & !<-- MPI task rank
                                               ,NUM_TRACERS_MET         & !<-- # of meteorological tracer variables
                                               ,NUM_TRACERS_CHEM          !<-- # of chemistry tracer variables
!
      INTEGER(KIND=KINT),INTENT(INOUT)       :: NTIMESTEP                 !<-- The current forecast timestep
!
      TYPE(ESMF_Time),INTENT(INOUT)          :: CURRTIME                  !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: HALFDFIINTVAL
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: TIMESTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: FILTER_METHOD
      
      END SUBROUTINE NMM_FWD_INTEGRATE
      
!-----------------------------------------------------------------------

      END MODULE MODULE_NMM_FWD_INTEGRATE
