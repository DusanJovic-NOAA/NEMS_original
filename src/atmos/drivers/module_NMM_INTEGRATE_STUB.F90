!-----------------------------------------------------------------------
!
      MODULE MODULE_NMM_INTEGRATE
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
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: ALARM_HISTORY,ALARM_RESTART,NMM_INTEGRATE
!
!-----------------------------------------------------------------------
      INCLUDE '../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      LOGICAL,SAVE :: RESTARTED_RUN_FIRST=.TRUE.
!
      TYPE(ESMF_Alarm),SAVE :: ALARM_HISTORY,ALARM_RESTART
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_INTEGRATE(CLOCK_DIRECTION                          &
                              ,ATM_GRID_COMP                            &
                              ,IMP_STATE_ATM                            &
                              ,EXP_STATE_ATM                            &
                              ,CLOCK_INTEGRATE                          &
                              ,CURRTIME                                 &
                              ,STARTTIME                                &
                              ,TIMESTEP                                 &
                              ,NTIMESTEP                                &
                              ,DT                                       &
                              ,INTERVAL_CLOCKTIME                       &
                              ,INTERVAL_HISTORY                         &
                              ,INTERVAL_RESTART                         &
                              ,FILTER_METHOD                            &
                              ,HALFDFIINTVAL                            &
                              ,HALFDFITIME                              &
                              ,NDFISTEP                                 &
                              ,NPE_PRINT                                &
                              ,RESTARTED_RUN                            &
                              ,RST_OUT_00                               &
                              ,I_AM_A_FCST_TASK                         &
                              ,NESTING                                  &
                              ,I_AM_A_NEST                              &
                              ,MY_DOMAIN_ID                             &
                              ,COMM_TO_MY_PARENT                        &
                              ,NUM_CHILDREN                             &
                              ,PARENT_CHILD_CPL                         &
                              ,IMP_STATE_CPL_NEST                       &
                              ,EXP_STATE_CPL_NEST                       &
                              ,PAR_CHI_TIME_RATIO                       &
                              ,MYPE)
!-----------------------------------------------------------------------
!
!-----------------
!*** Arguments IN
!-----------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_TO_MY_PARENT                &  !<-- MPI Communicator to parent of this domain
                                      ,FILTER_METHOD                    &  !<-- The type of digital filtering desired
                                      ,MYPE                             &  !<-- MPI task rank
                                      ,NPE_PRINT                        &  !<-- Task to print clocktimes
                                      ,NUM_CHILDREN                        !<-- # of children on this domain
!
      REAL(kind=KFPT),INTENT(IN) :: DT                                     !<-- Fundamental timestep of this domain (REAL) (s)
!
      LOGICAL(kind=KLOG),INTENT(IN) :: NESTING                          &  !<-- Are there any nested domains?
                                      ,RESTARTED_RUN                    &  !<-- Is this a restarted run?
                                      ,RST_OUT_00                          !<-- Shall we write 00h history in restarted run?

!
      CHARACTER(7),INTENT(IN) :: CLOCK_DIRECTION                           !<-- The direction of time in the Clock
!
      TYPE(ESMF_Logical),INTENT(IN) :: I_AM_A_FCST_TASK                 &  !<-- Am I in a forecast task?
                                      ,I_AM_A_NEST                         !<-- Am I in a nested domain?
!
      TYPE(ESMF_Time),INTENT(IN) :: STARTTIME                              !<-- The clock's start time
!
      TYPE(ESMF_TimeInterval),INTENT(IN)  :: TIMESTEP                      !<-- Fundamental timestep of this domain (ESMF) (s)
!
!--------------------
!*** Arguments INOUT
!--------------------
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NTIMESTEP                        !<-- The timestep count
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM component of this domain
!
      TYPE(ESMF_Time),INTENT(INOUT) :: CURRTIME                            !<-- The clock's current time
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_INTEGRATE                    !<-- This ATM Component's ESMF Clock
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_ATM                   &  !<-- Import state of this ATM component
                                       ,EXP_STATE_ATM                      !<-- Export state of this ATM component
!
      TYPE(ESMF_State),INTENT(INOUT),OPTIONAL:: IMP_STATE_CPL_NEST      &
                                               ,EXP_STATE_CPL_NEST
!
!------------------------
!***  Optional Arguments
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN),OPTIONAL :: MY_DOMAIN_ID            &  !<-- The ID of this domain (ATM component)
                                               ,NDFISTEP                &
                                               ,PAR_CHI_TIME_RATIO         !<-- Ratio of parent's timestep to this domain's
!
      TYPE(ESMF_Time),INTENT(IN),OPTIONAL :: HALFDFITIME
!
      TYPE(ESMF_TimeInterval),INTENT(IN),OPTIONAL :: HALFDFIINTVAL      &
                                                    ,INTERVAL_CLOCKTIME &  !<-- Time interval between clocktime prints
                                                    ,INTERVAL_HISTORY   &  !<-- Time interval between history output
                                                    ,INTERVAL_RESTART      !<-- Time interval between restart output
!
      TYPE(ESMF_CplComp),INTENT(IN),OPTIONAL :: PARENT_CHILD_CPL           !<-- Coupler component for parent-child/nest exchange
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      END SUBROUTINE NMM_INTEGRATE
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NMM_INTEGRATE
!
!-----------------------------------------------------------------------
