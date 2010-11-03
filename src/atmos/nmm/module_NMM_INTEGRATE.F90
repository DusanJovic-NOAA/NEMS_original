!-----------------------------------------------------------------------
!
      MODULE module_NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  This module holds the fundamental NMM integration runstream.
!***  It is called from subroutine NMM_RUN in module_NMM_GRID_COMP.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! PROGRAM HISTORY LOG:
!   2008-08     Colon - Moved NMM runstream from ATM_RUN into separate
!                       routines when adding digital filters.
!   2009-07-09  Black - Condense all three NMM integrate routines
!                       into one when when merging with nesting.
!   2010-03-24  Black - Revised for new structure.
!   2010-10-xx  Pyle  - Revised for digital filters.
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE module_CLOCKTIMES,ONLY: PRINT_CLOCKTIMES
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!!!   USE MODULE_DIGITAL_FILTER_NMM
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE
!
      USE module_WRITE_ROUTINES,ONLY: WRITE_ASYNC
!
      USE module_NESTING,ONLY: BOUNDARY_DATA_STATE_TO_STATE
!
      USE module_CONTROL,ONLY: TIMEF
!
      USE module_PARENT_CHILD_CPL_COMP, ONLY: NSTEP_CHILD_RECV
!
      USE module_INCLUDE
!
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: ALARM_HISTORY                                           &
               ,ALARM_RESTART                                           &
               ,NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!
      LOGICAL(kind=KLOG),SAVE :: RESTARTED_RUN_FIRST=.TRUE.
!
      CHARACTER(ESMF_MAXSTR) :: CWRT                                       !<-- Restart/History label
!
      TYPE(ESMF_Alarm),SAVE :: ALARM_HISTORY                            &  !<-- The ESMF Alarm for history output
                              ,ALARM_RESTART                               !<-- The ESMF Alarm for restart output
!
!-----------------------------------------------------------------------
!***  For determining clocktimes of various pieces of the Dynamics.
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0                                     &
                        ,atm_drv_run_1                                  &
                        ,atm_drv_run_2                                  &
                        ,atm_drv_run_3                                  &
                        ,atm_drv_run_cpl1                               &
                        ,atm_drv_run_cpl2                               &
                        ,cpl1_recv_tim                                  &
                        ,cpl2_send_tim                                  &
                        ,cpl2_comp_tim                                  &
                        ,cpl2_wait_tim                                  &
                        ,phase1_tim
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
                              ,DOMAIN_GRID_COMP                         &
                              ,IMP_STATE_DOMAIN                         &
                              ,EXP_STATE_DOMAIN                         &
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
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!*** The following USEs are needed for NMM-B time series output.
!-----------------------------------------------------------------------
!
      USE MODULE_DYNAMICS_INTERNAL_STATE, ONLY: DYNAMICS_INTERNAL_STATE &
                                               ,WRAP_DYN_INT_STATE
!
      USE MODULE_PHYSICS_INTERNAL_STATE, ONLY: PHYSICS_INTERNAL_STATE   &
                                             , WRAP_PHY_INT_STATE
!
      USE MODULE_TIMESERIES
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
      CHARACTER(8),INTENT(IN) :: CLOCK_DIRECTION                           !<-- The direction of time in the Clock
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
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(ESMF_Time),INTENT(INOUT) :: CURRTIME                            !<-- The clock's current time
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_INTEGRATE                    !<-- This DOMAIN Component's ESMF Clock
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_DOMAIN                &  !<-- Import state of this DOMAIN component 
                                       ,EXP_STATE_DOMAIN                   !<-- Export state of this DOMAIN component
!
      TYPE(ESMF_State),INTENT(INOUT),OPTIONAL:: IMP_STATE_CPL_NEST      &
                                               ,EXP_STATE_CPL_NEST
!
!------------------------
!***  Optional Arguments
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN),OPTIONAL :: MY_DOMAIN_ID            &  !<-- The ID of this domain 
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
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: YY,MM,DD,H,M,S
      INTEGER(kind=KINT) :: I,KOUNT_STEPS,N
      INTEGER(kind=KINT) :: IERR,RC,RC_INTEG
      INTEGER(kind=KINT), ALLOCATABLE :: LOC_PAR_CHILD_TIME_RATIO(:)
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF_INT)
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT
!
      TYPE(ESMF_Time) :: ALARM_HISTORY_RING                             &
                        ,ALARM_RESTART_RING                             &
                        ,ALARM_CLOCKTIME_RING
!
      TYPE(ESMF_Time) :: ADJTIME_HISTORY                                &
                        ,ADJTIME_RESTART                                &
                        ,ADJTIME_CLOCKTIME
!
      TYPE(ESMF_TimeInterval) :: TIMESTEP_FILTER                           !<-- Dynamics timestep during filter (s) (ESMF)
!
      TYPE(ESMF_Alarm) :: ALARM_CLOCKTIME                                  !<-- The ESMF Alarm for clocktime prints
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP
!
      TYPE(WRAP_DYN_INT_STATE) :: WRAP_DYN
!
      TYPE(WRAP_PHY_INT_STATE) :: WRAP_PHY
!
      TYPE(DYNAMICS_INTERNAL_STATE),POINTER :: DYN_INT_STATE
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: PHY_INT_STATE
!
      LOGICAL, SAVE :: TS_INITIALIZED = .FALSE.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      phase1_tim      =0
      atm_drv_run_1   =0.
      atm_drv_run_2   =0.
      atm_drv_run_3   =0.
      atm_drv_run_cpl1=0.
      atm_drv_run_cpl2=0.
!
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_INTEG=ESMF_SUCCESS
!
      FMT='(I2.2)'
      WRITE(INT_TO_CHAR,FMT)MY_DOMAIN_ID
!
      KOUNT_STEPS=0
!
!-----------------------------------------------------------------------
!***  For normal forecast integration set the Alarm ring times
!***  while accounting for restarts and digital filtering.
!-----------------------------------------------------------------------
!
      IF(FILTER_METHOD==0)THEN
!
        CALL RESET_ALARMS
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Before beginning the integration of the DOMAIN component,
!***  extract the internal state which will be needed for
!***  initial writing of history/restart files.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get DOMAIN Internal State in NMM_INTEGRATE"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Forecast tasks extract the Dynamics and Physics internal states
!***  needed for timeseries output.
!-----------------------------------------------------------------------
!
      IF(MYPE<domain_int_state%NUM_PES_FCST)THEN
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%DYN_GRID_COMP &  !<-- The Dynamics component
                                          ,WRAP_DYN                       &  !<-- The F90 wrap of the Dynamics internal state
                                          ,RC)
!
        DYN_INT_STATE => wrap_dyn%INT_STATE
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%PHY_GRID_COMP &  !<-- The Physics component
                                          ,WRAP_PHY                       &  !<-- The F90 wrap of the Physics internal state
                                          ,RC)
!
        PHY_INT_STATE => wrap_phy%INT_STATE
!
      END IF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE INTEGRATION TIME LOOP OF THE ATMOSPHERE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      timeloop_drv: DO WHILE (.NOT.ESMF_ClockIsStopTime(CLOCK_INTEGRATE &
                                                       ,RC) )

!
!-----------------------------------------------------------------------
!
        btim0=timef()
!
!-----------------------------------------------------------------------
!***  Call the 1st Phase of the Parent_Child coupler where children
!***  will recv from their parents.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_NEST==ESMF_TRUE.AND.I_AM_A_FCST_TASK==ESMF_TRUE)THEN
!
!!!       IF(ESMF_AlarmIsRinging(alarm=ALARM_RECV_FROM_PARENT           &  !<-- Alarm to alert child that it must recv from parent
!         IF((ESMF_AlarmIsRinging(alarm=ALARM_RECV_FROM_PARENT          &  !<-- Alarm to alert child that it must recv from parent
!                               ,rc   =RC)                              &
!            .or.ntimestep==0)  &        !<-- bandaid
          IF(MOD(KOUNT_STEPS,PAR_CHI_TIME_RATIO)==0                     &
                            .AND.                                       &
             COMM_TO_MY_PARENT/=-999)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Call Coupler: Children Recv from Parents"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL           &  !<-- The Nesting coupler component
                                ,importState=IMP_STATE_CPL_NEST         &  !<-- The Nesting coupler import state
                                ,exportState=EXP_STATE_CPL_NEST         &  !<-- The Nesting coupler export state
                                ,clock      =CLOCK_INTEGRATE            &  !<-- The DOMAIN Clock
                                ,phase      =1                          &  !<-- The phase (subroutine) of the Coupler to execute
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The new nest boundary data must be moved from the Parent-Child
!***  coupler into the nests' DOMAIN components.
!-----------------------------------------------------------------------
!

            CALL BOUNDARY_DATA_STATE_TO_STATE(state_in =EXP_STATE_CPL_NEST &  !<-- The Nesting coupler export state
                                             ,state_out=IMP_STATE_DOMAIN)     !<-- The DOMAIN import state

!
!-----------------------------------------------------------------------
!
          ENDIF
!
        ENDIF
!
        atm_drv_run_cpl1=atm_drv_run_cpl1+(timef()-btim0)
!
!-----------------------------------------------------------------------
!***  If filtering is not in effect and this is the start or restart
!***  of a forecast then write out a history file.
!-----------------------------------------------------------------------
!
        history_output_0_a: IF(NTIMESTEP==0                             &
                                 .AND.                                  &
                              .NOT.RESTARTED_RUN                        &
                                 .AND.                                  &
                               FILTER_METHOD==0                         &
                                 .AND.                                  &
                               domain_int_state%QUILTING)THEN
!
          CWRT='History'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_INTEGRATE                              &
                          ,MYPE                                         &
                          ,CWRT)
!
        ENDIF  history_output_0_a
!
        history_output_0_b: IF(RESTARTED_RUN                            &
                                 .AND.                                  &
                               RESTARTED_RUN_FIRST                      &
                                 .AND.                                  &
                               RST_OUT_00                               &
                                 .AND.                                  &
                               domain_int_state%QUILTING)THEN
!
          RESTARTED_RUN_FIRST=.FALSE.
          CWRT='History'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_INTEGRATE                              &
                          ,MYPE                                         &
                          ,CWRT)
!
        ENDIF  history_output_0_b
!
!-----------------------------------------------------------------------
!***  Initialize the timeseries output and write timestep 0 data
!***  for this domain.
!-----------------------------------------------------------------------
!
        time_series_0: IF(.NOT.TS_INITIALIZED) THEN
!
          IF(MYPE<domain_int_state%NUM_PES_FCST)THEN
!
            CALL TIMESERIES_INITIALIZE(DYN_INT_STATE                    &
                                      ,PHY_INT_STATE                    &
                                      ,MY_DOMAIN_ID                     &
                                      ,NTIMESTEP                        &
                                      ,IERR)
!
            IF (IERR /= 0) THEN
              CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
            END IF
!
            CALL TIMESERIES_RUN(DYN_INT_STATE                           &
                               ,PHY_INT_STATE                           &
                               ,MY_DOMAIN_ID                            &
                               ,NTIMESTEP                               &
                               ,IERR)
!
            IF (IERR /= 0) THEN
              CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
            END IF
!
          END IF
!
          TS_INITIALIZED = .TRUE.
!
        END IF time_series_0
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the DOMAIN components.
!-----------------------------------------------------------------------
!
        btim0=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Run DOMAIN Component "//INT_TO_CHAR
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP              &  !<-- The DOMAIN gridded component
                             ,importState=IMP_STATE_DOMAIN              &  !<-- The DOMAIN import state
                             ,exportState=EXP_STATE_DOMAIN              &  !<-- The DOMAIN export state
                             ,clock      =CLOCK_INTEGRATE               &  !<-- The ESMF DOMAIN Clock
                             ,phase      =1                             &  !<-- The phase (subroutine) of DOMAIN Run to execute
                             ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        phase1_tim = (timef()-btim0)
        atm_drv_run_1=atm_drv_run_1+phase1_tim
!
!-----------------------------------------------------------------------
!***  Call the 2nd Phase of the Parent_Child coupler where parents
!***  send data to their children (at the end of all parents'
!***  timesteps, but before any potential filter averaging).
!-----------------------------------------------------------------------
!
        btim0=timef()
!
        IF(NUM_CHILDREN>0)THEN                                             !<-- Call the coupler if there are children
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Call Coupler: Parents Send to Children"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL           &  !<-- The Nesting coupler component
                                ,importState=IMP_STATE_CPL_NEST         &  !<-- The Nesting coupler import state
                                ,exportState=EXP_STATE_CPL_NEST         &  !<-- The Nesting coupler export state
                                ,clock      =CLOCK_INTEGRATE            &  !<-- The DOMAIN Clock
                                ,phase      =2                          &  !<-- The phase (subroutine) of the Coupler to execute
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        atm_drv_run_cpl2=atm_drv_run_cpl2+(timef()-btim0)
!-----------------------------------------------------------------------
!***  If there is filtering, execute Phase 2 of the DOMAIN Run step.
!-----------------------------------------------------------------------
!
        btim0=timef()
!
        IF(FILTER_METHOD>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Run DOMAIN Filtering "
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP            &  !<-- The DOMAIN gridded component for this domain
                               ,importState=IMP_STATE_DOMAIN            &  !<-- The DOMAIN import state
                               ,exportState=EXP_STATE_DOMAIN            &  !<-- The DOMAIN export state
                               ,clock      =CLOCK_INTEGRATE             &  !<-- The ESMF DOMAIN Clock
                               ,phase      =2                           &  !<-- The phase (subroutine) of DOMAIN Run to execute
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
        atm_drv_run_2=atm_drv_run_2+(timef()-btim0)
!
!-----------------------------------------------------------------------
!***  Increment the timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_INTEGRATE                    &
                              ,rc   =RC)
        kount_steps=kount_steps+1
        IF(FILTER_METHOD > 0 .AND. MYPE == 0) THEN
          write(0,*) 'filter running, kount_steps: ', kount_steps
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the timestep from the DOMAIN clock and
!***  print forecast time.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK_INTEGRATE                 &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!
!
!-----------------------------------------------------------------------
!***  Write timeseries data for this timestep on this domain.
!-----------------------------------------------------------------------
!
        IF(FILTER_METHOD==0 .and. MYPE<domain_int_state%NUM_PES_FCST)THEN
          CALL TIMESERIES_RUN(DYN_INT_STATE                             &
                             ,PHY_INT_STATE                             &
                             ,MY_DOMAIN_ID                              &
                             ,NTIMESTEP                                 &
                             ,IERR)
          IF (IERR /= 0) THEN
            CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          END IF
!
        END IF
!
!-----------------------------------------------------------------------
!***  Now that the clock has been advanced, write the history output
!***  if it is time to do so.  This must be done through the DOMAIN
!***  component since its internal state contains the output data 
!***  so we call Phase 3 of the Run step for the DOMAIN components.
!-----------------------------------------------------------------------
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Call Run3 for History Output"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        btim0=timef()
        CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP              &  !<-- The DOMAIN gridded component
                             ,importState=IMP_STATE_DOMAIN              &  !<-- The DOMAIN import state
                             ,exportState=EXP_STATE_DOMAIN              &  !<-- The DOMAIN export state
                             ,clock      =CLOCK_INTEGRATE               &  !<-- The ESMF Clock for "mini" forecast
                             ,phase      =3                             &  !<-- The phase (subroutine) of DOMAIN Run to execute
                             ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        atm_drv_run_3=atm_drv_run_3+(timef()-btim0)
!
!-----------------------------------------------------------------------
!***  Lead forecast task prints timestep information in free forecast.
!-----------------------------------------------------------------------
!
        IF(MYPE==0.AND.FILTER_METHOD==0)THEN
!!!       IF(I_AM_A_FCST_TASK==ESMF_TRUE)THEN
          WRITE(0,25)NTIMESTEP-1,NTIMESTEP*DT/3600.,phase1_tim
   25     FORMAT(' Finished Timestep ',i5,' ending at ',f7.3,           &
                 ' hours: elapsed integration time ',g10.4)
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Print clocktimes of integration sections
!***  on MPI task of choice.
!-----------------------------------------------------------------------
!

        IF (FILTER_METHOD == 0) THEN
        IF(ESMF_AlarmIsRinging(alarm=ALARM_CLOCKTIME                    &  !<-- The alarm to print clocktimes used by model parts
                              ,rc   =RC))THEN
!
          CALL PRINT_CLOCKTIMES(NTIMESTEP,MYPE,NPE_PRINT)
!
        ENDIF
        ENDIF
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ENDDO timeloop_drv
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      driver_run_end: IF(FILTER_METHOD==0)THEN                             !<-- For standard integration, no filtering
!
!-----------------------------------------------------------------------
!***  Extract Clocktimes of the Parent-Child Coupler from that
!***  component's export state and print them.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_NEST==ESMF_TRUE.AND.I_AM_A_FCST_TASK==ESMF_TRUE)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl1 Recv Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl1_Recv_Time'                 &  !<-- Name of the attribute to extract
                                ,value=cpl1_recv_tim                    &  !<-- Clocktime for Recv in Phase 1 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
        IF(NUM_CHILDREN>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Wait Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Wait_Time'                 &  !<-- Name of the attribute to extract
                                ,value=cpl2_wait_tim                    &  !<-- Clocktime for Wait in Phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Comp Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Comp_Time'                 &  !<-- Name of the attribute to extract
                                ,value=cpl2_comp_tim                    &  !<-- Clocktime for Compute in Phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Send Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Send_Time'                 &  !<-- Name of the attribute to extract
                                ,value=cpl2_send_tim                    &  !<-- Clocktime for Send in Phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        WRITE(0,*)' '
        WRITE(0,*)' Clocktime NMM_RUN'
        WRITE(0,*)' '
        WRITE(0,*)'   Run Phase 1=',atm_drv_run_1
        WRITE(0,*)'   Run Phase 2=',atm_drv_run_2
        WRITE(0,*)'   Run Phase 3=',atm_drv_run_3
        WRITE(0,*)' '
!
        IF(NESTING.AND.I_AM_A_FCST_TASK==ESMF_TRUE)THEN
!
          IF(I_AM_A_NEST==ESMF_TRUE)THEN
            if (cpl1_recv_tim > 1.0) then
            WRITE(0,*)'   Recv in Cpl Phase 1=',cpl1_recv_tim
            WRITE(0,*)'   Total Cpl Phase 1=',atm_drv_run_cpl1
           endif
          ENDIF
!
!
          IF(NUM_CHILDREN>0)THEN
            if( cpl2_comp_tim >2.0) then
            WRITE(0,*)' '
            WRITE(0,*)'   Cpl Phase 2 Compute=',cpl2_comp_tim
            WRITE(0,*)'   Cpl Phase 2 Wait   =',cpl2_wait_tim
            WRITE(0,*)'   Cpl Phase 2 Send   =',cpl2_send_tim
            WRITE(0,*)'   Total Cpl Phase 2  =',atm_drv_run_cpl2
            WRITE(0,*)' '
           endif
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ELSE  driver_run_end                                                 !<-- Filtering is in effect
!
!-----------------------------------------------------------------------
!***  If we are completing the execution of digital filtering then reset
!***  the Clock and times.
!-----------------------------------------------------------------------
!
        IF(CLOCK_DIRECTION=='Bckward ')THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Get CurrTime and Timestep for Bckward"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_INTEGRATE                   &
                            ,currtime=CURRTIME                          &
                            ,timestep=TIMESTEP_FILTER                   &
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          TIMESTEP_FILTER=-TIMESTEP_FILTER                                 !<-- We must set the timestep back to positive.
!
!-----------------------------------------------------------------------
          filter_method_block : IF(FILTER_METHOD==3)THEN
!-----------------------------------------------------------------------
            ndfiloop: DO I=1,NDFISTEP
!
!-----------------------------------------------------------------------
!***  Now set the timestep of the children at which they receive
!***  data from their parent.  This must be known by the parents
!***  since it will provide the proper tag to the MPI data sent.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
              parents_only: IF(NUM_CHILDREN>0                           &
                                   .AND.                                &
                               I_AM_A_FCST_TASK==ESMF_TRUE) THEN
!-----------------------------------------------------------------------
!
                IF(.NOT.ALLOCATED(LOC_PAR_CHILD_TIME_RATIO)) THEN
                  ALLOCATE(LOC_PAR_CHILD_TIME_RATIO(1:NUM_CHILDREN))
                ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INTEGRATE: Parent/child DT Ratio for TDFI"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeGet(state    =IMP_STATE_CPL_NEST        &  !<-- The parent-child coupler import state
                                      ,name     ='Parent-Child Time Ratio' &  !<-- Name of the attribute to extract
                                      ,count    =NUM_CHILDREN              &  !<-- # of items in the Attribute
                                      ,valueList=LOC_PAR_CHILD_TIME_RATIO  &  !<-- Ratio of parent to child DTs
                                      ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                DO N=1,NUM_CHILDREN
                  NSTEP_CHILD_RECV(N)=NSTEP_CHILD_RECV(N) +  LOC_PAR_CHILD_TIME_RATIO(N)
                ENDDO
!
!-----------------------------------------------------------------------
!
              ENDIF parents_only
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="NMM_INTEGRATE: Advance Clock for TDFI"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_ClockAdvance(clock   =CLOCK_INTEGRATE           &
                                    ,timestep=TIMESTEP_FILTER           &  !<-- Advance the clock to the forward starttime
                                    ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            ENDDO ndfiloop
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INTEGRATE: Get Current Clock Time for TDFI"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock   =CLOCK_INTEGRATE                 &
                              ,currtime=CURRTIME                        &
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          ENDIF  filter_method_block
!
!-----------------------------------------------------------------------
!
        ELSEIF(CLOCK_DIRECTION=='Forward ')THEN
!
!-----------------------------------------------------------------------
!
          IF(FILTER_METHOD==1)THEN
!
            CURRTIME=HALFDFITIME-TIMESTEP
            NTIMESTEP=NTIMESTEP-(HALFDFIINTVAL/TIMESTEP)-1
            NTIMESTEP_ESMF=NTIMESTEP
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INTEGRATE: Set Time to Half Filter Interval"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockSet(clock       =CLOCK_INTEGRATE             &  !<-- Reset current time and timestep to the
                              ,currtime    =CURRTIME                    &  !    halfway point of the filter interval.
                              ,advanceCount=NTIMESTEP_ESMF              &
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF  driver_run_end
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_INTEG==ESMF_SUCCESS)THEN
!       WRITE(0,*)'NMM RUN step succeeded'
      ELSE
        WRITE(0,*)'NMM RUN step failed RC_INTEG=',RC_INTEG
      ENDIF
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      SUBROUTINE RESET_ALARMS           
!
!-----------------------------------------------------------------------
!***  For normal forecast integration set the Alarm ring times
!***  while accounting for restarts and digital filtering.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(CURRTIME==STARTTIME)THEN
        ALARM_HISTORY_RING  =CURRTIME
        ALARM_RESTART_RING  =CURRTIME
        ALARM_CLOCKTIME_RING=CURRTIME
      ELSE
        IF(RESTARTED_RUN)THEN
          ALARM_HISTORY_RING  =CURRTIME+INTERVAL_HISTORY
          ALARM_RESTART_RING  =CURRTIME+INTERVAL_RESTART
          ALARM_CLOCKTIME_RING=CURRTIME+INTERVAL_CLOCKTIME
        ELSE
          ALARM_HISTORY_RING  =STARTTIME+INTERVAL_HISTORY
          ALARM_RESTART_RING  =STARTTIME+INTERVAL_RESTART
          ALARM_CLOCKTIME_RING=STARTTIME+INTERVAL_CLOCKTIME
        ENDIF
      ENDIF
!
!-------------------------------------------------
!***  Adjust time of History Alarm if necessary.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_HISTORY_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_HISTORY_RING                         &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year 
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_HISTORY_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_HISTORY_RING                       &  !<-- Reset the time for initial history output
                         ,yy  =YY                                       &  !<-- Year 
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-------------------------------------------------
!***  Adjust time of Restart Alarm if necessary.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_RESTART_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_RESTART_RING                         &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year 
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_RESTART_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_RESTART_RING                       &  !<-- Reset the time for initial restart output
                         ,yy  =YY                                       &  !<-- Year 
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-------------------------------------------------------------
!***  Adjust time of Alarm for clocktime writes if necessary.
!-------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_CLOCKTIME_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_CLOCKTIME_RING                       &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year 
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_CLOCKTIME_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_CLOCKTIME_RING                     &  !<-- Reset the time for clocktime prints
                         ,yy  =YY                                       &  !<-- Year 
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Now create the three Alarms using the final ringtimes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Creating the Alarms."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALARM_HISTORY=ESMF_AlarmCreate(name             ='ALARM_HISTORY'     &
                                    ,clock            =CLOCK_INTEGRATE     &  !<-- DOMAIN Clock
                                    ,ringTime         =ALARM_HISTORY_RING  &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =INTERVAL_HISTORY    &  !<-- Time interval between history output
                                    ,ringTimeStepCount=1                   &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.             &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
!
      ALARM_RESTART=ESMF_AlarmCreate(name             ='ALARM_RESTART'     &
                                    ,clock            =CLOCK_INTEGRATE     &  !<-- DOMAIN Clock
                                    ,ringTime         =ALARM_RESTART_RING  &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =INTERVAL_RESTART    &  !<-- Time interval between  restart output (ESMF)
                                    ,ringTimeStepCount=1                   &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.             &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
!
      ALARM_CLOCKTIME=ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'     &
                                      ,clock            =CLOCK_INTEGRATE       &  !<-- DOMAIN Clock
                                      ,ringTime         =ALARM_CLOCKTIME_RING  &  !<-- Forecast start time (ESMF)
                                      ,ringInterval     =INTERVAL_CLOCKTIME    &  !<-- Time interval between clocktime prints (ESMF)
                                      ,ringTimeStepCount=1                     &  !<-- The Alarm rings for this many timesteps
                                      ,sticky           =.false.               &  !<-- Alarm does not ring until turned off
                                      ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RESET_ALARMS
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NMM_INTEGRATE
!
!-----------------------------------------------------------------------
