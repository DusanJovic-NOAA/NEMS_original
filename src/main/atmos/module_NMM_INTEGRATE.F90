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
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!
      PUBLIC :: NMM_INTEGRATE
!
!-----------------------------------------------------------------------
      INCLUDE '../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
!
!
!
!-----------------------------------------------------------------------
!***  FOR DETERMINING CLOCKTIMES OF VARIOUS PIECES OF THE DYNAMICS.
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT) :: btim,btim0
!
!-----------------------------------------------------------------------
!
      LOGICAL,SAVE :: RESTARTED_RUN_FIRST=.TRUE.
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_INTEGRATE(ATM_GRID_COMP                            &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
				  ,STARTTIME                            &
				  ,TIMEINTERVAL_CLOCKTIME               &
				  ,TIMEINTERVAL_HISTORY                 &
				  ,TIMEINTERVAL_RESTART                 &
				  ,MYPE                                 &
				  ,NUM_TRACERS_MET                      &
				  ,NUM_TRACERS_CHEM                     &
				  ,NTIMESTEP                            &
				  ,NPE_PRINT                            &
				  ,PHYSICS_ON                           &
				  ,RESTARTED_RUN)
!
!-----------------------------------------------------------------------
!
      USE MODULE_DYNAMICS_INTERNAL_STATE                                  !<-- Horizontal loop limits obtained here
      USE MODULE_PHYSICS_INTERNAL_STATE
!
      USE MODULE_CLOCKTIMES
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK

      USE MODULE_DIGITAL_FILTER_NMM
      USE MODULE_ATM_INTERNAL_STATE
      USE MODULE_CONTROL,ONLY: TIMEF
      USE MODULE_WRITE_ROUTINES,ONLY: WRITE_ASYNC                         !<-- These are routines used only when asynchronous
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
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME                  !<-- The current forecast time
      TYPE(ESMF_TimeInterval),INTENT(IN)     :: TIMEINTERVAL_CLOCKTIME
      TYPE(ESMF_TimeInterval),INTENT(IN)     :: TIMEINTERVAL_HISTORY
      TYPE(ESMF_TimeInterval),INTENT(IN)     :: TIMEINTERVAL_RESTART
      TYPE(ESMF_Time)                        ::  ALARM_HISTORY_RING,ALARM_RESTART_RING,ALARM_CLOCKTIME_RING
      TYPE(ESMF_Time)                        :: ADJTIME_HISTORY, ADJTIME_RESTART, ADJTIME_CLOCKTIME
!
      TYPE(ESMF_Alarm)         :: ALARM_CLOCKTIME           !<-- The ESMF Alarm for clocktime prints
      TYPE(ESMF_Alarm) 	     :: ALARM_HISTORY             !<-- The ESMF Alarm for history output
      TYPE(ESMF_Alarm) 	     :: ALARM_RESTART             !<-- The ESMF Alarm for restart output
!
      integer :: YY, MM, DD, H, M, S

      LOGICAL,INTENT(IN)                     :: PHYSICS_ON
      LOGICAL,INTENT(IN)                     :: RESTARTED_RUN
!
!
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT)         :: RC,RC_LOOP
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF                        !<-- The current forecast timestep (ESMF_INT)
!
      CHARACTER(ESMF_MAXSTR) :: CWRT                                      !<-- Restart/History label
      INTEGER(KIND=KINT)       :: HDIFF_ON = 1
!
!-----------------------------------------------------------------------
!***********************************************************************
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert HDIFF into Dynamics Export State"
       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =atm_int_state%EXP_STATE_DYN    &  !<-- Dynamics impor
                              ,name     ='HDIFF'                 &  !<-- The attribute'
                              ,value= HDIFF_ON     &  !<-- Insert this qu
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       MESSAGE_CHECK="ADUSTING ALARM RINGS"
       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       IF (CURRTIME .EQ. STARTTIME)THEN       
           ALARM_HISTORY_RING   = CURRTIME
           ALARM_RESTART_RING   = CURRTIME
           ALARM_CLOCKTIME_RING = CURRTIME
       ELSE
         IF(RESTARTED_RUN)THEN
           ALARM_HISTORY_RING   = CURRTIME + TIMEINTERVAL_HISTORY
           ALARM_RESTART_RING   = CURRTIME + TIMEINTERVAL_RESTART
           ALARM_CLOCKTIME_RING = CURRTIME + TIMEINTERVAL_CLOCKTIME
         ELSE
           ALARM_HISTORY_RING   = STARTTIME + TIMEINTERVAL_HISTORY
           ALARM_RESTART_RING   = STARTTIME + TIMEINTERVAL_RESTART
           ALARM_CLOCKTIME_RING = STARTTIME + TIMEINTERVAL_CLOCKTIME
         ENDIF
       ENDIF

       call ESMF_TimeGet(ALARM_HISTORY_RING, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC) 
        IF (M .ne. 0) THEN
         H=H+1
         M=0 
         CALL ESMF_TimeSet(ALARM_HISTORY_RING,yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC)
        ENDIF
       call ESMF_TimeGet(ALARM_RESTART_RING, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC)
        IF (M .ne. 0) THEN
         H=H+1
         M=0
         CALL ESMF_TimeSet(ALARM_RESTART_RING,yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC)
        ENDIF
       call ESMF_TimeGet(ALARM_CLOCKTIME_RING, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC)
        IF (M .ne. 0) THEN
         H=H+1
         M=0
         CALL ESMF_TimeSet(ALARM_CLOCKTIME_RING,yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=RC)
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       MESSAGE_CHECK="CREATING ALARMS "
       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


      ALARM_HISTORY=ESMF_AlarmCreate(name             ='ALARM_HISTORY'      &
                                    ,clock            =CLOCK_ATM            &  !<-- ATM Clock
                                    ,ringTime         =ALARM_HISTORY_RING             &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =TIMEINTERVAL_HISTORY &  !<-- Time interval between
                                    ,ringTimeStepCount=1                    &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.              &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
      ALARM_RESTART=ESMF_AlarmCreate(name             ='ALARM_RESTART'      &
                                    ,clock            =CLOCK_ATM            &  !<-- ATM Clock
                                    ,ringTime         =ALARM_RESTART_RING             &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =TIMEINTERVAL_RESTART &  !<-- Time interval between  restart output (ESMF)
                                    ,ringTimeStepCount=1                    &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.              &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
      ALARM_CLOCKTIME=ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'      &
                                    ,clock              =CLOCK_ATM              &  !<-- ATM Clock
                                    ,ringTime           =ALARM_CLOCKTIME_RING    &  !<-- Forecast start time (ESMF)
                                    ,ringInterval       =TIMEINTERVAL_CLOCKTIME &  !<-- Time interval between clocktime prints (ESMF)
                                    ,ringTimeStepCount  =1                      &  !<-- The Alarm rings for this many timesteps
                                    ,sticky             =.false.                &  !<-- Alarm does not ring until turned off
                                    ,rc                 =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!
      timeloop: DO WHILE (.NOT.ESMF_ClockIsStopTime(CLOCK_ATM,RC))
!
!-----------------------------------------------------------------------
!
        btim0=timef()
!
!-----------------------------------------------------------------------
!***  WRITE A HISTORY FILE AT THE START OF THE FIRST TIMESTEP
!***  OTHERWISE WRITE IT AT THE END OF THE APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!
       history_output_0: IF(NTIMESTEP==0)THEN
!
          IF(atm_int_state%QUILTING)THEN
            CWRT='History'
            CALL WRITE_ASYNC(ATM_GRID_COMP                              &
                            ,ATM_INT_STATE                              &
                            ,CLOCK_ATM                                  &
                            ,MYPE                                       &
                            ,CWRT)
          ENDIF
!
        ENDIF history_output_0
!
!-----------------------------------------------------------------------
!***  WRITE A HISTORY FILE AT THE BEGINNING OF THE RESTARTED RUN
!***  OTHERWISE WRITE IT AT THE END OF THE APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!
        restart_output_0: IF(RESTARTED_RUN.and.RESTARTED_RUN_FIRST) THEN
!
          IF(atm_int_state%QUILTING)THEN
            CWRT='History'
            CALL WRITE_ASYNC(ATM_GRID_COMP                              &
                            ,ATM_INT_STATE                              &
                            ,CLOCK_ATM                                  &
                            ,MYPE                                       &
                            ,CWRT)
          ENDIF
!
          RESTARTED_RUN_FIRST=.FALSE.
!
        ENDIF restart_output_0
!
!-----------------------------------------------------------------------
!***  THE FORECAST TASKS EXECUTE THE RUN STEP OF THE DYNAMICS.
!***  THIS IS THE RUN SUBROUTINE SPECIFIED IN
!***  THE DYNAMICS REGISTER ROUTINE CALLED IN
!***  ESMF_GridCompSetServices ABOVE.
!-----------------------------------------------------------------------
!
        fcst_pes: IF(MYPE<atm_int_state%NUM_PES_FCST)THEN                  !<-- Only the forecast tasks integrate
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run Step for Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompRun(gridcomp   =atm_int_state%DYN_GRID_COMP &  !<-- The dynamics component
                               ,importState=atm_int_state%IMP_STATE_DYN &  !<-- The dynamics import state
                               ,exportState=atm_int_state%EXP_STATE_DYN &  !<-- The dynamics export state
                               ,clock      =CLOCK_ATM                   &  !<-- The ATM clock
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  BRING EXPORT DATA FROM THE DYNAMICS INTO THE COUPLER
!***  AND EXPORT IT TO THE PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Coupler moves Information from Dynamics to Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The coupler import state = dynamics export state
                              ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The coupler export state = physics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM clock
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! - FIX later
 RC=ESMF_SUCCESS
 RC_LOOP=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  EXECUTE THE RUN STEP OF THE PHYSICS.
!-----------------------------------------------------------------------
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run Step for Physics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
         IF (PHYSICS_ON) THEN
!
          CALL ESMF_GridCompRun(gridcomp   =atm_int_state%PHY_GRID_COMP &  !<-- The physics component
                               ,importState=atm_int_state%IMP_STATE_PHY &  !<-- The physics import state
                               ,exportState=atm_int_state%EXP_STATE_PHY &  !<-- The physics export state
                               ,clock      =CLOCK_ATM                   &  !<-- The ATM Clock
                               ,rc         =RC)
         ELSE


          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%IMP_STATE_PHY        &  !<-- The coupler import state = physics export state
                              ,exportState=atm_int_state%EXP_STATE_PHY        &  !<-- The coupler export state = dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                              ,rc         =RC)

          ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  BRING EXPORT DATA FROM THE PHYSICS INTO THE COUPLER
!***  AND EXPORT IT TO THE DYNAMICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Coupler moves Information from Physics to Dynamics"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%EXP_STATE_PHY        &  !<-- The coupler import state = physics export state
                              ,exportState=atm_int_state%IMP_STATE_DYN        &  !<-- The coupler export state = dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                              ,rc         =RC)
!
!-----------------------------------------------------------------------
!
        ENDIF fcst_pes
        

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! - FIX later
 RC=ESMF_SUCCESS
 RC_LOOP=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_ATM                          &
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK AND PRINT FORECAST TIME.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &
                          ,currtime=CURRTIME                            &
                          ,rc      =RC)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  WRITE A HISTORY FILE IF THE ALARM INDICATES IT IS TIME TO DO SO
!***  AFTER TIMESTEP 0.
!-----------------------------------------------------------------------
!
        output: IF(ESMF_AlarmIsRinging(alarm=ALARM_HISTORY              &  !<-- The history output alarm
                                      ,rc   =RC))THEN
!
          IF(atm_int_state%QUILTING)THEN
            CWRT='History'
            CALL WRITE_ASYNC(ATM_GRID_COMP                              &
                            ,ATM_INT_STATE                              &
                            ,CLOCK_ATM                                  &
                            ,MYPE                                       &
                            ,CWRT)
          ENDIF
!
        ENDIF output
!
!-----------------------------------------------------------------------
!***  WRITE A RESTART FILE IF THE ALARM INDICATES IT IS TIME TO DO SO.
!-----------------------------------------------------------------------
!
        restart: IF(ESMF_AlarmIsRinging(alarm=ALARM_RESTART             &  !<-- The restart output alarm
                                       ,rc   =RC))THEN
!
          IF(atm_int_state%QUILTING)THEN
            CWRT='Restart'
            CALL WRITE_ASYNC(ATM_GRID_COMP                              &
                            ,ATM_INT_STATE                              &
                            ,CLOCK_ATM                                  &
                            ,MYPE                                       &
                            ,CWRT)
          ENDIF
!
        ENDIF restart
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        total_tim=total_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!***  PRINT CLOCKTIMES OF INTEGRATION SECTIONS
!***  ON MPI TASK OF CHOICE.
!-----------------------------------------------------------------------
!
        clocktimes: IF(ESMF_AlarmIsRinging(alarm=ALARM_CLOCKTIME        &  !<-- The alarm to print clocktimes
                                          ,rc   =RC))THEN
!
          CALL PRINT_CLOCKTIMES(NTIMESTEP,MYPE,NPE_PRINT)
!
        ENDIF clocktimes
!-----------------------------------------------------------------------
!
      ENDDO timeloop
!
      END SUBROUTINE NMM_INTEGRATE
!
      END MODULE MODULE_NMM_INTEGRATE
