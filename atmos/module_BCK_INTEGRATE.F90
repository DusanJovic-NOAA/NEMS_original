!-----------------------------------------------------------------------
!
      MODULE MODULE_BCK_INTEGRATE
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
      USE MODULE_DYNAMICS_INTERNAL_STATE                                  !<-- Horizontal loop limits obtained here
      USE MODULE_PHYSICS_INTERNAL_STATE
!
      USE MODULE_CLOCKTIMES
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
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
      PUBLIC :: NMM_BCK_INTEGRATE       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!-----------------------------------------------------------------------
      INCLUDE 'kind.inc'
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
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_BCK_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
				  ,STARTTIME                            &
				  ,HALFDFIINTVAL                        &
                                  ,TIMESTEP                             &
                                  ,FILTER_METHOD                        &
				  ,MYPE                                 &
				  ,NUM_TRACERS_MET                      &
				  ,NUM_TRACERS_CHEM                     &
				  ,NTIMESTEP                            &
				  ,NPE_PRINT)                         
!
!-----------------------------------------------------------------------
!
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
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      TYPE(ESMF_Time)          :: DFITIME
      TYPE(ESMF_Time)          :: HALFDFITIME
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: HALFDFIINTVAL
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: TIMESTEP
      INTEGER(KIND=KINT),INTENT(INOUT)       :: FILTER_METHOD
!
!
      INTEGER(KIND=KINT)       :: NDFISTEP, I
      INTEGER(KIND=KINT)       :: MEAN_ON
      INTEGER(KIND=KINT)       :: HDIFF_ON = 0

!
!
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT)         :: RC,RC_LOOP
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF                        !<-- The current forecast timestep (ESMF_INT)
!
      CHARACTER(ESMF_MAXSTR) :: CWRT                                      !<-- Restart/History label
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
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

          NDFISTEP = HALFDFIINTVAL / TIMESTEP
          HALFDFITIME = STARTTIME - HALFDFIINTVAL
          DFITIME = HALFDFITIME - HALFDFIINTVAL

         IF (FILTER_METHOD .EQ. 2) THEN
         MEAN_ON=1
         ELSE IF (FILTER_METHOD .EQ. 3) THEN
         MEAN_ON=0
         ENDIF

          CALL ESMF_ClockSet(CLOCK_ATM   &
         , direction=ESMF_MODE_REVERSE  &
         ,starttime=STARTTIME          &
         , rc=rc)


         timeloop: DO WHILE (CURRTIME .GE. DFITIME) 
!
!-----------------------------------------------------------------------
!
        btim0=timef()
!
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

          CALL ESMF_CplCompRun(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dynamics-physics coupler component
                              ,importState=atm_int_state%IMP_STATE_PHY        &  !<-- The coupler import state = physics export state
                              ,exportState=atm_int_state%EXP_STATE_PHY        &  !<-- The coupler export state = dynamics import state
                              ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                              ,rc         =RC)

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
!-----------------------------------------------------------------------
!
!
          IF (MYPE==0) THEN
          IF(CURRTIME==STARTTIME)THEN
            CALL DIGITAL_FILTER_DYN_INIT_NMM(atm_int_state%IMP_STATE_DYN &
                                            ,NDFISTEP               &
                                            ,NUM_TRACERS_MET             &
                                            ,NUM_TRACERS_CHEM)

! -------------------- inital summation  ------------------------------
!
          ENDIF
!
! -------------------- summation stage ---------------------------------
!
          IF (CURRTIME .LE. STARTTIME)THEN
          CALL DIGITAL_FILTER_DYN_SUM_NMM(atm_int_state%IMP_STATE_DYN   &
                                         ,MEAN_ON                       &
                                         ,NUM_TRACERS_MET               &
                                         ,NUM_TRACERS_CHEM)
          ENDIF
!
! ----------------------------------------------------------------------
!
!
! ----------------------------------------------------------------------
!
! --------------------- final stage ------------------------------------
!
          IF(CURRTIME==DFITIME)THEN
                print *,' dfi at finaldfitime '
                CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(atm_int_state%IMP_STATE_DYN &
                                                 ,NUM_TRACERS_MET             &
                                                 ,NUM_TRACERS_CHEM)
!
! ----------------------------------------------------------------------
!
                CALL ESMF_ClockPrint(clock  =CLOCK_ATM                    &
                                  ,options="currtime string"            &
                                  ,rc     =RC)
          ENDIF
        ENDIF
!
!
!        IF( CURRTIME .NE. DFITIME) THEN
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_ATM                          &
                              ,rc   =RC)

        CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &
                          ,currtime=CURRTIME                            &
                          ,rc      =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!        ENDIF
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!

      ENDDO timeloop

          CALL ESMF_ClockSet(CLOCK_ATM                  &
              , direction=ESMF_MODE_FORWARD                 &
              ,starttime=CURRTIME                           &
              , rc=rc)

          IF (FILTER_METHOD .EQ. 2) THEN
            CALL ESMF_ClockAdvance(CLOCK_ATM              &
              , rc=rc)
          ELSE IF (FILTER_METHOD .EQ. 3) THEN
            DO I=1,NDFISTEP+1
            CALL ESMF_ClockAdvance(CLOCK_ATM              &
              , rc=rc)
            ENDDO
          ENDIF

          CALL ESMF_ClockGet(CLOCK_ATM              &
               ,currtime=CURRTIME                     &
               ,rc= RC)



      END SUBROUTINE NMM_BCK_INTEGRATE
      
!-----------------------------------------------------------------------
!
      END MODULE MODULE_BCK_INTEGRATE
