!-----------------------------------------------------------------------
!
      MODULE MODULE_FWD_INTEGRATE
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
      PUBLIC :: NMM_FWD_INTEGRATE,GFS_FWD_INTEGRATE       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
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
      SUBROUTINE NMM_FWD_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
				  ,HALFDFITIME                          &
                                  ,HALFDFIINTVAL                        &
				  ,SDFITIME                             &
				  ,DFITIME                              &
				  ,STARTTIME                            &
				  ,ALARM_CLOCKTIME                      &
				  ,ALARM_HISTORY                        &
				  ,ALARM_RESTART                        &
				  ,MYPE                                 &
				  ,NUM_TRACERS_MET                      &
				  ,NUM_TRACERS_CHEM                     &
				  ,NTIMESTEP                            &
				  ,NDFISTEP                             &
				  ,NPE_PRINT                            &
				  ,DFIHR                                &
				  ,PHYSICS_ON                           &
				  ,RESTARTED_RUN)
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
      INTEGER(KIND=KINT),INTENT(IN)          :: NPE_PRINT,NDFISTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE                    & !<-- MPI task rank
                                               ,NUM_TRACERS_MET         & !<-- # of meteorological tracer variables
                                               ,NUM_TRACERS_CHEM          !<-- # of chemistry tracer variables
!
      INTEGER(KIND=KINT),INTENT(INOUT)       :: NTIMESTEP                 !<-- The current forecast timestep
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR
!
      TYPE(ESMF_Time),INTENT(INOUT)          :: CURRTIME                  !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)          :: HALFDFITIME
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: HALFDFIINTVAL
      TYPE(ESMF_Time),INTENT(INOUT)          :: SDFITIME
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      TYPE(ESMF_Time),INTENT(INOUT)          :: DFITIME
!
      TYPE(ESMF_Alarm),INTENT(INOUT)         :: ALARM_CLOCKTIME           !<-- The ESMF Alarm for clocktime prints
      TYPE(ESMF_Alarm),INTENT(INOUT) 	     :: ALARM_HISTORY             !<-- The ESMF Alarm for history output
      TYPE(ESMF_Alarm),INTENT(INOUT) 	     :: ALARM_RESTART             !<-- The ESMF Alarm for restart output
!
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
!
!-----------------------------------------------------------------------
!***********************************************************************
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
        dfihr_gt0: IF(DFIHR>0)THEN
!
!-----------------------------------------------------------------------
!
          IF (MYPE == 0) THEN
          IF(CURRTIME==SDFITIME)THEN
            CALL DIGITAL_FILTER_DYN_INIT_NMM(atm_int_state%IMP_STATE_DYN &
                                            ,NDFISTEP                    &
                                            ,NUM_TRACERS_MET             &
                                            ,NUM_TRACERS_CHEM)
!
! -------------------- inital summation  ------------------------------
!
            IF(PHYSICS_ON)THEN
              CALL DIGITAL_FILTER_PHY_INIT_NMM(atm_int_state%IMP_STATE_PHY)
            ENDIF
          ENDIF
!
! -------------------- summation stage ---------------------------------
!
          CALL DIGITAL_FILTER_DYN_SUM_NMM(atm_int_state%IMP_STATE_DYN   &
                                         ,NUM_TRACERS_MET               &
                                         ,NUM_TRACERS_CHEM)
!
! ----------------------------------------------------------------------
!
          IF(PHYSICS_ON)THEN
!
            IF(CURRTIME==HALFDFITIME)THEN
                CALL DIGITAL_FILTER_PHY_SAVE_NMM(atm_int_state%IMP_STATE_PHY)
            ENDIF
          ENDIF
!
! ----------------------------------------------------------------------
!
! --------------------- final stage ------------------------------------
!
          IF(CURRTIME==DFITIME)THEN
!              IF (ESMF_AlarmIsRinging(alarm=ALARM_FILTER(3), rc=RC)) then
              print *,' dfi at finaldfitime '
              CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(atm_int_state%IMP_STATE_DYN &
                                                 ,NUM_TRACERS_MET             &
                                                 ,NUM_TRACERS_CHEM)
              IF(PHYSICS_ON)THEN
                CALL DIGITAL_FILTER_PHY_RESTORE_NMM(atm_int_state%IMP_STATE_PHY)
              ENDIF
!
! ----------------------------------------------------------------------
!
!             CALL ESMF_ClockSet(clock   =CLOCK_ATM                     &
!                               ,currtime=HALFDFITIME                   &
!                               ,rc      =RC)
              DFITIME = STARTTIME
              DFIHR = 0
              CALL ESMF_ClockPrint(clock  =CLOCK_ATM                    &
                                  ,options="currtime string"            &
                                  ,rc     =RC)
!           ENDIF
          ENDIF
        ENDIF
!
        ENDIF dfihr_gt0
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
      END SUBROUTINE NMM_FWD_INTEGRATE
      
!-----------------------------------------------------------------------

      SUBROUTINE GFS_FWD_INTEGRATE(gc_gfs_dyn,gc_gfs_phy,gc_atm_cpl,imp_gfs_dyn,exp_gfs_dyn &
                                   ,imp_gfs_phy,exp_gfs_phy,CLOCK_MAIN,ALARM_HISTORY,CURRTIME,HALFDFITIME &
                                   ,HALFDFIINTVAL,SDFITIME,DFITIME,STARTTIME,NDFISTEP,DFIHR,MYPE,PHYSICS_ON)

      USE MODULE_DIGITAL_FILTER_GFS
      USE MODULE_CONTROL,ONLY: TIMEF
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: gc_gfs_phy
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: gc_atm_cpl
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_phy,exp_gfs_phy
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_MAIN                         !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)             :: CURRTIME                           !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)             :: HALFDFITIME
      TYPE(ESMF_TimeInterval),INTENT(INOUT)     :: HALFDFIINTVAL
      TYPE(ESMF_Alarm),INTENT(INOUT)            :: ALARM_HISTORY
      TYPE(ESMF_Time),INTENT(INOUT)             :: SDFITIME
      TYPE(ESMF_Time),INTENT(INOUT)             :: STARTTIME
      TYPE(ESMF_Time),INTENT(INOUT)             :: DFITIME
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR
      INTEGER(KIND=KINT),INTENT(IN)          :: NDFISTEP, MYPE
      LOGICAL,INTENT(IN)                     :: PHYSICS_ON
      INTEGER(KIND=KINT)                     :: RC,RC_LOOP

          do while (.not.esmf_clockisstoptime(CLOCK_MAIN,rc))
          call esmf_logwrite("execute dynamics",esmf_log_info,rc=rc)
          call esmf_gridcomprun(gridcomp   =gc_gfs_dyn          &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_MAIN           &
                               ,rc         =RC)
!
          call err_msg(RC,'execute dynamics',RC_LOOP)
          call err_msg(RC,'execute dynamics',RC_LOOP)
          call esmf_logwrite("couple dyn_exp-to-phy_imp",      &
                             esmf_log_info,rc=rc)
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl          &
                              ,importstate=exp_gfs_dyn         &
                              ,exportstate=imp_gfs_phy         &
                              ,clock      =CLOCK_MAIN           &
                              ,rc         =RC)
!
          call err_msg(RC,'couple dyn-to-phy',RC_LOOP)
!
          IF (PHYSICS_ON) THEN
            call esmf_logwrite("execute physics",esmf_log_info,rc=rc)
            call esmf_gridcomprun(gridcomp   =gc_gfs_phy            &
                                 ,importstate=imp_gfs_phy           &
                                 ,exportstate=exp_gfs_phy           &
                                 ,clock      =CLOCK_MAIN             &
                                 ,rc         =RC)
           call err_msg(RC,'execute physics',RC_LOOP)
          ELSE
           call esmf_logwrite("pass phy_imp to phy_exp ",       &
                               esmf_log_info,rc=rc)
!
            call esmf_cplcomprun(            gc_atm_cpl          &
                                ,importstate=imp_gfs_phy         &
                                ,exportstate=exp_gfs_phy         &
                                ,clock      =CLOCK_MAIN          &
                                ,rc         =RC)
!
            call err_msg(RC,'pass phy_imp-to-phy_exp',RC_LOOP)
          ENDIF

 !
!
!----------------
!-----------------------------------------------------------------------
!***  bring export data from the physics into the coupler
!***  and export it to the dynamics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple phy_exp-to-dyn_imp",         &
                             esmf_log_info,rc=RC)
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl             &
                              ,importstate=exp_gfs_phy            &
                              ,exportstate=imp_gfs_dyn            &
                              ,clock      =CLOCK_MAIN             &
                              ,rc         =RC)
!
          call err_msg(RC,'couple phy_exp-to-dyn_imp',RC_LOOP)

!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_MAIN                          &
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
        CALL ESMF_ClockGet(clock       =CLOCK_MAIN                       &
                          ,currtime    =CURRTIME                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
!
        IF (MYPE==0) THEN
        IF (DFIHR .GT. 0) THEN
! --------------------- first stage -----------------------------------
           IF ( CURRTIME .eq. SDFITIME ) THEN
           CALL DIGITAL_FILTER_DYN_INIT_GFS(imp_gfs_dyn,NDFISTEP)

! -------------------- inital summation  ------------------------------
!
                IF ( PHYSICS_ON ) THEN
                        CALL DIGITAL_FILTER_PHY_INIT_GFS(imp_gfs_phy)
                ENDIF
             ENDIF

!
! -------------------- summation stage ---------------------------------
!
              CALL DIGITAL_FILTER_DYN_SUM_GFS(imp_gfs_dyn,mype)
!
! ----------------------------------------------------------------------
!
              IF( PHYSICS_ON ) then
                IF( CURRTIME .eq. HALFDFITIME ) THEN
                    CALL DIGITAL_FILTER_PHY_SAVE_GFS(imp_gfs_phy)
                ENDIF
              ENDIF
!
! ----------------------------------------------------------------------
!
! --------------------- final stage ------------------------------------
!
            IF( CURRTIME .eq. DFITIME ) then
                 print *,' dfi at finaldfitime '
                CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(imp_gfs_dyn,mype)
                IF( PHYSICS_ON ) then
                  CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
                ENDIF
!
! ----------------------------------------------------------------------
!
!             CALL ESMF_ClockSet(clock   =CLOCK_MAIN                     &
!                               ,currtime=HALFDFITIME                   &
!                               ,rc      =RC)

              DFITIME = STARTTIME
              DFIHR = 0
              CALL ESMF_ClockPrint(clock=CLOCK_MAIN                      &
                              ,options="currtime string"                &
                              ,rc=RC)
             ENDIF
        ENDIF !< -- DFIHR
        ENDIF

        total_tim=total_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!

       ENDDO

        MESSAGE_CHECK="last step dynamics"
                call esmf_gridcomprun(gridcomp=gc_gfs_dyn       &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_MAIN          &
                               ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)


      END SUBROUTINE GFS_FWD_INTEGRATE
!
      END MODULE MODULE_FWD_INTEGRATE
