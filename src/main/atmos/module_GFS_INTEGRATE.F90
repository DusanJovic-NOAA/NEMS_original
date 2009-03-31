!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_INTEGRATE
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
      USE MODULE_ERR_MSG
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
      PUBLIC :: GFS_INTEGRATE
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

      SUBROUTINE GFS_INTEGRATE(gc_gfs_dyn                               &
                                 ,gc_gfs_phy                            &
                                 ,gc_atm_cpl                            &
                                 ,imp_gfs_dyn                           &
                                 ,exp_gfs_dyn                           &
                                 ,imp_gfs_phy                           &
                                 ,exp_gfs_phy                           &
                                 ,CLOCK_MAIN                            &
                                 ,CURRTIME                              &
                                 ,STARTTIME                             &
                                 ,NTIMESTEP                             &
                                 ,TIMESTEP                              &
                                 ,DFIHR                                 &
                                 ,MYPE                                  &
                                 ,PHYSICS_ON)

!
!-----------------------------------------------------------------------
!


      USE MODULE_DIGITAL_FILTER_GFS

      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: gc_gfs_phy
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: gc_atm_cpl
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_phy,exp_gfs_phy
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_MAIN                         !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)             :: CURRTIME                           !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)             :: STARTTIME
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR, NTIMESTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(IN)          :: TIMESTEP                      !<-- The ESMF timestep (s)
      LOGICAL,INTENT(IN)                     :: PHYSICS_ON
      INTEGER(KIND=KINT)                     :: RC,RC_LOOP,I
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF                        !<-- The current forecast timestep (ESMF_INT)
      INTEGER(KIND=KINT)                     :: NDFISTEP
      TYPE(ESMF_Time)                        :: HALFDFITIME
      TYPE(ESMF_Time)                        :: DFITIME
      TYPE(ESMF_TimeInterval)                :: HALFDFIINTVAL


      IF (DFIHR .GT. 0) THEN
            CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,h=DFIHR,rc=RC)
            NDFISTEP = HALFDFIINTVAL / TIMESTEP
            HALFDFITIME = STARTTIME + HALFDFIINTVAL
            DFITIME = HALFDFITIME + HALFDFIINTVAL
      ENDIF



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
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK AND PRINT FORECAST TIME.
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------

!        IF (MYPE==0) THEN
        IF (DFIHR .GT. 0) THEN

! --------------------- first stage -----------------------------------
           IF ( CURRTIME .eq. STARTTIME ) THEN
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
              CALL DIGITAL_FILTER_DYN_SUM_GFS(imp_gfs_dyn)
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
                CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(imp_gfs_dyn)
                IF( PHYSICS_ON ) then
                  CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
                ENDIF
!
! ----------------------------------------------------------------------
!
              call ESMF_ClockSet(CLOCK_MAIN                             &
                                ,currtime=HALFDFITIME                   &
                                ,rc      =rc)

              DFITIME = STARTTIME
              DFIHR = 0
              print *,' dfi reset time to '
              CALL ESMF_ClockPrint(clock=CLOCK_MAIN                      &
                              ,options="currtime string"                &
                              ,rc=RC)
             ENDIF
        ENDIF !< -- DFIHR
!        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_MAIN                          &
                              ,rc   =RC)
        CALL ESMF_ClockGet(clock       =CLOCK_MAIN                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
        NTIMESTEP=NTIMESTEP_ESMF
        CALL ESMF_ClockGet(clock       =CLOCK_MAIN                       &
                          ,currTime=CURRTIME                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 
       ENDDO

                call esmf_gridcomprun(gridcomp=gc_gfs_dyn       &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_MAIN          &
                               ,rc         =RC)


      END SUBROUTINE GFS_INTEGRATE
!
      END MODULE MODULE_GFS_INTEGRATE
