!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PRIMARY INTEGRATION RUNSTREAM OF THE GFS
!***  WITHIN SUBROUTINE GFS_INTEGRATE.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2009-12-23  Lu    - GFS_INTEGRATE modified to loop thru dyn, phy, &
!                       chem gridded component
!   2010-02-04  Lu    - GOCART_INTEGRATE added
!   2010-02-05  WANG  - change alarm set up for restart option of GFS
!   2010-03-09  Lu    - Add CHEM2PHY CPL
!-----------------------------------------------------------------------

      USE ESMF_MOD
      USE MODULE_ERR_MSG

      USE MODULE_DIGITAL_FILTER_GFS
      USE MODULE_GFS_WRITE,        ONLY: WRITE_ASYNC_GFS
      USE MODULE_GOCART_ROUTINES,  ONLY: GOCART_INTEGRATE
      USE MODULE_CONTROL,          ONLY: TIMEF
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
      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GFS_INTEGRATE(gc_gfs_dyn                               &
                                 ,gc_gfs_phy                            &
                                 ,GC_GFS_CHEM                           &
                                 ,gc_atm_cpl                            &
                                 ,GC_PHY2CHEM_CPL                       &
                                 ,GC_CHEM2PHY_CPL                       &
                                 ,wrt_comps                             &
                                 ,imp_gfs_dyn                           &
                                 ,exp_gfs_dyn                           &
                                 ,imp_gfs_phy                           &
                                 ,exp_gfs_phy                           &
                                 ,IMP_GFS_CHEM                          &
                                 ,EXP_GFS_CHEM                          &
                                 ,imp_gfs_wrt                           &
                                 ,exp_gfs_wrt                           &
                                 ,CLOCK_ATM                             &
                                 ,OUTPUT_INTERVAL                       &
                                 ,quilting                              &
                                 ,WRITE_GROUP_READY_TO_GO               &
                                 ,CURRTIME                              &
                                 ,STARTTIME                             &
                                 ,NTIMESTEP                             &
                                 ,TIMESTEP                              &
                                 ,DFIHR                                 &
                                 ,MYPE                                  &
                                 ,PHYSICS_ON                            &
                                 ,CHEMISTRY_ON)

!
!-----------------------------------------------------------------------
!

      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: gc_gfs_phy &
                                          ,GC_GFS_CHEM                     !<-- The Chemistry component
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: gc_atm_cpl
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_PHY2CHEM_CPL                  !<-- The Phy-to-Chem coupler component
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_CHEM2PHY_CPL                  !<-- The Chem-to-Phy coupler component

!jw
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: wrt_comps(:)
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_phy,exp_gfs_phy  &
                                       ,IMP_GFS_CHEM,EXP_GFS_CHEM          !<-- The import/export states for Chemistry component

!jw
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_wrt,exp_gfs_wrt
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                          !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)          :: CURRTIME                           !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR, NTIMESTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(IN)     :: TIMESTEP                           !<-- The ESMF timestep (s)
      TYPE(ESMF_Logical),INTENT(IN)          :: PHYSICS_ON                         !<-- Is physics on (true) or off (false)?
      TYPE(ESMF_Logical),INTENT(IN)          :: CHEMISTRY_ON                       !<-- Is chemistry on (true) or off (false)?
!jw
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: output_interval
      LOGICAL,INTENT(IN)                     :: QUILTING
      INTEGER(KIND=KINT),INTENT(INOUT)       :: WRITE_GROUP_READY_TO_GO
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT)                     :: RC,RC_LOOP,I
      INTEGER(kind=ESMF_KIND_I8)             :: NTIMESTEP_ESMF                   & !<-- The current forecast timestep (ESMF_INT)
                                               ,NTIMESTEPH                         !<-- The timestep at fdfi

      INTEGER(KIND=KINT)                     :: NDFISTEP

      TYPE(ESMF_Time)                        :: HALFDFITIME
      TYPE(ESMF_Time)                        :: DFITIME
      TYPE(ESMF_Time)                        :: CurrentTime
      TYPE(ESMF_TimeInterval)                :: HALFDFIINTVAL
!jw
      TYPE(ESMF_Time)                        :: ALARM_OUTPUT_RING
      TYPE(ESMF_Alarm), SAVE                 :: ALARM_OUTPUT
      TYPE(ESMF_LOGICAL)                     :: Cpl_flag
      LOGICAL, SAVE                          :: write_flag = .true.
      LOGICAL, SAVE                          :: first      = .true.
      LOGICAL, SAVE                          :: first_dfi  = .true.
      INTEGER                                :: YY, MM, DD, H, M, S
!
!-----------------------------------------------------------------------
!***  Set up alarm for output,alarm starts from current time
!-----------------------------------------------------------------------
!
       IF(first) THEN
           ALARM_OUTPUT_RING=CURRTIME+OUTPUT_INTERVAL
           CALL ESMF_TimeGet(time=ALARM_OUTPUT_RING                          &
                            ,yy  =YY                                         &
                            ,mm  =MM                                         &
                            ,dd  =DD                                         &
                            ,h   =H                                          &
                            ,m   =M                                          &
                            ,s   =S                                          &
                            ,rc  =RC)
!
           IF(M /= 0) THEN
               H = H + 1
               M = 0
               CALL ESMF_TimeSet(time=ALARM_OUTPUT_RING                        &
                                ,yy  =YY                                       &
                                ,mm  =MM                                       &
                                ,dd  =DD                                       &
                                ,h   =H                                        &
                                ,m   =M                                        &
                                ,s   =S                                        &
                                ,rc  =RC)
           END IF

           write(0,*)'alarm_output_ring,H=',H,'m=',m,'s=',s
!
           ALARM_OUTPUT =ESMF_AlarmCreate(name             ='ALARM_OUTPUT'     &
                                         ,clock            =CLOCK_ATM          &  !<-- ATM Clock
                                         ,ringTime         =ALARM_OUTPUT_RING  &  !<-- Forecast/Restart start time (ESMF)
                                         ,ringInterval     =OUTPUT_INTERVAL    &  !<-- Time interval between
                                         ,ringTimeStepCount=1                  &  !<-- The Alarm rings for this many timesteps
                                         ,sticky           =.false.            &  !<-- Alarm does not ring until turned off
                                         ,rc               =RC)
            first = .false.
        END IF
!
      IF(DFIHR > 0) THEN
            CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL               &
                                     ,h           =DFIHR                       &
                                     ,rc          =RC)
            NDFISTEP    = HALFDFIINTVAL / TIMESTEP
            HALFDFITIME = STARTTIME     + HALFDFIINTVAL
            DFITIME     = HALFDFITIME   + HALFDFIINTVAL
      END IF
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Dynamics component
!-----------------------------------------------------------------------
!
      integrate: DO WHILE(.NOT.ESMF_ClockIsStopTime(CLOCK_ATM, rc = RC))

          CALL ESMF_LogWrite("Execute GFS Dynamics",ESMF_LOG_INFO,rc=RC)
!
          CALL ESMF_GridCompRun(gridcomp   =GC_GFS_DYN                    &
                               ,importstate=IMP_GFS_DYN                   &
                               ,exportstate=EXP_GFS_DYN                   &
                               ,clock      =CLOCK_ATM                     &
                               ,rc         =RC)
!
          CALL ERR_MSG(RC,'execute dynamics',RC_LOOP)
!
          CALL ESMF_LogWrite("couple dyn_exp-to-phy_imp"                  &
                             ,esmf_log_info,rc=rc)
!
          CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                            ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                            ,currTime    =currentTime                     &  !<-- # of times the clock has advanced
                            ,rc          =RC)

          NTIMESTEP=NTIMESTEP_ESMF
          IF(DFIHR>0 .and. currentTime>HALFDFITIME .and.                 &
            currentTime<=DFITIME .and. first_dfi) THEN
            call ESMF_AlarmDisable(ALARM_OUTPUT,rc=rc)
            first_dfi=.false.
          ENDIF
          print *,'in gfsintg, NTIMESTEP=',NTIMESTEP,'alarmenable=',   &
           ESMF_AlarmISEnabled(ALARM_OUTPUT,rc=rc)
!
!-----------------------------------------------------------------------
!***  Call the Write component if it is time.
!-----------------------------------------------------------------------
!
 outputdyn: IF(((ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC) .AND. &
               ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)) .OR.    & !<-- The history output alarm
               NTIMESTEP == 1)  .AND. write_flag) THEN
                 CALL WRITE_ASYNC_GFS(WRT_COMPs,exp_gfs_dyn               &
                                  ,imp_gfs_wrt,exp_gfs_wrt                &
                                  ,CLOCK_ATM                              &
                                  ,MYPE                                   &
                                  ,WRITE_GROUP_READY_TO_GO)
            ELSE
                write_flag = .true.
            END IF outputdyn
!
!-----------------------------------------------------------------------
!***  Bring export data from the Dynamics into the coupler
!***  and export it to the Physics.
!-----------------------------------------------------------------------
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl          &
                              ,importstate=exp_gfs_dyn         &
                              ,exportstate=imp_gfs_phy         &
                              ,clock      =CLOCK_ATM            &
                              ,rc         =RC)
!
          call err_msg(RC,'couple dyn-to-phy',RC_LOOP)
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Physics Component
!-----------------------------------------------------------------------
!
          IF (PHYSICS_ON==ESMF_True) THEN
            call esmf_logwrite("execute physics",esmf_log_info,rc=rc)
            call esmf_gridcomprun(gridcomp   =gc_gfs_phy            &
                                 ,importstate=imp_gfs_phy           &
                                 ,exportstate=exp_gfs_phy           &
                                 ,clock      =CLOCK_ATM              &
                                 ,rc         =RC)
            call err_msg(RC,'execute physics',RC_LOOP)
!check time step
            CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                           ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                           ,rc          =RC)
!
            NTIMESTEP=NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***  Invoke GOCART
!-----------------------------------------------------------------------
          IF (CHEMISTRY_ON==ESMF_True) THEN

              MESSAGE_CHECK="Execute GOCART module"

              CALL GOCART_INTEGRATE(                                    &
                                   GC_GFS_CHEM,                         &
                                   GC_PHY2CHEM_CPL,                     &
                                   GC_CHEM2PHY_CPL,                     &
                                   EXP_GFS_PHY,                         &
                                   IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                   CLOCK_ATM, RC                    )

             CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

          ENDIF
!
!-----------------------------------------------------------------------
!***  Skip the Physics is the user has turned it off.
!-----------------------------------------------------------------------
          ELSE
           call esmf_logwrite("pass phy_imp to phy_exp ",       &
                               esmf_log_info,rc=rc)
!
            call esmf_cplcomprun(            gc_atm_cpl          &
                                ,importstate=imp_gfs_phy         &
                                ,exportstate=exp_gfs_phy         &
                                ,clock      =CLOCK_ATM           &
                                ,rc         =RC)
!
            call err_msg(RC,'pass phy_imp-to-phy_exp',RC_LOOP)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Bring export data from the Physics into the coupler
!***  and export it to the Dynamics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple phy_exp-to-dyn_imp",         &
                             esmf_log_info,rc=RC)
!
          call esmf_cplcomprun(cplcomp    =gc_atm_cpl             &
                              ,importstate=exp_gfs_phy            &
                              ,exportstate=imp_gfs_dyn            &
                              ,clock      =CLOCK_ATM              &
                              ,rc         =RC)
!
          call err_msg(RC,'couple phy_exp-to-dyn_imp',RC_LOOP)
!
!-----------------------------------------------------------------------
!***  Digital filter
!-----------------------------------------------------------------------
!
        filter_block: IF(DFIHR>0)THEN
!
!--------------------------
!***  Filter's first stage
!--------------------------
!
          IF(CURRTIME==STARTTIME)THEN
            CALL DIGITAL_FILTER_DYN_INIT_GFS(EXP_GFS_DYN,NDFISTEP)
!
!---------------------------
!***  The initial summation
!---------------------------
!
            IF(PHYSICS_ON==ESMF_True)THEN
              CALL DIGITAL_FILTER_PHY_INIT_GFS(imp_gfs_phy)
            ENDIF
!
          ENDIF
!
!-------------------------
!***  The summation stage
!-------------------------
!
          CALL DIGITAL_FILTER_DYN_SUM_GFS(EXP_GFS_DYN)
!
          IF(PHYSICS_ON==ESMF_True)THEN
            IF(CURRTIME==HALFDFITIME)THEN
              CALL DIGITAL_FILTER_PHY_SAVE_GFS(IMP_GFS_PHY)
              NTIMESTEPH=NTIMESTEP_ESMF
            ENDIF
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          IF(CURRTIME==DFITIME)THEN
            write(0,*)' DFI at final DFITIME '
            CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(exp_gfs_dyn)
            IF(PHYSICS_ON==ESMF_True)THEN
              CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
            ENDIF
!
            CALL ESMF_ClockSet(clock       =CLOCK_ATM                   &
                              ,currtime    =HALFDFITIME                 &
                              ,advanceCount=NTIMESTEPH                  &
                              ,rc          =RC)
!
            DFITIME = STARTTIME
            DFIHR = 0
            write(0,*)' DFI reset time to:'
            CALL ESMF_ClockPrint(clock  =CLOCK_ATM                      &
                                ,options="currtime string"              &
                                ,rc     =RC)

!
            CALL ESMF_AlarmEnable(alarm=ALARM_OUTPUT                    &
                                ,rc     =RC)

          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF  filter_block
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(imp_gfs_dyn, 'Cpl_flag', Cpl_flag, rc = rc)
        IF(Cpl_flag == ESMF_FALSE) THEN
            CALL ESMF_ClockAdvance(clock = CLOCK_ATM, rc = RC)
        END IF

        CALL ESMF_ClockGet(clock       =CLOCK_ATM                        &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
        NTIMESTEP=NTIMESTEP_ESMF
        CALL ESMF_ClockGet(clock       =CLOCK_ATM                        &
                          ,currTime=CURRTIME                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 
       ENDDO integrate

                call esmf_gridcomprun(gridcomp=gc_gfs_dyn       &
                               ,importstate=imp_gfs_dyn         &
                               ,exportstate=exp_gfs_dyn         &
                               ,clock      =CLOCK_ATM           &
                               ,rc         =RC)
!jws
    output2: IF(ESMF_AlarmIsRinging(alarm=ALARM_OUTPUT, rc = RC)) THEN    !<-- The history output alarm
                 CALL WRITE_ASYNC_GFS(WRT_COMPs,exp_gfs_dyn            &
                         ,imp_gfs_wrt,exp_gfs_wrt                      &
                         ,CLOCK_ATM                                    &
                         ,MYPE                                         &
                         ,WRITE_GROUP_READY_TO_GO)
                 write_flag = .false.
             ELSE
                 write_flag = .true.
             END IF output2
!jwe
!
      END SUBROUTINE GFS_INTEGRATE
!
      END MODULE MODULE_GFS_INTEGRATE
