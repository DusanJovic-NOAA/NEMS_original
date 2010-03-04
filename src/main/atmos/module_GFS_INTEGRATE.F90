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
!-----------------------------------------------------------------------

!
      USE ESMF_MOD
      USE MODULE_ERR_MSG
!
      USE MODULE_DIGITAL_FILTER_GFS
      USE MODULE_WRITE_ROUTINES_GFS,ONLY: WRITE_ASYNC_GFS
      USE MODULE_GOCART_INTEGRATE, ONLY: GOCART_INTEGRATE
      USE MODULE_CONTROL,ONLY: TIMEF
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
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GFS_INTEGRATE(GC_GFS_DYN                               &
                              ,GC_GFS_PHY                               &
                              ,GC_GFS_CHEM                              &
                              ,GC_ATM_CPL                               &
                              ,GC_PHY2CHEM_CPL                          &
                              ,WRT_COMPS                                &
                              ,IMP_GFS_DYN                              &
                              ,EXP_GFS_DYN                              &
                              ,IMP_GFS_PHY                              &
                              ,EXP_GFS_PHY                              &
                              ,IMP_GFS_CHEM                             &
                              ,EXP_GFS_CHEM                             &
                              ,IMP_GFS_WRT                              &
                              ,EXP_GFS_WRT                              &
                              ,CLOCK_ATM                                &
                              ,OUTPUT_INTERVAL                          &
                              ,QUILTING                                 &
                              ,WRITE_GROUP_READY_TO_GO                  &
                              ,CURRTIME                                 &
                              ,STARTTIME                                &
                              ,NTIMESTEP                                &
                              ,TIMESTEP                                 &
                              ,DFIHR                                    &
                              ,MYPE                                     &
                              ,PHYSICS_ON                               &
                              ,CHEMISTRY_ON)

!
!-----------------------------------------------------------------------
!
!------------------
!***  Arguments IN
!------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: MYPE                                !<-- The local task ID
!
      LOGICAL,INTENT(IN) :: QUILTING                                       !<-- Are separate quilt tasks being used for output?
!
      TYPE(ESMF_Logical),INTENT(IN) :: PHYSICS_ON                          !<-- Is physics on (true) or off (false)?
      TYPE(ESMF_Logical),INTENT(IN) :: CHEMISTRY_ON                        !<-- Is chemistry on (true) or off (false)?
!
      TYPE(ESMF_TimeInterval),INTENT(IN) :: TIMESTEP                       !<-- The ESMF timestep (s)
!
!---------------------
!***  Arguments INOUT
!---------------------
!
      INTEGER(kind=KINT),INTENT(INOUT) :: DFIHR                         &  !<-- Filter duration
                                         ,NTIMESTEP                     &  !<-- The forecast timestep
                                         ,WRITE_GROUP_READY_TO_GO          !<-- The number of the current Write group
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_DYN                   &  !<-- The Dynamics component
                                          ,GC_GFS_PHY                   &  !<-- The Physics component
                                          ,GC_GFS_CHEM                     !<-- The Chemistry component
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: WRT_COMPS(:)                    !<-- The Write components for output
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_ATM_CPL                       !<-- The Dynamics-Physics coupler component
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_PHY2CHEM_CPL                  !<-- The Phy-to-Chem coupler component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_GFS_DYN,EXP_GFS_DYN         &  !<-- The import/export states for Dynamics component
                                       ,IMP_GFS_PHY,EXP_GFS_PHY         &  !<-- The import/export states for Physics component
                                       ,IMP_GFS_CHEM,EXP_GFS_CHEM       &  !<-- The import/export states for Chemistry component
                                       ,IMP_GFS_WRT,EXP_GFS_WRT            !<-- The import/export states for Write components
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_ATM                          !<-- The ATM Component's ESMF Clock
!
      TYPE(ESMF_Time),INTENT(INOUT) :: CURRTIME                         &  !<-- The current forecast time
                                      ,STARTTIME                           !<-- The forecast start time
!
      TYPE(ESMF_TimeInterval),INTENT(INOUT) :: OUTPUT_INTERVAL             !<-- Frequency of history output
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,N_GROUP,NDFISTEP,RC,RC_LOOP
      INTEGER(kind=KINT) :: YY,MM,DD,H,M,S
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                       & !<-- The current forecast timestep (ESMF_INT)
                                   ,NTIMESTEPH                             !<-- The timestep at fdfi
!
      TYPE(ESMF_Time) :: ALARM_OUTPUT_RING,DFITIME,HALFDFITIME,CurrentTime
!
      TYPE(ESMF_TimeInterval) :: HALFDFIINTVAL
!
      TYPE(ESMF_Alarm) :: ALARM_OUTPUT
!
      LOGICAL :: first=.true.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Set up alarm for output,alarm starts from current time
!-----------------------------------------------------------------------
!
      ALARM_OUTPUT_RING=CURRTIME+OUTPUT_INTERVAL
!
      ALARM_OUTPUT=ESMF_AlarmCreate(name             ='ALARM_OUTPUT'    &
                                   ,clock            =CLOCK_ATM         &  !<-- ATM Clock
                                   ,ringTime         =ALARM_OUTPUT_RING &  !<-- Forecast/Restart start time (ESMF)
                                   ,ringInterval     =OUTPUT_INTERVAL   &  !<-- Time interval between
                                   ,ringTimeStepCount=1                 &  !<-- The Alarm rings for this many timesteps
                                   ,sticky           =.false.           &  !<-- Alarm does not ring until turned off
                                   ,rc               =RC)
!
      IF(DFIHR>0)THEN
!
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL            &
                                 ,h           =DFIHR                    &
                                 ,rc          =RC)
!
        NDFISTEP = HALFDFIINTVAL / TIMESTEP
        HALFDFITIME = STARTTIME + HALFDFIINTVAL
        DFITIME = HALFDFITIME + HALFDFIINTVAL
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Dynamics component
!-----------------------------------------------------------------------
!
      integrate: DO WHILE(.NOT.ESMF_ClockIsStopTime(CLOCK_ATM,RC))
!
!-----------------------------------------------------------------------
!
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
!
         NTIMESTEP=NTIMESTEP_ESMF
         IF(DFIHR>0 .and. currentTime>HALFDFITIME .and.                 &
           currentTime<=DFITIME .and. first) THEN
           call ESMF_AlarmDisable(ALARM_OUTPUT,rc=rc)
           first=.false.
         ENDIF
         print *,'in gfsintg, NTIMESTEP=',NTIMESTEP,'alarmenable=',   &
          ESMF_AlarmISEnabled(ALARM_OUTPUT,rc=rc)
!
!-----------------------------------------------------------------------
!***  Call the Write component if it is time.
!-----------------------------------------------------------------------
!
        outputdyn: IF((ESMF_AlarmIsEnabled(alarm=ALARM_OUTPUT,rc=RC).and.&
                      ESMF_AlarmIsRinging(alarm=ALARM_OUTPUT,rc=rc)) .or.& !<-- The history output alarm
                      NTIMESTEP==1) THEN
!
          CALL WRITE_ASYNC_GFS(WRT_COMPS                                &
                              ,EXP_GFS_DYN                              &
                              ,IMP_GFS_WRT                              &
                              ,EXP_GFS_WRT                              &
                              ,CLOCK_ATM                                &
                              ,MYPE                                     &
                              ,WRITE_GROUP_READY_TO_GO)
        ENDIF outputdyn
!
!-----------------------------------------------------------------------
!***  Bring export data from the Dynamics into the coupler
!***  and export it to the Physics.
!-----------------------------------------------------------------------
!
        CALL ESMF_CplCompRun(cplcomp    =GC_ATM_CPL                     &
                            ,importstate=EXP_GFS_DYN                    &
                            ,exportstate=IMP_GFS_PHY                    &
                            ,clock      =CLOCK_ATM                      &
                            ,rc         =RC)
!
        CALL ERR_MSG(RC,'couple dyn-to-phy',RC_LOOP)

!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Physics Component
!-----------------------------------------------------------------------
!
        IF (PHYSICS_ON==ESMF_True) THEN
          CALL ESMF_LogWrite("execute physics",ESMF_LOG_INFO,rc=RC)
!
          CALL ESMF_GridCompRun(gridcomp   =GC_GFS_PHY                  &
                               ,importstate=IMP_GFS_PHY                 &
                               ,exportstate=EXP_GFS_PHY                 &
                               ,clock      =CLOCK_ATM                   &
                               ,rc         =RC)
!
          call err_msg(RC,'execute physics',RC_LOOP)
!check time step
          CALL ESMF_ClockGet(clock       =CLOCK_ATM                     &
                            ,advanceCount=NTIMESTEP_ESMF                &  !<-- # of times the clock has advanced
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
                                   EXP_GFS_PHY,                         &
                                   IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                   CLOCK_ATM, RC                    )

             CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

          ENDIF
!
!
!-----------------------------------------------------------------------
!***  Skip the Physics is the user has turned it off. 
!-----------------------------------------------------------------------
!
        ELSE
          CALL ESMF_LogWrite("pass phy_imp to phy_exp "                 &
                             ,ESMF_LOG_INFO,rc=rc)
!
          CALL ESMF_CplCompRun(cplcomp    =GC_ATM_CPL                   &
                              ,importstate=IMP_GFS_PHY                  &
                              ,exportstate=EXP_GFS_PHY                  &
                              ,clock      =CLOCK_ATM                    &
                              ,rc         =RC)
!
          CALL ERR_MSG(RC,'pass phy_imp-to-phy_exp',RC_LOOP)
!
        ENDIF

!
!-----------------------------------------------------------------------
!***  Bring export data from the Physics into the coupler
!***  and export it to the Dynamics.
!-----------------------------------------------------------------------
!
        CALL ESMF_LogWrite("couple phy_exp-to-dyn_imp"                  &
                           ,ESMF_LOG_INFO,rc=RC)
!
        CALL ESMF_CplCompRun(cplcomp    =GC_ATM_CPL                     &
                            ,importstate=EXP_GFS_PHY                    &
                            ,exportstate=IMP_GFS_DYN                    &
                            ,clock      =CLOCK_ATM                      &
                            ,rc         =RC)
!
        CALL ERR_MSG(RC,'couple phy_exp-to-dyn_imp',RC_LOOP)
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
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockAdvance(clock=CLOCK_ATM                          &
                              ,rc   =RC)
!
        CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &
                          ,currTime=CURRTIME                            &  !<-- The current forecast time
                          ,rc      =RC)
!
!-----------------------------------------------------------------------
! 
      ENDDO  integrate
!
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompRun(gridcomp   =GC_GFS_DYN                      &
                           ,importstate=IMP_GFS_DYN                     &
                           ,exportstate=EXP_GFS_DYN                     &
                           ,clock      =CLOCK_ATM                       &
                           ,rc         =RC)
!
      outputdyn2: IF(ESMF_AlarmIsRinging(alarm=ALARM_OUTPUT             &  !<-- The history output alarm
                                        ,rc   =RC))THEN
!
        CALL WRITE_ASYNC_GFS(WRT_COMPS                                  &
                            ,EXP_GFS_DYN                                &
                            ,IMP_GFS_WRT                                &
                            ,EXP_GFS_WRT                                &
                            ,CLOCK_ATM                                  &
                            ,MYPE                                       &
                            ,WRITE_GROUP_READY_TO_GO)
      ENDIF outputdyn2
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_INTEGRATE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_INTEGRATE
!
!-----------------------------------------------------------------------
