#include "../../ESMFVersionDefine.h"

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
!   2010-08-18  WANG  - output filtered fields at HALFDFITIME
!   2010-11-10  WANG  - reset integration loop for dfi
!   2011-02   W Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                       ESMF 5 series library and the the
!                       ESMF 3.1.0rp2 library.
!   2011-03   W Yang  - Modified the digiter filter code for turning off
!                       the digiter filter case.
!   2011-10-01  Wang/Lu - MYPE added to GOCART_INTEGRATE argument
!-----------------------------------------------------------------------

      USE ESMF_MOD
      USE MODULE_ERR_MSG
      USE MODULE_INCLUDE

      USE MODULE_DIGITAL_FILTER_GFS
      USE MODULE_GFS_WRITE,        ONLY: WRITE_ASYNC_GFS
      USE MODULE_GOCART_ROUTINES,  ONLY: GOCART_INTEGRATE
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
!
      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GFS_INTEGRATE(gc_gfs_dyn                               &
                                 ,gc_gfs_phy                            &
                                 ,GC_GFS_CHEM                           &
                                 ,gc_gfs_cpl                            &
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
                                 ,CLOCK_GFS                             &
                                 ,OUTPUT_INTERVAL                       &
                                 ,quilting                              &
                                 ,WRITE_GROUP_READY_TO_GO               &
                                 ,CURRTIME                              &
                                 ,STARTTIME                             &
                                 ,NTIMESTEP                             &
                                 ,TIMESTEP                              &
                                 ,DFIHR                                 &
                                 ,LDFIFLTO                              &
                                 ,MYPE                                  &
                                 ,PHYSICS_ON                            &
                                 ,CHEMISTRY_ON)

!
!-----------------------------------------------------------------------
!

      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)	     :: gc_gfs_phy &
                                               ,GC_GFS_CHEM        !<-- The Chemistry component
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: gc_gfs_cpl
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: GC_PHY2CHEM_CPL    !<-- The Phy-to-Chem coupler component
      TYPE(ESMF_CplComp),INTENT(INOUT)       :: GC_CHEM2PHY_CPL    !<-- The Chem-to-Phy coupler component

!jw
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: wrt_comps(:)
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_phy,exp_gfs_phy  &
                                               ,IMP_GFS_CHEM,EXP_GFS_CHEM  !<-- The import/export states for Chemistry component

!jw
      TYPE(ESMF_State),INTENT(INOUT)         :: imp_gfs_wrt,exp_gfs_wrt
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_GFS                          !<-- The GFS Component's ESMF Clock
      TYPE(ESMF_Time),INTENT(INOUT)          :: CURRTIME                           !<-- The current forecast time
      TYPE(ESMF_Time),INTENT(INOUT)          :: STARTTIME
      INTEGER(KIND=KINT),INTENT(INOUT)       :: DFIHR, NTIMESTEP
      INTEGER(KIND=KINT),INTENT(IN)          :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(INout)  :: TIMESTEP                 !<-- The ESMF timestep (s)
      TYPE(ESMF_Logical),INTENT(IN)          :: PHYSICS_ON               !<-- Is physics on (true) or off (false)?
      TYPE(ESMF_Logical),INTENT(IN)          :: CHEMISTRY_ON                       !<-- Is chemistry on (true) or off (false)?
!jw
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: output_interval
      LOGICAL,INTENT(IN)                     :: QUILTING
      LOGICAL,INTENT(IN)                     :: LDFIFLTO
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
      TYPE(ESMF_TimeInterval)                :: HALFDFIINTVAL
!jw
      TYPE(ESMF_Time)                        :: ALARM_OUTPUT_RING
      TYPE(ESMF_Alarm), SAVE                 :: ALARM_OUTPUT
      LOGICAL                                :: Cpl_flag
      TYPE(ESMF_LOGICAL)                     :: Cpl_flag_ESMF
      LOGICAL, SAVE                          :: write_flag = .true.
      LOGICAL, SAVE                          :: first      = .true.
      LOGICAL, SAVE                          :: first_dfi  = .true.
      LOGICAL, SAVE                          :: end_dfi  = .false.
      LOGICAL, SAVE                          :: LDFI  = .false.
      LOGICAL                                :: LSPC 
      LOGICAL                                :: lalm1,lalm2,lskip
      INTEGER                                :: YY, MM, DD, H, M, S
!
      TYPE(ESMF_Field)             :: FIELD
      real(ESMF_KIND_R8), dimension(:,:), pointer :: tmp_ptr2d
!-----------------------------------------------------------------------
!***  Set up alarm for output,alarm starts from current time
!-----------------------------------------------------------------------
!
       IF(first) THEN
           ALARM_OUTPUT_RING = CURRTIME + OUTPUT_INTERVAL
           CALL ESMF_TimeGet(time = ALARM_OUTPUT_RING                         &
                            ,yy   = YY                                        &
                            ,mm   = MM                                        &
                            ,dd   = DD                                        &
                            ,h    = H                                         &
                            ,m    = M                                         &
                            ,s    = S                                         &
                            ,rc   = RC)

           write(0,*)'set up alarm_output_ring,H=',H,'m=',m,'s=',s
!
           ALARM_OUTPUT =ESMF_AlarmCreate(name             ='ALARM_OUTPUT'     &
                                         ,clock            =CLOCK_GFS          &  !<-- GFS Clock
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
            LDFI=.true.
            NTIMESTEPH  = NDFISTEP
!
            CALL ESMF_ClockGet(clock       =CLOCK_GFS                     &
                            ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has adva nced
                            ,currTime    =currTime                        &  !<-- # of times the clock has advance d
                            ,rc          =RC)
            NTIMESTEP = NTIMESTEP_ESMF
            IF(NTIMESTEP>2*NTIMESTEPH) LDFI=.false.
      ELSE
          HALFDFITIME = STARTTIME
          DFITIME     = STARTTIME
      END IF
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Dynamics component
!-----------------------------------------------------------------------
!
      integrate: DO WHILE(.NOT.ESMF_ClockIsStopTime(CLOCK_GFS, rc = RC)   &
                          .OR. (LDFI.and.DFIHR>0) )

          CALL ESMF_LogWrite("Execute GFS Dynamics",ESMF_LOG_INFO,rc=RC)
!
          CALL ESMF_GridCompRun(gridcomp   =GC_GFS_DYN                    &
                               ,importstate=IMP_GFS_DYN                   &
                               ,exportstate=EXP_GFS_DYN                   &
                               ,clock      =CLOCK_GFS                     &
                               ,rc         =RC)
!
          CALL ERR_MSG(RC,'execute dynamics',RC_LOOP)
!
          CALL ESMF_LogWrite("couple dyn_exp-to-phy_imp"                  &
                             ,esmf_log_info,rc=rc)
!
          CALL ESMF_ClockGet(clock       =CLOCK_GFS                       &
                            ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                            ,currTime    =currTime                     &  !<-- # of times the clock has advanced
                            ,rc          =RC)

          NTIMESTEP = NTIMESTEP_ESMF
!
!*** decide when to disable alarm

          IF(DFIHR > 0) THEN
              LALM1=.not.LDFIFLTO.and.CURRTIME>HALFDFITIME     &
               .and.CURRTIME<=DFITIME .and. first_dfi

              LALM2=LDFIFLTO.and.CURRTIME>=HALFDFITIME         &
               .and.CURRTIME<=DFITIME .and. first_dfi

              IF(LALM1.or.LALM2) THEN
                call ESMF_AlarmDisable(ALARM_OUTPUT,rc=rc)
                first_dfi=.false.
              ENDIF
          END IF
!
!-----------------------------------------------------------------------
!***  Call the Write component if it is time.
!-----------------------------------------------------------------------
!
           LSPC= NTIMESTEP==1.OR.LDFI.and.LDFIFLTO.and.NTIMESTEP==NTIMESTEPH 
!           LSPC= NTIMESTEP==1.OR.LDFI.and.LDFIFLTO.and.NTIMESTEP==NTIMESTEPH &
!                 .or.LDFI.and..not.LDFIFLTO.and.NTIMESTEP==NTIMESTEPH+1

 outputdyn: IF(((ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC) .AND. &
               ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)) .OR.    & !<-- The history output alarm
               LSPC)  .AND. write_flag) THEN
               CALL WRITE_ASYNC_GFS(WRT_COMPs,exp_gfs_dyn               &
                                ,imp_gfs_wrt,exp_gfs_wrt                &
                                ,CLOCK_GFS                              &
                                ,MYPE                                   &
                                ,WRITE_GROUP_READY_TO_GO)
            ELSE
                write_flag = .true.
            END IF outputdyn
!
!*** skip phys when output filtered fields at halfdfitime
         IF(DFIHR > 0) THEN
             LSKIP=CURRTIME==DFITIME
         ELSE
             LSKIP = .false.
         END IF
!         write(0,*)'in gfs intg,aft wrt, NTIMESTEP=',NTIMESTEP,'LSKIP=',LSKIP
!
         lskipif: if(.not.LSKIP) then
!
!-----------------------------------------------------------------------
!***  Bring export data from the Dynamics into the coupler
!***  and export it to the Physics.
!-----------------------------------------------------------------------

!
         call esmf_cplcomprun(cplcomp     = gc_gfs_cpl          &
                             ,importstate = exp_gfs_dyn         &
                             ,exportstate = imp_gfs_phy         &
                             ,clock       = CLOCK_GFS           &
                             ,rc          = RC)
!
         call err_msg(RC,'couple dyn-to-phy',RC_LOOP)
!
         IF (PHYSICS_ON == ESMF_True) THEN
!
!-----------------------------------------------------------------------
!***              Execute the Run step of the Physics Component
!-----------------------------------------------------------------------
!
           call esmf_logwrite("execute physics",esmf_log_info,rc=rc)
           call esmf_gridcomprun(gridcomp    = gc_gfs_phy            &
                                ,importstate = imp_gfs_phy           &
                                ,exportstate = exp_gfs_phy           &
                                ,clock       = CLOCK_GFS             &
                                ,rc          = RC)
           call err_msg(RC,'execute physics',RC_LOOP)
!check time step
           CALL ESMF_ClockGet(clock       =CLOCK_GFS                       &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
           NTIMESTEP = NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***              Invoke GOCART
!-----------------------------------------------------------------------
           IF (CHEMISTRY_ON == ESMF_True) THEN

              MESSAGE_CHECK = "Execute GOCART module"

              CALL GOCART_INTEGRATE(                                    &
                                   GC_GFS_CHEM,                         &
                                   GC_PHY2CHEM_CPL,                     &
                                   GC_CHEM2PHY_CPL,                     &
                                   EXP_GFS_PHY,                         &
                                   IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                   CLOCK_GFS, MYPE, RC                  )

              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

           ENDIF
!
         ELSE
!-----------------------------------------------------------------------
!***              Skip the Physics if the user has turned it off.
!-----------------------------------------------------------------------
            call esmf_logwrite("pass phy_imp to phy_exp ",       &
                                esmf_log_info,rc=rc)
!
            call esmf_cplcomprun(              gc_gfs_cpl          &
                                ,importstate = imp_gfs_phy         &
                                ,exportstate = exp_gfs_phy         &
                                ,clock       = CLOCK_GFS           &
                                ,rc          = RC)
!
           call err_msg(RC,'pass phy_imp-to-phy_exp',RC_LOOP)
         ENDIF
!
!-----------------------------------------------------------------------
!***  Bring export data from the Physics into the coupler
!***  and export it to the Dynamics.
!-----------------------------------------------------------------------
!
         call esmf_logwrite("couple phy_exp-to-dyn_imp",           &
                             esmf_log_info,rc=RC)
!
         call esmf_cplcomprun(cplcomp      = gc_gfs_cpl            &
                             ,importstate = exp_gfs_phy            &
                             ,exportstate = imp_gfs_dyn            &
                             ,clock       = CLOCK_GFS              &
                             ,rc          = RC)
!
         call err_msg(RC,'couple phy_exp-to-dyn_imp',RC_LOOP)
!
!
         endif lskipif
!-----------------------------------------------------------------------
!***  Digital filter
!-----------------------------------------------------------------------
!
         filter_block: IF(DFIHR>0.and.LDFI)THEN
!
!--------------------------
!***  Filter's first stage
!--------------------------
!
           IF(CURRTIME == STARTTIME) THEN
             CALL DIGITAL_FILTER_DYN_INIT_GFS(EXP_GFS_DYN,NDFISTEP)
!
!---------------------------
!***  The initial summation
!---------------------------
!
             IF(PHYSICS_ON == ESMF_True) THEN
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
           IF(PHYSICS_ON == ESMF_True) THEN
             IF(CURRTIME == HALFDFITIME) THEN
               CALL DIGITAL_FILTER_PHY_SAVE_GFS(IMP_GFS_PHY)
               NTIMESTEPH = NTIMESTEP_ESMF
             ENDIF
           ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
           IF(CURRTIME == DFITIME)THEN
             CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(exp_gfs_dyn)
!
             IF(PHYSICS_ON == ESMF_True) THEN
               CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
             ENDIF
!
             CALL ESMF_ClockSet(clock       =CLOCK_GFS             &
                               ,currtime    =HALFDFITIME           &
                               ,advanceCount=NTIMESTEPH            &
                               ,rc          =RC)
!
             DFITIME = STARTTIME
             DFIHR   = 0
             end_dfi=.true.
!

!            CALL ESMF_ClockPrint(clock  =CLOCK_GFS                      &
!                                ,options="currtime string"              &
!                                 ,rc     =RC)

           ENDIF
!
!-----------------------------------------------------------------------
!
         ENDIF  filter_block
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
         CALL ESMF_AttributeGet(imp_gfs_dyn, 'Cpl_flag', Cpl_flag_ESMF, rc = rc)

         IF(Cpl_flag_ESMF == ESMF_FALSE .AND. .NOT. LSKIP) THEN
             CALL ESMF_ClockAdvance(clock = CLOCK_GFS, rc = RC)
         END IF
#else
         CALL ESMF_AttributeGet(imp_gfs_dyn, 'Cpl_flag', Cpl_flag, rc = rc)

         IF(.NOT. Cpl_flag .AND. .NOT. LSKIP) THEN
             CALL ESMF_ClockAdvance(clock = CLOCK_GFS, rc = RC)
         END IF
#endif

         CALL ESMF_ClockGet(clock       =CLOCK_GFS             &
                          ,advanceCount=NTIMESTEP_ESMF        &  !<-- # of times the clock has advanced
                          ,rc          =RC)
         NTIMESTEP = NTIMESTEP_ESMF
         CALL ESMF_ClockGet(clock       =CLOCK_GFS             &
                          ,currTime=CURRTIME                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
!-----------------------------------------------------------------------
!*** enable alarm
!
         LALM1=end_dfi.and..not.LDFIFLTO.and.currTime>HALFDFITIME
         LALM2=end_dfi.and.LDFIFLTO.and.currTime==HALFDFITIME
!
         IF (lalm1 .or.lalm2 ) then
             CALL ESMF_AlarmEnable(alarm=ALARM_OUTPUT                    &
                                ,rc     =RC)
             end_dfi=.false.
         endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 
       ENDDO integrate

       call esmf_gridcomprun(gridcomp=gc_gfs_dyn             &
                            ,importstate=imp_gfs_dyn         &
                            ,exportstate=exp_gfs_dyn         &
                            ,clock      =CLOCK_GFS           &
                            ,rc         =RC)
    output2: IF(ESMF_AlarmIsRinging(alarm=ALARM_OUTPUT, rc = RC)) THEN    !<-- The history output alarm
                 CALL WRITE_ASYNC_GFS(WRT_COMPs,exp_gfs_dyn            &
                         ,imp_gfs_wrt,exp_gfs_wrt                      &
                         ,CLOCK_GFS                                    &
                         ,MYPE                                         &
                         ,WRITE_GROUP_READY_TO_GO)
                 write_flag = .false.
             ELSE
                 write_flag = .true.
             END IF output2
!
      END SUBROUTINE GFS_INTEGRATE
!
      END MODULE MODULE_GFS_INTEGRATE
