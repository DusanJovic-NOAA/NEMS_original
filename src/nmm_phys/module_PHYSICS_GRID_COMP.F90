!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PHYSICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE MAIN GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS PHYSICS INITIALIZE, ETC.) 
!***  IN MODULE_ATM_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-28  Vasic - Removed counters (computed in SET_INTERNAL_STATE_PHY)
!   2008-10-    Vasic - Restart capability
!   2008-08-23  Janjic - Removed uz0h, vz0h
!   2008-08-23  Janjic - General hybrid coordinate
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_PHYSICS_INTERNAL_STATE
      USE MODULE_PHYSICS_FIELDS,ONLY : ARRAY_T                          &
                                      ,ARRAY_U,ARRAY_V                  &
                                      ,ARRAY_Q2,ARRAY_PD                &
                                      ,ARRAY_TRACERS                    &
                                      ,ALLOC_FIELDS_PHY
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,MYPE_SHARE                          &
                                   ,DSTRB,IDSTRB
!
      USE MODULE_CONTROL,ONLY : CAPPA,TIMEF
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_FLTBNDS,ONLY : POLEHN,POLEWN,SWAPHN,SWAPWN
!
      USE MODULE_RADIATION    ,ONLY : GFDL_INIT,RADIATION,RDTEMP        &
                                     ,RRTMINIT,SWINIT,TIME_MEASURE
      USE MODULE_TURBULENCE   ,ONLY : MYJPBL_INIT,TURBL
      USE MODULE_SURFACE_LAYER,ONLY : MYJSFC_INIT
      USE MODULE_LANDSURFACE  ,ONLY : DZSOIL,NOAH_LSM_INIT              &
                                     ,NUM_SOIL_LAYERS,SLDPTH
      USE MODULE_CONVECTION   ,ONLY : BMJ_INIT,CUCNVC
      USE MODULE_MICROPHYSICS_NMM ,ONLY : FERRIER_INIT,GSMDRIVE             &
                                     ,WSM3INIT,MICRO_RESTART
      USE MODULE_H_TO_V       ,ONLY : H_TO_V,H_TO_V_TEND
      USE MODULE_GWD          ,ONLY : GWD_INIT
      USE MODULE_PRECIP_ADJUST
!
      USE MODULE_EXCHANGE
      USE MODULE_DIAGNOSE,ONLY: TWR,VWR
!
      USE MODULE_PHYSICS_OUTPUT,ONLY: POINT_PHYSICS_OUTPUT
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
      PUBLIC :: PHY_REGISTER
!
!-----------------------------------------------------------------------
!!!!  INCLUDE 'kind.inc'
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PUBLIC :: IM,JM,LM
!
      INTEGER(KIND=KINT) :: MYPE,NUM_PES
      INTEGER(KIND=KINT) :: START_YEAR,START_MONTH,START_DAY,START_HOUR &
                           ,START_MINUTE,START_SECOND
!
      INTEGER(KIND=KINT),SAVE :: JC,NSTEPS_HIST,NSTEPS_PER_HOUR
!
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),SAVE :: DT
!
      REAL(KIND=KFPT) :: btim,btim0
!
!-----------------------------------------------------------------------
!
!***  FOR NOW, SET THE DOMAIN'S TOP PRESSURE HERE.
!***  THIS WILL BE STANDARDIZED WITH THE DYNAMICS SOON.
!***  IN THE DYNAMICS COMPONENT IT IS SET IN MODULE_CONTROL
!***  IN SUBROUTINE CONSTS.
!
!!!!  REAL :: PT=5000.   !<--  This is read in from input file in subroutine PHYSICS_INITIALIZE
      REAL :: PT
!
!-----------------------------------------------------------------------
!
      TYPE(INTERNAL_STATE),POINTER :: INT_STATE    ! The internal state pointer.
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Physics Gridded Component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_REG=ESMF_SUCCESS
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,PHY_INITIALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type
                                     ,PHY_RUN                           &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type
                                     ,PHY_FINALIZE                      &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CHECK THE ERROR SIGNAL VARIABLE.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' PHY_REGISTER SUCCEEDED'
      ELSE
        WRITE(0,*)' PHY_REGISTER FAILED RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_INITIALIZE(GRID_COMP,IMP_STATE,EXP_STATE,CLOCK     &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SET UP THE MODEL PHYSICS.
!-----------------------------------------------------------------------
!
!     USE MODULE_ESMF_State
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!***  WRAP_INTERNAL_STATE IS DEFINED IN THE INTERNAL STATE MODULE.
!-----------------------------------------------------------------------
!
      INTEGER                      :: L,N,RC
!
      TYPE(WRAP_INTERNAL_STATE)    :: WRAP                               !<-- This wrap is a derived type which contains
                                                                         !    only a pointer to the internal state.  It is needed
                                                                         !    for using different architectures or compilers.
!
      TYPE(ESMF_State)        :: IMP_STATE_WRITE         
      TYPE(ESMF_Grid)         :: GRID
      TYPE(ESMF_VM)           :: VM                                      !<-- The virtual machine
      TYPE(ESMF_TimeInterval) :: DT_ESMF                                 !<-- The timestep from the ATM Clock
!
      INTEGER :: IDENOMINATOR_DT,INTEGER_DT,NUMERATOR_DT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim=timef()
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE PHYSICS TIMERS.
!-----------------------------------------------------------------------
!
      phy_init_tim=0.
      phy_run_tim=0.
      phy_sum_tim=0.
      update_phy_int_state_tim=0.
      pole_swap_phy_tim=0.
      exch_phy_tim=0.
      cucnvc_tim=0.
      gsmdrive_tim=0.
      h_to_v_tim=0.
      radiation_tim=0.
      rdtemp_tim=0.
      turbl_tim=0.
      adjppt_tim=0. 
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE PHYSICS INTERNAL STATE POINTER.
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE,STAT=RC)
!
!-----------------------------------------------------------------------
!***  ATTACH THE INTERNAL STATE TO THE PHYSICS GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      WRAP%INT_STATE=>INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach Physics Internal State to the Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GRID_COMP                      &  !<-- Physics gridcomp
                                        ,WRAP                           &  !<-- Data pointer to internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE IMPORT STATE OF THE WRITE GRIDDED COMPONENT
!***  FROM THE PHYSICS EXPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Write Import State from Physics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =EXP_STATE                          &  !<-- The Physics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
                        ,nestedState=IMP_STATE_WRITE                    &  !<-- Extract Write component import state from Physics export
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE LOCAL DOMAIN STARTING LIMITS AND THE HALO WIDTH INTO
!***  THE PHYSICS INTERNAL STATE.
!-----------------------------------------------------------------------
!
!     IF(IHALO==JHALO)THEN
!       int_state%NHALO=IHALO
!     ELSE
!       RC_INIT=ESMF_FAILURE
!       WRITE(0,*)'Error due to ihalo /= jhalo'
!     ENDIF
!
      int_state%ITS=ITS
      int_state%ITE=ITE
      int_state%JTS=JTS
      int_state%JTE=JTE
!
!-----------------------------------------------------------------------
!***  USE ESMF UTILITIES TO GET INFORMATION FROM THE CONFIGURATION FILE.
!***  THE FUNCTION IS SIMILAR TO READING A NAMELIST.  THE GET_CONFIG
!***  ROUTINE IS THE USER'S.  IT EXTRACTS VALUES FRON THE CONFIG FILE
!***  AND PLACES THEM IN THE NAMELIST COMPONENTS OF THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_PHY(GRID_COMP,INT_STATE,RC)
      IM=int_state%IM
      JM=int_state%JM
      LM=int_state%LM
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM TO OBTAIN THE TASK ID AND TOTAL NUMBER OF TASKS
!***  FOR THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get VM from the Physics Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &  !<-- The Physics gridded component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Task IDs and Number of MPI Tasks from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &
                     ,localpet=int_state%MYPE                           &  !<-- local task rank
                     ,petcount=int_state%NUM_PES                        &  !<-- total # of MPI tasks
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  int_state%NUM_PES TAKEN FROM VM IS THE TOTAL NUMBER OF TASKS
!***  IN THE RUN INCLUDING QUILT TASKS.  ACTUALLY WE WANT JUST THE
!***  NUMBER OF FORECAST TASKS.
!-----------------------------------------------------------------------
!
      int_state%NUM_PES=int_state%INPES*int_state%JNPES
!
      NUM_PES=int_state%NUM_PES  ! The number of forecast tasks
      MYPE=int_state%MYPE        ! The local PE
!
!-----------------------------------------------------------------------
!***  ONLY FORECAST TASKS ARE NEEDED FOR THE REMAINING
!***  INITIALIZATION PROCESS.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<NUM_PES)THEN                                    !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  SET UP THE PHYSICS INTERNAL STATE VARIABLES.
!-----------------------------------------------------------------------
!
        CALL SET_INTERNAL_STATE_PHY(GRID_COMP,INT_STATE)
!
!-----------------------------------------------------------------------
!***  ASSIGN THE FUNDAMENTAL TIMESTEP RETRIEVED FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Fundamental Timestep ffrom ATM Clock" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockGet(clock   =CLOCK                               &  !<-- The ATM Clock
                          ,timeStep=DT_ESMF                             &  !<-- Fundamental timestep (s) (ESMF)
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Real Timestep from ESMF Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                  &  !<-- the ESMF timestep
                                 ,s           =INTEGER_DT               &  !<-- the integer part of the timestep in seconds
                                 ,sN          =NUMERATOR_DT             &  !<-- the numerator of the fractional second
                                 ,sD          =IDENOMINATOR_DT          &  !<-- the denominator of the fractional second
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        int_state%DT=REAL(INTEGER_DT)+REAL(NUMERATOR_DT)                &  !<-- Fundamental tiemstep (s) (REAL)
                                     /REAL(IDENOMINATOR_DT)
        DT=int_state%DT
!
        NSTEPS_PER_HOUR=NINT(3600./DT)
        NSTEPS_HIST=int_state%NHOURS_HISTORY*NSTEPS_PER_HOUR
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE PHYSICS SCHEMES. 
!-----------------------------------------------------------------------
!
        CALL PHYSICS_INITIALIZE(int_state%GFS                           &
                               ,int_state%SHORTWAVE                     &
                               ,int_state%LONGWAVE                      &
                               ,int_state%CONVECTION                    &
                               ,int_state%MICROPHYSICS                  &
                               ,int_state%SFC_LAYER                     &
                               ,int_state%TURBULENCE                    &
                               ,int_state%LAND_SURFACE                  &
                               ,int_state%CO2TF                         &
                               ,IDS,IDE,JDS,JDE,LM                      &
                               ,IMS,IME,JMS,JME                         &
                               ,ITS,ITE,JTS,JTE)
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF Arrays FOR THE IMPORT/EXPORT STATES.
!***  FOR NOW SEND ALLOC_FIELDS_PHY THE ENTIRE INTERNAL STATE
!***  FROM WHICH THE DESIRED VARIABLES WILL BE EXTRACTED FOR
!***  INSERTION INTO THE IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the ESMF Grid from the Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGet(gridcomp=GRID_COMP                        &  !<-- The Physics gridded component
                             ,grid    =GRID                             &  !<-- The ESMF Grid
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ALLOC_FIELDS_PHY(GRID,INT_STATE)
!
!-----------------------------------------------------------------------
!***  ADD THE DESIRED ESMF Arrays TO THE PHYSICS EXPORT STATE.
!***  THE POINTERS INSIDE THE Arrays ARE POINTING TO THE APPROPRIATE
!***  VARIABLES INSIDE THE INTERNAL STATE (see ALLOC_FIELDS_PHY
!***  IN module_PHYSICS_FIELDS.F).
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ADD THE 3D QUANTITIES TO THE EXPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add 3-D Arrays to Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_T                                &  !<-- Temperature
                          ,rc   =RC)
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_U                                &  !<-- U wind component
                          ,rc   =RC)
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_V                                &  !<-- V wind component
                          ,rc   =RC)
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_Q2                               &  !<-- TKE
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ADD THE 2D QUANTITIES TO THE EXPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add 2-D Arrays to Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_PD                               &  !<-- Vertical pressure difference, sigma range
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ADD THE 4D TRACERS ARRAY TO THE EXPORT STATE.
!***  THE NUMBER OF 3D CONSTITUENTS IS GIVEN BY NUM_TRACERS_TOTAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add 4-D Tracer Data to Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=ARRAY_TRACERS                          &  !<-- Tracer variables
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ALSO INSERT THE VALUE OF NUM_TRACERS_TOTAL INTO THE EXPORT STATE.
!***  THIS WILL TELL THE Dyn-Phy Coupler HOW MANY CONSTITUENTS
!***  THERE ARE TO TRANSFER IN THE 4-D TRACERS ARRAY.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert NUM_TRACERS_TOTAL into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='NUM_TRACERS_TOTAL'                &  !<-- The inserted quantity will have this name
                              ,value=int_state%NUM_TRACERS_TOTAL        &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  EXTRACT ALL FORECAST TASKS' HORIZONTAL SUBDOMAIN LIMITS
!***  FROM THE PHYSICS IMPORT STATE AND GIVE THEM TO THE 
!***  PHYSICS INTERNAL STATE.
!***  THIS IS NECESSARY IF QUILTING IS SELECTED BECAUSE THESE
!***  LIMITS WILL BE TAKEN FROM THE DYNAMICS/PHYSICS INTERNAL
!***  STATES, PLACED INTO THE WRITE COMPONENTS' IMPORT STATES
!***  AND USED FOR THE COMBINING OF LOCAL DOMAIN DATA ONTO THE
!***  GLOBAL DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Local Domain Limits to Physics Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The write component import state
                              ,name     ='LOCAL_ISTART'                 &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_ISTART         &  !<-- Extract this attribute from import state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The write component import state
                              ,name     ='LOCAL_IEND'                   &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_IEND           &  !<-- Extract this attribute from import state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The write component import state
                              ,name     ='LOCAL_JSTART'                 &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_JSTART         &  !<-- Extract this attribute from import state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The write component import state
                              ,name     ='LOCAL_JEND'                   &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_JEND           &  !<-- Extract this attribute from import state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        CALL POINT_PHYSICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!
      RC=0
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PHY INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY INITIALIZE STEP FAILED  RC_INIT=',RC_INIT
      ENDIF
!
      phy_init_tim=timef()-btim
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_RUN(GRID_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_RUN)
!
!-----------------------------------------------------------------------
!***  THE INTEGRATION OF THE MODEL PHYSICS IS DONE
!***  THROUGH THIS ROUTINE.
!-----------------------------------------------------------------------
!
!!!   USE MODULE_INTEGRATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: I,J,IRTN,ISTAT,JULDAY,JULYR,L               &
                           ,N,NPRECIP,NSTEPS_PREC,NTIMESTEP,RC          &
			   ,NTIMESTEP_RAD
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      REAL :: JULIAN,PDTOP,SECONDS_TOTAL,XTIME
!
      REAL,DIMENSION(LM) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(LM+1) :: PSG1,SG2
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1) :: RQVBLTEN,RTHBLTEN  ! For WRF physics

!
      LOGICAL :: CALL_LONGWAVE,CALL_PRECIP,CALL_SHORTWAVE,CALL_TURBULENCE
!
      TYPE(ESMF_Field)    :: HOLD_FIELD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS 
      MYPE=MYPE_SHARE
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from ATM Clock in Physics Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK                             &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,rc          =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF
      int_state%NTSD=NTIMESTEP
!
!-- Call radiation so that updated fields are written to the history files after 0 h.
!
      IF (NTIMESTEP == 0) THEN
         NTIMESTEP_RAD=NTIMESTEP
      ELSE
         NTIMESTEP_RAD=NTIMESTEP+1
      ENDIF
!
!
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  DEREFERENCE SOME INTERNAL STATE COMPONENTS FOR CONVENIENCE.
!-----------------------------------------------------------------------
!
      NPRECIP=int_state%NPRECIP
      PDTOP=int_state%PDTOP
      PT=int_state%PT
!
      DO L=1,LM
        DSG2(L)=int_state%DSG2(L)
        PDSG1(L)=int_state%PDSG1(L)
        PSGML1(L)=int_state%PSGML1(L)
        SGML2(L)=int_state%SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        SG2(L)=INT_STATE%SG2(L)
        PSG1(L)=INT_STATE%PSG1(L)
      ENDDO
!
!-----------------------------------------------------------------------
!***  UPDATE THE PHYSICS INTERNAL STATE WITH DATA FROM
!***  THE IMPORT STATE.  THIS MUST BE DONE EVERY TIME STEP
!***  SINCE THE TEMPERATURE IS UPDATED EVERY TIMESTEP.
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL UPDATE_INTERNAL_STATE_PHY(IMP_STATE,INT_STATE)
      update_phy_int_state_tim=update_phy_int_state_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  SET LOGICAL SWITCHES FOR CALLING EACH OF THE PHYSICS SCHEMES.
!-----------------------------------------------------------------------
!
      CALL_SHORTWAVE=MOD(NTIMESTEP_RAD,int_state%NRADS)==0
      CALL_LONGWAVE=MOD(NTIMESTEP_RAD,int_state%NRADL)==0
      CALL_TURBULENCE=MOD(NTIMESTEP,int_state%NPHS)==0
      CALL_PRECIP=MOD(NTIMESTEP,NPRECIP)==0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CALL THE INDIVIDUAL PHYSICAL PROCESSES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***      Call READPCP to
!***            1) READ IN PRECIPITATION FOR HOURS 1, 2 and 3;
!***            2) Initialize DDATA to 999. (this is the amount
!***               of input precip allocated to each physics time step
!***               in ADJPPT; TURBL/SURFCE, which uses DDATA, is called
!***               before ADJPPT)
!***            3) Initialize LSPA to zero
!***
!-----------------------------------------------------------------------
      IF (int_state%NTSD==0) THEN
        IF (int_state%PCPFLG) THEN
          CALL READPCP(MYPE,int_state%PPTDAT,int_state%DDATA,int_state%LSPA  &
     &      ,IDS,IDE,JDS,JDE,LM                                    &
     &      ,IMS,IME,JMS,JME                                    &
     &      ,ITS,ITE,JTS,JTE,int_state%PCPHR)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  RADIATION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  RADIATION NEEDS SOME SPECIFIC TIME QUANTITIES.
!
      CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR     &
                       ,START_MINUTE,START_SECOND                       &
                       ,NTIMESTEP,int_state%DT                          &
                       ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
      radiatn: IF(CALL_SHORTWAVE.OR.CALL_LONGWAVE)THEN
!-----------------------------------------------------------------------
!
        btim=timef()
!
!-----------------------------------------------------------------------
!***  EMPTY THE ACFRST AND ACFRCV ARRAYS IF IT IS TIME.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,NSTEPS_HIST)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACFRST(I,J)=0.
            int_state%ACFRCV(I,J)=0.
            int_state%NCFRST(I,J)=0
            int_state%NCFRCV(I,J)=0
          ENDDO
          ENDDO
        ENDIF

        CALL RADIATION(NTIMESTEP_RAD,int_state%DT,JULDAY,JULYR,XTIME,JULIAN &
                      ,START_HOUR,int_state%NPHS                        &
                      ,int_state%GLAT,int_state%GLON                    &
                      ,int_state%NRADS,int_state%NRADL                  &
                      ,DSG2,SGML2,PDSG1,PSGML1                          &
                      ,int_state%PT,int_state%PD                        &
                      ,int_state%T,int_state%Q,int_state%CW             &
                      ,int_state%THS,int_state%ALBEDO,int_state%EPSR    &
                      ,int_state%F_ICE,int_state%F_RAIN                 &
                      ,int_state%P_QV,int_state%P_QC,int_state%P_QR     &
                      ,int_state%P_QI,int_state%P_QS,int_state%P_QG     &
                      ,int_state%F_QV,int_state%F_QC,int_state%F_QR     &
                      ,int_state%F_QI,int_state%F_QS,int_state%F_QG     &
                      ,int_state%SM,int_state%CLDFRA                    &
                      ,int_state%NUM_WATER,int_state%WATER              &
                      ,int_state%RLWTT,int_state%RSWTT                  &
                      ,int_state%RLWIN,int_state%RSWIN                  &
                      ,int_state%RSWINC,int_state%RSWOUT                &
                      ,int_state%RLWTOA,int_state%RSWTOA                &
                      ,int_state%CZMEAN,int_state%SIGT4                 &
                      ,int_state%CFRACL,int_state%CFRACM                &
                      ,int_state%CFRACH                                 &
                      ,int_state%ACFRST,int_state%NCFRST                &
                      ,int_state%ACFRCV,int_state%NCFRCV                &
                      ,int_state%CUPPT,int_state%VEGFRC,int_state%SNO   &
                      ,int_state%HTOP,int_state%HBOT                    &
                      ,int_state%SHORTWAVE,int_state%LONGWAVE           &
                      ,LM)
!
        radiation_tim=radiation_tim+timef()-btim
!
      ENDIF radiatn
!
!-----------------------------------------------------------------------
!***  UPDATE THE TEMPERATURE WITH THE RADIATIVE TENDENCY.
!-----------------------------------------------------------------------
!
      btim=timef()
!
	if (ITS .le. 50 .and. ITE .ge. 50 .and. JTS .le. 50 .and. JTE .ge. 50) then
!	write(0,*) 'call RDTEMP with these arguments'
!	write(0,*) 'int_state%GLAT,int_state%GLON: ', int_state%GLAT(50,50),int_state%GLON(50,50) 
!	write(0,*) 'CZEN, CZMEAN, T, RSWTT, RLWTT: ', int_state%CZEN(50,50),int_state%CZMEAN(50,50),int_state%T(50,50,10), int_state%RSWTT(50,50,10),int_state%RLWTT(50,50,10)
	endif
      CALL RDTEMP(NTIMESTEP,int_state%DT,JULDAY,JULYR,START_HOUR        &
                 ,int_state%GLAT,int_state%GLON                         &
                 ,int_state%CZEN,int_state%CZMEAN,int_state%T           &
                 ,int_state%RSWTT,int_state%RLWTT                       &
                 ,LM)
	if (ITS .le. 50 .and. ITE .ge. 50 .and. JTS .le. 50 .and. JTE .ge. 50) then
!	write(0,*) 'int_state%T(50,50,10) after RDTEMP: ', int_state%T(50,50,10)
	endif

 973 continue
!
      rdtemp_tim=rdtemp_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  POLES AND EAST-WEST BOUNDARY.
!-----------------------------------------------------------------------
!
      IF(int_state%GLOBAL)THEN
        btim=timef()
!
        CALL SWAPHN(int_state%RSWIN,IMS,IME,JMS,JME,1,int_state%INPES)
        CALL POLEHN(int_state%RSWIN,IMS,IME,JMS,JME,1                   &
                   ,int_state%INPES,int_state%JNPES)
!
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                      &
                   ,int_state%INPES,int_state%JNPES)
!
        pole_swap_phy_tim=pole_swap_phy_tim+timef()-btim
      ENDIF
!
!-----------------------------------------------------------------------
!***  EMPTY THE RADIATION,FLUX ARRAYS IF IT IS TIME.
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP,int_state%NRDLW)==0)THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%ALWIN(I,J)=0.
          int_state%ALWOUT(I,J)=0.
          int_state%ALWTOA(I,J)=0.
        ENDDO
        ENDDO
!
!       int_state%ARDLW=0.   ! Precomputed
!
      ENDIF
!
      IF(MOD(NTIMESTEP,int_state%NRDSW)==0)THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%ASWIN(I,J)=0.
          int_state%ASWOUT(I,J)=0.
          int_state%ASWTOA(I,J)=0.
        ENDDO
        ENDDO
!
!       int_state%ARDSW=0.   ! Precomputed
!
      ENDIF
!
      IF(MOD(NTIMESTEP,int_state%NSRFC)==0)THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%SNOPCX(I,J)=0.
          int_state%POTEVP(I,J)=0.
          int_state%SFCEVP(I,J)=0.
          int_state%SFCLHX(I,J)=0.
          int_state%SFCSHX(I,J)=0.
          int_state%SUBSHX(I,J)=0.
          int_state%BGROFF(I,J)=0.
          int_state%SSROFF(I,J)=0.
        ENDDO
        ENDDO
!
!       int_state%ASRFC=0.   ! Precomputed
!
      ENDIF
!
!     IF(MOD(NTIMESTEP,int_state%NPHS)==0)THEN
!       int_state%APHTIM=0.  ! Precomputed
!     ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  TURBULENCE, SFC LAYER, AND LAND SURFACE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      turbulence: IF(CALL_TURBULENCE)THEN
!
        btim=timef()
!
        DO L=1,NUM_SOIL_LAYERS
          DZSOIL(L)=SLDPTH(L)
        ENDDO
!
        CALL TURBL(NTIMESTEP,int_state%DT,int_state%NPHS               &
                  ,int_state%NUM_WATER,NUM_SOIL_LAYERS,SLDPTH,DZSOIL   &
                  ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1,PT                 &
                  ,int_state%SM,int_state%CZEN,int_state%CZMEAN        &
                  ,int_state%SIGT4,int_state%RLWIN,int_state%RSWIN     &
                  ,int_state%RADOT                                     &
                  ,int_state%PD,int_state%T                            &
                  ,int_state%Q,int_state%CW                            &
                  ,int_state%F_ICE,int_state%F_RAIN,int_state%SR       &
                  ,int_state%Q2,int_state%U,int_state%V                &
                  ,int_state%DUDT,int_state%DVDT                       &
                  ,int_state%THS,int_state%TSKIN,int_state%SST         &
                  ,int_state%PREC,int_state%SNO                        &
                  ,int_state%WATER                                     &
                  ,int_state%P_QV,int_state%P_QC,int_state%P_QR        &
                  ,int_state%P_QI,int_state%P_QS,int_state%P_QG        &
                  ,int_state%F_QV,int_state%F_QC,int_state%F_QR        &
                  ,int_state%F_QI,int_state%F_QS,int_state%F_QG        &
                  ,int_state%FIS,int_state%Z0,int_state%Z0BASE         &
                  ,int_state%USTAR,int_state%PBLH,int_state%LPBL       &
                  ,int_state%XLEN_MIX,int_state%RMOL                   &
                  ,int_state%EXCH_H,int_state%AKHS,int_state%AKMS      &
                  ,int_state%AKHS_OUT,int_state%AKMS_OUT               &
                  ,int_state%THZ0,int_state%QZ0                        &
                  ,int_state%UZ0,int_state%VZ0                         &
                  ,int_state%UZ0H,int_state%VZ0H                       &
                  ,int_state%QSH,int_state%MAVAIL                      &
                  ,int_state%STC,int_state%SMC,int_state%CMC           &
                  ,int_state%SMSTAV,int_state%SMSTOT                   &
                  ,int_state%SSROFF,int_state%BGROFF                   &
                  ,int_state%IVGTYP,int_state%ISLTYP,int_state%VEGFRC  &
                  ,int_state%SHDMIN,int_state%SHDMAX,int_state%GRNFLX  &
                  ,int_state%SFCEXC,int_state%ACSNOW,int_state%ACSNOM  &
                  ,int_state%SNOPCX,int_state%SICE                     &
                  ,int_state%TG,int_state%SOILTB                       &
                  ,int_state%ALBASE,int_state%MXSNAL,int_state%ALBEDO  &
                  ,int_state%SH2O,int_state%SI,int_state%EPSR          &
                  ,int_state%U10,int_state%V10                         &
                  ,int_state%TH10,int_state%Q10                        &
                  ,int_state%TSHLTR,int_state%QSHLTR,int_state%PSHLTR  &
                  ,int_state%PSFC,int_state%T2                         &
                  ,int_state%QSG,int_state%QVG,int_state%QCG           &
                  ,int_state%SOILT1,int_state%TSNAV                    &
                  ,int_state%TWBS,int_state%QWBS                       &
                  ,int_state%SFCSHX,int_state%SFCLHX,int_state%SFCEVP  &
                  ,int_state%POTEVP,int_state%POTFLX,int_state%SUBSHX  &
                  ,int_state%APHTIM                                    &
                  ,int_state%ARDSW,int_state%ARDLW                     &
                  ,int_state%ASRFC                                     &
                  ,int_state%CROT,int_state%SROT                       &
                  ,int_state%HSTDV,int_state%HCNVX,int_state%HASYW     &
                  ,int_state%HASYS,int_state%HASYSW,int_state%HASYNW   &
                  ,int_state%HLENW,int_state%HLENS,int_state%HLENSW    &
                  ,int_state%HLENNW,int_state%HANGL,int_state%HANIS    &
                  ,int_state%HSLOP,int_state%HZMAX                     &
                  ,int_state%RSWOUT,int_state%RSWTOA,int_state%RLWTOA  &
                  ,int_state%ASWIN,int_state%ASWOUT,int_state%ASWTOA   &
                  ,int_state%ALWIN,int_state%ALWOUT,int_state%ALWTOA   &
                  ,int_state%RTHBLTEN,int_state%RQVBLTEN               &
                  ,int_state%GWDFLG,int_state%PCPFLG                   &
                  ,int_state%DDATA,int_state%UCMCALL                   &
                  ,int_state%TURBULENCE,int_state%SFC_LAYER            &
                  ,int_state%LAND_SURFACE,int_state%LONGWAVE           &
                  ,int_state%MICROPHYSICS                              &
                  ,IDS,IDE,JDS,JDE,LM                                  &
                  ,IMS,IME,JMS,JME                                     &
                  ,ITS,ITE,JTS,JTE)


!
        turbl_tim=turbl_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  EXCHANGE WIND TENDENCIES.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%UZ0H,1,int_state%VZ0H,1,1,1)
        CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,1,1)
!
        exch_phy_tim=exch_phy_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  NOW INTERPOLATE WIND TENDENCIES
!***  FROM H TO V POINTS.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL H_TO_V(int_state%UZ0H,int_state%UZ0)
        CALL H_TO_V(int_state%VZ0H,int_state%VZ0)
        CALL H_TO_V_TEND(int_state%DUDT,int_state%DT,int_state%NPHS,LM  &
                        ,int_state%U)
        CALL H_TO_V_TEND(int_state%DVDT,int_state%DT,int_state%NPHS,LM  &
                        ,int_state%V)
!
        h_to_v_tim=h_to_v_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  POLES AND EAST-WEST BOUNDARY.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                   &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM                   &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM        &
                     ,int_state%INPES,int_state%JNPES)
!
          pole_swap_phy_tim=pole_swap_phy_tim+timef()-btim
        ENDIF
!
!-----------------------------------------------------------------------
!***  EXCHANGE WIND COMPONENTS AND TKE.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                    &
                      ,2,2)
!
        CALL HALO_EXCH(int_state%UZ0,1,int_state%VZ0,1                  &
                      ,int_state%Q2,LM                                  &
                      ,1,1)
!
        exch_phy_tim=exch_phy_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF turbulence
!
!----------------------------------------------------------------------- 
!***  EMPTY THE PRECIPITATION ARRAYS IF IT IS TIME.                      
!-----------------------------------------------------------------------
!
      NSTEPS_PREC=int_state%NHRS_PREC*NSTEPS_PER_HOUR
      IF(MOD(NTIMESTEP+1-int_state%NPRECIP,NSTEPS_PREC)==0)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%ACPREC(I,J)=0.
          int_state%CUPREC(I,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CONVECTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      convection: IF(CALL_PRECIP.AND.int_state%CONVECTION/='none')THEN
!
        btim=timef()
!
        CALL CUCNVC(NTIMESTEP,int_state%DT,int_state%NPRECIP           &
                   ,int_state%NRADS,int_state%NRADL                    &
                   ,int_state%NHOURS_HISTORY                           &
                   ,int_state%DYH,int_state%RESTART,int_state%HYDRO    &
                   ,int_state%CLDEFI,int_state%NUM_WATER               &
                   ,int_state%F_ICE,int_state%F_RAIN                   &
                   ,int_state%P_QV,int_state%P_QC,int_state%P_QR       &
                   ,int_state%P_QI,int_state%P_QS,int_state%P_QG       &
                   ,int_state%F_QV,int_state%F_QC,int_state%F_QR       &
                   ,int_state%F_QI,int_state%F_QS,int_state%F_QG       &
                   ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1                   &
                   ,int_state%PT,int_state%PD                          &
                   ,int_state%T,int_state%Q                            &
                   ,int_state%CW,int_state%TCUCN,int_state%WATER       &
                   ,int_state%OMGALF                                   &
                   ,int_state%U,int_state%V                            &
                   ,int_state%FIS,int_state%W0AVG                      &
                   ,int_state%PREC,int_state%ACPREC,int_state%CUPREC   &
                   ,int_state%CUPPT,int_state%CPRATE                   &
                   ,int_state%CNVBOT,int_state%CNVTOP                  &
                   ,int_state%SM,int_state%LPBL                        &
                   ,int_state%HTOP,int_state%HTOPD,int_state%HTOPS     &
                   ,int_state%HBOT,int_state%HBOTD,int_state%HBOTS     &
                   ,int_state%AVCNVC,int_state%ACUTIM                  &
                   ,int_state%RSWIN,int_state%RSWOUT                   &
                   ,int_state%CONVECTION                               &
!mep               ,IDS,IDE-1,JDS,JDE-1,LM                             &
                   ,IDS,IDE,JDS,JDE,LM                                 &
                   ,IMS,IME,JMS,JME                                    &
                   ,ITS,ITE,JTS,JTE)
!
        cucnvc_tim=cucnvc_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  POLES AND EAST-WEST BOUNDARY.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          pole_swap_phy_tim=pole_swap_phy_tim+timef()-btim
        ENDIF
!
      ENDIF convection
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  MICROPHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      microphysics: IF(CALL_PRECIP)THEN

!
        btim=timef()
!
        CALL GSMDRIVE(NTIMESTEP,int_state%DT                            &
                     ,NPRECIP,int_state%NUM_WATER                       &
                     ,int_state%DXH(JC),int_state%DYH                   &
                     ,int_state%SM,int_state%FIS                        &
                     ,DSG2,SGML2,PDSG1,PSGML1                           &
                     ,int_state%PT,int_state%PD                         &
                     ,int_state%T,int_state%Q                           &
                     ,int_state%CW,int_state%OMGALF                     &
                     ,int_state%WATER                                   &
                     ,int_state%TRAIN,int_state%SR                      &
                     ,int_state%F_ICE,int_state%F_RAIN,int_state%F_RIMEF &
                     ,int_state%P_QV,int_state%P_QC,int_state%P_QR      &
                     ,int_state%P_QI,int_state%P_QS,int_state%P_QG      &
                     ,int_state%F_QV,int_state%F_QC,int_state%F_QR      &
                     ,int_state%F_QI,int_state%F_QS,int_state%F_QG      &
                     ,int_state%PREC,int_state%ACPREC,int_state%AVRAIN  &
                     ,int_state%MP_RESTART_STATE                        &
                     ,int_state%TBPVS_STATE,int_state%TBPVS0_STATE      &
                     ,int_state%SPECIFIED,int_state%NESTED              &
                     ,int_state%MICROPHYSICS                            &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
!
        gsmdrive_tim=gsmdrive_tim+timef()-btim

!-----------------------------------------------------------------------
!---------PRECIPITATION ASSIMILATION------------------------------------
!-----------------------------------------------------------------------
!
        IF (int_state%PCPFLG) THEN
!
        btim=timef()
          CALL CHKSNOW(MYPE,int_state%NTSD,int_state%DT,int_state%NPHS,int_state%SR,int_state%PPTDAT                 &
     &      ,IDS,IDE,JDS,JDE,LM                                    &
     &      ,IMS,IME,JMS,JME                                    &
     &      ,ITS,ITE,JTS,JTE,int_state%PCPHR)
          CALL ADJPPT(MYPE,int_state%NTSD,int_state%DT,int_state%NPHS,int_state%PREC,int_state%LSPA,int_state%PPTDAT,int_state%DDATA     &
     &      ,IDS,IDE,JDS,JDE,LM                                    &
     &      ,IMS,IME,JMS,JME                                    &
     &      ,ITS,ITE,JTS,JTE,int_state%PCPHR)
!
          adjppt_tim=adjppt_tim+timef()-btim
        ENDIF
!

!
!-----------------------------------------------------------------------
!***  POLES AND EAST-WEST BOUNDARY.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                    &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                   &
                     ,int_state%INPES,int_state%JNPES)
!
          pole_swap_phy_tim=pole_swap_phy_tim+timef()-btim
        ENDIF
!
!-----------------------------------------------------------------------
!***  EXCHANGE Q AND CW.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL HALO_EXCH(int_state%Q,LM,int_state%CW,LM                   &
                      ,1,1)
!
        exch_phy_tim=exch_phy_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF microphysics
!
!-----------------------------------------------------------------------
!***  ALWAYS EXCHANGE TEMPERATURE ARRAY SINCE RADIATIVE UPDATES
!***  ARE DONE EVERY TIMESTEP.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL HALO_EXCH(int_state%T,LM                                     &
                    ,1,1)
!
      exch_phy_tim=exch_phy_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  NOTE:  THE PHYSICS EXPORT STATE IS FULLY UPDATED NOW
!***         BECAUSE SUBROUTINE PHY_INITIALIZE INSERTED THE
!***         APPROPRIATE ESMF Fields INTO IT.  THOSE FIELDS
!***         CONTAIN POINTERS TO THE ACTUAL DATA AND THOSE
!***         POINTERS ARE NEVER RE-DIRECTED.
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PHY RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      phy_run_tim=phy_run_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_FINALIZE(GRID_COMP                                 &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK_ATM                                 &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE PHYSICS COMPONENT.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Physics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE                 !<-- The Physics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE_WRITE                 !<-- The Physics export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM component's ESMF Clock.
!
      INTEGER            ,INTENT(OUT)   :: RCFINAL
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      MYPE=MYPE_SHARE
!
      IF(MYPE==0)THEN
        WRITE(0,*)' Physics Completed Normally.'
      ENDIF
!
!-----------------------------------------------------------------------
!***  DO NOT DEALLOCATE THE PHYSICS INTERNAL STATE POINTER
!***  WITHOUT DEALLOCATING ITS CONTENTS.
!-----------------------------------------------------------------------
!
!!!   DEALLOCATE(INT_STATE,stat=RC)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_FINALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_INITIALIZE(GFS                                 &
                                   ,SHORTWAVE                           &
                                   ,LONGWAVE                            &
                                   ,CONVECTION                          &
                                   ,MICROPHYSICS                        &
                                   ,SFC_LAYER                           &
                                   ,TURBULENCE                          &
                                   ,LAND_SURFACE                        &
                                   ,CO2TF                               &
                                   ,IDS,IDE,JDS,JDE,LM                  &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE)
!
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,ONLY : A,CLIQ,CV,DTR,PI                      &
                                 ,RHOAIR0,RHOWATER,RHOSNOW
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: CO2TF
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE
!
      LOGICAL,INTENT(IN) :: GFS
!
      CHARACTER(99),INTENT(IN) :: CONVECTION,LONGWAVE,MICROPHYSICS      &
                                 ,SFC_LAYER,SHORTWAVE,TURBULENCE        &
                                 ,LAND_SURFACE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,I_HI,I_LO,IHRST,II,IRTN,J,J_HI,J_LO,JJ,JULDAY,JULYR  &
                ,K,KFLIP,L,LPT2,N,NFCST,NRECS_SKIP_FOR_PT               &
                ,NSOIL,NSTEPS_PER_HOUR,NTIMESTEP
!
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
!
      INTEGER :: IYEAR_FCST,IMONTH_FCST,IDAY_FCST,IHOUR_FCST            &
                ,IMINUTE_FCST
!
      INTEGER,DIMENSION(3) :: IDAT
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ITEMP,LOWLYR
!
      REAL :: SECOND_FCST
!
      REAL :: SWRAD_SCAT=1.
!
      REAL :: ALM,ANUM,APH,AVE,CTLM,CTPH,CTPH0,DELX,DELY,DENOM          &
             ,DLM,DLMD,DPH,DPHD,DT,DT_MICRO,DTPHS                       &
             ,GMT,JULIAN,PDBOT,PDTOP,PDTOT,PT_CB,RELM,RPDTOT            &
             ,SB,SPH,STLM,STPH,STPH0,THETA_HALF                         &
             ,TLM,TLM_BASE,TPH,TPH_BASE,TPH0,TPV,WB,XTIME
!
      REAL,DIMENSION(LM) :: DSG1,DSG2,PDSG1,PSGML1,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM                         &
                             ,SFULL,SFULL_FLIP,SMID,SMID_FLIP
!
!zj      REAL,DIMENSION(:),ALLOCATABLE,TARGET :: DXH,DXV
      REAL,DIMENSION(:),ALLOCATABLE,TARGET :: DXH,DXV,RDXH,RDXV !zj
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: EMISS
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP1,TEMP_GWD
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: TEMPSOIL
      REAL,DIMENSION(NUM_SOIL_LAYERS)   :: SOIL1DIN
!
      CHARACTER(ESMF_MAXSTR) :: INFILE
!
      LOGICAL,SAVE :: ALLOWED_TO_READ=.TRUE.
      LOGICAL :: OPENED,RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NSOIL=NUM_SOIL_LAYERS                                              !<-- From Landsurface module
!
!-----------------------------------------------------------------------
!***  DEREFERENCE THE START TIME.
!-----------------------------------------------------------------------
!
      START_YEAR=int_state%START_YEAR
      START_MONTH=int_state%START_MONTH
      START_DAY=int_state%START_DAY
      START_HOUR=int_state%START_HOUR
      START_MINUTE=int_state%START_MINUTE
      START_SECOND=int_state%START_SECOND
      DT=int_state%DT
!
!***  RADIATION NEEDS SOME SPECIFIC TIME QUANTITIES.
!
      CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR     &
                       ,START_MINUTE,START_SECOND                       &
                       ,NTIMESTEP,DT                                    &
                       ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
! *** OPEN AND READ GWD DATA FILE (14 OROGRAPHY FIELDS)
!-----------------------------------------------------------------------
!
      gwd_read: IF(int_state%GWDFLG) THEN
!
        select_GWD_unit: DO N=51,59
          INQUIRE(N,OPENED=OPENED)
          IF(.NOT.OPENED)THEN
            NFCST=N
            EXIT select_GWD_unit
          ENDIF
        ENDDO select_GWD_unit
!
        ALLOCATE(TEMP_GWD(IDS:IDE,JDS:JDE))
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          OPEN(unit=NFCST,file='GWD.bin',status='old',form='unformatted')
        ENDIF
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HSTDV,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HCNVX,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HASYW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HASYS,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HASYSW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HASYNW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HLENW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HLENS,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HLENSW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HLENNW,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HANGL,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HANIS,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HSLOP,1,1,1,1,1)
!-----------------------------------------------------------------------
        IF(MYPE==0)THEN
          READ(NFCST)TEMP_GWD
        ENDIF
!
        CALL DSTRB(TEMP_GWD,int_state%HZMAX,1,1,1,1,1)
!-----------------------------------------------------------------------
!
        IF(MYPE==0)THEN
          CLOSE(NFCST)
        ENDIF
!
        DEALLOCATE(TEMP_GWD)
!-----------------------------------------------------------------------
!
      ENDIF gwd_read
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  OPEN THE INPUT DATA FILE
!-----------------------------------------------------------------------
!
      select_unit: DO N=51,59
        INQUIRE(N,OPENED=OPENED)
        IF(.NOT.OPENED)THEN
          NFCST=N
          EXIT select_unit
        ENDIF
      ENDDO select_unit

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
         IF(.NOT.int_state%RESTART) THEN           ! COLD START

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      INFILE='main_input_filename'
!
!-----------------------------------------------------------------------
!***  FIRST WE NEED THE VALUE OF PT (PRESSURE AT TOP OF DOMAIN)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted')
!        NRECS_SKIP_FOR_PT=6+5*LM+21 !<-- For current WPS input
        NRECS_SKIP_FOR_PT=6+5*LM+23 !zj +21
!
        DO N=1,NRECS_SKIP_FOR_PT
          READ(NFCST)
        ENDDO
!
        READ(NFCST)PT
        WRITE(0,*)' PT FROM INPUT FILE EQUALS ',PT
        int_state%PT=PT
        REWIND NFCST
      ENDIF
!
      CALL MPI_BCAST(int_state%PT,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      PT=int_state%PT
      PT_CB=PT*1.E-3   !<-- Convert pascals to centibars for GFDL initialization
!
!-----------------------------------------------------------------------
!***  VERTICAL LAYER INFORMATION IS NEEDED IN ORDER TO SEND IT TO
!***  SOME SPECIFIC SCHEMES' INITIALIZATION ROUTINES THAT FOLLOW
!***  BELOW.
!-----------------------------------------------------------------------
!

      IF(MYPE==0)THEN
        READ(NFCST)RUN,IDAT,IHRST
        READ(NFCST)PT,PDTOP,LPT2,SGM,SG1,DSG1,SGML1,SG2,DSG2,SGML2

	write(0,*) 'PDTOP, SG1: ', PDTOP, SG1
	write(0,*) 'SG2: ', SG2
!
!***  CHECK TO SEE IF THE STARTING DATE/TIME IN THE INPUT DATA FILE
!***  AGREES WITH THAT IN THE CONFGURE FILE.
!
        IF(IDAT(2)/=START_MONTH.OR.                                     &
           IDAT(1)/=START_DAY.OR.                                       &
           IDAT(3)/=START_YEAR.OR.                                      &
           IHRST  /=START_HOUR)THEN
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' DATES IN INPUT FILE AND CONFIGURE FILE DISAGREE!!'
          WRITE(0,*)' INPUT: HOUR=',IHRST,' DAY=',IDAT(1)               &
                    ,' MONTH=',IDAT(2),' YEAR=',IDAT(3)
          WRITE(0,*)' CONFIG: HOUR=',START_HOUR,' DAY=',START_DAY       &
                    ,' MONTH=',START_MONTH,' YEAR=',START_YEAR
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
        ENDIF
      ENDIF
!
      CALL MPI_BCAST(SGM(1),LM+1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG1(1),LM+1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG1(1),LM,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML1(1),LM,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG2(1),LM+1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG2(1),LM,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML2(1),LM,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(PDTOP,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(LPT2,1,MPI_INTEGER,0,MPI_COMM_COMP,IRTN)
!
      CALL MPI_BCAST(IDAT         ,3,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHRST        ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!
      IYEAR_FCST  = IDAT(3)
      IMONTH_FCST = IDAT(2)
      IDAY_FCST   = IDAT(1)
      IHOUR_FCST  = IHRST
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
!       WRITE(0,*)' PHYSICS_INITIALIZE l=',l,' pdsg1=',pdsg1(l),' dsg1=',dsg1(l) &
!               ,' psgml1=',psgml1(l),' sgml1=',sgml1(l),' pdtop=',pdtop,' pt=',pt
      ENDDO
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
!
!-----------------------------------------------------------------------
!***  BEFORE MOVING ON, TRANSFER VALUES TO THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
        int_state%SGML2(L)=SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
      ALLOCATE(TEMP1(IDS:IDE,JDS:JDE),STAT=I)
!
!-----------------------------------------------------------------------
!***  PROCEED WITH GETTING FIELDS FROM INPUT FILE.
!***  NOTE: TWO RECORDS WERE ALREADY READ AT THE TOP OF THIS ROUTINE.
!-----------------------------------------------------------------------
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3) !zj
!
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, CW
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! U
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! V
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! T
	write(0,*) 'min max of T read in: ', K, minval(TEMP1), maxval(TEMP1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!
      DO K=1,LM
        JJ=LDIM2-1
        DO J=JMS,JME
          JJ=JJ+1
          II=LDIM1-1
          DO I=IMS,IME
            II=II+1
            int_state%WATER(II,JJ,K,int_state%P_QV)=                      & ! WRF water array uses mixing ratio for vapor
                      int_state%Q(II,JJ,K)/(1.-int_state%Q(II,JJ,K))     
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CWM
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO 
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!
! **** EPSR
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max of EPSR: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!*** SNOW ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max of MXSNAL: ', minval(TEMP1), maxval(TEMP1)
      ENDIF

      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SST/TSK
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1        ! actually NMM_TSK from WRF
	write(0,*) 'min, max of TSKIN: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1        ! actually NMM_TSK from WRF
	write(0,*) 'min, max of SST: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SNO, SICE, STC, SMC, ISLTYP, IVGTYP, VEGFRC
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max of SNO: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max for SR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max for USTAR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
	write(0,*) 'min, max for Z0: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
!zj      go to 11111  !<-- For current WPS input
!
      IF(MYPE==0)THEN !zj
        READ(NFCST)TEMP1 !zj
        write(0,*) 'min, max for Z0BASE: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF !zj
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1) !zj
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3) !zj
!
      IF(MYPE==0)THEN !zj
        READ(NFCST)TEMP1 !zj
        write(0,*) 'min, max for STDH: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF !zj
      CALL DSTRB(TEMP1,int_state%STDH,1,1,1,1,1) !zj
      CALL HALO_EXCH(int_state%STDH,1,3,3) !zj
!
11111 continue
!
      ALLOCATE(TEMPSOIL(NSOIL,IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
	write(0,*) 'min, max for STC: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,4),1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
        write(0,*) 'min, max for SMC: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SMC(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SMC(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SMC(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SMC(:,:,4),1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
        write(0,*) 'min, max for SH2O: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SH2O(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SH2O(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SH2O(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SH2O(:,:,4),1,1,1,1,1)
!
      DEALLOCATE(TEMPSOIL)
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
        write(0,*) 'min, max for ISLTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
        write(0,*) 'min, max for IVGTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
        write(0,*) 'min, max for VEGFRC: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST) SOIL1DIN
	write(0,*) 'SOIL1DIN: ', SOIL1DIN
      ENDIF
!
      IF(MYPE==0)THEN
        READ(NFCST) SOIL1DIN
	write(0,*) 'SOIL1DIN: ', SOIL1DIN
      ENDIF
!
      DO N=1,NSOIL
        int_state%SLDPTH(N)=SLDPTH(N)
      ENDDO
!
      IF(MYPE==0)THEN
        READ(NFCST) PT
        write(0,*) 'read in ptop in PHYS: ', PT
      ENDIF
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
      CLOSE(NFCST)
!----------------------------------------------------------------------
!
      DEALLOCATE(ITEMP)
      DEALLOCATE(TEMP1)
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
         ELSE                                      ! RESTART READ INIT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      INFILE='restart_file'
!
!-----------------------------------------------------------------------
!***  FIRST WE NEED THE VALUE OF PT (PRESSURE AT TOP OF DOMAIN)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted')
      ENDIF
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: INTEGER SCALARS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) IYEAR_FCST
        READ(NFCST) IMONTH_FCST
        READ(NFCST) IDAY_FCST
        READ(NFCST) IHOUR_FCST
        READ(NFCST) IMINUTE_FCST
        READ(NFCST) SECOND_FCST
        READ(NFCST) ! NTSD
        READ(NFCST) ! IM
        READ(NFCST) ! JM
        READ(NFCST) ! LM
        READ(NFCST) IHRST
        READ(NFCST) LPT2
      ENDIF
!
      IF(MYPE==0)THEN
        write(0,*)'**** read in physics ****************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        write(0,*)'*************************************'
      ENDIF
!
      CALL MPI_BCAST(IYEAR_FCST   ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IMONTH_FCST  ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IDAY_FCST    ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHOUR_FCST   ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IMINUTE_FCST ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SECOND_FCST  ,1,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHRST        ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: INTEGER 1D ARRAYS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) IDAT
      ENDIF
      CALL MPI_BCAST(IDAT,3,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: INTEGER SCALARS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! MP_PHYSICS
        READ(NFCST) ! SF_SURFACE_PHYSICS
        READ(NFCST) NSOIL
        READ(NFCST) ! NPHS
        READ(NFCST) ! NCLOD
        READ(NFCST) ! NHEAT
        READ(NFCST) ! NPREC
        READ(NFCST) ! NRDLW
        READ(NFCST) ! NRDSW
        READ(NFCST) ! NSRFC
      ENDIF
      CALL MPI_BCAST(NSOIL,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL SCALARS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! DT
        READ(NFCST) ! DYH
        READ(NFCST) PDTOP
      ENDIF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)PT
        int_state%PT=PT
        READ(NFCST) ! TLM0D
        READ(NFCST) ! TPH0D
        READ(NFCST) ! TSTART
        READ(NFCST) ! DPHD
        READ(NFCST) ! DLMD
      ENDIF
!
      CALL MPI_BCAST(int_state%PT,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      PT=int_state%PT
      PT_CB=PT*1.E-3   !<-- Convert pascals to centibars for GFDL initialization
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL 1D ARRAYS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! DXH
        READ(NFCST) SG1
        READ(NFCST) SG2
        READ(NFCST) DSG1
        READ(NFCST) DSG2
        READ(NFCST) SGML1
        READ(NFCST) SGML2
        READ(NFCST) SGM
        READ(NFCST) ! APHTIM
        READ(NFCST) ! ARDLW
        READ(NFCST) ! ARDSW
        READ(NFCST) ! ASRFC
        READ(NFCST) ! AVCNVC
        READ(NFCST) ! AVRAIN
        READ(NFCST) SLDPTH
        READ(NFCST) int_state%MP_RESTART_STATE
        READ(NFCST) int_state%TBPVS_STATE
        READ(NFCST) int_state%TBPVS0_STATE
      ENDIF
!
      CALL MPI_BCAST(SGM(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG1(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG1(1)  ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML1(1) ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG2(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG2(1)  ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML2(1) ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(PDTOP    ,1     ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(LPT2     ,1     ,MPI_INTEGER,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SLDPTH   ,NSOIL ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)


      CALL MPI_BCAST(int_state%MP_RESTART_STATE(1) ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TBPVS_STATE(1)      ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TBPVS0_STATE(1)     ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
!
      DO N=1,NSOIL
        int_state%SLDPTH(N)=SLDPTH(N)
      ENDDO
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: LOGICAL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! GLOBAL
        READ(NFCST) RUN
        READ(NFCST) ! ADIABATIC
      ENDIF
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
!
!-----------------------------------------------------------------------
!***  BEFORE MOVING ON, TRANSFER VALUES TO THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
        int_state%SGML2(L)=SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: INTEGER 2D ARRAYS
!-----------------------------------------------------------------------
!
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRCV)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRST)
!
      DEALLOCATE(ITEMP)
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL 2D ARRAYS
!-----------------------------------------------------------------------
!
      ALLOCATE(TEMP1(IDS:IDE,JDS:JDE),STAT=I)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) !GLAT
        READ(NFCST) !GLON
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) !VLAT
        READ(NFCST) !VLON
        READ(NFCST)TEMP1    ! PDO
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!              SKIP FROM RESTART FILE: REAL 3D ARRAYS (only from DYN)
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! W
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! DWDT
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! PINT
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! OMGALF
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RRW
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! DIV
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RTOP
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCU
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCV
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCT
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TP
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! UP
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! VP
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! E2
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM-1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! PSGDT
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Z
        ENDIF
      ENDDO 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL 2D ARRAYS (contd.)
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRCV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKHS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKMS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%BGROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CLDEFI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CPRATE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPPT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZEN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZMEAN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%GRNFLX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PBLH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Q10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RADOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWINC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEXC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCLHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SIGT4,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTAV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOPCX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SOILTB,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SSROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SUBSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TH10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%U10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%UZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%V10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWTOA,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTFLX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RMOL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%T2,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL 3D ARRAYS
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, Q2, CW, F_ICE, F_RIMEF, F_RAIN, RTHBLTEN, RQVBLTEN
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CLDFRA
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CLDFRA,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CWM
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
!      DO K=1,LM
!        IF(MYPE==0)THEN
!          READ(NFCST)TEMP1   ! EXCH_H
!        ENDIF
!      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q2
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q2,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RLWTT
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RLWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RSWTT
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RSWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! T
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! TCUCN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TCUCN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! TRAIN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! U
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! V
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! XLEN_MIX
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%XLEN_MIX,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_ICE
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_ICE,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RIMEF
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RIMEF,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RAIN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RAIN,1,1,1,LM,K)
      ENDDO 
!-----------------------------------------------------------------------
!
      DO K=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RTHBLTEN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RTHBLTEN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RTHBLTEN,1,1,1,LM+1,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RQVBLTEN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RQVBLTEN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RQVBLTEN,1,1,1,LM+1,K)
      ENDDO
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SH2O: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SH2O,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SMC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SMC,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for STC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%STC,1,1,1,NSOIL,K)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  TRACERS
!-----------------------------------------------------------------------
      DO N=1,int_state%NUM_TRACERS_TOTAL
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TRACERS
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRACERS(I,J,K,N)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRACERS(:,:,:,N),1,1,1,LM,K)
      ENDDO 
      ENDDO 
!rv do not use HALO_EXCH
!-----------------------------------------------------------------------
!
      DEALLOCATE(TEMP1)
!
!-----------------------------------------------------------------------
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
         CLOSE(NFCST)
      ENDIF
!-----------------------------------------------------------------------
!
         ENDIF                                     ! COLD START /RESTART
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  MAKE UP A POTENTIAL SKIN TEMPERATURE.
!-----------------------------------------------------------------------
!
      IF(.NOT.int_state%RESTART) THEN
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        int_state%THS(I,J)=int_state%TSKIN(I,J)                        &
                          *(100000./(SG2 (LM+1)*int_state%PD(I,J)      &
                                    +PSG1(LM+1)))**CAPPA
      ENDDO
      ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  RECREATE SIGMA VALUES AT LAYER INTERFACES FOR THE FULL VERTICAL
!***  DOMAIN. 
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        SFULL(L)=SGM(L)
      ENDDO
!
      DO L=1,LM
        SMID(L)=(SFULL(L)+SFULL(L+1))*0.5
      ENDDO
      SMID(LM+1)=-9999999.
!
!-----------------------------------------------------------------------
!***  THE RADIATIVE EMISSIVITY
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        EMISS(I,J)=1.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      ALLOCATE(DXH(JDS:JDE),STAT=I)
      ALLOCATE(DXV(JDS:JDE),STAT=I)
!
      ALLOCATE(RDXH(JDS:JDE),STAT=I) !zj
      ALLOCATE(RDXV(JDS:JDE),STAT=I) !zj
!
!----------------------------------------------------------------------
!***  GEOGRAPHIC LATITUDE/LONGITUDE
!----------------------------------------------------------------------
!
      SB=int_state%SBD*DTR
      WB=int_state%WBD*DTR
      TPH0=int_state%TPH0D*DTR
      STPH0=SIN(TPH0)
      CTPH0=COS(TPH0)
!
      IF(int_state%GLOBAL)THEN
!
        I_LO=MAX(IMS,IDS)
        I_HI=MIN(IME,IDE)
        J_LO=MAX(JMS,JDS)
        J_HI=MIN(JME,JDE)
!
        DPHD=-int_state%SBD*2./REAL(JDE-3)
        DLMD=-int_state%WBD*2./REAL(IDE-3)
        DPH=DPHD*DTR
        DLM=DLMD*DTR
!        TPH_BASE=SB-DPH-DPH
!
!        DO J=J_LO,J_HI
!          TPH=TPH_BASE+(J-JDS+1)*DPH
!          STPH=SIN(TPH)
!          CTPH=COS(TPH)
!
!          TLM_BASE=WB-DLM
!          DO I=I_LO,I_HI
!            TLM=TLM_BASE+(I-IDS+1)*DLM
!            STLM=SIN(TLM)
!            CTLM=COS(TLM)
!            SPH=CTPH0*STPH+STPH0*CTPH*CTLM
!            APH=ASIN(SPH)
!            int_state%GLAT(I,J)=APH
!            ANUM=CTPH*STLM
!            DENOM=(CTLM*CTPH-STPH0*SPH)/CTPH0
!            RELM=ATAN2(ANUM,DENOM)
!            ALM=RELM+int_state%TLM0D*DTR
!            IF(ALM>PI)ALM=ALM-PI-PI
!            IF(ALM<-PI)ALM=ALM+PI+PI
!            int_state%GLON(I,J)=ALM
         tph_base=sb-dph-dph
         do j=j_lo,j_hi
           tlm_base=wb-dlm-dlm
           aph=tph_base+(j-jds+1)*dph
           do i=i_lo,i_hi
             alm=tlm_base+(i-ids+1)*dlm
             if(alm> pi) alm=alm-pi-pi
             if(alm<-pi) alm=alm+pi+pi
             int_state%GLAT(I,J)=aph
             int_state%GLON(I,J)=alm
           enddo
         enddo
!
      ELSE  ! regional

        DPHD=-int_state%SBD*2./REAL(JDE-1)
        DLMD=-int_state%WBD*2./REAL(IDE-1)
        DPH=DPHD*DTR
        DLM=DLMD*DTR
        TPH_BASE=SB-DPH
!
        DO J=JTS,JTE
          TPH=TPH_BASE+(J-JDS+1)*DPH
          STPH=SIN(TPH)
          CTPH=COS(TPH)
!
          TLM_BASE=WB-DLM
          DO I=ITS,ITE
            TLM=TLM_BASE+(I-IDS+1)*DLM
            STLM=SIN(TLM)
            CTLM=COS(TLM)
            SPH=CTPH0*STPH+STPH0*CTPH*CTLM
            APH=ASIN(SPH)
            int_state%GLAT(I,J)=APH
            ANUM=CTPH*STLM
            DENOM=(CTLM*CTPH-STPH0*SPH)/CTPH0
            RELM=ATAN2(ANUM,DENOM)
            ALM=RELM+int_state%TLM0D*DTR
            IF(ALM>PI)ALM=ALM-PI-PI
            IF(ALM<-PI)ALM=ALM+PI+PI
            int_state%GLON(I,J)=ALM
          ENDDO
        ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------
!***  DELTA X AND Y
!----------------------------------------------------------------------
!
      int_state%DYH=A*DPH
      int_state%DYV=A*DPH
!
!----------------------------------------------------------------------
      global_regional_dx: IF(int_state%GLOBAL)THEN
!----------------------------------------------------------------------
        TPH=SB
        TPV=TPH+DPH*0.5
!
!----------------------------------------------------------------------
!***  SOUTH POLE
!----------------------------------------------------------------------
!
        DXH(JDS+1)=0.
        RDXH(JDS+1)=0.
        DXV(JDS+1)=A*DLM*COS(TPV)
!
!----------------------------------------------------------------------
!***  BETWEEN THE POLES
!
!----------------------------------------------------------------------
!
        DO J=JDS+2,JDE-2
          TPH=SB+(J-JDS-1)*DPH
          TPV=TPH+DPH*0.5
          DXH(J)=A*DLM*COS(TPH)
          DXV(J)=A*DLM*COS(TPV)
          RDXH(J)=1./DXH(J)
          RDXV(J)=1./DXV(J)
        ENDDO
!
!-----------------------------------------------------------------------
!***  GHOST LINE BEYOND SOUTH POLE
!-----------------------------------------------------------------------
!
        DXH(JDS)=DXH(JDS+2)
        DXV(JDS)=DXV(JDS+1)
        RDXH(JDS)=RDXH(JDS+2)
        RDXV(JDS)=RDXV(JDS+1)
!
!-----------------------------------------------------------------------
!***  NORTH POLE
!-----------------------------------------------------------------------
!
        DXH(JDE-1)=0.
        RDXH(JDE-1)=0.
!
!-----------------------------------------------------------------------
!***  GHOST LINE BEYOND NORTH POLE
!-----------------------------------------------------------------------
!
        DXH(JDE)=DXH(JDE-2)
        DXV(JDE-1)=DXV(JDE-2)
        DXV(JDE)=DXV(JDE-2)
        RDXH(JDE)=RDXH(JDE-2)
        RDXV(JDE-1)=RDXV(JDE-2)
        RDXV(JDE)=RDXV(JDE-2)
!
!-----------------------------------------------------------------------
!***  AVERAGE OVER HEIGHT LATITUDES FOR ACCURACY.
!-----------------------------------------------------------------------
!
        DO J=JDS,JDE/2
          AVE=(DXH(J)+DXH(JDE+1-J))*0.5
          DXH(J)=AVE
          DXH(JDE+1-J)=AVE
          AVE=(RDXH(J)+RDXH(JDE+1-J))*0.5
          RDXH(J)=AVE
          RDXH(JDE+1-J)=AVE
        ENDDO
!
!-----------------------------------------------------------------------
!***  AVERAGE OVER WIND LATITUDES FOR ACCURACY.
!-----------------------------------------------------------------------
!
        DO J=JDS,(JDE-1)/2
          AVE=(DXV(J)+DXV(JDE-J))*0.5
          DXV(J)=AVE
          DXV(JDE-J)=AVE
          AVE=(RDXV(J)+RDXV(JDE-J))*0.5
          RDXV(J)=AVE
          RDXV(JDE-J)=AVE
        ENDDO
!
!-----------------------------------------------------------------------
      ELSE global_regional_dx  ! Regional
!-----------------------------------------------------------------------
!
        DO J=JDS,JDE
          TPH=SB+(J-JDS)*DPH
          TPV=TPH+DPH*0.5
          DXH(J)=A*DLM*COS(TPH)
          DXV(J)=A*DLM*COS(TPV)
          RDXH(J)=1./DXH(J)
          RDXV(J)=1./DXV(J)
        ENDDO
!
!-----------------------------------------------------------------------
      ENDIF global_regional_dx
!-----------------------------------------------------------------------
!
      DO J=JDS,JDE
        int_state%DXH(J)=DXH(J)
        int_state%DXV(J)=DXV(J)
        int_state%RDXH(J)=RDXH(J)
        int_state%RDXV(J)=RDXV(J)
      ENDDO
!
      DEALLOCATE(DXH)
      DEALLOCATE(DXV)
      DEALLOCATE(RDXH)
      DEALLOCATE(RDXV)
!
!-----------------------------------------------------------------------
!***  CHOOSE A J INDEX FOR AN "AVERAGE" DX.
!***  SELECT THE J THAT DIVIDES THE DOMAINS AREA IN HALF.
!-----------------------------------------------------------------------
!
      THETA_HALF=ASIN(0.5*SIN(-SB))
      JC=NINT(0.5*(JDE-JDS+1)+THETA_HALF/DPH)
!
!-----------------------------------------------------------------------
!***  SET TIME VARIABLES NEEDED FOR HISTORY OUTPUT.
!-----------------------------------------------------------------------
!
      NSTEPS_PER_HOUR=3600./int_state%DT
      int_state%NPREC=NSTEPS_PER_HOUR*int_state%NHRS_PREC
      int_state%NCLOD=NSTEPS_PER_HOUR*int_state%NHRS_CLOD
      int_state%NHEAT=NSTEPS_PER_HOUR*int_state%NHRS_HEAT
      int_state%NRDLW=NSTEPS_PER_HOUR*int_state%NHRS_RDLW
      int_state%NRDSW=NSTEPS_PER_HOUR*int_state%NHRS_RDSW
      int_state%NSRFC=NSTEPS_PER_HOUR*int_state%NHRS_SRFC
!
!-----------------------------------------------------------------------
!***  FINALLY INITIALIZE INDIVIDUAL SCHEMES.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE GFS PHYSICS SUITE IS CONSIDERED A SINGLE PACKAGE HERE.
!----------------------------------------------------------------------
!
      package: IF(GFS)THEN
!!!!!!  CALL GFS_PHYSICS_INIT
!
!----------------------------------------------------------------------
!***  IF NOT SELECTING THE GFS SUITE, EACH OF THE PHYSICS GROUPS IS
!***  TREATED INDIVIDUALLY.
!----------------------------------------------------------------------
!
      ELSE
!
!----------------------------------------------------------------------
!***  LONGWAVE RADIATION
!----------------------------------------------------------------------
!
        SELECT CASE (longwave)  
          CASE ('gfdl')
!
!***  WE ARE CALLING A WRF ROUTINE THUS FLIP THE VERTICAL.
!
            DO K=1,LM
              KFLIP=LM+1-K
              SFULL_FLIP(KFLIP)=SFULL(K+1)
              SMID_FLIP(KFLIP)=SMID(K)
            ENDDO
            SFULL_FLIP(LM+1)=SFULL(1)
!
            GMT=REAL(IHOUR_FCST)
            CALL GFDL_INIT(EMISS,SFULL_FLIP,SMID_FLIP,PT_CB            &
                          ,JULYR,START_MONTH,START_DAY,GMT             &
                          ,CO2TF                                       &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM)
          CASE ('rrtm')
!!!         CALL RRTMINIT(RTHRATEN,RTHRATENLW                          &
!!!                      ,int_state%CLDFRA,RESTART                     &
            CALL RRTMINIT(int_state%RESTART                            &
                         ,ALLOWED_TO_READ                              &
                         ,IDS,IDE,JDS,JDE,1,LM+1                       &
                         ,IMS,IME,JMS,JME,1,LM+1                       &
                         ,ITS,ITE,JTS,JTE,1,LM)
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF LONGWAVE SCHEME: INIT '
        END SELECT
!
!----------------------------------------------------------------------
!***  SHORTWAVE RADIATION
!----------------------------------------------------------------------
!
        SELECT CASE (shortwave)
          CASE ('gfdl')
!           WRITE(0,*)' Already called GFDL_INIT from LONGWAVE'
!!!       CASE ('gsfc')
!!!         CALL GSFC_INIT
          CASE ('dudh')
            CALL SWINIT(SWRAD_SCAT,int_state%RESTART                   &
                       ,ALLOWED_TO_READ                                &
                       ,IDS,IDE,JDS,JDE,1,LM+1                         &
                       ,IMS,IME,JMS,JME,1,LM+1                         &
                       ,ITS,ITE,JTS,JTE,1,LM)
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF SHORTWAVE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  SURFACE LAYER
!----------------------------------------------------------------------
!
        ALLOCATE(LOWLYR(IMS:IME,JMS:JME),STAT=I)
!
        SELECT CASE (sfc_layer)
          CASE ('myj')
            CALL MYJSFC_INIT(LOWLYR                                    &  !<-- Placeholder (computed in TURBULENCE)
                            ,int_state%USTAR,int_state%Z0              &
                            ,int_state%SM,int_state%SICE               &
                            ,int_state%IVGTYP,int_state%RESTART        &            
                            ,ALLOWED_TO_READ                           &
                            ,IDS,IDE,JDS,JDE,1,LM+1                    &
                            ,IMS,IME,JMS,JME,1,LM+1                    &
                            ,ITS,ITE,JTS,JTE,1,LM)   
!!!       CASE ('mm5')
!!!         CALL SFCLYR_INIT
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF SURFACE LAYER SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  TURBULENCE
!----------------------------------------------------------------------
!
        SELECT CASE (turbulence)
          CASE ('myj')
            CALL MYJPBL_INIT(int_state%EXCH_H,int_state%RESTART        &
                            ,IDS,IDE,JDS,JDE,LM                        &
                            ,IMS,IME,JMS,JME                           &
                            ,ITS,ITE,JTS,JTE)
!!!       CASE ('ysu')
!!!         CALL YSU_INIT
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF TURBULENCE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  LAND SURFACE
!----------------------------------------------------------------------
!
        SELECT CASE (land_surface)
          CASE ('noah')

          CALL NOAH_LSM_INIT(int_state%CMC,     int_state%ISLTYP       &
                            ,int_state%STC,     int_state%SMC          &
                            ,int_state%SH2O,    NUM_SOIL_LAYERS        &
                            ,int_state%RESTART, ALLOWED_TO_READ        &
                            ,IDS,IDE, JDS,JDE                          &
                            ,IMS,IME, JMS,JME                          &
                            ,ITS,ITE, JTS,JTE                         )

          CASE ('liss')

!!!         CALL LSM_INIT

          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF LAND SURFACE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!****  CONVECTION
!----------------------------------------------------------------------
!
        SELECT CASE (convection)
          CASE ('bmj')
            CALL BMJ_INIT(int_state%CLDEFI,int_state%RESTART           &
                         ,IDS,IDE,JDS,JDE,1,LM+1                       &
                         ,IMS,IME,JMS,JME,1,LM+1                       &
                         ,ITS,ITE,JTS,JTE,1,LM)
                    
!!!       CASE('kf')
!!!         CALL KF_INIT
!!!       CASE ('sas')
!!!         CALL SAS_INIT
!!!       CASE ('gd')
!!!         CALL GD_INIT
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF CONVECTION SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  MICROPHYSICS
!----------------------------------------------------------------------
!
        SELECT CASE (microphysics)
!
          CASE ('fer')
            DT_MICRO=int_state%NPRECIP*DT
            DELX=-2.*int_state%WBD*111.3/REAL(int_state%IM) !DX at rotated equator (km)
            DELY=-2.*int_state%SBD*111.3/REAL(int_state%JM) !DY at rotated equator (km)
!
            CALL FERRIER_INIT(DT_MICRO,DT,DELX,DELY,int_state%RESTART  &
                             ,int_state%F_ICE                          &
                             ,int_state%F_RAIN                         &
                             ,int_state%F_RIMEF                        &
                             ,int_state%MP_RESTART_STATE               &
                             ,int_state%TBPVS_STATE                    &
                             ,int_state%TBPVS0_STATE                   &
                             ,int_state%AVRAIN                         &
                             ,ALLOWED_TO_READ                          &
                             ,IDS,IDE,JDS,JDE,1,LM+1                   &
                             ,IMS,IME,JMS,JME,1,LM+1                   &
                             ,ITS,ITE,JTS,JTE,1,LM)
!
          CASE ('wsm3')
            CALL WSM3INIT(RHOAIR0,RHOWATER,RHOSNOW,CLIQ,CV             &
                         ,ALLOWED_TO_READ )
!!!       CASE ('kes')
!!!         CALL KESSLER_INIT
!!!       CASE ('tho')
!!!         CALL THOMPSON_INIT
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF MICROPHYSICS SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!****  GRAVITY WAVE DRAG (GWD) & MOUNTAIN BLOCKING (MB) INIT
!----------------------------------------------------------------------
!
        DTPHS=int_state%DT*int_state%NPHS
!
        CALL GWD_init(DTPHS,int_state%GLOBAL,int_state%RESTART          &
                      ,int_state%TPH0D,int_state%TLM0D                  &
                      ,int_state%GLAT,int_state%GLON                    &
                      ,int_state%CROT,int_state%SROT,int_state%HANGL    &
                      ,IDS,IDE,JDS,JDE,1,LM                             &
                      ,IMS,IME,JMS,JME,1,LM                             &
                      ,ITS,ITE,JTS,JTE,1,LM )
!
! uncomment this for output in future
!
!       IF(.NOT.int_state%RESTART)THEN
!         DO J=JMS,JME
!         DO I=IMS,IME
!           UGWDsfc(I,J)=0.
!           VGWDsfc(I,J)=0.
!         ENDDO
!         ENDDO
!       ENDIF
!
!
!----------------------------------------------------------------------
!
      ENDIF package
!
!----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_INITIALIZE
!
!-----------------------------------------------------------------------
      SUBROUTINE EXIT(NAME,T,Q,U,V,Q2,NTSD)
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      INCLUDE "../../inc/mpif.h"
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: NTSD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: T,Q,U,V,Q2
      CHARACTER(*),INTENT(IN) :: NAME
!
      INTEGER :: I,J,K,IEND,IERR,IRET
!----------------------------------------------------------------------
      IRET=0
  100 FORMAT(' EXIT ',A,' AT NTSD=',I5)
      IEND=ITE
!
      DO J=JTS,JTE
      IEND=ITE
      DO K=1,LM
!
      DO I=ITS,IEND
        IF(T(I,J,K)>330..OR.T(I,J,K)<180..OR.T(I,J,K)/=T(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,200)I,J,K,T(I,J,K),MYPE,NTSD
  200     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' T=',E12.5      &
      ,          ' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
  205     FORMAT(' EXIT ',A,' TEMPERATURE=',E12.5                      &
      ,          ' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ELSEIF(Q(I,J,K)<-1.E-4.OR.Q(I,J,K)>30.E-3                      &
               .OR.Q(I,J,K)/=Q(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,300)I,J,K,Q(I,J,K),MYPE,NTSD
  300     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' Q=',E12.5      &
      ,          ' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
  305     FORMAT(' EXIT ',A,' SPEC HUMIDITY=',E12.5                    &
      ,          ' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      IEND=ITE
      DO K=1,LM
      DO I=ITS,IEND
        IF(ABS(U(I,J,K))>125..OR.ABS(V(I,J,K))>125.                    &
     &         .OR.U(I,J,K)/=U(I,J,K).OR.V(I,J,K)/=V(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,400)I,J,K,U(I,J,K),V(I,J,K),MYPE,NTSD
  400     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' U=',E12.5      &
     &,          ' V=',E12.5,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
  405     FORMAT(' EXIT ',A,' U=',E12.5,' V=',E12.5                    &
     &,          ' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!----------------------------------------------------------------------
      END SUBROUTINE EXIT
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
