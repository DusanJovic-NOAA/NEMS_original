!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PHYSICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS PHYSICS INITIALIZE, ETC.) 
!***  IN MODULE_ATM_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-28  Vasic    - Removed counters (computed in SET_INTERNAL_STATE_PHY)
!   2008-10-    Vasic    - Restart capability
!   2008-08-23  Janjic   - Removed uz0h, vz0h
!   2008-08-23  Janjic   - General hybrid coordinate
!   2009-07-01  Vasic    - Added GFS physics package
!   2009-08-10  Black    - Merge with nest code
!   2009-11-03  W. Wang  - Remove WRF driver and flipping (for Ferrier)
!   2009-11-24  Sarah Lu - RAIN and RAINC added to GBPHYS arguments
!   2009-12-10  Sarah Lu - GRRAD output instant cloud cover
!   2009-12-11  Sarah Lu - GRRAD calling argument modified: remove ldiag3d;
!                           reverse flxur/cldcov sequence 
!   2009-12-15  Sarah Lu - GBPHYS calling argument modified: add dqdt
!   2009-11     Jovic    - Modified for ownership/import/export specification
!   2010-02-08  W. Wang  - Add wsm6 microphysics 
!   2010-03-10  W. Wang  - ADD species advection option
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_VARS_STATE
      USE MODULE_PHYSICS_INTERNAL_STATE
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &
                                   ,MPI_COMM_COMP                       &
                                   ,MYPE_SHARE                          &
                                   ,DSTRB,IDSTRB
!
      USE MODULE_CONTROL,ONLY : CAPPA,TIMEF
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_FLTBNDS,ONLY : POLEHN,POLEWN,SWAPHN,SWAPWN
!
      USE MODULE_RADIATION    ,ONLY : RADIATION
      USE MODULE_RA_GFDL      ,ONLY : GFDL_INIT,RDTEMP,TIME_MEASURE
      USE MODULE_RA_RRTM      ,ONLY : RRTM_INIT
      USE MODULE_TURBULENCE   ,ONLY : TURBL
      USE MODULE_SF_JSFC      ,ONLY : JSFC_INIT
      USE MODULE_BL_MYJPBL    ,ONLY : MYJPBL_INIT
      USE MODULE_LS_NOAHLSM   ,ONLY : DZSOIL,NOAH_LSM_INIT              &
                                     ,NUM_SOIL_LAYERS,SLDPTH
      USE MODULE_CU_BMJ       ,ONLY : BMJ_INIT
      USE MODULE_CONVECTION   ,ONLY : CUCNVC
!rv
      USE MODULE_CU_BMJ_DEV       ,ONLY : BMJ_INIT_DEV
      USE MODULE_CONVECTION_DEV   ,ONLY : CUCNVC_DEV
!rv
!      USE MODULE_MICROPHYSICS_NMM ,ONLY : FERRIER_INIT,GSMDRIVE         &
!                                         ,WSM3INIT,MICRO_RESTART
      USE MODULE_MICROPHYSICS_NMM ,ONLY : GSMDRIVE                      &
                                         ,MICRO_RESTART
      USE MODULE_MP_ETANEW, ONLY : FERRIER_INIT
      USE MODULE_MP_WSM6,   ONLY : WSM6INIT

      USE MODULE_H_TO_V       ,ONLY : H_TO_V,H_TO_V_TEND
      USE MODULE_GWD          ,ONLY : GWD_INIT
      USE MODULE_PRECIP_ADJUST
!
      USE MODULE_EXCHANGE
      USE MODULE_DIAGNOSE,ONLY: TWR,VWR,EXIT,EXIT_PHY,WRT_2D
!
      USE MODULE_PHYSICS_OUTPUT,ONLY: POINT_PHYSICS_OUTPUT
!
      USE MODULE_CLOCKTIMES,ONLY : cucnvc_tim,exch_phy_tim              &
                                  ,gsmdrive_tim,h_to_v_tim              &
                                  ,phy_init_tim,phy_run_tim,phy_sum_tim &
                                  ,pole_swap_phy_tim                    &
                                  ,radiation_tim,rdtemp_tim             &
                                  ,turbl_tim,update_phy_int_state_tim   &
                                  ,adjppt_tim,gfs_phy_tim
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE MODULE_PHYSICS_INIT_READ_BIN
      USE MODULE_PHYSICS_INIT_READ_NEMSIO
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
      INTEGER(KIND=KINT) :: MY_DOMAIN_ID,MYPE,NUM_PES
      INTEGER(KIND=KINT) :: START_YEAR,START_MONTH,START_DAY,START_HOUR &
                           ,START_MINUTE,START_SECOND
!
      INTEGER(KIND=KINT),SAVE :: JC,NSTEPS_PER_HOUR
!
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),SAVE :: DT,PT
!
      REAL(KIND=KDBL) :: btim,btim0
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                   !<-- The Physics internal state pointer
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
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Physics Gridded Component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Phy Register
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_REG=ESMF_SUCCESS
                                                                                                                                              
!-----------------------------------------------------------------------
!***  Register the Physics Initialize subroutine.  Since it is just one
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN,
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Initialize phase 1"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,PHY_INITIALIZE_1                  &  !<-- User's subroutine name
                                     ,1                                 &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Initialize phase 2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,PHY_INITIALIZE_2                  &  !<-- User's subroutine name
                                     ,2                                 &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Physics Run subroutine.
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
!***  Register the Physics Finalize subroutine.
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
!***  Check the error signal variable.
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
      SUBROUTINE PHY_INITIALIZE_1(GRID_COMP                             &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_ATM                             &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  CARRY OUT ALL NECESSARY SETUPS FOR THE MODEL PHYSICS.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP			   !<-- The Physics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE			   !<-- The Physics Initialize step's import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE			   !<-- The Physics Initialize step's export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM			   !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT) :: L,N,RC
!
      TYPE(WRAP_PHY_INT_STATE) :: WRAP                                     !<-- This wrap is a derived type which contains
                                                                           !    only a pointer to the internal state.  It is needed
                                                                           !    for using different architectures or compilers.
!
      TYPE(ESMF_Grid)         :: GRID
      TYPE(ESMF_VM)           :: VM                                        !<-- The virtual machine
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize the Physics timers.
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
      gfs_phy_tim=0.
!
!-----------------------------------------------------------------------
!***  Allocate the Physics internal state pointer.
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE,STAT=RC)
!
!-----------------------------------------------------------------------
!***  Attach the internal state to the Physics gridded component.
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

!-----------------------------------------------------------------------
!***  Insert the local domain starting limits and the halo width into
!***  the Physics internal state.
!-----------------------------------------------------------------------
!
!     IF(IHALO==JHALO)THEN
!       int_state%NHALO=IHALO
!     ELSE
!        RC_INIT=ESMF_FAILURE
!        WRITE(0,*)'Error due to ihalo /= jhalo'
!      ENDIF
!
      int_state%ITS=ITS
      int_state%ITE=ITE
      int_state%JTS=JTS
      int_state%JTE=JTE
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values fron the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_PHY_DIMS(GRID_COMP                                &
                              ,int_state%INPES,int_state%JNPES          &
                              ,LM                                       &
                              ,int_state%NUM_TRACERS_MET                &
                              ,int_state%NUM_TRACERS_CHEM               &
                              ,int_state%PCPHR                          &
                              ,int_state%GFS                            &
                              ,int_state%MICROPHYSICS                   &
                              ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM to obtain the task ID and total number of tasks
!***  for the internal state.
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
!***  int_state%NUM_PES taken from VM is the total number of tasks
!***  in the run including quilt tasks.  Actually we want just the
!***  number of forecast tasks.
!-----------------------------------------------------------------------
!
      int_state%NUM_PES=int_state%INPES*int_state%JNPES
!
      NUM_PES=int_state%NUM_PES                                            !<-- The number of forecast tasks
      MYPE=int_state%MYPE                                                  !<-- The local task ID
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are needed for the remaining
!***  initialization process.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<NUM_PES)THEN                                     !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  Allocate all necessary internal state variables.  Those that
!***  are owned/exported are pointed into allocated memory within
!***  the Physics' composite VARS array.
!-----------------------------------------------------------------------
!
        CALL SET_INTERNAL_STATE_PHY_1(INT_STATE,LM)
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
!-----------------------------------------------------------------------
!***  Put the allocated pointers of all export variables (they must be
!***  owned) into the Physics export state.
!-----------------------------------------------------------------------
!
        CALL PUT_VARS_IN_STATE(int_state%VARS,int_state%NUM_VARS,'X',GRID,EXP_STATE)
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
!-----------------------------------------------------------------------
!
      phy_init_tim=(timef()-btim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_INITIALIZE_1
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_INITIALIZE_2(GRID_COMP                             &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK                                 &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Point all unowned internal state variables at allocated memory,
!***  get time information, read input data.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Physics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Physics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Physics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK                           !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL,    INTENT(OUT)  :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT) :: L,N,RC                                      &
                           ,IDENOMINATOR_DT,INTEGER_DT,NUMERATOR_DT
!
      TYPE(ESMF_State)        :: IMP_STATE_WRITE
      TYPE(ESMF_Grid)         :: GRID
      TYPE(ESMF_TimeInterval) :: DT_ESMF                                   !<-- The timestep from the ATM Clock
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For unowned (unallocated) internal state variables, point their
!***  appropriate corresponding locations in the composite VARS array
!***  at allocated memory from the owned Dynamics variable seen in the
!***  Physics import state.
!***  For owned (allocated) variables that are also owned/exported by
!***  the Dynamics, transfer the data into them from the import state.
!-----------------------------------------------------------------------
!
      CALL GET_VARS_FROM_STATE(int_state%VARS, int_state%NUM_VARS, IMP_STATE)
!
!-----------------------------------------------------------------------
!***  Now re-point the internal state variables associated with the
!***  VARS array into that array.  Unowned variables will now point
!***  to allocated memory in the Physics.
!-----------------------------------------------------------------------
!
      CALL SET_INTERNAL_STATE_PHY_2(INT_STATE,LM)
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values fron the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_PHY(GRID_COMP,INT_STATE,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IM=int_state%IM
      JM=int_state%JM
!d      LM=int_state%LM
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are needed for the remaining
!***  initialization process.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(int_state%MYPE<int_state%NUM_PES)THEN                 !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  Assign the fundamental timestep retrieved from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Fundamental Timestep from ATM's Clock"
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
!
!-----------------------------------------------------------------------
!***  Retrieve the domain ID from the Physics import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Domain ID from Physics Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                           &  !<-- The Physics import state
                              ,name ='DOMAIN_ID'                         &  !<-- Name of variable to get from Physics import state
                              ,value=MY_DOMAIN_ID                        &  !<-- Put extracted value here
                              ,rc =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
                               ,int_state%SBD                           &
                               ,int_state%WBD                           &
                               ,int_state%DPHD                          &
                               ,int_state%DLMD                          &
                               ,int_state%TPH0D                         &
                               ,int_state%TLM0D                         &
                               ,MY_DOMAIN_ID                            &
                               ,IDS,IDE,JDS,JDE,LM                      &
                               ,IMS,IME,JMS,JME                         &
                               ,ITS,ITE,JTS,JTE)
!
!-----------------------------------------------------------------------
!***  Create the ESMF Fields for the import/export states.
!***  For now send ALLOC_FIELDS_PHY the entire internal state
!***  from which the desired variables will be extracted for
!***  insertion into the import/export states.
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
!
!-----------------------------------------------------------------------
!***  Also insert the value of NUM_TRACERS_TOTAL into the export state.
!***  This will tell the Dyn-Phy Coupler how many constituents
!***  there are to transfer in the 4-D Tracers Field.
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
!***  Also insert the index values of the 4-D Tracers array where
!***  Q and CW reside.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_Q into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_Q'                           &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_Q                   &  !<-- The location of Q in TRACERS
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_CW into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_CW'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_CW                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract all forecast tasks' horizontal subdomain limits
!***  from the Physics import state and give them to the 
!***  Physics internal state.
!***  This is necessary if quilting is selected because these
!***  limits will be taken from the Dynamics/Physics internal
!***  states, placed into the Write components' import states
!***  and used for the combining of local domain data onto the
!***  global domain.
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
!***  Retrieve the import state of the Write gridded component
!***  from the Physics export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Retrieve Write Import State from Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE                        &  !<-- The Physics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Physics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write component import state from Physics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL POINT_PHYSICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!
        if (mype==0)        CALL ESMF_StatePrint(EXP_STATE)
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks

      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PHY INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_INITIALIZE_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_RUN(GRID_COMP                                      &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK                                          &
                        ,RC_RUN )
!
!-----------------------------------------------------------------------
!***  THE INTEGRATION OF THE MODEL PHYSICS IS DONE
!***  THROUGH THIS ROUTINE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Only for GFS physics:
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,            ONLY: CP,R,RHOWATER,STBOLT,XLV
!
      USE NAMELIST_PHYSICS_DEF,        ONLY: FHSWR                      &
                                            ,IAER,IALB,ICO2,IEMS        &
                                            ,IOVR_LW,IOVR_SW,ISOL       &
                                            ,LDIAG3D,LSCCA              &
                                            ,LSLWR,LSM,LSSAV,LSSWR      &
                                            ,PRE_RAD,RAS,SASHAL
!
      USE LAYOUT1,                    ONLY : IPT_LATS_NODE_R            &
                                            ,LATS_NODE_R
!
      USE DATE_DEF,                   ONLY : FHOUR
      USE MODULE_RADIATION_DRIVER,    ONLY : GRRAD,RADINIT
      USE MODULE_RADIATION_ASTRONOMY, ONLY : ASTRONOMY
      USE MERSENNE_TWISTER
      USE RESOL_DEF,                  ONLY : LATR,LONR                  &
                                            ,NCLD,NFXR,NMTVR            &
                                            ,NTCW,NTOZ                  &
                                            ,NUM_P2D,NUM_P3D
!
      USE OZNE_DEF,                   ONLY : LEVOZP,PL_COEFF,PL_PRES
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Physics component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The Physics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                       !<-- The Physics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK                           !<-- The ATM Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!---------------------------------
!***  GFS physics local variables
!---------------------------------
!
      TYPE(ESMF_Time)                              :: CURRTIME
      LOGICAL,SAVE                                 :: FIRST=.true.
      LOGICAL                                      :: LPRNT=.false.
      LOGICAL                                      :: LSFWD,OPENED,FLIPV,CHANGE
      INTEGER,PARAMETER                            :: IFLIP=1,NTRAC=3            !!!!!! later ntrac read form namelist
      INTEGER                                      :: II, JJ, ICWP,IMJM, IDATE(4), JDAT(8)
      INTEGER                                      :: KFLIP,ISEED
      INTEGER ,SAVE                                :: ID,IDAY,IMON,MIDMON,MIDM,MIDP,K1OZ,K2OZ,SEED0
      INTEGER ,DIMENSION(13)                       :: DAYS
      INTEGER ,DIMENSION(JTS:JTE)                  :: LONSPERLAR, GLOBAL_LATS_R
!
      REAL (KIND=KDBL)                             :: T850,FACOZ,DTLW,DTSW,DTLWI,DTSWI,RTvR,CLSTP,DTP,DTF,SOLHR
      REAL (KIND=KDBL)                             :: XLVRW,XLVRWI,DTPHS,DTPHSI,RoCP,MINDT
      REAL (KIND=KDBL) ,DIMENSION(1)               :: RCS2_V,XLAT,FLGMIN_L,CV,CVB,CVT  ! (cv, cvb, cvt not in use when ntcw-1 > 0)
      REAL (KIND=KDBL) ,DIMENSION(1)               :: TSEA,TISFC,ZORL,SLMSK,SNWDPH,SHELEG,SNCOVR,SNOALB
      REAL (KIND=KDBL) ,DIMENSION(1)               :: XSIHFCS,XSICFCS,XSLPFCS,XTG3FCS,XVEGFCS,XVETFCS,XSOTFCS
      REAL (KIND=KDBL) ,DIMENSION(1)               :: ALVSF,ALNSF,ALVWF,ALNWF,FACSF,FACWF
      REAL (KIND=KDBL) ,DIMENSION(1)               :: WRK, DPSHC, GQ, RANNUM_V
      REAL (KIND=KDBL) ,DIMENSION(1)               :: ORO, EVAP, HFLX, CDQ, QSS
      REAL (KIND=KDBL) ,DIMENSION(LM)              :: CLDCOV_V,PRSL,PRSLK,GU,GV,GT,GR,VVEL,F_ICE,F_RAIN,R_RIME
      REAL (KIND=KDBL) ,DIMENSION(LM)              :: ADT,ADU,ADV,PHIL
      REAL (KIND=KDBL) ,DIMENSION(LM,NTRAC)        :: GR3,ADR
      REAL (KIND=KDBL) ,DIMENSION(LM+1)            :: PRSI,PRSIK,RSGM,PHII
      REAL (KIND=KDBL) ,DIMENSION(JTS:JTE)         :: SINLAT_R,COSLAT_R
      REAL (KIND=KDBL) ,DIMENSION(ITS:ITE,JTS:JTE) :: XLON,COSZEN,COSZDG,XKT2
      REAL (KIND=KDBL) ,DIMENSION((ITE-ITS+1)*(JTE-JTS+1)) :: RANNUM
!
      REAL (KIND=KDBL) ,DIMENSION(27)              :: FLUXR_V
      REAL (KIND=KDBL) ,DIMENSION(1,LM,NTRAC-1)    :: GR1
!
      REAL (KIND=KDBL) ,DIMENSION(1)               :: SFCNSW, SFCDSW, SFALB, SFCDLW, TSFLW
      REAL (KIND=KDBL) ,DIMENSION(LM)              :: SWH, HLW
!--- gbphys ---
      LOGICAL                                      :: OLD_MONIN, CNVGWD, NEWSAS
      INTEGER ,DIMENSION(2)                        :: NCW
      REAL (KIND=KDBL)                             :: CCWF
      REAL (KIND=KDBL) ,DIMENSION(1)               :: BENGSH, GESHEM, TPRCP, SRFLAG, SHDMIN, SHDMAX, CANOPY
      REAL (KIND=KDBL) ,DIMENSION(1)               :: RAIN, RAINC
      REAL (KIND=KDBL) ,DIMENSION(1)               :: ACV, ACVB, ACVT
      REAL (KIND=KDBL) ,DIMENSION(2)               :: FLGMIN
      REAL (KIND=KDBL) ,DIMENSION(3)               :: CRTRH
      REAL (KIND=KDBL) ,DIMENSION(NUM_SOIL_LAYERS) :: SMC_V, STC_V, SLC_V
      REAL (KIND=KDBL) ,DIMENSION(14)              :: HPRIME
      REAL (KIND=KDBL) ,DIMENSION(LM)              :: UPD_MF, DWN_MF, DET_MF   !!!!!!!!!!! not in use
      REAL (KIND=KDBL) ,DIMENSION(LM)              :: DQDT                     !!!!!!!!!!! not in use
      REAL (KIND=KDBL) ,DIMENSION(LM,9)            :: DQ3DT                    !!!!!!!!!!!  (9=5+pl_coeff)
      REAL (KIND=KDBL) ,DIMENSION(LM,6)            :: DT3DT                    !!!!!!!!!!! while
      REAL (KIND=KDBL) ,DIMENSION(LM,4)            :: DU3DT, DV3DT             !!!!!!!!!!! LDIAG3D =.FALSE.

      REAL (KIND=KDBL) ,DIMENSION(:,:)  ,ALLOCATABLE :: OZPLOUT_V
      REAL (KIND=KDBL) ,DIMENSION(:,:,:),ALLOCATABLE :: OZPLOUT

      REAL (KIND=KDBL) ,DIMENSION(3)               :: PHY_F2DV   ! NUM_P2D for Zhao =3, Ferr=1 (fix later)
      REAL (KIND=KDBL) ,DIMENSION(LM,4)            :: PHY_F3DV   ! NUM_P3D for Zhao =4, Ferr=3 (fix later)
!--- gbphys output
      REAL (KIND=KDBL) ,DIMENSION(1)               :: EVBSA, EVCWA, TRANSA, SBSNOA, SNOWCA, CLDWRK, PSMEAN
      REAL (KIND=KDBL) ,DIMENSION(1)               :: CHH, CMM, EP, EPI, DLWSFCI, ULWSFCI, USWSFCI, DSWSFCI
      REAL (KIND=KDBL) ,DIMENSION(1)               :: DLWSFC, ULWSFC, DTSFC, DQSFC, DUSFC, DVSFC, GFLUX
      REAL (KIND=KDBL) ,DIMENSION(1)               :: DTSFCI, DQSFCI, GFLUXI, T1, Q1, U1, V1
      REAL (KIND=KDBL) ,DIMENSION(1)               :: ZLVL, SOILM, RUNOFF, SRUNOFF
      REAL (KIND=KDBL) ,DIMENSION(1)               :: F10M, UUSTAR, FFMM, FFHH
      REAL (KIND=KDBL) ,DIMENSION(1)               :: PSURF, U10M, V10M, T2M, Q2M, HPBL, PWAT
!
!---------------------------
!***  Other local variables
!---------------------------
!
      INTEGER(KIND=KINT) :: I,J,IRTN,ISTAT,JULDAY,JULYR,L               &
                           ,N,NPRECIP,NTIMESTEP,RC,NTIMESTEP_RAD,IMICRO,&
                            SPECADV
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      REAL(KIND=KFPT) :: JULIAN,PDTOP,SECONDS_TOTAL,XTIME
!
      REAL(KIND=KFPT),DIMENSION(LM) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL(KIND=KFPT),DIMENSION(LM+1) :: PSG1,SG2
!
      LOGICAL(KIND=KLOG) :: CALL_LONGWAVE                               &
                           ,CALL_SHORTWAVE                              &
                           ,CALL_TURBULENCE                             &
                           ,CALL_PRECIP                                 &
                           ,CALL_GFS_PHY
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      DATA DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS 
      MYPE=MYPE_SHARE
!
!-----------------------------------------------------------------------
!***  Retrieve the timestep from the Clock.
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
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF
      int_state%NTSD=NTIMESTEP
!
!-----------------------------------------------------------------------
!***  Call radiation so that updated fields are written to the
!***  history files after 0 hours.
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP==0)THEN
         NTIMESTEP_RAD=NTIMESTEP
      ELSE
         NTIMESTEP_RAD=NTIMESTEP+1
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Dereference some internal state components for convenience.
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
!***  Update the Physics internal state with data from
!***  the import state.  This must be done every time step
!***  since the temperature is updated every timestep.
!-----------------------------------------------------------------------
!
!d      btim=timef()
      CALL GET_VARS_FROM_STATE(int_state%VARS, int_state%NUM_VARS, IMP_STATE)
!d      CALL UPDATE_INTERNAL_STATE_PHY(IMP_STATE,INT_STATE)
!d      update_phy_int_state_tim=update_phy_int_state_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update the max/min values of the Temperature in the lowest
!***  model layer.  Reset those values at the start of each hour.
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP*int_state%DT,3600.)==0)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%TLMAX(I,J)=-999.
          int_state%TLMIN(I,J)=999.
        ENDDO
        ENDDO
      ENDIF
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        int_state%TLMAX(I,J)=MAX(int_state%TLMAX(I,J),int_state%T(I,J,LM))  !<--- Hourly max lowest layer T
        int_state%TLMIN(I,J)=MIN(int_state%TLMIN(I,J),int_state%T(I,J,LM))  !<--- Hourly min lowest layer T
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      gfs_phys_test: IF(.NOT.int_state%GFS)THEN                            !<-- NMM-B physics is NOT the GFS package
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Set logical switches for calling each of the Physics schemes.
!-----------------------------------------------------------------------
!
        CALL_SHORTWAVE=MOD(NTIMESTEP_RAD,int_state%NRADS)==0
        CALL_LONGWAVE=MOD(NTIMESTEP_RAD,int_state%NRADL)==0
        CALL_TURBULENCE=MOD(NTIMESTEP,int_state%NPHS)==0
        CALL_PRECIP=MOD(NTIMESTEP,NPRECIP)==0
!
!-----------------------------------------------------------------------
!***  Update WATER array from CWM, F_ICE, F_RAIN for Ferrier 
!***  microphysics but only if any of the Physics subroutines 
!***  are called (subroutine UPDATE_WATER is after subroutine
!***  PHYSICS_INITIALIZE in this module).
!
!***  Expanded to also update CWM, F_ICE, F_RAIN, F_RIMEF for non-Ferrier
!     microphysics
!-----------------------------------------------------------------------
!
           SPECADV = 0
!           IF(int_state%SPEC_ADV) SPECADV=1
           IF(int_state%SPEC_ADV .OR. int_state%MICROPHYSICS=='wsm6') SPECADV=1
!        update_wtr: IF (int_state%MICROPHYSICS=='fer' .AND.             &
        update_wtr: IF ( (int_state%MICROPHYSICS=='fer' .OR.            &
                          int_state%MICROPHYSICS=='wsm6') .AND.         &!for wsm6, 5/28/2010,
                                                                         !maybe remove int_state%microphys.
                          (CALL_SHORTWAVE .OR. CALL_LONGWAVE .OR.       &
                           CALL_TURBULENCE .OR. CALL_PRECIP) ) THEN
           SELECT CASE (trim(int_state%MICROPHYSICS))
           CASE ('wsm6')
              IMICRO=1
           CASE ('wsm3')
              IMICRO=2
           CASE DEFAULT
              IMICRO=0
           END SELECT
           CALL UPDATE_WATER(int_state%CW                               &
                            ,int_state%F_ICE                            &
                            ,int_state%F_RAIN                           &
                            ,int_state%F_RIMEF                          &
                            ,int_state%NUM_WATER                        &
                            ,int_state%WATER                            &
                            ,int_state%T                                &
                            ,int_state%P_QC                             &
                            ,int_state%P_QR                             &
                            ,int_state%P_QS                             &
                            ,int_state%P_QI                             &
                            ,int_state%P_QG                             &
                            ,IMICRO                                     &
                            ,SPECADV                                    &
                            ,NTIMESTEP                                  &
                            ,IDS,IDE,JDS,JDE,LM                         &
                            ,IMS,IME,JMS,JME                            &
                            ,ITS,ITE,JTS,JTE)
        ENDIF update_wtr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CALL THE INDIVIDUAL PHYSICAL PROCESSES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
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
!
        IF(int_state%NTSD==0)THEN
          IF(int_state%PCPFLG)THEN
            CALL READPCP(MYPE                                           &
                        ,int_state%PPTDAT                               &
                        ,int_state%DDATA                                &
                        ,int_state%LSPA                                 &
                        ,int_state%PCPHR                                &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE)
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Call the individual physical processes.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Radiation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  Radiation needs some specific time quantities.
!
        CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR   &
                         ,START_MINUTE,START_SECOND                     &
                         ,NTIMESTEP,int_state%DT                        &
                         ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
        radiatn: IF(CALL_SHORTWAVE.OR.CALL_LONGWAVE)THEN
!-----------------------------------------------------------------------
!
          btim=timef()
!
!-----------------------------------------------------------------------
!***  Empty the ACFRST and ACFRCV accumulation arrays if it is time
!***  to do so prior to their being updated by the radiation.
!-----------------------------------------------------------------------
!
          IF(MOD(NTIMESTEP,int_state%NCLOD)==0)THEN
            DO J=JTS,JTE
            DO I=ITS,ITE
              int_state%ACFRST(I,J)=0.
              int_state%ACFRCV(I,J)=0.
              int_state%NCFRST(I,J)=0
              int_state%NCFRCV(I,J)=0
            ENDDO
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!***  Temporary switch between radiation schemes placed in PHY_RUN
!***  rather than inside RADIATION_DRIVER (will be done later)
!-----------------------------------------------------------------------
!
          CALL ESMF_ClockGet(clock       =CLOCK                         &  !<-- The ESMF Clock
                            ,currTime    =CURRTIME                      &  !<-- The current time (ESMF) on the clock
                            ,rc          =RC)
!
          CALL ESMF_TimeGet(time=CURRTIME                               &  !<-- The cuurent forecast time (ESMF)
                           ,yy  =JDAT(1)                                &  !<-- The current forecast year (integer)
                           ,mm  =JDAT(2)                                &  !<-- The current forecast month (integer)
                           ,dd  =JDAT(3)                                &  !<-- The current forecast day (integer)
                           ,h   =JDAT(5)                                &  !<-- The current forecast hour (integer)
                           ,m   =JDAT(6)                                &  !<-- The current forecast minute (integer)
                           ,s   =JDAT(7)                                &  !<-- The current forecast second (integer)
                           ,rc  =RC)
          JDAT(4)=0
          JDAT(8)=0
!
          CALL RADIATION(NTIMESTEP_RAD                                  &
                        ,int_state%DT,JULDAY,JULYR,XTIME,JULIAN         &
                        ,START_HOUR,int_state%NPHS                      &
                        ,int_state%GLAT,int_state%GLON                  &
                        ,int_state%NRADS,int_state%NRADL                &
                        ,DSG2,SGML2,PDSG1,PSGML1                        &
                        ,int_state%PT,int_state%PD                      &
                        ,int_state%T,int_state%Q                        &
                        ,int_state%THS,int_state%ALBEDO                 &
                        ,int_state%P_QV,int_state%P_QC,int_state%P_QR   &
                        ,int_state%P_QI,int_state%P_QS,int_state%P_QG   &
                        ,int_state%F_QV,int_state%F_QC,int_state%F_QR   &
                        ,int_state%F_QI,int_state%F_QS,int_state%F_QG   &
                        ,int_state%SM,int_state%CLDFRA                  &
                        ,int_state%NUM_WATER,int_state%WATER            &
                        ,int_state%RLWTT,int_state%RSWTT                &
                        ,int_state%RLWIN,int_state%RSWIN                &
                        ,int_state%RSWINC,int_state%RSWOUT              &
                        ,int_state%RLWTOA,int_state%RSWTOA              &
                        ,int_state%CZMEAN,int_state%SIGT4               &
                        ,int_state%CFRACL,int_state%CFRACM              &
                        ,int_state%CFRACH                               &
                        ,int_state%ACFRST,int_state%NCFRST              &
                        ,int_state%ACFRCV,int_state%NCFRCV              &
                        ,int_state%CUPPT,int_state%SNO                  &
                        ,int_state%HTOP,int_state%HBOT                  &
                        ,int_state%SHORTWAVE,int_state%LONGWAVE         &
!---- RRTM part ---------------------------------------------------------
                        ,int_state%DT_INT,JDAT                          &
                        ,int_state%CW,int_state%O3                      &
                        ,int_state%F_ICE,int_state%F_RAIN               &
                        ,int_state%F_RIMEF                              &
                        ,int_state%SI,int_state%TSKIN                   &
                        ,int_state%Z0,int_state%SICE                    &
                        ,int_state%MXSNAL,int_state%SGM                 &
                        ,int_state%STDH,int_state%OMGALF                &
!------------------------------------------------------------------------
                        ,LM)
!
          radiation_tim=radiation_tim+(timef()-btim)
!
        ENDIF radiatn
!
!-----------------------------------------------------------------------
!***  Update the temperature with the radiative tendency.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL RDTEMP(NTIMESTEP,int_state%DT,JULDAY,JULYR,START_HOUR      &
                   ,int_state%GLAT,int_state%GLON                       &
                   ,int_state%CZEN,int_state%CZMEAN,int_state%T         &
                   ,int_state%RSWTT,int_state%RLWTT                     &
                   ,LM)
!
        rdtemp_tim=rdtemp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL SWAPHN(int_state%RSWIN,IMS,IME,JMS,JME,1,int_state%INPES)
          CALL POLEHN(int_state%RSWIN,IMS,IME,JMS,JME,1                  &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                     &
                     ,int_state%INPES,int_state%JNPES)
!
          pole_swap_phy_tim=pole_swap_phy_tim+(timef()-btim)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Empty the accumulators of sfc energy flux and sfc hydrology if
!***  it is time to do so prior to their being updated by turbulence.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NRDLW)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ALWIN(I,J) =0.
            int_state%ALWOUT(I,J)=0.
            int_state%ALWTOA(I,J)=0.
            int_state%ARDLW(I,J) =0.                                       !<-- An artificial 2-D array (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NRDSW)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ASWIN(I,J)=0.
            int_state%ASWOUT(I,J)=0.
            int_state%ASWTOA(I,J)=0.
            int_state%ARDSW(I,J) =0.                                       !<-- An artificial 2-D array (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NSRFC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%SFCSHX(I,J)=0.
            int_state%SFCLHX(I,J)=0.
            int_state%SUBSHX(I,J)=0.
            int_state%SNOPCX(I,J)=0.
            int_state%POTFLX(I,J)=0.
            int_state%ASRFC(I,J) =0.                                       !<-- An artificial 2-D array (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACSNOW(I,J)=0.
              int_state%ACSNOM(I,J)=0.
            int_state%SSROFF(I,J)=0.
            int_state%BGROFF(I,J)=0.
            int_state%SFCEVP(I,J)=0.
            int_state%POTEVP(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
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
          CALL TURBL(NTIMESTEP,int_state%DT,int_state%NPHS              &
                    ,int_state%NUM_WATER,NUM_SOIL_LAYERS,SLDPTH,DZSOIL  &
                    ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1,PT                &
                    ,int_state%SM,int_state%CZEN,int_state%CZMEAN       &
                    ,int_state%SIGT4,int_state%RLWIN,int_state%RSWIN    &
                    ,int_state%RADOT                                    &
                    ,int_state%PD,int_state%T                           &
                    ,int_state%Q,int_state%CW                           &
                    ,int_state%F_ICE,int_state%F_RAIN,int_state%SR      &
                    ,int_state%Q2,int_state%U,int_state%V               &
                    ,int_state%DUDT,int_state%DVDT                      &
                    ,int_state%THS,int_state%TSKIN,int_state%SST        &
                    ,int_state%PREC,int_state%SNO                       &
                    ,int_state%WATER                                    &
                    ,int_state%P_QV,int_state%P_QC,int_state%P_QR       &
                    ,int_state%P_QI,int_state%P_QS,int_state%P_QG       &
                    ,int_state%F_QV,int_state%F_QC,int_state%F_QR       &
                    ,int_state%F_QI,int_state%F_QS,int_state%F_QG       &
                    ,int_state%FIS,int_state%Z0,int_state%Z0BASE        &
                    ,int_state%USTAR,int_state%PBLH,int_state%LPBL      &
                    ,int_state%XLEN_MIX,int_state%RMOL                  &
                    ,int_state%EXCH_H,int_state%AKHS,int_state%AKMS     &
                    ,int_state%AKHS_OUT,int_state%AKMS_OUT              &
                    ,int_state%THZ0,int_state%QZ0                       &
                    ,int_state%UZ0,int_state%VZ0                        &
                    ,int_state%QSH,int_state%MAVAIL                     &
                    ,int_state%STC,int_state%SMC,int_state%CMC          &
                    ,int_state%SMSTAV,int_state%SMSTOT                  &
                    ,int_state%SSROFF,int_state%BGROFF                  &
                    ,int_state%IVGTYP,int_state%ISLTYP,int_state%VEGFRC &
                    ,int_state%SHDMIN,int_state%SHDMAX,int_state%GRNFLX &
                    ,int_state%SFCEXC,int_state%ACSNOW,int_state%ACSNOM &
                    ,int_state%SNOPCX,int_state%SICE                    &
                    ,int_state%TG,int_state%SOILTB                      &
                    ,int_state%ALBASE,int_state%MXSNAL,int_state%ALBEDO &
                    ,int_state%SH2O,int_state%SI,int_state%EPSR         &
                    ,int_state%U10,int_state%V10                        &
                    ,int_state%TH10,int_state%Q10                       &
                    ,int_state%TSHLTR,int_state%QSHLTR,int_state%PSHLTR &
                    ,int_state%PSFC,int_state%T2                        &
                    ,int_state%QSG,int_state%QVG,int_state%QCG          &
                    ,int_state%SOILT1,int_state%TSNAV                   &
                    ,int_state%TWBS,int_state%QWBS                      &
                    ,int_state%SFCSHX,int_state%SFCLHX,int_state%SFCEVP &
                    ,int_state%POTEVP,int_state%POTFLX,int_state%SUBSHX &
                    ,int_state%APHTIM                                   &
                    ,int_state%ARDSW,int_state%ARDLW                    &
                    ,int_state%ASRFC                                    &
                    ,int_state%CROT,int_state%SROT                      &
                    ,int_state%HSTDV,int_state%HCNVX,int_state%HASYW    &
                    ,int_state%HASYS,int_state%HASYSW,int_state%HASYNW  &
                    ,int_state%HLENW,int_state%HLENS,int_state%HLENSW   &
                    ,int_state%HLENNW,int_state%HANGL,int_state%HANIS   &
                    ,int_state%HSLOP,int_state%HZMAX                    &
                    ,int_state%RSWOUT,int_state%RSWTOA,int_state%RLWTOA &
                    ,int_state%ASWIN,int_state%ASWOUT,int_state%ASWTOA  &
                    ,int_state%ALWIN,int_state%ALWOUT,int_state%ALWTOA  &
                    ,int_state%RTHBLTEN,int_state%RQVBLTEN              &
                    ,int_state%GWDFLG,int_state%PCPFLG                  &
                    ,int_state%DDATA,int_state%UCMCALL                  &
                    ,int_state%TURBULENCE,int_state%SFC_LAYER           &
                    ,int_state%LAND_SURFACE                             &
                    ,int_state%MICROPHYSICS                             &
                    ,IDS,IDE,JDS,JDE,LM                                 &
                    ,IMS,IME,JMS,JME                                    &
                    ,ITS,ITE,JTS,JTE)
!
          turbl_tim=turbl_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Exchange wind tendencies.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,1,1)
!
          exch_phy_tim=exch_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Now interpolate wind tendencies from H to V points.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL H_TO_V_TEND(int_state%DUDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%U)
          CALL H_TO_V_TEND(int_state%DVDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%V)
!
          h_to_v_tim=h_to_v_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM      &
                       ,int_state%INPES,int_state%JNPES)
!
            pole_swap_phy_tim=pole_swap_phy_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange wind components and TKE.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                  &
                        ,2,2)
!
          CALL HALO_EXCH(int_state%UZ0,1,int_state%VZ0,1                &
                        ,int_state%Q2,LM                                &
                        ,1,1)
!
          exch_phy_tim=exch_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        ENDIF turbulence
!
!----------------------------------------------------------------------- 
!***  Empty the accumulators of precipitation and latent heating if is
!***  is time prior to their being updated by convection/microphysics.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACPREC(I,J)=0.
            int_state%CUPREC(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
!
      IF(MOD(NTIMESTEP,int_state%NHEAT)==0)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%AVCNVC(I,J)=0.   !- was a scalar, now 2D for ESMF
          int_state%AVRAIN(I,J)=0.   !- was a scalar, now 2D for ESMF
        ENDDO
        ENDDO
        DO L=1,LM
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%TRAIN(I,J,L)=0.
          int_state%TCUCN(I,J,L)=0.
        ENDDO
        ENDDO
        ENDDO
      ENDIF    !-- IF(MOD(NTSD_BUCKET,NHEAT)==0)THEN
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Convection
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        convection: IF(CALL_PRECIP.AND.int_state%CONVECTION/='none')THEN
!
          btim=timef()
!
!-----------------------------------------------------------------------
!***  Temporary switch between two convection schemes (bmj & bmj_dev)
!***  placed here in PHY_RUN
!-----------------------------------------------------------------------
          IF(int_state%CONVECTION=='bmj')THEN
!
          CALL CUCNVC(NTIMESTEP,int_state%DT,int_state%NPRECIP          &
                     ,int_state%NRADS,int_state%NRADL                   &
                     ,int_state%MINUTES_HISTORY                         &
                     ,int_state%DYH,int_state%RESTART,int_state%HYDRO   &
                     ,int_state%CLDEFI,int_state%NUM_WATER              &
                     ,int_state%F_ICE,int_state%F_RAIN                  &
                     ,int_state%P_QV,int_state%P_QC,int_state%P_QR      &
                     ,int_state%P_QI,int_state%P_QS,int_state%P_QG      &
                     ,int_state%F_QV,int_state%F_QC,int_state%F_QR      &
                     ,int_state%F_QI,int_state%F_QS,int_state%F_QG      &
                     ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1                  &
                     ,int_state%PT,int_state%PD                         &
                     ,int_state%T,int_state%Q                           &
                     ,int_state%CW,int_state%TCUCN,int_state%WATER      &
                     ,int_state%OMGALF                                  &
                     ,int_state%U,int_state%V                           &
                     ,int_state%FIS,int_state%W0AVG                     &
                     ,int_state%PREC,int_state%ACPREC,int_state%CUPREC  &
                     ,int_state%CUPPT,int_state%CPRATE                  &
                     ,int_state%CNVBOT,int_state%CNVTOP                 &
                     ,int_state%SM,int_state%LPBL                       &
                     ,int_state%HTOP,int_state%HTOPD,int_state%HTOPS    &
                     ,int_state%HBOT,int_state%HBOTD,int_state%HBOTS    &
                     ,int_state%AVCNVC,int_state%ACUTIM                 &
                     ,int_state%RSWIN,int_state%RSWOUT                  &
                     ,int_state%CONVECTION                              &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
!
          ELSEIF(int_state%CONVECTION=='bmj_dev')THEN
!
          CALL CUCNVC_DEV(NTIMESTEP,int_state%DT,int_state%NPRECIP      &
                     ,int_state%NRADS,int_state%NRADL                   &
                     ,int_state%MINUTES_HISTORY                         &
                     ,int_state%ENTRAIN,int_state%NEWALL                &
                     ,int_state%NEWSWAP,int_state%NEWUPUP               &
                     ,int_state%NODEEP                                  &
                     ,int_state%DYH,int_state%RESTART,int_state%HYDRO   &
                     ,int_state%CLDEFI,int_state%NUM_WATER              &
                     ,int_state%F_ICE,int_state%F_RAIN                  &
                     ,int_state%P_QV,int_state%P_QC,int_state%P_QR      &
                     ,int_state%P_QI,int_state%P_QS,int_state%P_QG      &
                     ,int_state%F_QV,int_state%F_QC,int_state%F_QR      &
                     ,int_state%F_QI,int_state%F_QS,int_state%F_QG      &
                     ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1                  &
                     ,int_state%PT,int_state%PD                         &
                     ,int_state%T,int_state%Q                           &
                     ,int_state%CW,int_state%TCUCN,int_state%WATER      &
                     ,int_state%OMGALF                                  &
                     ,int_state%U,int_state%V                           &
                     ,int_state%FIS,int_state%W0AVG                     &
                     ,int_state%PREC,int_state%ACPREC,int_state%CUPREC  &
                     ,int_state%CUPPT,int_state%CPRATE                  &
                     ,int_state%CNVBOT,int_state%CNVTOP                 &
                     ,int_state%SM,int_state%LPBL                       &
                     ,int_state%HTOP,int_state%HTOPD,int_state%HTOPS    &
                     ,int_state%HBOT,int_state%HBOTD,int_state%HBOTS    &
                     ,int_state%AVCNVC,int_state%ACUTIM                 &
                     ,int_state%RSWIN,int_state%RSWOUT                  &
                     ,int_state%CONVECTION                              &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
!
          ELSE
!
          write(0,*)'Wrong convection scheme choice'
          STOP
!
          ENDIF
!
          cucnvc_tim=cucnvc_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***    Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            pole_swap_phy_tim=pole_swap_phy_tim+(timef()-btim)
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
          CALL GSMDRIVE(NTIMESTEP,int_state%DT                             &
                       ,NPRECIP,int_state%NUM_WATER                        &
                       ,int_state%DXH(JC),int_state%DYH                    &
                       ,int_state%SM,int_state%FIS                         &
                       ,DSG2,SGML2,PDSG1,PSGML1                            &
                       ,int_state%PT,int_state%PD                          &
                       ,int_state%T,int_state%Q                            &
                       ,int_state%CW,int_state%OMGALF                      &
                       ,int_state%WATER                                    &
                       ,int_state%TRAIN,int_state%SR                       &
                       ,int_state%F_ICE,int_state%F_RAIN,int_state%F_RIMEF &
                       ,int_state%P_QV,int_state%P_QC,int_state%P_QR       &
                       ,int_state%P_QI,int_state%P_QS,int_state%P_QG       &
                       ,int_state%F_QV,int_state%F_QC,int_state%F_QR       &
                       ,int_state%F_QI,int_state%F_QS,int_state%F_QG       &
                       ,int_state%PREC,int_state%ACPREC,int_state%AVRAIN   &
                       ,int_state%MP_RESTART_STATE                         &
                       ,int_state%TBPVS_STATE,int_state%TBPVS0_STATE       &
                       ,int_state%SPECIFIED,int_state%NESTED               &
                       ,int_state%MICROPHYSICS                             &
                       ,IDS,IDE,JDS,JDE,LM                                 &
                       ,IMS,IME,JMS,JME                                    &
                       ,ITS,ITE,JTS,JTE)
!
          gsmdrive_tim=gsmdrive_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Precipitation Assimilation
!-----------------------------------------------------------------------
!
          IF (int_state%PCPFLG) THEN
!
            btim=timef()
            CALL CHKSNOW(MYPE                                           &
                        ,int_state%NTSD                                 &
                        ,int_state%DT                                   &
                        ,int_state%NPHS                                 &
                        ,int_state%SR                                   &
                        ,int_state%PPTDAT                               &
                        ,int_state%PCPHR                                &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE)
!
            CALL ADJPPT(MYPE                                            &
                       ,int_state%NTSD                                  &
                       ,int_state%DT                                    &
                       ,int_state%NPHS                                  &
                       ,int_state%PREC                                  &
                       ,int_state%LSPA                                  &
                       ,int_state%PPTDAT                                &
                       ,int_state%DDATA                                 &
                       ,int_state%PCPHR                                 &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE)
!
            adjppt_tim=adjppt_tim+(timef()-btim)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            pole_swap_phy_tim=pole_swap_phy_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***    Exchange Q and CW.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%Q,LM,int_state%CW,LM                 &
                        ,1,1)
!
          exch_phy_tim=exch_phy_tim+(timef()-btim)

!
!-----------------------------------------------------------------------
!
        ENDIF microphysics
!
!-----------------------------------------------------------------------
!***  Always exchange Temperature array since radiative updates
!***  are done every timestep.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,1,1)
!
        exch_phy_tim=exch_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  NOTE:  The Physics export state is fully updated now
!***         because subroutine PHY_INITIALIZE inserted the
!***         appropriate ESMF Fields into it.  Those Fields
!***         contain pointers to the actual data and those
!***         pointers are never re-directed.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ELSE gfs_phys_test                                                   !<-- Use GFS physics package
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!#######################################################################
!#######################################################################
!############ G F S   P H Y S I C S   D R I V E R ######################
!#######################################################################
!#######################################################################
!
        btim=timef()
!
        CALL ESMF_ClockGet(clock       =CLOCK                           &  !<-- The ESMF Clock
                          ,currTime    =CURRTIME                        &  !<-- The current time (ESMF) on the clock
                          ,rc          =RC)
!
        CALL ESMF_TimeGet(time=CURRTIME                                 &  !<-- The cuurent forecast time (ESMF)
                         ,yy  =JDAT(1)                                  &  !<-- The current forecast year (integer)
                         ,mm  =JDAT(2)                                  &  !<-- The current forecast month (integer)
                         ,dd  =JDAT(3)                                  &  !<-- The current forecast day (integer)
                         ,h   =JDAT(5)                                  &  !<-- The current forecast hour (integer)
                         ,m   =JDAT(6)                                  &  !<-- The current forecast minute (integer)
                         ,s   =JDAT(7)                                  &  !<-- The current forecast second (integer)
                         ,rc  =RC)
        JDAT(4)=0
        JDAT(8)=0
!
        DO J=JTS,JTE
          GLOBAL_LATS_R(J) = J-JTS+1
          LONSPERLAR(J)    = ITE-ITS+1
          SINLAT_R(J)      = SIN(int_state%GLAT( (ITS+ITE)/2 ,J))
          COSLAT_R(J)      = SQRT( 1.d0 - SINLAT_R(J)*SINLAT_R(J) )
          DO I=ITS,ITE
            XLON(I,J)        = int_state%GLON(I,J)
            IF(int_state%GLON(I,J)<0) &
             XLON(I,J)        = 2.0d0*3.14159d0+XLON(I,J)
            COSZEN(I,J)      = int_state%CZEN(I,J)
            COSZDG(I,J)      = int_state%CZMEAN(I,J)
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  GFS Radiation
!-----------------------------------------------------------------------
!
        CALL_GFS_PHY = MOD(NTIMESTEP,int_state%NPHS)==0

        FHSWR        = FLOAT(int_state%NRADS)*int_state%DT/3600.   ! [h]
        LSCCA        = MOD(NTIMESTEP+1,int_state%NRADS)==0         ! logical true during a step for which convective clouds
                                                                   ! are calculated from convective precipitation rates
        LSSWR        = MOD(NTIMESTEP,int_state%NRADS)==0
        LSLWR        = MOD(NTIMESTEP,int_state%NRADL)==0
!
!-----------------------------------------------------------------------
        lw_or_sw: IF (LSSWR .OR. LSLWR ) THEN
!-----------------------------------------------------------------------
!
          DO L=1,LM+1
            KFLIP=LM-L+2
            RSGM(KFLIP)=int_state%SGM(L)
          ENDDO
! 
          ICWP=0                  ! control flag for cloud generation schemes
          IF (NTCW > 0) ICWP = 1  ! 0: use diagnostic cloud scheme
                                  ! 1: use prognostic cloud scheme (default)
!
! ----
          CALL RADINIT ( RSGM, LM, IFLIP, NUM_P3D, ISOL, ICO2,  &
                         ICWP, IALB, IEMS, IAER, JDAT, MYPE )
! ----

          IF (NTOZ .LE. 0) THEN                ! Climatological Ozone
!
            IDAY   = JDAT(3)
            IMON   = JDAT(2)
            MIDMON = DAYS(IMON)/2 + 1
            CHANGE = FIRST .OR. ( (IDAY .EQ. MIDMON) .AND. (JDAT(5).EQ.0) )
!
            IF (CHANGE) THEN
              IF (IDAY .LT. MIDMON) THEN
                 K1OZ = MOD(IMON+10,12) + 1
                 MIDM = DAYS(K1OZ)/2 + 1
                 K2OZ = IMON
                 MIDP = DAYS(K1OZ) + MIDMON
              ELSE
                 K1OZ = IMON
                 MIDM = MIDMON
                 K2OZ = MOD(IMON,12) + 1
                 MIDP = DAYS(K2OZ)/2 + 1 + DAYS(K1OZ)
              ENDIF
            ENDIF
!
            IF (IDAY .LT. MIDMON) THEN
              ID = IDAY + DAYS(K1OZ)
            ELSE
              ID = IDAY
            ENDIF
!
            FACOZ = REAL (ID-MIDM) / REAL (MIDP-MIDM)
!
          ELSE
!
            K1OZ = 0
            K2OZ = 0
            FACOZ = 1.0D0
!
          ENDIF
!
        FLGMIN_L(1)     = 0.2D0      ! --- for ferrier (for now, any number)

! ----
          CALL ASTRONOMY                                                &
!  ---  inputs:
             ( LONSPERLAR, GLOBAL_LATS_R, SINLAT_R, COSLAT_R, XLON,     &
               FHSWR, JDAT,                                             &
               LONR, LATS_NODE_R, LATR, IPT_LATS_NODE_R, LSSWR, MYPE,   &
!  ---  outputs:
               int_state%SOLCON, int_state%SLAG, int_state%SDEC,        &
               int_state%CDEC, COSZEN, COSZDG )
!
!-----------------------------------------------------------------------
!
        ENDIF  lw_or_sw
!
!-----------------------------------------------------------------------
!
!---
        IF (FIRST) THEN
!
          SEED0 = JDAT(4) + JDAT(3) + JDAT(2) + JDAT(1)
          CALL RANDOM_SETSEED(SEED0)
          CALL RANDOM_NUMBER(WRK)
          SEED0 = SEED0 + NINT(WRK(1)*1000.0)
          FIRST = .FALSE.
!
        ENDIF
!---
        FHOUR=NTIMESTEP*int_state%DT/3600.d0
        ISEED = MOD(100.0*SQRT(FHOUR*3600),1.0d9) + 1 + SEED0
        CALL RANDOM_SETSEED(ISEED)
        CALL RANDOM_NUMBER(RANNUM)
        N=0
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          N=N+1
          XKT2(I,J) = RANNUM(N)
        ENDDO
        ENDDO
!---
        LSFWD=NTIMESTEP==0
        DTP=2.*int_state%DT
        DTF=int_state%DT
        IF(LSFWD) DTP=DTF
!---
        SOLHR=MOD(FHOUR+START_HOUR,24.d0)
!---
!...  set switch for saving convective clouds
        IF(LSCCA.AND.LSSWR) THEN
          CLSTP=1100+MIN(FHSWR,FHOUR,99.d0)  !initialize,accumulate,convert
        ELSEIF(LSCCA) THEN
          CLSTP=0100+MIN(FHSWR,FHOUR,99.d0)  !accumulate,convert
        ELSEIF(LSSWR) THEN
          CLSTP=1100                         !initialize,accumulate
        ELSE
          CLSTP=0100                         !accumulate
        ENDIF
!---
!---- OZONE ------------------------------------------------------------
!
        IF(.NOT.ALLOCATED(OZPLOUT_V)) &
                           ALLOCATE (OZPLOUT_V(LEVOZP,        PL_COEFF))
        IF(.NOT.ALLOCATED(OZPLOUT  )) &
                           ALLOCATE (OZPLOUT  (LEVOZP,JTS:JTE,PL_COEFF))
!
        IDATE(1)=JDAT(5)
        IDATE(2)=JDAT(2)
        IDATE(3)=JDAT(3)
        IDATE(4)=JDAT(1)
!
        IF (NTOZ .GT. 0) THEN
          CALL OZINTERPOL(MYPE,LATS_NODE_R,LATS_NODE_R,IDATE,FHOUR,     &
                          int_state%JINDX1,int_state%JINDX2,            &
                          int_state%OZPLIN,OZPLOUT,int_state%DDY)
        ENDIF
!
!---- OZONE ------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Set diagnostics to 0.
!-----------------------------------------------------------------------
!
        DT3DT=0.0d0
        DU3DT=0.0d0
        DV3DT=0.0d0
        DQ3DT=0.0d0
!
        CV (1) = 0.d0       !!!!! not in use if ntcw-1 > 0
        CVB(1) = 0.d0       !!!!! not in use if ntcw-1 > 0
        CVT(1) = 0.d0       !!!!! not in use if ntcw-1 > 0
!
!-----------------------------------------------------------------------
!***  Empty the radiation flux and precipitation arrays if it is time.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NRDLW)==0)THEN
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ALWIN(I,J) =0.
            int_state%ALWOUT(I,J)=0.
            int_state%ALWTOA(I,J)=0.
            int_state%ARDLW (I,J)=0.                                       !<-- An artificial 2-D array (ESMF cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
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
            int_state%ARDSW (I,J)=0.                                       !<-- An artificial 2-D array (ESMF cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
!
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NSRFC)==0)THEN
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACSNOW(I,J)=0.
            int_state%POTEVP(I,J)=0.
            int_state%SFCEVP(I,J)=0.
            int_state%SFCLHX(I,J)=0.
            int_state%SFCSHX(I,J)=0.
            int_state%SUBSHX(I,J)=0.
            int_state%BGROFF(I,J)=0.
            int_state%SSROFF(I,J)=0.
            int_state%ASRFC (I,J)=0.   !<-- An artificial 2-D array (ESMF cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
!
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACPREC(I,J)=0.
            int_state%CUPREC(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
        gfs_physics: IF(CALL_GFS_PHY)THEN
!-----------------------------------------------------------------------
!
             DTLW    = FLOAT(int_state%NRADL)*int_state%DT   ! [s]
             DTSW    = FLOAT(int_state%NRADS)*int_state%DT   ! [s]
             DTLWI   = 1./DTLW
             DTSWI   = 1./DTSW
             MINDT   = 1./MIN(DTLW,DTSW)
             DTPHS   = int_state%NPHS*int_state%DT
             DTPHSI  = 1./DTPHS
             XLVRW   = XLV*RHOWATER
             XLVRWI  = 1./XLVRW
             RoCP    = R/CP
!
!-----------------------------------------------------------------------
! ***  MAIN GFS-PHYS DOMAIN LOOP
!-----------------------------------------------------------------------
!
          j_loop: DO J=JTS,JTE
            JJ=J-JTS+1+JHALO   ! Pointer indices (Q & CW)
!
            i_loop: DO I=ITS,ITE
              II=I-ITS+1+IHALO ! Pointer indices (Q & CW)
!
              int_state%ACUTIM(I,J) = int_state%ACUTIM(I,J) + 1.     ! advance counters
              int_state%APHTIM(I,J) = int_state%APHTIM(I,J) + 1.
              int_state%ARDLW(I,J)  = int_state%ARDLW(I,J)  + 1.
              int_state%ARDSW(I,J)  = int_state%ARDSW(I,J)  + 1.
              int_state%ASRFC(I,J)  = int_state%ASRFC(I,J)  + 1.
              int_state%AVRAIN(I,J) = int_state%AVRAIN(I,J) + 1.
              int_state%AVCNVC(I,J) = int_state%AVCNVC(I,J) + 1.
!
              T1(1)           = 0.0D0        ! initialize all local variables
              Q1(1)           = 0.0D0        ! used in gfs_physics
              U1(1)           = 0.0D0
              V1(1)           = 0.0D0
              QSS(1)          = 0.0D0
              CDQ(1)          = 0.0D0
              HFLX(1)         = 0.0D0
              EVAP(1)         = 0.0D0
              DTSFC(1)        = 0.0D0
              DQSFC(1)        = 0.0D0
              DUSFC(1)        = 0.0D0
              DVSFC(1)        = 0.0D0
              PSMEAN(1)       = 0.0D0
              EPI(1)          = 0.0D0
              EVBSA(1)        = 0.0D0
              EVCWA(1)        = 0.0D0
              TRANSA(1)       = 0.0D0
              SBSNOA(1)       = 0.0D0
              SOILM(1)        = 0.0D0
              SNOWCA(1)       = 0.0D0
              CLDWRK(1)       = 0.0D0
              ZLVL(1)         = 0.0D0
              PHII            = 0.0D0
              PHIL            = 0.0D0
              CHH(1)          = 0.0D0
              HPBL(1)         = 0.0D0
              PSURF(1)        = 100.0D0
              T2M(1)          = 273.0D0
              Q2M(1)          = 0.0D0
              U10M(1)         = 0.0D0
              V10M(1)         = 0.0D0
              ADR             = 0.0D0
              ADT             = 0.0D0
              ADU             = 0.0D0
              ADV             = 0.0D0

         IF(int_state%TSKIN(I,J) .LT. 50. ) THEN
             TSEA(1)         = int_state%SST(I,J)
             TISFC(1)        = int_state%SST(I,J)
         ELSE
             TSEA(1)         = int_state%TSKIN(I,J)
             TISFC(1)        = int_state%TSKIN(I,J)
         ENDIF

         IF(int_state%SICE(I,J) > 0.5 ) THEN                                ! slmsk - ocean  - 0
             SLMSK(1)        = 2.0D0                                        !         land   - 1
         ELSE                                                               !         seaice - 2
             SLMSK(1)        = 1.0D0-int_state%SM(I,J)                      !
         ENDIF

         DO L=1,LM
            KFLIP=LM+1-L
!            CLDCOV_V(KFLIP) = 0.0D0                      ! GRRAD now returns instant cloud cover (Sarah Lu)
             F_ICE(KFLIP)    = int_state%F_ICE(I,J,L)                       ! for ferrier phy, do init first
             F_RAIN(KFLIP)   = int_state%F_RAIN(I,J,L)
             R_RIME(KFLIP)   = int_state%F_RIMEF(I,J,L)
         ENDDO

             XLAT(1)         = int_state%GLAT(I,J)
             ZORL(1)         = int_state%ZORFCS(I,J)
             SNCOVR(1)       = int_state%SNO(I,J)/(int_state%SNO(I,J)+70.)  ! FORMULATION OF MARSHALL ET AL. 1994
                                                                            ! change this later only initially, add new int_state
             SNWDPH(1)       = int_state%SI(I,J)                            ! snwdph[mm]
             SHELEG(1)       = int_state%SNO(I,J)                           ! snow water eq.[mm]
             SNOALB(1)       = int_state%MXSNAL(I,J)
             ALVSF(1)        = int_state%ALBFC1(I,J,1)                      ! VIS, direct
             ALVWF(1)        = int_state%ALBFC1(I,J,2)                      ! VIS, diffuse
             ALNSF(1)        = int_state%ALBFC1(I,J,3)                      ! NIR, direct
             ALNWF(1)        = int_state%ALBFC1(I,J,4)                      ! NIR, diffuse
             FACSF(1)        = int_state%ALFFC1(I,J,1)                      ! direct
             FACWF(1)        = int_state%ALFFC1(I,J,2)                      ! diffuse
!
             PRSI (LM+1)     = int_state%PT/1000.                           ! [kPa]
             PRSIK(LM+1)     = (PRSI(LM+1)*0.01d0)**RoCP
         DO L=1,LM
            KFLIP=LM+1-L
             PRSI (KFLIP)    = PRSI(KFLIP+1) + &
                                (DSG2(L)*int_state%PD(I,J)+PDSG1(L))/1000.d0 ! (pressure on interface) [kPa]
             PRSIK(KFLIP)    = (PRSI(KFLIP)*0.01d0)**RoCP

             PRSL (KFLIP)    = (PRSI(KFLIP)+PRSI(KFLIP+1))*0.5d0             ! (pressure on mid-layer) [kPa]
             PRSLK(KFLIP)    = (PRSL(KFLIP)*0.01d0)**RoCP
!
             RTvR = 1. / ( R * (int_state%Q(II,JJ,L)*0.608+1.-int_state%CW(II,JJ,L) ) * int_state%T(I,J,L) )
             VVEL(KFLIP)     = int_state%OMGALF(I,J,L) * 1000.d0*PRSL(KFLIP) * RTvR
!
             GU(KFLIP)       = (int_state%U(I,J  ,L) + int_state%U(I-1,J  ,L) +                    &
                                int_state%U(I,J-1,L) + int_state%U(I-1,J-1,L))*COSLAT_R(J)*0.25d0
             GV(KFLIP)       = (int_state%V(I,J  ,L) + int_state%V(I-1,J  ,L) +                    &
                                int_state%V(I,J-1,L) + int_state%V(I-1,J-1,L))*COSLAT_R(J)*0.25d0
             GT(KFLIP)       = int_state%T(I,J,L)
             GR(KFLIP)       = int_state%Q(II,JJ,L)
             GR3(KFLIP,1)    = int_state%Q(II,JJ,L)
           IF (NTIMESTEP == 0 ) THEN
             GR3(KFLIP,2)    = 0.0d0
             GR3(KFLIP,3)    = 0.0d0
           ELSE
             GR3(KFLIP,2)    = int_state%O3(II,JJ,L)
             GR3(KFLIP,3)    = int_state%CW(II,JJ,L)
           ENDIF
             GR1(1,KFLIP,1)  = GR3(KFLIP,2)
             GR1(1,KFLIP,2)  = int_state%CW(II,JJ,L)
         ENDDO
!---
             DLWSFC(1)       = int_state%ALWIN(I,J)
             ULWSFC(1)       = int_state%ALWOUT(I,J)
             DLWSFCI(1)      = int_state%RLWIN(I,J)
             ULWSFCI(1)      = int_state%RADOT(I,J)
             DSWSFCI(1)      = int_state%RSWIN(I,J)
             USWSFCI(1)      = int_state%RSWOUT(I,J)
!---
             GFLUX(1)        = 0.0D0
             DTSFCI(1)       = -int_state%TWBS(I,J)
             DQSFCI(1)       = -int_state%QWBS(I,J)
             GFLUXI(1)       = 0.0D0
             EP(1)           = int_state%POTEVP(I,J)*XLVRW
!---
             XSIHFCS(1)      = int_state%SIHFCS(I,J)
             XSICFCS(1)      = int_state%SICFCS(I,J)
             XSLPFCS(1)      = int_state%SLPFCS(I,J)
             XTG3FCS(1)      = int_state%TG3FCS(I,J)
             XVEGFCS(1)      = int_state%VEGFCS(I,J)
             XVETFCS(1)      = int_state%VETFCS(I,J)
             XSOTFCS(1)      = int_state%SOTFCS(I,J)
!---
             FLUXR_V         = 0.0D0
             IF(.NOT.LSLWR) FLUXR_V(1) = int_state%RLWTOA(I,J)*DTLW
             IF(.NOT.LSSWR) FLUXR_V(2) = int_state%RSWTOA(I,J)*DTSW
!---
             HPRIME (1)      = int_state%HSTDV(I,J)
             HPRIME (2)      = int_state%HCNVX(I,J)
             HPRIME (3)      = int_state%HASYW(I,J)
             HPRIME (4)      = int_state%HASYS(I,J)
             HPRIME (5)      = int_state%HASYSW(I,J)
             HPRIME (6)      = int_state%HASYNW(I,J)
             HPRIME (7)      = int_state%HLENW(I,J)
             HPRIME (8)      = int_state%HLENS(I,J)
             HPRIME (9)      = int_state%HLENSW(I,J)
             HPRIME(10)      = int_state%HLENNW(I,J)
             HPRIME(11)      = int_state%HANGL(I,J)*180.D0/3.14159D0
             HPRIME(12)      = int_state%HANIS(I,J)
             HPRIME(13)      = int_state%HSLOP(I,J)*0.2d0   ! calculated from different terrain file
             HPRIME(14)      = int_state%HZMAX(I,J)
!---
             RUNOFF(1)       = int_state%BGROFF(I,J)*0.001D0
             SRUNOFF(1)      = int_state%SSROFF(I,J)*0.001D0
!---
             SFCNSW(1)       = int_state%SFCNSW(I,J)
             SFCDSW(1)       = int_state%SFCDSW(I,J)
             SFALB(1)        = int_state%SFALB(I,J)
             SFCDLW(1)       = int_state%SFCDLW(I,J)
             TSFLW(1)        = int_state%TSFLW(I,J)
           DO L=1,LM
             SWH(L)          = int_state%SWH(I,J,L)
             HLW(L)          = int_state%HLW(I,J,L)
           ENDDO
           DO N=1,3                                    ! for Zhao =3, Ferr=1
             PHY_F2DV(N)     = int_state%PHY_F2DV (I,J,N)
           ENDDO
           DO N=1,4                                    ! for Zhao =4, Ferr=3
           DO L=1,LM
            PHY_F3DV(L,N)    = int_state%PHY_F3DV (I,J,L,N)
           ENDDO
           ENDDO
!
!-----------------------------------------------------------------------
          CALL GRRAD                                                 &
!-----------------------------------------------------------------------
!  ---  inputs:
           ( PRSI, PRSL, PRSLK, GT, GR, GR1, VVEL, SLMSK,            &
             XLON(I,J), XLAT, TSEA,                                  &
             SNWDPH, SNCOVR, SNOALB,                                 &
             ZORL, HPRIME(1),                                        &
             ALVSF, ALNSF, ALVWF,                                    &
             ALNWF, FACSF, FACWF,                                    &
             XSICFCS, TISFC,                                         &
             int_state%SOLCON,                                       &
             COSZEN(I,J), COSZDG(I,J), K1OZ, K2OZ, FACOZ,            &
             CV, CVT, CVB,                                           &
             IOVR_SW, IOVR_LW, F_ICE, F_RAIN, R_RIME, FLGMIN_L,      &
             NUM_P3D, NTCW-1, NCLD, NTOZ-1, NTRAC-1, NFXR,           &
             DTLW, DTSW, LSSWR, LSLWR, LSSAV, SASHAL,                &
             1, 1, LM, IFLIP, MYPE, LPRNT,                           &
!  ---  outputs:
             SWH, SFCNSW, SFCDSW,                                    &
             SFALB,                                                  &
             HLW, SFCDLW, TSFLW,                                     &
             CLDCOV_V,                                               & 
!  ---  input/output:
             FLUXR_V                                                 &
           )
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  GBPHYS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---
          DPSHC(1)    = 0.3 * PRSI(1)
          GQ(1)       = PRSI(1)
!---
          RANNUM_V(1) = XKT2(I,J)
!---
          RCS2_V(1)   = 1.0d0/(0.0001d0+COS(XLAT(1))*COS(XLAT(1))) ! fixed XLAT=+/-90
!---
          TPRCP(1)    = int_state%PREC(I,J)
          BENGSH(1)   = int_state%CUPREC(I,J)
          GESHEM(1)   = int_state%ACPREC(I,J)

         DO L=1,NUM_SOIL_LAYERS
          SMC_V(L)    = int_state%SMC(I,J,L)
          STC_V(L)    = int_state%STC(I,J,L)
          SLC_V(L)    = int_state%SH2O(I,J,L)
         ENDDO

         SHDMIN(1)   = int_state%SHDMIN(I,J)
         SHDMAX(1)   = int_state%SHDMAX(I,J)

         UUSTAR(1)   = int_state%USTAR(I,J)
         CANOPY(1)   = int_state%CMC(I,J)*1000.
         FLGMIN      = 0.2d0 !!! put this in ferrier init
         CRTRH       = 0.85d0
         FLIPV       = .FALSE.
         NCW(1)      = 50
         NCW(2)      = 150
         OLD_MONIN   = .TRUE.
         CNVGWD      = .FALSE.
         NEWSAS      = .FALSE.
         CCWF        = 0.5d0  ! only for RAS scheme
!
! ---- ESTIMATE T850 FOR RAIN-SNOW DECISION ----------------------------
!
          T850 = GT(1)
      DO L = 1, LM - 1
        IF(PRSL(L) .GT. 85.d0 .AND. PRSL(L+1) .LE. 85.d0) THEN
          T850 = GT(L) - (PRSL(L)-85.) / (PRSL(L)-PRSL(L+1)) * (GT(L)-GT(L+1))
        ENDIF
      ENDDO
!
      SRFLAG(1) = 0.0d0
      IF(T850 .LE. 273.16d0) SRFLAG(1) = 1.0d0
!
!---- OZONE ------------------------------------------------------------
      IF (NTOZ .GT. 0) THEN
        DO N=1,PL_COEFF
          DO L=1,LEVOZP
              OZPLOUT_V(L,N) = OZPLOUT(L,j,N)
          ENDDO
        ENDDO
      ELSE
              OZPLOUT_V      = 0.0d0
      ENDIF
!-----------------------------------------------------------------------
!
      CALL GBPHYS(1, 1, LM, NUM_SOIL_LAYERS, LSM, NTRAC, NCLD,                 &
           NTOZ, NTCW, NMTVR, LONR, LATR, 62, RAS, LONR, RANNUM_V, 1, PRE_RAD, &
           GU, GV, GQ, GT, GR3, VVEL,                                          &
           ADT, ADR, ADU, ADV, SINLAT_R(J), COSLAT_R(J), RCS2_V,               &
           PRSI, PRSL, PRSLK, PRSIK, PHII, PHIL, DPSHC, FHOUR, LSSAV, SOLHR,   &
           LSFWD, CLSTP, DTP, DTF, PL_PRES, OZPLOUT_V, LEVOZP, PL_COEFF,       &
           XSIHFCS, XSICFCS, TISFC, SFCDSW,                                    &
           TPRCP, SRFLAG,                                                      &
           SLC_V, SNWDPH, XSLPFCS, SHDMIN, SHDMAX, SNOALB, SFALB,              &
           CHH, CMM, EPI, DLWSFCI, ULWSFCI, USWSFCI, DSWSFCI, DTSFCI,          &
           DQSFCI, GFLUXI,           SRUNOFF     , T1, Q1, U1, V1, ZLVL,       &
           EVBSA, EVCWA,                                                       &
           TRANSA, SBSNOA, SNOWCA,                                             &
!          SOILM, TSEA, SHELEG, SNCOVR, XTG3FCS, ZORL, CV, CVB, CVT,           &
           SOILM, RAIN, RAINC, TSEA, SHELEG, SNCOVR, XTG3FCS, ZORL, CV, CVB, CVT,&   
           SLMSK, XVEGFCS, CANOPY, F10M, XVETFCS, XSOTFCS, UUSTAR, FFMM, FFHH, &
           int_state%TMPMIN(I,J), int_state%TMPMAX(I,J), GESHEM,               &
           DUSFC, DVSFC, DTSFC, DQSFC, DLWSFC, ULWSFC, GFLUX,RUNOFF, EP,       &
           CLDWRK, int_state%DUGWD(I,J), int_state%DVGWD(I,J),                 &
           PSMEAN, BENGSH, XLON(I,J),                                          &
           COSZEN(I,J), SFCNSW, XLAT,                                          &
           SFCDLW, TSFLW, PSURF, U10M, V10M, T2M, Q2M,                         &
           HPBL, PWAT, SWH, HLW, SMC_V, STC_V, HPRIME, int_state%SLAG,         &
           int_state%SDEC, int_state%CDEC,                                     &
           ACV, ACVB, ACVT,                                                    &
           PHY_F3DV, PHY_F2DV, NUM_P3D, NUM_P2D, FLGMIN,                       &
!          DT3DT, DQ3DT, DU3DT, DV3DT, UPD_MF, DWN_MF, DET_MF, LDIAG3D,        &  
           DT3DT, DQ3DT, DU3DT, DV3DT, DQDT, UPD_MF, DWN_MF, DET_MF, LDIAG3D,  &  
           FLIPV, MYPE, NTIMESTEP, J-JTS+1, ORO,                               &
           CRTRH, NCW, OLD_MONIN, CNVGWD, CCWF, SASHAL, NEWSAS)
!
!-----------------------------------------------------------------------
! ***     UPDATE AFTER PHYSICS
!-----------------------------------------------------------------------
             int_state%SIHFCS(I,J)        = XSIHFCS(1)
             int_state%SICFCS(I,J)        = XSICFCS(1)

             int_state%CZEN(I,J)          = COSZEN(I,J)
             int_state%CZMEAN(I,J)        = COSZDG(I,J)

             int_state%SFCNSW(I,J)        = SFCNSW(1)
             int_state%SFCDSW(I,J)        = SFCDSW(1)
             int_state%SFALB (I,J)        = SFALB(1)
             int_state%SFCDLW(I,J)        = SFCDLW(1)
             int_state%TSFLW (I,J)        = TSFLW(1)

             int_state%SI(I,J)            = SNWDPH(1)
             int_state%SNO(I,J)           = SHELEG(1)
             int_state%MXSNAL(I,J)        = SNOALB(1)

         DO L=1,NUM_SOIL_LAYERS
             int_state%SMC(I,J,L)         = SMC_V(L)
             int_state%STC(I,J,L)         = STC_V(L)
             int_state%SH2O(I,J,L)        = SLC_V(L)
         ENDDO
             int_state%CMC(I,J)           = CANOPY(1)*0.001
             int_state%USTAR(I,J)         = UUSTAR(1)
             int_state%SMSTOT(I,J)        = SOILM(1)*1000.

         DO L=1,LM
             int_state%SWH(I,J,L)         = SWH(L)
             int_state%HLW(I,J,L)         = HLW(L)
         ENDDO

         DO N=1,3                                    ! for Zhao =3, Ferr=1
             int_state%PHY_F2DV (I,J,N)   = PHY_F2DV(N)
         ENDDO
         DO N=1,4                                    ! for Zhao =4, Ferr=3
         DO L=1,LM
             int_state%PHY_F3DV (I,J,L,N) = PHY_F3DV(L,N)
         ENDDO
         ENDDO

             int_state%ALWIN(I,J)         = int_state%ALWIN(I,J)  + int_state%RLWIN(I,J)
             int_state%ALWOUT(I,J)        = int_state%ALWOUT(I,J) - int_state%RADOT(I,J)
             int_state%ASWIN(I,J)         = int_state%ASWIN(I,J)  + int_state%RSWIN(I,J)
             int_state%ASWOUT(I,J)        = int_state%ASWOUT(I,J) - int_state%RSWOUT(I,J)
             int_state%RLWIN(I,J)         = DLWSFCI(1)
             int_state%RADOT(I,J)         = ULWSFCI(1)
             int_state%RSWIN(I,J)         = DSWSFCI(1)
             int_state%RSWOUT(I,J)        = USWSFCI(1)
             int_state%RSWINC(I,J)        = int_state%RSWIN(I,J)/(1.-int_state%ALBEDO(I,J))

             int_state%RLWTOA(I,J)        = FLUXR_V(1)*DTLWI
             int_state%RSWTOA(I,J)        = FLUXR_V(2)*DTSWI
             int_state%ALWTOA(I,J)        = int_state%ALWTOA(I,J) + FLUXR_V(1)*DTLWI
             int_state%ASWTOA(I,J)        = int_state%ASWTOA(I,J) + FLUXR_V(2)*DTSWI

             int_state%SFCSHX(I,J)        = int_state%SFCSHX(I,J) - HFLX(1)*1000.*PRSL(1)*RTvR*CP     !need HFLX from GFS
             int_state%SFCLHX(I,J)        = int_state%SFCLHX(I,J) - EVAP(1)*1000.*PRSL(1)*RTvR*XLV    !need EVAP from GFS
             int_state%SUBSHX(I,J)        = int_state%SUBSHX(I,J) + GFLUXI(1)
             int_state%TWBS(I,J)          = -DTSFCI(1)
             int_state%QWBS(I,J)          = -DQSFCI(1)
             int_state%GRNFLX(I,J)        = GFLUXI(1)
             int_state%POTEVP(I,J)        = EP(1)*XLVRWI
             int_state%POTFLX(I,J)        = -EP(1)*DTPHSI
             int_state%SFCEVP(I,J)        = int_state%SFCEVP(I,J) + DQSFCI(1)*DTPHS*XLVRWI

             int_state%SFCEXC(I,J)        = CDQ(1)                                                    !need CDQ  from GFS
             int_state%PBLH(I,J)          = HPBL(1)
             int_state%PSFC(I,J)          = PSURF(1)*1000.
             int_state%PREC(I,J)          = TPRCP(1)
             int_state%CUPPT(I,J)         = BENGSH(1)-int_state%CUPREC(I,J)
             int_state%CPRATE(I,J)        = BENGSH(1)-int_state%CUPREC(I,J)
             int_state%CUPREC(I,J)        = BENGSH(1)
             int_state%ACPREC(I,J)        = GESHEM(1)

             int_state%BGROFF(I,J)        = RUNOFF(1)*1000.
             int_state%SSROFF(I,J)        = SRUNOFF(1)*1000.

             int_state%TSKIN(I,J)         = TISFC(1)
             int_state%SST(I,J)           = TSEA(1)
             int_state%SOILTB(I,J)        = XTG3FCS(1)
             IF( SRFLAG(1) >= 0.5 .AND. SLMSK(1) >= 0.5 ) &
             int_state%ACSNOW(I,J)        = int_state%ACSNOW(I,J) + int_state%ACPREC(I,J)

             int_state%PSHLTR(I,J)        = PSURF(1)*EXP(0.06823/T2M(1))*1000.
             int_state%TSHLTR(I,J)        = T2M(1)
             int_state%QSHLTR(I,J)        = Q2M(1)
             int_state%QSH(I,J)           = QSS(1)                                                    !need QSS  from GFS
             int_state%T2(I,J)            = T2M(1)
             int_state%TH02(I,J)          = T2M(1)*(100./PSURF(1))**RoCP
             int_state%Q02(I,J)           = Q2M(1)
             int_state%U10(I,J)           = U10M(1)
             int_state%V10(I,J)           = V10M(1)
             int_state%THS(I,J)           = TSFLW(1)*(100./PSURF(1))**RoCP
             int_state%SIGT4(I,J)         = int_state%T(I,J,LM)*int_state%T(I,J,LM) * &
                                            int_state%T(I,J,LM)*int_state%T(I,J,LM) * STBOLT
 
         DO L=1,LM
            KFLIP=LM+1-L
             int_state%T(I,J,L)           = ADT(KFLIP)
             int_state%DUDT(I,J,L)        = (ADU(KFLIP) - GU(KFLIP)) / (COSLAT_R(J) + 0.0001) / DTP
             int_state%DVDT(I,J,L)        = (ADV(KFLIP) - GV(KFLIP)) / (COSLAT_R(J) + 0.0001) / DTP
!*           int_state%CLDFRA(I,J,L)      = CLDCOV_V(KFLIP) * MINDT  ! scaling not needed (Sarah Lu)
             int_state%CLDFRA(I,J,L)      = CLDCOV_V(KFLIP) 
             int_state%Q (II,JJ,L)        = ADR(KFLIP,1)
             int_state%O3(II,JJ,L)        = ADR(KFLIP,2)
             int_state%CW(II,JJ,L)        = ADR(KFLIP,3)
         ENDDO
!
!-----------------------------------------------------------------------
! ***     END UPDATE AFTER PHYSICS
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! ***  END GFS-PHYS DOMAIN LOOP
!-----------------------------------------------------------------------
!
            ENDDO  i_loop
!
          ENDDO    j_loop
!
!-----------------------------------------------------------------------
!
          gfs_phy_tim=gfs_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  EXCHANGE WIND TENDENCIES
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,3,3)
!
          exch_phy_tim=exch_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  NOW INTERPOLATE WIND TENDENCIES FROM H TO V POINTS.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL H_TO_V_TEND(int_state%DUDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%U)
          CALL H_TO_V_TEND(int_state%DVDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%V)
!
          h_to_v_tim=h_to_v_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  POLES AND EAST-WEST BOUNDARY.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM      &
                       ,int_state%INPES,int_state%JNPES)
!
            pole_swap_phy_tim=pole_swap_phy_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***  EXCHANGE U, V, T, Q and CW
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%T,LM                                 &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%Q,LM,int_state%CW,LM                 &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%O3,LM                                &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                  &
                        ,3,3)
!
          exch_phy_tim=exch_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        ENDIF gfs_physics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!#######################################################################
!#######################################################################
!######### E N D   O F   G F S   P H Y S I C S   D R I V E R ###########
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!
      ENDIF  gfs_phys_test 
!
!-----------------------------------------------------------------------
!***  Write precipitation files for ADJPPT regression test
!-----------------------------------------------------------------------
!
      IF( int_state%WRITE_PREC_ADJ   .AND.                              &
          MOD(XTIME,60.) <= 0.001    .AND.                              &
          INT(XTIME/60.) <= int_state%PCPHR ) THEN
        CALL WRT_2D(int_state%PREC,MYPE,NUM_PES,MPI_COMM_COMP           &
                ,INT(XTIME/60.)+1                                       &
                ,IDS,IDE,JDS,JDE                                        &
                ,IMS,IME,JMS,JME                                        &
                ,ITS,ITE,JTS,JTE)
      ENDIF
!
!-----------------------------------------------------------------------
!     if(ntimestep<=5)then
!       call twr(int_state%t,lm,'t_phy',ntimestep,mype,num_pes,mpi_comm_comp &
!               ,ids,ide,jds,jde &
!               ,ims,ime,jms,jme &
!               ,its,ite,jts,jte)
!     endif
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
      phy_run_tim=phy_run_tim+(timef()-btim0)
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
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
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
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The Physics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Physics export state
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
                                   ,SBD,WBD                             &
                                   ,DPHD,DLMD                           &
                                   ,TPH0D,TLM0D                         &
                                   ,MY_DOMAIN_ID                        &
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
!***  Only for GFS physics
!-----------------------------------------------------------------------
!
      USE FUNCPHYS
      USE MERSENNE_TWISTER
      USE LAYOUT1,          ONLY : LATS_NODE_R,IPT_LATS_NODE_R
      USE TRACER_CONST,     ONLY : SET_TRACER_CONST
      USE DATE_DEF,         ONLY : FHOUR
      USE RESOL_DEF,        ONLY : LSOIL,LEVR,NXPT,JCAP,LEVS,NYPT       &
                                  ,JINTMX,THERMODYN_ID,SFCPRESS_ID      &
                                  ,NUM_P3D,NUM_P2D,NTOZ,NTCW,NCLD       &
                                  ,NMTVR,NFXR,LONR,LATR

      USE OZNE_DEF,         ONLY: LEVOZC,LATSOZP,BLATC,TIMEOZC,TIMEOZ   &
                                 ,KOZPL,LEVOZP,PL_TIME,PL_LAT,PL_PRES   &
                                 ,KOZC,DPHIOZC,LATSOZC,PL_COEFF
!
      USE NAMELIST_PHYSICS_DEF, ONLY: ISOL,ICO2,IALB,IEMS,IAER,IOVR_SW  &
                                     ,IOVR_LW,LSSAV,LDIAG3D,FHCYC       &
                                     ,SASHAL,PRE_RAD,RAS,LSM
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(KIND=KINT),INTENT(IN) :: CO2TF,MY_DOMAIN_ID
!
      INTEGER(KIND=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE,LM               &
                                      ,IMS,IME,JMS,JME                  &
                                      ,ITS,ITE,JTS,JTE
!
      REAL(KIND=KFPT),INTENT(INOUT) :: DLMD,DPHD                        &
                                      ,TPH0D,TLM0D                      &
                                      ,SBD,WBD
!
      LOGICAL,INTENT(IN) :: GFS
!
      CHARACTER(99),INTENT(IN) :: CONVECTION,LONGWAVE,MICROPHYSICS      &
                                 ,SFC_LAYER,SHORTWAVE,TURBULENCE        &
                                 ,LAND_SURFACE
!
!---------------------
!***  Local variables
!---------------------
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
             ,DLM,DPH,DT,DT_MICRO,DTPHS                                 &
             ,GMT,JULIAN,PDBOT,PDTOP,PDTOT,PT_CB,RELM,RPDTOT            &
             ,SB,SPH,STLM,STPH,STPH0,THETA_HALF                         &
             ,TLM,TLM_BASE,TPH,TPH_BASE,TPH0,TPV,WB,XTIME
!
      REAL,DIMENSION(LM) :: DSG1,DSG2,PDSG1,PSGML1,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM                          &
                             ,SFULL,SFULL_FLIP,SMID,SMID_FLIP
!
      REAL,DIMENSION(:),ALLOCATABLE,TARGET :: DXH,DXV,RDXH,RDXV
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: EMISS
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP1,TEMP_GWD
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: TEMPSOIL
      REAL,DIMENSION(NUM_SOIL_LAYERS)   :: SOIL1DIN
!
      CHARACTER(ESMF_MAXSTR) :: INFILE
!
      LOGICAL,SAVE :: ALLOWED_TO_READ=.TRUE.
      LOGICAL :: OPENED
!
!---------------------------------
!***  GFS physics local variables
!---------------------------------
!
      CHARACTER(80)   :: GFS_PHY_NAMELIST
      INTEGER         :: JDAT(8),RC,NLUNIT,NTRAC,IRET,IMJM,NIJ
      REAL(KIND=KDBL) :: DELTIM,GAUL

      REAL *4              :: BLATC4
      REAL *4, ALLOCATABLE :: PL_LAT4(:), PL_PRES4(:), PL_TIME4(:), TEMPIN(:)

      REAL(KIND=KDBL),DIMENSION(:),ALLOCATABLE ::                             &
                     SIG1T, RLA, RLO, SLMASK, OROG, AISFCS,                   &
                     SIHFCS, SICFCS, SITFCS, SWDFCS, VMNFCS, VMXFCS, SLPFCS,  &
                     ABSFCS, TSFFCS, SNOFCS, ZORFCS, TG3FCS, CNPFCS, SLIFCS,  &
                     F10MFCS, VEGFCS, VETFCS, SOTFCS, CVFCS, CVBFCS, CVTFCS
!
      REAL(KIND=KDBL),DIMENSION(:,:),ALLOCATABLE :: ALFFC1              &
                                                   ,SMCFC1              &
                                                   ,STCFC1              &
                                                   ,SLCFC1              &
                                                   ,ALBFC1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NSOIL=NUM_SOIL_LAYERS                                              !<-- From Landsurface module
!
!-----------------------------------------------------------------------
!***  Dereference the start time.
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
!-----------------------------------------------------------------------
!***  Radiation needs some specific time quantities.
!-----------------------------------------------------------------------
!
      CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR     &
                       ,START_MINUTE,START_SECOND                       &
                       ,NTIMESTEP,DT                                    &
                       ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
!***  Open and read GWD data file (14 orography fields)
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
        INFILE='GWD.bin'
!
!-----------------------------------------------------------------------
!
        CALL PHYSICS_READ_GWD(INFILE,NFCST,INT_STATE                    &
                             ,MYPE,MPI_COMM_COMP                        &
                             ,IDS,IDE,JDS,JDE)
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
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(.NOT.int_state%RESTART) THEN           ! COLD START
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        IF(.NOT.int_state%NEMSIO_INPUT)THEN    
!!!       INFILE='main_input_filename'
          WRITE(INFILE,'(A,I2.2)')'input_domain_',MY_DOMAIN_ID
!
          CALL PHYSICS_READ_INPUT_BINARY(INFILE,NFCST                   &
                                        ,MYPE,MPI_COMM_COMP             &
                                        ,IDAT,IHRST,PT                  &
                                        ,INT_STATE                      &
                                        ,NSOIL,LM                       &
                                        ,IDS,IDE,JDS,JDE                &
                                        ,IMS,IME,JMS,JME                &
                                        ,IRTN )
!
        ELSE
          WRITE(INFILE,'(A,I2.2,A)')'input_domain_',MY_DOMAIN_ID,'_nemsio'
          CALL PHYSICS_READ_INPUT_NEMSIO(INFILE                         &
                                        ,MYPE,MPI_COMM_COMP             &
                                        ,IDAT,IHRST,PT                  &
                                        ,INT_STATE                      &
                                        ,NSOIL,LM                       &
                                        ,IDS,IDE,JDS,JDE                &
                                        ,IMS,IME,JMS,JME                &
                                        ,IRTN )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        PT_CB=PT*1.0E-3      !<-- Convert pascals to centibars for GFDL initialization
!
!-----------------------------------------------------------------------
!***  CHECK TO SEE IF THE STARTING DATE/TIME IN THE INPUT DATA FILE
!***  AGREES WITH THAT IN THE CONFGURE FILE.
!-----------------------------------------------------------------------
!
        IF(MYPE==0)THEN
          IF(IDAT(2)/=START_MONTH.OR.                                   &
             IDAT(1)/=START_DAY.OR.                                     &
             IDAT(3)/=START_YEAR.OR.                                    &
             IHRST  /=START_HOUR)THEN
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' DATES IN INPUT FILE AND CONFIGURE FILE DISAGREE!!'
            WRITE(0,*)' INPUT: HOUR=',IHRST,' DAY=',IDAT(3)               &
                      ,' MONTH=',IDAT(2),' YEAR=',IDAT(1)
            WRITE(0,*)' CONFIG: HOUR=',START_HOUR,' DAY=',START_DAY       &
                      ,' MONTH=',START_MONTH,' YEAR=',START_YEAR
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          ENDIF
        ENDIF
!
        IYEAR_FCST =IDAT(3)
        IMONTH_FCST=IDAT(2)
        IDAY_FCST  =IDAT(1)
        IHOUR_FCST =IHRST
!
!-----------------------------------------------------------------------
!
      ELSE                                      ! RESTART READ INIT
!
!-----------------------------------------------------------------------
!
        IF(.NOT.int_state%NEMSIO_INPUT)THEN  
!
          WRITE(INFILE,'(A,I2.2)')'restart_file_',MY_DOMAIN_ID
!
          CALL PHYSICS_READ_RESTT_BINARY(INFILE,NFCST                   &
                                        ,MYPE,MPI_COMM_COMP             &
                                        ,IYEAR_FCST                     &
                                        ,IMONTH_FCST                    &
                                        ,IDAY_FCST                      &
                                        ,IHOUR_FCST                     &
                                        ,IMINUTE_FCST                   &
                                        ,SECOND_FCST                    &
                                        ,IHRST,IDAT,PT                  &
                                        ,INT_STATE                      &
                                        ,NSOIL,LM                       &
                                        ,IDS,IDE,JDS,JDE                &
                                        ,IMS,IME,JMS,JME                &
                                        ,IRTN )
!
        ELSE
!
          WRITE(INFILE,'(A,I2.2,A)')'restart_file_',MY_DOMAIN_ID,'_nemsio'
!
          CALL PHYSICS_READ_RESTT_NEMSIO(INFILE                         &
                                        ,MYPE,MPI_COMM_COMP             &
                                        ,IYEAR_FCST,IMONTH_FCST         &
                                        ,IDAY_FCST,IHOUR_FCST           &
                                        ,IMINUTE_FCST,SECOND_FCST       &
                                        ,IHRST,IDAT,PT                  &
                                        ,INT_STATE                      &
                                        ,NSOIL,LM                       &
                                        ,IDS,IDE,JDS,JDE                &
                                        ,IMS,IME,JMS,JME                &
                                        ,IRTN )
!
        ENDIF
!
        PT_CB=PT*1.0E-3   !<-- Convert pascals to centibars for GFDL initialization
!
!-----------------------------------------------------------------------
!
      ENDIF                                     ! COLD START /RESTART
!
!-----------------------------------------------------------------------
!***  MAKE UP A POTENTIAL SKIN TEMPERATURE.
!-----------------------------------------------------------------------
!
      IF(.NOT.int_state%RESTART) THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%THS(I,J)=int_state%TSKIN(I,J)                       &
                       *(100000./(int_state%SG2(LM+1)*int_state%PD(I,J) &
                                 +int_state%PSG1(LM+1)))**CAPPA
        ENDDO
        ENDDO
!
      ENDIF

!-----------------------------------------------------------------------
!*** INITIALIZING TLMAX, TLMIN
!-----------------------------------------------------------------------

      DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%TLMAX(I,J)=int_state%T(I,J,1)
          int_state%TLMIN(I,J)=int_state%T(I,J,1)
       ENDDO
     ENDDO

!
!-----------------------------------------------------------------------
!***  RECREATE SIGMA VALUES AT LAYER INTERFACES FOR THE FULL VERTICAL
!***  DOMAIN. 
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        SFULL(L)=int_state%SGM(L)
      ENDDO
!
      DO L=1,LM
        SMID(L)=(SFULL(L)+SFULL(L+1))*0.5
      ENDDO
!
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
!     SB=SBD*DTR
!     WB=WBD*DTR
      TPH0=int_state%TPH0D*DTR
!     TPH0=TPH0D*DTR
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
!!!!!   DPHD=-int_state%SBD*2./REAL(JDE-3)
!!!!!   DLMD=-int_state%WBD*2./REAL(IDE-3)
        DPHD=int_state%DPHD
        DLMD=int_state%DLMD
        DPH=DPHD*DTR
        DLM=DLMD*DTR
        TPH_BASE=SB-DPH-DPH
!
        DO J=J_LO,J_HI
          TLM_BASE=WB-DLM-DLM
          APH=TPH_BASE+(J-JDS+1)*DPH
          DO I=I_LO,I_HI
            ALM=TLM_BASE+(I-IDS+1)*DLM
            IF(ALM> PI) ALM=ALM-PI-PI
            IF(ALM<-PI) ALM=ALM+PI+PI
            int_state%GLAT(I,J)=APH
            int_state%GLON(I,J)=ALM
          ENDDO
        ENDDO
!
      ELSE  ! regional

!!!     DPHD=-int_state%SBD*2./REAL(JDE-1)
!!!     DLMD=-int_state%WBD*2./REAL(IDE-1)
        DPHD=int_state%DPHD
        DLMD=int_state%DLMD
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
!!!   THETA_HALF=ASIN(0.5*SIN(-SB))
      theta_half=0.
      JC=NINT(0.5*(JDE-JDS+1)+THETA_HALF/DPH)
!
!-----------------------------------------------------------------------
!***  SET TIME VARIABLES NEEDED FOR HISTORY OUTPUT.
!-----------------------------------------------------------------------
!
      NSTEPS_PER_HOUR=NINT(3600./int_state%DT)
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
!
!-----------------------------------------------------------------------
!***  THE GFS PHYSICS SUITE IS CONSIDERED A SINGLE PACKAGE HERE.
!----------------------------------------------------------------------
!
      namelist_unit: DO N=101,151
        INQUIRE(N,OPENED=OPENED)
        IF(.NOT.OPENED)THEN
          NLUNIT=N
          EXIT namelist_unit
        ENDIF
      ENDDO namelist_unit

      GFS_PHY_NAMELIST = 'atm_namelist'
      DELTIM           = int_state%DT
      LEVS             = LM
      CALL COMPNS_PHYSICS(DELTIM,    IRET,                   &
                   NTRAC,   NXPT,    NYPT,  JINTMX,          &
                   JCAP,    LEVS,    LEVR,    LONR,   LATR,  &
                   NTOZ,    NTCW,    NCLD,                   &
                   LSOIL,   NMTVR,   NUM_P3D, NUM_P2D,       &
                   THERMODYN_ID,     SFCPRESS_ID,            &
                   NLUNIT,  MYPE,    GFS_PHY_NAMELIST)
!--------------
        IALB            = 0
        DELTIM          = int_state%DT
        LONR            = ITE-ITS+1  ! this is changed in compns_physics from
        LATR            = JTE-JTS+1  ! atm_namelist (restore it back)
        LATS_NODE_R     = JTE-JTS+1
        IPT_LATS_NODE_R = 1
        NFXR            = 27
        LSSAV           = .TRUE.   ! logical flag for store 3-d cloud field
        LDIAG3D         = .FALSE.  ! logical flag for store 3-d diagnostic fields

!--------------
      CALL SET_SOILVEG(MYPE,NLUNIT)
      CALL SET_TRACER_CONST(NTRAC,MYPE,NLUNIT)
      CALL GFUNCPHYS
!
!-----------------------------------------------------------------------
!***  INIT OZONE
!-----------------------------------------------------------------------
     IF( NTOZ .LE. 0 ) THEN        ! Diagnostic ozone
!
!!        rewind (kozc)            ! there is no header in global_o3clim.txt file
!!        read   (kozc,end=101) latsozc, levozc, timeozc, blatc4
!!101     if (levozc .lt. 10 .or. levozc .gt. 100) then
!!          rewind (kozc)
            levozc  = 17
            latsozc = 18
            blatc   = -85.0
!!        else
!!          blatc   = blatc4
!!        endif
          latsozp   = 2
          levozp    = 1
          timeoz    = 1
          pl_coeff  = 1  !!!  0 (MUST set > 0, used in GBPHYS for allocation)
          timeozc   = 12 !!!  this is not in header
!
     ELSE                          ! Prognostic Ozone
!
         rewind (kozpl)
         read   (kozpl) pl_coeff, latsozp, levozp, timeoz
       IF(.NOT.ALLOCATED(pl_lat))THEN
         allocate (pl_lat (latsozp), pl_pres (levozp),pl_time (timeoz+1))
       ENDIF
       IF(.NOT.ALLOCATED(pl_lat4))THEN
         allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1))
       ENDIF
       IF(.NOT.ALLOCATED(tempin)) allocate (tempin(latsozp))
         rewind (kozpl)
         read (kozpl) pl_coeff, latsozp, levozp, timeoz, pl_lat4, pl_pres4, pl_time4
         pl_pres(:) = pl_pres4(:)
         pl_lat(:)  = pl_lat4(:)
         pl_time(:) = pl_time4(:)
!
         DO J=JTS,JTE

           gaul=int_state%GLAT( (ITS+ITE)/2 ,J)*180.0d0/3.14159d0

           int_state%jindx2(j) = latsozp + 1
           do i=1,latsozp
             if (gaul.lt. pl_lat(i)) then
               int_state%jindx2(j) = i
               exit
             endif
           enddo
           int_state%jindx1(j) = max(int_state%jindx2(j)-1,1)
           int_state%jindx2(j) = min(int_state%jindx2(j),latsozp)

           if (int_state%jindx2(j) .ne. int_state%jindx1(j)) then
             int_state%ddy(j) = (gaul                        - pl_lat(int_state%jindx1(j)))  &
                              / (pl_lat(int_state%jindx2(j)) - pl_lat(int_state%jindx1(j)))
           else
             int_state%ddy(j) = 1.0
           endif

         ENDDO
!
         DO I=1,TIMEOZ
           DO N=1,PL_COEFF
             DO K=1,LEVOZP
               READ(KOZPL) TEMPIN
               int_state%OZPLIN(:,K,N,I) = TEMPIN(:)
             ENDDO
           ENDDO
         ENDDO
!
     ENDIF                          ! Diagnostic/Prognostic Ozone
!
        dphiozc = -(blatc+blatc)/(latsozc-1)
!-----------------------------------------------------------------------
!***  END INIT OZONE
!-----------------------------------------------------------------------
!
       IMJM=(ITE-ITS+1)*(JTE-JTS+1)
       ALLOCATE(SIG1T(IMJM),RLA(IMJM),RLO(IMJM),SLMASK(IMJM),OROG(IMJM)   &
               ,AISFCS(IMJM),SIHFCS(IMJM),SICFCS(IMJM),SITFCS(IMJM)       &
               ,SWDFCS(IMJM),VMNFCS(IMJM),VMXFCS(IMJM),SLPFCS(IMJM)       &
               ,ABSFCS(IMJM),TSFFCS(IMJM),SNOFCS(IMJM),ZORFCS(IMJM)       &
               ,TG3FCS(IMJM),CNPFCS(IMJM),SLIFCS(IMJM),F10MFCS(IMJM)      &
               ,VEGFCS(IMJM),VETFCS(IMJM),SOTFCS(IMJM),CVFCS(IMJM)        &
               ,CVBFCS(IMJM),CVTFCS(IMJM),ALFFC1(IMJM,2),ALBFC1(IMJM,4)   &
               ,SMCFC1(IMJM,NSOIL),STCFC1(IMJM,NSOIL),SLCFC1(IMJM,NSOIL) )

       SIHFCS  = 0.0d0
       SICFCS  = 0.0d0
       SITFCS  = 0.0d0
       SWDFCS  = 0.0d0
       VMNFCS  = 0.0d0
       VMXFCS  = 0.0d0
       SLPFCS  = 0.0d0
       ABSFCS  = 0.0d0
       TSFFCS  = 0.0d0
       SNOFCS  = 0.0d0
       ZORFCS  = 0.0d0
       TG3FCS  = 0.0d0
       CNPFCS  = 0.0d0
       SLIFCS  = 0.0d0
       F10MFCS = 0.0d0
       VEGFCS  = 0.0d0
       VETFCS  = 0.0d0
       SOTFCS  = 0.0d0
       CVFCS   = 0.0d0
       CVBFCS  = 0.0d0
       CVTFCS  = 0.0d0
       ALFFC1  = 0.0d0
       ALBFC1  = 0.0d0
       SMCFC1  = 0.0d0
       STCFC1  = 0.0d0
       SLCFC1  = 0.0d0


       JDAT(1)=IDAT(3)
       JDAT(2)=IDAT(2)
       JDAT(3)=IDAT(1)
       JDAT(4)=0
       JDAT(5)=IHRST
       JDAT(6)=0
       JDAT(7)=0
       JDAT(8)=0
       FHOUR=FLOAT(IHRST)
!
    NIJ=0
    DO J=JTS,JTE
      DO I=ITS,ITE

       NIJ=NIJ+1

       SIG1T(NIJ)    = 0.d0
       RLA(NIJ)      = int_state%GLAT(I,J)*180.0d0/3.14159d0
       RLO(NIJ)      = int_state%GLON(I,J)*180.0d0/3.14159d0
       SLMASK(NIJ)   = 1.0d0-int_state%SM(I,J)
       OROG(NIJ)     = int_state%FIS(I,J)/9.81d0
       AISFCS(NIJ)   = int_state%SICE(I,J)

      ENDDO
    ENDDO
!
      CALL SFCCYCLE(204,(ITE-ITS+1)*(JTE-JTS+1),4,SIG1T,FHCYC              &
     &,             JDAT(1), JDAT(2), JDAT(3), JDAT(5), FHOUR              &
     &,             RLA, RLO, SLMASK, OROG                                 &
     &,             SIHFCS,   SICFCS, SITFCS                               &
     &,             SWDFCS,   SLCFC1                                       &
     &,             VMNFCS,   VMXFCS, SLPFCS, ABSFCS                       &
     &,             TSFFCS,   SNOFCS, ZORFCS, ALBFC1, TG3FCS               &
     &,             CNPFCS,   SMCFC1, STCFC1, SLIFCS, AISFCS, F10MFCS      &
     &,             VEGFCS,   VETFCS, SOTFCS, ALFFC1                       &
     &,             CVFCS,    CVBFCS, CVTFCS, MYPE, NLUNIT, IALB)
!
     NIJ=0
     DO J=JTS,JTE
       DO I=ITS,ITE
        NIJ=NIJ+1

        SIHFCS(NIJ) = int_state%SICE(I,J) * 1.0d0  ! initialize like this
        SICFCS(NIJ) = int_state%SICE(I,J) * 0.9d0  ! initialize like this

        int_state%ZORFCS(I,J)   = ZORFCS(NIJ)
        int_state%SIHFCS(I,J)   = SIHFCS(NIJ)
        int_state%SICFCS(I,J)   = SICFCS(NIJ)
        int_state%SLPFCS(I,J)   = SLPFCS(NIJ)
        int_state%TG3FCS(I,J)   = TG3FCS(NIJ)
        int_state%VEGFCS(I,J)   = VEGFCS(NIJ)
        int_state%VETFCS(I,J)   = VETFCS(NIJ)
        int_state%SOTFCS(I,J)   = SOTFCS(NIJ)

        int_state%ALBFC1(I,J,1) = ALBFC1(NIJ,1)
        int_state%ALBFC1(I,J,2) = ALBFC1(NIJ,2)
        int_state%ALBFC1(I,J,3) = ALBFC1(NIJ,3)
        int_state%ALBFC1(I,J,4) = ALBFC1(NIJ,4)

        int_state%ALFFC1(I,J,1) = ALFFC1(NIJ,1)
        int_state%ALFFC1(I,J,2) = ALFFC1(NIJ,2)

       ENDDO
     ENDDO

     DEALLOCATE(SIG1T,RLA,RLO,SLMASK,OROG     &
               ,AISFCS,SIHFCS,SICFCS,SITFCS   &
               ,SWDFCS,VMNFCS,VMXFCS,SLPFCS   &
               ,ABSFCS,TSFFCS,SNOFCS,ZORFCS   &
               ,TG3FCS,CNPFCS,SLIFCS,F10MFCS  &
               ,VEGFCS,VETFCS,SOTFCS,CVFCS    &
               ,CVBFCS,CVTFCS,ALFFC1,ALBFC1   &
               ,SMCFC1,STCFC1,SLCFC1)
!
!----------------------------------------------------------------------
!***  Set fluxes to zero
!----------------------------------------------------------------------
!
             int_state%ALWIN(I,J)         = 0.
             int_state%ALWOUT(I,J)        = 0.
             int_state%ASWIN(I,J)         = 0.
             int_state%ASWOUT(I,J)        = 0.
             int_state%RLWIN(I,J)         = 0.
             int_state%RADOT(I,J)         = 0.
             int_state%RSWIN(I,J)         = 0.
             int_state%RSWOUT(I,J)        = 0.

             int_state%ALWTOA(I,J)        = 0.
             int_state%ASWTOA(I,J)        = 0.
             int_state%RLWTOA(I,J)        = 0.
             int_state%RSWTOA(I,J)        = 0.

             int_state%SFCSHX(I,J)        = 0.
             int_state%SFCLHX(I,J)        = 0.
             int_state%TWBS(I,J)          = 0.
             int_state%QWBS(I,J)          = 0.

             int_state%BGROFF(I,J)        = 0.
             int_state%SSROFF(I,J)        = 0.
             int_state%ACSNOW(I,J)        = 0.

             int_state%CUPPT(I,J)         = 0.
!
!----------------------------------------------------------------------
!***  END GFS PACKAGE INIT
!----------------------------------------------------------------------
!
      ELSE
!
!----------------------------------------------------------------------
!***  IF NOT SELECTING THE GFS SUITE, EACH OF THE PHYSICS GROUPS IS
!***  TREATED INDIVIDUALLY.
!----------------------------------------------------------------------
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
!           CALL RRTMINIT(int_state%RESTART                            &
!                        ,ALLOWED_TO_READ                              &
!                        ,IDS,IDE,JDS,JDE,1,LM+1                       &
!                        ,IMS,IME,JMS,JME,1,LM+1                       &
!                        ,ITS,ITE,JTS,JTE,1,LM)
            DO K=1,LM
              KFLIP=LM+1-K
              SFULL_FLIP(KFLIP)=SFULL(K+1)
              SMID_FLIP(KFLIP)=SMID(K)
            ENDDO
            SFULL_FLIP(LM+1)=SFULL(1)
!
            GMT=REAL(IHOUR_FCST)

            CALL RRTM_INIT(EMISS,SFULL_FLIP,SMID_FLIP,PT_CB            &
                          ,JULYR,START_MONTH,START_DAY,GMT             &
                          ,CO2TF                                       &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM)
!
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
          CASE ('rrtm')
!           WRITE(0,*)' Already called RRTM_INIT from LONGWAVE'
!!!       CASE ('gsfc')
!!!         CALL GSFC_INIT
          CASE ('dudh')
!!!         CALL SWINIT(SWRAD_SCAT,int_state%RESTART                   &
!!!                    ,ALLOWED_TO_READ                                &
!!!                    ,IDS,IDE,JDS,JDE,1,LM+1                         &
!!!                    ,IMS,IME,JMS,JME,1,LM+1                         &
!!!                    ,ITS,ITE,JTS,JTE,1,LM)
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
            CALL JSFC_INIT(LOWLYR                                      &  !<-- Placeholder (computed in TURBULENCE)
                          ,int_state%USTAR,int_state%Z0                &
                          ,int_state%SM,int_state%SICE                 &
                          ,int_state%IVGTYP,int_state%RESTART          &            
                          ,ALLOWED_TO_READ                             &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
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
!rv
          CASE ('bmj_dev')
            CALL BMJ_INIT_DEV(int_state%CLDEFI,int_state%RESTART       &
                         ,IDS,IDE,JDS,JDE,1,LM+1                       &
                         ,IMS,IME,JMS,JME,1,LM+1                       &
                         ,ITS,ITE,JTS,JTE,1,LM)
!rv
!
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
                             ,ALLOWED_TO_READ                          &
                             ,IDS,IDE,JDS,JDE,1,LM+1                   &
                             ,IMS,IME,JMS,JME,1,LM                     &
                             ,ITS,ITE,JTS,JTE,1,LM)
!
          CASE ('wsm6')
             CALL WSM6INIT(RHOAIR0,RHOWATER,RHOSNOW,CLIQ,CV             &
                          ,ALLOWED_TO_READ )
          CASE ('wsm3')
!            CALL WSM3INIT(RHOAIR0,RHOWATER,RHOSNOW,CLIQ,CV             &
!                         ,ALLOWED_TO_READ )
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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_WATER(CWM,F_ICE,F_RAIN,F_RIMEF                  &
                             ,NUM_WATER,WATER,T                         &
                             ,P_QC,P_QR,P_QS,P_QI,P_QG, IMICRO,SPECADV  &
                             ,NTIMESTEP                                 &
                             ,IDS,IDE,JDS,JDE,LM                        &
                             ,IMS,IME,JMS,JME                           &
                             ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    UPDATE_WATER          UPDATE WATER ARRAY
!   PRGRMMR: FERRIER         ORG: NP22     DATE: 3 AUG 2009
!
! ABSTRACT:
!     UPDATE WATER ARRAY FOR FERRIER MICROPHYSICS
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   2009-08     FERRIER     - Synchronize WATER array with CWM, F_rain, F_ice arrays
!
! USAGE: CALL UPDATE_WATER FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!-----------------------------------------------------------------------
      USE MODULE_CONSTANTS,ONLY : EPSQ,TIW
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!-- INPUT FIELDS:
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE                             &
                           ,NUM_WATER                                   &
                           ,P_QC,P_QR,P_QS,P_QI,P_QG,IMICRO,SPECADV     &
                           ,NTIMESTEP 
     ! REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: CWM,T          &
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM,T          &
                                                  ,F_ICE,F_RAIN,F_RIMEF
     ! REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,NUM_WATER),INTENT(OUT) :: WATER
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,NUM_WATER),INTENT(INOUT) :: WATER
!-----------------------------------------------------------------------
!--  LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER :: I,J,K
      REAL :: LIQW
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      micro_update:  IF (IMICRO <= 0) THEN
!-----------------------------------------------------------------------
!***  UPDATE CONDENSATE FIELDS IN WATER ARRAY FOR FERRIER MICROPHYSICS ONLY
!-----------------------------------------------------------------------
         DO K=1,LM
           DO J=JMS,JME
             DO I=IMS,IME
              !IF (specadv == 0) THEN       ! specadv
              IF (specadv == 0 .or. NTIMESTEP .LE. 1) THEN       ! specadv, 
                                                                 ! assign cloud water to WATER in the initial time step
                IF (CWM(I,J,K)>EPSQ) THEN
                   LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                   WATER(I,J,K,P_QC)=(1.-F_rain(I,J,K))*LIQW
                   WATER(I,J,K,P_QR)=F_rain(I,J,K)*LIQW
                   WATER(I,J,K,P_QS)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                   WATER(I,J,K,P_QC)=0.
                   WATER(I,J,K,P_QR)=0.
                   WATER(I,J,K,P_QS)=0.
                ENDIF
              ELSE                         !specadv
               ! IF species are advected, CWM(i,j,k) is changed.
                  CWM(I,J,K) = WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)+WATER(I,J,K,P_QS)
                     F_ICE(I,J,K)=0.0
                !  IF (CWM(I,J,K) > EPSQ) THEN
                  IF (WATER(I,J,K,P_QS) > EPSQ) THEN
                     F_ICE(I,J,K)=WATER(I,J,K,P_QS)/CWM(I,J,K)
                  ENDIF
                     F_RAIN(I,J,K)=0.
                  IF (WATER(I,J,K,P_QR) > EPSQ) THEN
                     LIQW=WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)
                     F_RAIN(I,J,K)=WATER(I,J,K,P_QR)/LIQW
                  ENDIF

              ENDIF                        !sepcadv
             ENDDO
           ENDDO
         ENDDO
      ELSE IF (IMICRO == 1) then micro_update
!-----------------------------------------------------------------------
!***  UPDATE CWM, F_rain, F_ice, F_RIMEF from WATER array for WSM6 microphysics
!-----------------------------------------------------------------------
         DO K=1,LM
           DO J=JMS,JME
             DO I=IMS,IME
               IF (NTIMESTEP .LE. 1) THEN         !! 5-28-2010 assign WATER at initital time 
                IF (CWM(I,J,K)>EPSQ) THEN
                   LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                   WATER(I,J,K,P_QC)=(1.-F_rain(I,J,K))*LIQW
                   WATER(I,J,K,P_QR)=F_rain(I,J,K)*LIQW
                   WATER(I,J,K,P_QI)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                   WATER(I,J,K,P_QC)=0.
                   WATER(I,J,K,P_QR)=0.
                   WATER(I,J,K,P_QI)=0.
                ENDIF
               ENDIF                             !!  END of WATER 
                CWM(I,J,K)=WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)+WATER(I,J,K,P_QI)  &
                          +WATER(I,J,K,P_QS)+WATER(I,J,K,P_QG)
                IF (CWM(I,J,K) > EPSQ) THEN
                   LIQW=WATER(I,J,K,P_QI)+WATER(I,J,K,P_QS)+WATER(I,J,K,P_QG)
                   F_ICE(I,J,K)=LIQW/CWM(I,J,K)
                ELSE
                   F_ICE(I,J,K)=0.
                ENDIF
                IF (WATER(I,J,K,P_QR) > EPSQ) THEN
                   LIQW=WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)
                   F_RAIN(I,J,K)=WATER(I,J,K,P_QR)/LIQW
                ELSE
                   F_RAIN(I,J,K)=0.
                ENDIF
                IF (WATER(I,J,K,P_QG) > EPSQ) THEN
!-- Crudely update F_RIMEF based on a factor of 5 higher density graupel (500 kg/m**3)
!   over unrimed snow (100 kg/m**3)
                   LIQW=5.*WATER(I,J,K,P_QG)+WATER(I,J,K,P_QS)
                   F_RIMEF(I,J,K)=LIQW/(WATER(I,J,K,P_QS)+WATER(I,J,K,P_QG))
                ELSE
                   F_RIMEF(I,J,K)=1.
                ENDIF
             ENDDO
           ENDDO
         ENDDO
      ELSE IF (IMICRO == 2) then micro_update
!-----------------------------------------------------------------------
!***  UPDATE CWM, F_rain, F_ice, F_RIMEF from WATER array for WSM3 microphysics
!-----------------------------------------------------------------------
         DO K=1,LM
           DO J=JMS,JME
             DO I=IMS,IME
                CWM(I,J,K)=WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)
                F_RIMEF(I,J,K)=1.
                IF (T(I,J,K) >= TIW) THEN
                   F_ICE(I,J,K)=0.
                   IF (WATER(I,J,K,P_QR) > EPSQ) THEN
                      LIQW=WATER(I,J,K,P_QC)+WATER(I,J,K,P_QR)
                      F_RAIN(I,J,K)=WATER(I,J,K,P_QR)/LIQW
                   ELSE
                      F_RAIN(I,J,K)=0.
                   ENDIF
                ELSE
                   F_ICE(I,J,K)=1.
                   F_RAIN(I,J,K)=0.
                ENDIF
             ENDDO
           ENDDO
         ENDDO
      ENDIF  micro_update
!
!----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_WATER
!
!----------------------------------------------------------------------
!######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
