!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.) 
!***  IN MODULE_ATM_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-30  Janjic - Add CONVECTION='none' to OPERATIONAL_PHYSICS.
!               Janjic - Fix lower J limit in FFTFHN(WATER).
!   2008-08-23  Janjic - General pressure-sigma hybrid
!               Janjic - Consistent nonhydrostatic correction in the
!                        first term of the pressure gradient force
!   2008-09-03  Black  - Added initialization of boundary arrays
!                        for nests.
!   2009-03-12  Black  - Changes for general hybrid coordinate.
!   2009-11     Jovic  - Modified for ownership/import/export specification
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_VARS_STATE
      USE MODULE_DYNAMICS_INTERNAL_STATE                                  !<-- Horizontal loop limits obtained here
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &  
                                   ,MPI_COMM_COMP                       &
                                   ,MYPE_SHARE

      USE MODULE_GET_CONFIG_DYN
!
      USE MODULE_CONTROL,ONLY : TIMEF
!
      USE MODULE_DIAGNOSE,ONLY : FIELD_STATS,TWR,VWR,EXIT
!
      USE MODULE_DYNAMICS_OUTPUT,ONLY: POINT_DYNAMICS_OUTPUT
!
      USE MODULE_CLOCKTIMES,ONLY : adv1_tim,adv2_tim                    &
                                  ,bocoh_tim,bocov_tim                  &
                                  ,cdwdt_tim,cdzdt_tim,consts_tim       &
                                  ,ddamp_tim,dht_tim                    &
                                  ,dyn_init_tim,dyn_run_tim             &
                                  ,exch_dyn_tim                         &
                                  ,fftfhn_tim,fftfwn_tim,hadv2_tim      &
                                  ,hdiff_tim,init_tim,mono_tim          &
                                  ,pdtsdt_tim,pgforce_tim,poavhn_tim    &
                                  ,polehn_tim,polewn_tim                &
                                  ,prefft_tim,presmud_tim               &
                                  ,swaphn_tim,swapwn_tim                &
                                  ,update_dyn_int_state_tim,updatet_tim &
                                  ,vadv2_tim,vsound_tim,vtoa_tim
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
      PUBLIC :: DYN_REGISTER                                            
!
!-----------------------------------------------------------------------
      INCLUDE '../../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PUBLIC :: IM,JM,LM  
      INTEGER(kind=KINT) :: MY_DOMAIN_ID,MYPE,NUM_PES
!
      LOGICAL(kind=KLOG) :: ADVECT_TRACERS                              &  !<-- Flag for advecting tracers
                           ,I_AM_A_NEST                                 &  !<-- Flag indicating if ATM Component is a nest
                           ,OLD_PASSIVE                                 &  !<-- Flag for old passive advection
                           ,OPERATIONAL_PHYSICS                            !<-- Flag to designate use of operational physics suite
!
      CHARACTER(2) :: DOMAIN_CHAR
      CHARACTER(6) :: FMT='(I2.2)'
!
      TYPE(DYNAMICS_INTERNAL_STATE),POINTER :: INT_STATE                   !<-- The Dynamics component internal state pointer.
!
      TYPE(ESMF_Logical) :: INPUT_READY                                 &  !<-- Flag indicating if nest has its own input data
                           ,NEST_FLAG                                      !<-- Flag indicating if ATM Component is a nest
!
!-----------------------------------------------------------------------
!***  FOR DETERMINING CLOCKTIMES OF VARIOUS PIECES OF THE DYNAMICS.
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Dynamics Gridded Component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Dyn register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_REG=ESMF_SUCCESS                                                 !<-- Initialize error signal variable
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Initialize phase 1"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- The gridded component
                                     ,ESMF_SETINIT                      &  !<-- Predefined subroutine type
                                     ,DYN_INITIALIZE_1                  &  !<-- User's subroutineName
                                     ,1                                 &  !<-- phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Initialize phase 2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- The gridded component
                                     ,ESMF_SETINIT                      &  !<-- Predefined subroutine type
                                     ,DYN_INITIALIZE_2                  &  !<-- User's subroutineName
                                     ,2                                 &  !<-- phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- gridcomp
                                     ,ESMF_SETRUN                       &  !<-- subroutineType
                                     ,DYN_RUN                           &  !<-- user's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- gridcomp
                                     ,ESMF_SETFINAL                     &  !<-- subroutineType
                                     ,DYN_FINALIZE                      &  !<-- user's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
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
!       WRITE(0,*)" DYN_REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" DYN_REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_INITIALIZE_1(GRID_COMP                             &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_ATM                             &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  CARRY OUT ALL NECESSARY SETUPS FOR THE MODEL DYNAMICS.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_ATM                       !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL,    INTENT(OUT)  :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: RC
!
      TYPE(WRAP_DYN_INT_STATE) :: WRAP                                     ! <-- This wrap is a derived type which contains
                                                                           !     only a pointer to the internal state.  It is needed
                                                                           !     for using different architectures or compilers.
!
      TYPE(ESMF_Grid) :: GRID                                              !<-- The ESMF Grid
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF Virtual Machine
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
!***  Initialize the Dynamics timers.
!-----------------------------------------------------------------------
!
      adv1_tim=0.
      adv2_tim=0.
      bocoh_tim=0.
      bocov_tim=0.
      cdwdt_tim=0.
      cdzdt_tim=0.
      consts_tim=0.
      ddamp_tim=0.
      dht_tim=0.
      dyn_run_tim=0.
      exch_dyn_tim=0.
      fftfhn_tim=0.
      fftfwn_tim=0.
      hadv2_tim=0.
      hdiff_tim=0.
      init_tim=0.
      mono_tim=0.
      pdtsdt_tim=0.
      pgforce_tim=0.
      poavhn_tim=0.
      polehn_tim=0.
      polewn_tim=0.
      prefft_tim=0.
      presmud_tim=0.
      swaphn_tim=0.
      swapwn_tim=0.
      updatet_tim=0.
      vadv2_tim=0.
      vsound_tim=0.
      vtoa_tim=0.
!
!-----------------------------------------------------------------------
!***  Allocate the Dynamics internal state pointer.
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE,STAT=RC)
!
!-----------------------------------------------------------------------
!***  Attach the internal state to the Dynamics gridded component.
!-----------------------------------------------------------------------
!
      WRAP%INT_STATE=>INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach Dynamics Internal State to the Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GRID_COMP                      &  !<-- The Dynamics gridded component
                                        ,WRAP                           &  !<-- Pointer to the Dynamics internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!***  Insert the local domain starting limits and the halo width into
!***  the Dynamics internal state.
!-----------------------------------------------------------------------
!
      IF(IHALO==JHALO)THEN
        int_state%NHALO=IHALO
      ELSE
        RC_INIT=ESMF_FAILURE
        WRITE(0,*)'Error due to ihalo /= jhalo'
      ENDIF
!
      int_state%ITS=ITS
      int_state%ITE=ITE
      int_state%JTS=JTS
      int_state%JTE=JTE
!
      int_state%IMS=IMS
      int_state%IME=IME
      int_state%JMS=JMS
      int_state%JME=JME
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values from the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_DYN_DIMS(GRID_COMP                                &
                              ,int_state%INPES,int_state%JNPES          &
                              ,LM                                       &
                              ,int_state%NUM_TRACERS_MET                &
                              ,int_state%NUM_TRACERS_CHEM               &
                              ,int_state%MICROPHYSICS                   &
                              ,int_state%LNSH, int_state%LNSV           &
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
      MESSAGE_CHECK="Get VM from the Dynamics Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &  !<-- The Dynamics gridded component
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
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The ESMF virtual machine
                     ,localpet=int_state%MYPE                           &  !<-- My task's rank
                     ,petcount=int_state%NUM_PES                        &  !<-- Total number of MPI tasks
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  int_state%NUM_PES taken from VM is the total number of tasks 
!***  on this domain including Write/Quilt tasks.  We want only the
!***  number of forecast tasks.
!-----------------------------------------------------------------------
!
      int_state%NUM_PES=int_state%INPES*int_state%JNPES
!
      NUM_PES=int_state%NUM_PES                                            !<-- The number of forecast tasks
      MYPE=int_state%MYPE                                                  !<-- The local task ID
!
      MYPE_SHARE=int_state%MYPE  ! This statement passes MYPE to
                                 ! module_DM_PARALLEL using
                                 ! MYPE_SHARE.
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
!***  the Dynamics' composite VARS array.  
!-----------------------------------------------------------------------
!
        CALL SET_INTERNAL_STATE_DYN_1(INT_STATE,LM)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the ESMF Grid from the Dynamics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGet(gridcomp=GRID_COMP                        &  !<-- The Dynamics gridded component
                             ,grid    =GRID                             &  !<-- The ESMF Grid
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Put the allocated pointers of all export variables (they must be
!***  owned) into the Dynamics export state.  
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
!       WRITE(0,*)'DYN INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      dyn_init_tim=(timef()-btim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_INITIALIZE_1
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_INITIALIZE_2(GRID_COMP                             &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK                                 &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Point all unowned internal state variables at allocated memory,
!***  get time information, read input data, compute variety of 
!***  constant quantities.
!-----------------------------------------------------------------------
!
      USE MODULE_CONTROL,ONLY : DT,HYDRO                                &  !  <--
                               ,ICYCLE                                  &  !  <-- Variables
                               ,NPES                                    &  !  <--
!                                                                            
                               ,BOUNDARY_INIT,CONSTS                       !  <-- Subroutines
!
      USE MODULE_DYNAMICS_INIT_READ_BIN,ONLY : DYNAMICS_READ_BINARY
      USE MODULE_DYNAMICS_INIT_READ_NEMSIO,ONLY : DYNAMICS_READ_NEMSIO
!
#ifdef IBM
      USE MODULE_FLTBNDS,ONLY : PREFFT
#else
      USE MODULE_FLTBNDS,ONLY : PREFFT, PRESMUD
#endif

!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK                           !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL,    INTENT(OUT)  :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: IDENOMINATOR_DT,IEND,IERR,INTEGER_DT,JEND   &
                           ,KSE,KSS,L,N,NUMERATOR_DT,RC
!
      LOGICAL(kind=KLOG) :: RUN_LOCAL
!
      CHARACTER(20) :: FIELD_NAME
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE                                  !<-- The Dynamics import state  
      TYPE(ESMF_Grid) :: GRID                                              !<-- The ESMF Grid
      TYPE(ESMF_TimeInterval) :: DT_ESMF                                   !<-- The ESMF fundamental timestep (s)
!
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
!***  at allocated memory from the owned Physics variable seen in the
!***  Dynamics import state.
!***  For owned (allocated) variables that are also owned/exported by
!***  the Physics, transfer the data into them from the import state.
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
      CALL SET_INTERNAL_STATE_DYN_2(INT_STATE,LM)
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values from the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_DYN(GRID_COMP,INT_STATE,RC)                          !<-- User's routine to extract config file information
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
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are needed for the remaining
!***  initialization process.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(int_state%MYPE<int_state%NUM_PES)THEN                  !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  Assign the fundamental timestep retrieved from the clock.
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
!-----------------------------------------------------------------------
!***  Retrieve the domain ID from the Dynamics import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Domain ID from Dynamics Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Dynamics import state
                              ,name ='DOMAIN_ID'                        &  !<-- Name of variable to get from Dynamics import state
                              ,value=MY_DOMAIN_ID                       &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the import state of the Write gridded component
!***  from the Dynamics export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Import State from Dynamics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE                        &  !<-- The Dynamics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Dynamics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract write component import state from Dynamics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Read the input file.
!-----------------------------------------------------------------------
!
        KSS=1        
        KSE=int_state%NUM_TRACERS_MET
!
        btim=timef()
!
        IF(.NOT.int_state%NEMSIO_INPUT)THEN
!
          CALL DYNAMICS_READ_BINARY(int_state%GLOBAL                    &
                   ,KSS,KSE                                             &
                   ,int_state%PDTOP,int_state%PT,int_state%LPT2         &
                   ,int_state%SG1,int_state%DSG1                        &
                   ,int_state%PSG1,int_state%PDSG1                      &
                   ,int_state%SG2,int_state%DSG2,int_state%SGM          &
                   ,int_state%SGML1,int_state%PSGML1,int_state%SGML2    &
                   ,int_state%SBD,int_state%WBD                         &
                   ,int_state%DPHD,int_state%DLMD                       &
                   ,int_state%TPH0D,int_state%TLM0D                     &
                   ,int_state%FIS,int_state%SM,int_state%SICE           &
                   ,int_state%PD,int_state%PDO,int_state%PINT           &
                   ,int_state%U,int_state%V,int_state%Q2,int_state%E2   &
                   ,int_state%T,int_state%Q,int_state%CW                &
                   ,int_state%TP,int_state%UP,int_state%VP              &
                   ,int_state%O3,int_state%DWDT,int_state%W             &
                   ,int_state%OMGALF,int_state%DIV,int_state%Z          &
                   ,int_state%RTOP                                      &
                   ,int_state%TCU,int_state%TCV,int_state%TCT           &
                   ,int_state%TRACERS_PREV                              &
                   ,int_state%INDX_Q,int_state%INDX_CW                  &
                   ,int_state%INDX_O3,int_state%INDX_Q2                 &
                   ,int_state%NTSTI,int_state%NTSTM                     &
                   ,int_state%IHR,int_state%IHRST,int_state%IDAT        &
                   ,RUN_LOCAL,int_state%RESTART                         &
                   ,int_state%NUM_WATER,int_state%WATER                 &
                   ,int_state%NUM_TRACERS_TOTAL,int_state%TRACERS       &
                   ,int_state%P_QV,int_state%P_QC,int_state%P_QR        &
                   ,int_state%P_QI,int_state%P_QS,int_state%P_QG        &
                   ,DT,int_state%NHOURS_FCST                            &
                   ,int_state%LNSV                                      &
                   ,int_state%UBS,int_state%UBN                         &
                   ,int_state%UBW,int_state%UBE                         &
                   ,int_state%VBS,int_state%VBN                         &
                   ,int_state%VBW,int_state%VBE                         &
                   ,MY_DOMAIN_ID                                        &
                   ,int_state%I_PAR_STA,int_state%J_PAR_STA )
!
        ELSE
!
          CALL DYNAMICS_READ_NEMSIO(int_state%GLOBAL                    &
                  ,KSS,KSE                                              &
                  ,int_state%PDTOP,int_state%PT,int_state%LPT2          &
                  ,int_state%SG1,int_state%DSG1                         &
                  ,int_state%PSG1,int_state%PDSG1                       &
                  ,int_state%SG2,int_state%DSG2,int_state%SGM           &
                  ,int_state%SGML1,int_state%PSGML1,int_state%SGML2     &
                  ,int_state%SBD,int_state%WBD                          &
                  ,int_state%DPHD,int_state%DLMD                        &
                  ,int_state%TPH0D,int_state%TLM0D                      &
                  ,int_state%FIS,int_state%SM,int_state%SICE            &
                  ,int_state%PD,int_state%PDO,int_state%PINT            &
                  ,int_state%U,int_state%V,int_state%Q2,int_state%E2    &
                  ,int_state%T,int_state%Q,int_state%CW                 &
                  ,int_state%TP,int_state%UP,int_state%VP               &
                  ,int_state%O3,int_state%DWDT,int_state%W              &
                  ,int_state%OMGALF,int_state%DIV,int_state%Z           &
                  ,int_state%RTOP                                       &
                  ,int_state%TCU,int_state%TCV,int_state%TCT            &
                  ,int_state%TRACERS_PREV                               &
                  ,int_state%INDX_Q,int_state%INDX_CW                   &
                  ,int_state%INDX_O3,int_state%INDX_Q2                  &
                  ,int_state%NTSTI,int_state%NTSTM                      &
                  ,int_state%IHR,int_state%IHRST,int_state%IDAT         &
                  ,RUN_LOCAL,int_state%RESTART                          &
                  ,int_state%NUM_WATER,int_state%WATER                  &
                  ,int_state%NUM_TRACERS_TOTAL,int_state%TRACERS        &
                  ,int_state%P_QV,int_state%P_QC,int_state%P_QR         &
                  ,int_state%P_QI,int_state%P_QS,int_state%P_QG         &
                  ,DT,int_state%NHOURS_FCST                             &
                  ,int_state%LNSV                                       &
                  ,int_state%UBS,int_state%UBN                          &
                  ,int_state%UBW,int_state%UBE                          &
                  ,int_state%VBS,int_state%VBN                          &
                  ,int_state%VBW,int_state%VBE                          &
                  ,MY_DOMAIN_ID                                         &
                  ,int_state%I_PAR_STA,int_state%J_PAR_STA )
!
        ENDIF
!
      if (mype==0) then
        write(0,*)'dynamics'
        write(0,*)'ihr,ihrst,lpt2,ntsti,ntstm=',int_state%ihr,int_state%ihrst,int_state%lpt2,int_state%ntsti,int_state%ntstm
        write(0,*)'idat=',int_state%idat(1),int_state%idat(2),int_state%idat(3)
        write(0,*)'dsg1=',minval(int_state%dsg1),maxval(int_state%dsg1)
        write(0,*)'pdsg1=',minval(int_state%pdsg1),maxval(int_state%pdsg1)
        write(0,*)'psgml1=',minval(int_state%psgml1),maxval(int_state%psgml1)
        write(0,*)'sgml1=',minval(int_state%sgml1),maxval(int_state%sgml1)
        write(0,*)'sgml2=',minval(int_state%sgml2),maxval(int_state%sgml2)
        write(0,*)'psg1=',minval(int_state%psg1),maxval(int_state%psg1)
        write(0,*)'sg1=',minval(int_state%sg1),maxval(int_state%sg1)
        write(0,*)'sg2=',minval(int_state%sg2),maxval(int_state%sg2)
        write(0,*)'fis=',minval(int_state%fis),maxval(int_state%fis)
        write(0,*)'pd=',minval(int_state%pd),maxval(int_state%pd)
        write(0,*)'pdo=',minval(int_state%pdo),maxval(int_state%pdo)
        write(0,*)'sice=',minval(int_state%sice),maxval(int_state%sice)
        write(0,*)'sm=',minval(int_state%sm),maxval(int_state%sm)
        write(0,*)'cw=',minval(int_state%cw),maxval(int_state%cw)
        write(0,*)'dwdt=',minval(int_state%dwdt),maxval(int_state%dwdt)
        write(0,*)'q=',minval(int_state%q),maxval(int_state%q)
        write(0,*)'q2=',minval(int_state%q2),maxval(int_state%q2)
        write(0,*)'o3=',minval(int_state%o3),maxval(int_state%o3)
        write(0,*)'omgalf=',minval(int_state%omgalf),maxval(int_state%omgalf)
        write(0,*)'div=',minval(int_state%div),maxval(int_state%div)
        write(0,*)'z=',minval(int_state%z),maxval(int_state%z)
        write(0,*)'rtop=',minval(int_state%rtop),maxval(int_state%rtop)
        write(0,*)'tcu=',minval(int_state%tcu),maxval(int_state%tcu)
        write(0,*)'tcv=',minval(int_state%tcv),maxval(int_state%tcv)
        write(0,*)'tct=',minval(int_state%tct),maxval(int_state%tct)
        write(0,*)'t=',minval(int_state%t),maxval(int_state%t)
        write(0,*)'tp=',minval(int_state%tp),maxval(int_state%tp)
        write(0,*)'u=',minval(int_state%u),maxval(int_state%u)
        write(0,*)'up=',minval(int_state%up),maxval(int_state%up)
        write(0,*)'v=',minval(int_state%v),maxval(int_state%v)
        write(0,*)'vp=',minval(int_state%vp),maxval(int_state%vp)
        write(0,*)'e2=',minval(int_state%e2),maxval(int_state%e2)
        write(0,*)'w=',minval(int_state%w),maxval(int_state%w)
        write(0,*)'pint=',minval(int_state%pint),maxval(int_state%pint)
        write(0,*)'water=',minval(int_state%water),minval(int_state%water)
        write(0,*)'tracers=',minval(int_state%tracers),maxval(int_state%tracers)
!       write(0,*)'sp=',minval(int_state%sp),maxval(int_state%sp)
        write(0,*)'run=',int_state%run 
      endif
!
!-----------------------------------------------------------------------
!***  Check if starting Date/Time in input data file agrees with
!***  the configure file.
!-----------------------------------------------------------------------
!
        IF(.NOT.int_state%RESTART.AND.MYPE==0)THEN
          IF(int_state%START_HOUR /=int_state%IHRST.OR.                 &
             int_state%START_DAY  /=int_state%IDAT(1).OR.               &
             int_state%START_MONTH/=int_state%IDAT(2).OR.               &
             int_state%START_YEAR /=int_state%IDAT(3))THEN
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' DATES IN INPUT AND CONFIGURE FILES DISAGREE!!'
            WRITE(0,*)' INPUT: HOUR=',int_state%IHRST                   &
                      ,       ' DAY=',int_state%IDAT(1)                 &
                      ,     ' MONTH=',int_state%IDAT(2)                 &
                      ,      ' YEAR=',int_state%IDAT(3)
            WRITE(0,*)' CONFIG: HOUR=',int_state%START_HOUR             &
                      ,        ' DAY=',int_state%START_DAY              &
                      ,      ' MONTH=',int_state%START_MONTH            &
                      ,       ' YEAR=',int_state%START_YEAR
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(RUN_LOCAL)THEN
          int_state%RUN=ESMF_TRUE
        ELSE
          int_state%RUN=ESMF_FALSE
        ENDIF
!
        init_tim=init_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Nested domains do not have boundary condition files since the
!***  boundary values come from their parents.  However the boundary
!***  variable arrays need to contain initial values before tendencies
!***  from the parent can be added.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  Retrieve the Nest/Not_A_Nest flag from the Dynamics import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Nest/Not-a-Nest Flag from Dynamics Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Dynamics import state
                              ,name ='I-Am-A-Nest Flag'                 &  !<-- Name of variable to get from Dynamics import state
                              ,value=NEST_FLAG                          &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        I_AM_A_NEST=NEST_FLAG                                              !<-- Convert from ESMF_Logical to LOGICAL

        IF(I_AM_A_NEST)THEN
!
!-----------------------------------------------------------------------
!
          CALL BOUNDARY_INIT(ITS,ITE,JTS,JTE,LM                         &
                            ,IMS,IME,JMS,JME                            &
                            ,IDS,IDE,JDS,JDE                            &
                            ,int_state%LNSH,int_state%LNSV              &
                            ,int_state%PD                               &
                            ,int_state%PDBS,int_state%PDBN              &
                            ,int_state%PDBW,int_state%PDBE              &
                            ,int_state%T                                &
                            ,int_state%TBS,int_state%TBN                &
                            ,int_state%TBW,int_state%TBE                &
                            ,int_state%Q                                &
                            ,int_state%QBS,int_state%QBN                &
                            ,int_state%QBW,int_state%QBE                &
                            ,int_state%CW                               &
                            ,int_state%WBS,int_state%WBN                &
                            ,int_state%WBW,int_state%WBE                &
                            ,int_state%U                                &
                            ,int_state%UBS,int_state%UBN                &
                            ,int_state%UBW,int_state%UBE                &
                            ,int_state%V                                &
                            ,int_state%VBS,int_state%VBN                &
                            ,int_state%VBW,int_state%VBE                &
                            ,int_state%RESTART                          &
                              )
!-----------------------------------------------------------------------
!***  Also we need to retrieve the Parent-Child timestep ratio in order
!***  to know how often to update the boundary tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Parent-Child Time Ratio from Dynamics Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Dynamics import state
                                ,name ='Parent-Child Time Ratio'          &  !<-- Name of variable to get from Dynamics import state
                                ,value=int_state%PARENT_CHILD_TIME_RATIO  &  !<-- Put extracted value here
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Assign grid-related constants after dereferencing needed variables.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL CONSTS(int_state%GLOBAL                                    &
                   ,int_state%SMAG2                                     &
                   ,int_state%CODAMP,int_state%WCOR                     &
                   ,int_state%PT                                        &
                   ,int_state%TPH0D,int_state%TLM0D                     &
                   ,int_state%SBD,int_state%WBD                         &
                   ,int_state%DPHD,int_state%DLMD                       &
                   ,int_state%DXH,int_state%RDXH                        &
                   ,int_state%DXV,int_state%RDXV                        &
                   ,int_state%DYH,int_state%RDYH                        &
                   ,int_state%DYV,int_state%RDYV                        &
                   ,int_state%DDV,int_state%RDDV                        &
                   ,int_state%DDMPU,int_state%DDMPV                     &
                   ,int_state%EF4T,int_state%WPDAR                      &
                   ,int_state%FCP,int_state%FDIV                        &
                   ,int_state%CURV,int_state%F                          &
                   ,int_state%FAD,int_state%FAH                         &
                   ,int_state%DARE,int_state%RARE                       &
                   ,int_state%GLAT,int_state%GLON                       &
                   ,int_state%VLAT,int_state%VLON                       &
                   ,int_state%HDACX,int_state%HDACY                     &
                   ,int_state%HDACVX,int_state%HDACVY                   &
                   ,int_state%LNSH,int_state%LNSAD                      &
                   ,int_state%NBOCO,int_state%TBOCO                     &
                   ,MY_DOMAIN_ID)
!
        consts_tim=consts_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Initialize the FFT filters.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL PREFFT(int_state%DLMD,int_state%DPHD,int_state%SBD,LM      &
                     ,int_state%KHFILT,int_state%KVFILT                   &
                     ,int_state%HFILT,int_state%VFILT                     &
#ifdef IBM
                     ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3  &
                     ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3  &
#else
                     ,int_state%WFFTRH,int_state%NFFTRH                   &
                     ,int_state%WFFTRW,int_state%NFFTRW                   &
#endif
                     ,int_state%INPES,int_state%JNPES)
!
          prefft_tim=prefft_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
#ifndef IBM
          btim=timef()
!
!-----------------------------------------------------------------------
!***  Initialize the polar filter for unfiltered variables.
!-----------------------------------------------------------------------
!
          CALL PRESMUD(int_state%DLMD,int_state%DPHD,int_state%SBD      &
                      ,int_state%NHSMUD)
!
          presmud_tim=presmud_tim+(timef()-btim)
#endif
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve the ESMF Grid then create the ESMF Fields on that Grid
!***  for the Dynamics import/export states.
!***  Send subroutine ALLOC_FIELDS_DYN the entire Dynamics internal
!***  state from which the desired arrays will be extracted for
!***  insertion into the import/export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Retrieve ESMF Grid in Dynamics Initialize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGet(gridcomp=GRID_COMP                        &  !<-- The Dynamics gridded component
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
        MESSAGE_CHECK="Insert NUM_TRACERS into Dynamics Export State"
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
!***  INSERT THIS TASK'S INTEGRATION INDEX LIMITS INTO THE
!***  EXPORT STATE ALONG WITH THE FULL DOMAIN LIMITS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Task Integration Limits to Dynamics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='ITS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%ITS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='ITE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%ITE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='JTS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JTS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='JTE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JTE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='LM'                               &  !<-- The inserted quantity will have this name
                              ,value=int_state%LM                       &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='NHALO'                            &  !<-- The inserted quantity will have this name
                              ,value=int_state%NHALO                    &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='IDS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%IDS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='IDE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%IDE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='JDS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JDS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='JDE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JDE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE DOMAIN'S TOP PRESSURE, THE PRESSURE THICKNESS OF THE
!***  PRESSURE DOMAIN, THE MID-LAYER PRESSURES IN THE PRESSURE DOMAIN
!***  AND THE MID-LAYER SIGMAS IN THE SIGMA DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert PT into Dynamics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='PT'                               &  !<-- The inserted quantity will have this name
                              ,value=int_state%PT                       &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='PDTOP'                            &  !<-- The inserted quantity will have this name
                              ,value=int_state%PDTOP                    &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Dynamics export state
                              ,name     ='PSGML1'                       &  !<-- The inserted quantity will have this name
                              ,count    =LM                             &  !<-- The data has this many items
                              ,valueList=int_state%PSGML1               &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Dynamics export state
                              ,name     ='SGML2'                        &  !<-- The inserted quantity will have this name
                              ,count    =LM                             &  !<-- The data has this many items
                              ,valueList=int_state%SGML2                &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Dynamics export state
                              ,name     ='SG1'                          &  !<-- The inserted quantity will have this name
                              ,count    =LM+1                           &  !<-- The data has this many items
                              ,valueList=int_state%SG1                  &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Dynamics export state
                              ,name     ='SG2'                          &  !<-- The inserted quantity will have this name
                              ,count    =LM+1                           &  !<-- The data has this many items
                              ,valueList=int_state%SG2                  &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FINALLY INSERT THE VALUE OF LNSH AND LNSV (THE WIDTH OF THE
!***  BLENDING REGION ALONG THE BOUNDARIES FOR H AND V POINTS).
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert LNSH, LNSV into Dynamics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='LNSH'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%LNSH                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Dynamics export state
                              ,name ='LNSV'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%LNSV                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  LET DYN_RUN KNOW THAT THE FIRST TIMESTEP IS SPECIAL.
!-----------------------------------------------------------------------
!
        int_state%FIRST=.TRUE.
!
!-----------------------------------------------------------------------
!***  EXTRACT ALL FORECAST TASKS' HORIZONTAL SUBDOMAIN LIMITS
!***  FROM THE DYNAMICS IMPORT STATE AND GIVE THEM TO THE
!***  DYNAMICS INTERNAL STATE.
!***  THIS IS NECESSARY IF QUILTING IS SELECTED BECAUSE THESE
!***  LIMITS WILL BE TAKEN FROM THE DYNAMICS/PHYSICS INTERNAL
!***  STATES, PLACED INTO THE WRITE COMPONENTS' IMPORT STATES
!***  AND USED FOR THE COMBINING OF LOCAL DOMAIN DATA ONTO THE
!***  GLOBAL DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Local Domain Limits to Dynamics Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Dynamics import state
                              ,name     ='LOCAL_ISTART'                 &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_ISTART         &  !<-- Insert Attribute into Dynamics internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Dynamics import state
                              ,name     ='LOCAL_IEND'                   &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_IEND           &  !<-- Insert Attribute into Dynamics internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Dynamics import state
                              ,name     ='LOCAL_JSTART'                 &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_JSTART         &  !<-- Insert Attribute into Dynamics internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Dynamics import state
                              ,name     ='LOCAL_JEND'                   &  !<-- Name of the attribute to extract
                              ,count    =NUM_PES                        &  !<-- # of items in attribute
                              ,valueList=int_state%LOCAL_JEND           &  !<-- Insert Attribute into Dynamics internal state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE RESTART OUTPUT FILE MUST CONTAIN THE WINDS FROM THE
!***  BOUNDARY ARRAYS UBS,UBN,.... BUT:
!***   (1) THEY MUST BE PASSED FROM THE DYNAMICS TO THE WRITE
!***       COMPONENT AND SINCE THEY ARE NOT ON THE ESMF GRID
!***       THEY MUST BE PASSED AS 1-D Attributes.
!***   (2) WE DO NOT WANT TO WASTE CLOCKTIME INSERTING THESE
!***       BC WINDS INTO THE 1-D ARRAYS EVERY TIMESTEP WHEN
!***       THEY ARE ONLY NEEDED AT RESTART OUTPUT TIMES
!***       SO WE MUST INFORM THE DYNAMICS WHEN TO FILL THOSE
!***       ARRAYS.
!
!***  THE 1-D ARRAYS ARE PLACED INTO THE Write Component's
!***  IMPORT STATE IN POINT_DYNAMICS_OUTPUT.  THEY ARE UNLOADED
!***  IN WRT_RUN AND SENT TO THE LEAD FORECAST TASK TO ASSEMBLE
!***  INTO A FULL-DOMAIN 1-D DATASTRING THAT CAN BE SENT TO THE
!***  LEAD WRITE TASK FOR INSERTION INTO THE RESTART FILE.
!-----------------------------------------------------------------------
!
        int_state%NSTEPS_BC_RESTART=NINT((int_state%NHOURS_RESTART*3600)  &  !<-- Timestep frequency for BC data insertion into
                                         /int_state%DT)                      !    1-D local datastrings
!
        IEND=MIN(ITE,IDE-1)
        JEND=MIN(JTE,JDE-1)
!
!       IF(JTS==1)THEN                                                       !<-- South boundary tasks
          int_state%NUM_WORDS_BC_SOUTH=2*2*int_state%LNSV*LM*(IEND-ITS+1)
          ALLOCATE(int_state%RST_BC_DATA_SOUTH(1:int_state%NUM_WORDS_BC_SOUTH))
          DO N=1,int_state%NUM_WORDS_BC_SOUTH
            int_state%RST_BC_DATA_SOUTH(N)=0.
          ENDDO
!       ENDIF
!
!       IF(JTE==JM)THEN                                                      !<-- North boundary tasks
          int_state%NUM_WORDS_BC_NORTH=2*2*int_state%LNSV*LM*(IEND-ITS+1)
          ALLOCATE(int_state%RST_BC_DATA_NORTH(1:int_state%NUM_WORDS_BC_NORTH))
          DO N=1,int_state%NUM_WORDS_BC_NORTH
            int_state%RST_BC_DATA_NORTH(N)=0.
          ENDDO
!       ENDIF
!
!       IF(ITS==1)THEN                                                       !<-- West boundary tasks
          int_state%NUM_WORDS_BC_WEST=2*2*int_state%LNSV*LM*(JEND-JTS+1)
          ALLOCATE(int_state%RST_BC_DATA_WEST(1:int_state%NUM_WORDS_BC_WEST))
          DO N=1,int_state%NUM_WORDS_BC_WEST
            int_state%RST_BC_DATA_WEST(N)=0.
          ENDDO
!       ENDIF
!
!       IF(ITE==IM)THEN                                                      !<-- East boundary tasks
          int_state%NUM_WORDS_BC_EAST=2*2*int_state%LNSV*LM*(JEND-JTS+1)
          ALLOCATE(int_state%RST_BC_DATA_EAST(1:int_state%NUM_WORDS_BC_EAST))
          DO N=1,int_state%NUM_WORDS_BC_EAST
            int_state%RST_BC_DATA_EAST(N)=0.
          ENDDO
!       ENDIF
!
!-----------------------------------------------------------------------
!***  INSERT HISTORY DATA POINTERS INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
        CALL POINT_DYNAMICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!***  SET FLAG FOR THE OPERATIONAL PHYSICS SUITE.
!***  THIS WILL BE USED TO SAVE CLOCKTIME BY SKIPPING
!***  FREQUENT UPDATES OF THE MOIST ARRAY AND INSTEAD
!***  UPDATE IT ONLY WHEN IT IS NEEDED FOR PHYSICS.
!-----------------------------------------------------------------------
!
        OPERATIONAL_PHYSICS=.FALSE.
!
        IF(int_state%SHORTWAVE   =='gfdl' .AND.                         &
           int_state%LONGWAVE    =='gfdl' .AND.                         &
           int_state%SFC_LAYER   =='myj'  .AND.                         &
           int_state%TURBULENCE  =='myj'  .AND.                         &
          (int_state%CONVECTION  =='bmj'  .OR.                          &
           int_state%CONVECTION  =='none').AND.                         &
           int_state%MICROPHYSICS=='fer' ) THEN
!
          OPERATIONAL_PHYSICS=.TRUE.
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  WILL THIS RUN ADVECT TRACERS?
!-----------------------------------------------------------------------
!
        ADVECT_TRACERS=int_state%ADVECT_TRACERS
!
        OLD_PASSIVE   =.NOT.ADVECT_TRACERS                                 !<-- The old scheme and new scheme are mutually exclusive
!
!-----------------------------------------------------------------------
!

        if (mype==0)        CALL ESMF_StatePrint(EXP_STATE)

      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYN INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_INITIALIZE_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_RUN(GRID_COMP                                      &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  THE INTEGRATION OF THE MODEL DYNAMICS IS DONE
!***  THROUGH THIS ROUTINE.
!-----------------------------------------------------------------------
!
      USE MODULE_CONTROL  ,ONLY : E_BDY,N_BDY,S_BDY,W_BDY
      USE MODULE_CONSTANTS,ONLY : CP,G
!
      USE MODULE_DYNAMICS_ROUTINES,ONLY: ADV1,ADV2,AVEQ2                &
                                        ,CDWDT,CDZDT,DDAMP,DHT          &
                                        ,HADV2,HADV2_SCAL,HDIFF         &
                                        ,IUNIT_ADVEC_SUMS               &
                                        ,IUNIT_POLE_SUMS                &
                                        ,MONO,PDTSDT,PGFORCE            &
                                        ,UPDATES,UPDATET,UPDATEUV       &
                                        ,VADV2,VADV2_SCAL,VSOUND,VTOA
!
      USE MODULE_EXCHANGE,ONLY: HALO_EXCH
!
      USE MODULE_FLTBNDS,ONLY: BOCOH,BOCOV,FFTFHN,FFTFUVN               &
                              ,POAVHN,POLEHN,POLEWN,READ_BC             &
                              ,SWAPHN,SWAPWN
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM                       !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IER,INPES,IRTN,ISTAT,J,JNPES              &
                           ,K,KFLIP,KS,L,N,NSTEPS_HISTORY,NTIMESTEP,RC
!
      INTEGER(kind=KINT),SAVE :: HDIFF_ON                               &
                                ,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG          &
                                ,PARENT_CHILD_TIME_RATIO
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      LOGICAL(kind=KLOG) :: READBC
! wang 3-09-2010
      INTEGER(kind=KINT):: KSE1,specadv

!
!------------------------------------------------------------------
!***  The following SAVEs are for dereferenced constant variables.
!------------------------------------------------------------------
!
      INTEGER(kind=KINT),SAVE :: IDTAD,IDTADT,IHRSTBC,KSE,KSS           &
                                ,LNSAD,LNSH,LNSV,LPT2,NBOCO
!
      INTEGER(kind=KINT),DIMENSION(3),SAVE :: IDATBC
!
      REAL(kind=KFPT) :: FICE,FRAIN,QI,QR,QW,SECONDS_TOTAL,WC
!
      REAL(kind=KFPT),SAVE :: DDMPV,DT,DYH,DYV,EF4T,PDTOP,PT            &
                             ,RDYH,RDYV,TBOCO
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: DSG2             &
                                                      ,PDSG1,PSGML1     &
                                                      ,SGML2
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: SG1,SG2
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: CURV             &
                                                      ,DARE,DDMPU,DXV   &
                                                      ,FAD,FAH          &
                                                      ,FCP,FDIV         &
                                                      ,RARE,RDXH,RDXV   &
                                                      ,WPDAR
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE,SAVE :: F,FIS          &
                                                        ,HDACX,HDACY    &
                                                        ,HDACVX,HDACVY  &
                                                        ,SICE,SM
!
      LOGICAL(kind=KLOG),SAVE :: FIRST_PASS=.TRUE.
      LOGICAL(kind=KLOG),SAVE :: GLOBAL,HYDRO,RUNBC,SECADV
      LOGICAL(kind=KLOG)      :: COMPUTE_BC
!
      INTEGER(kind=KINT),SAVE :: N_PRINT_STATS                            !<--- Timesteps between statistics prints
      LOGICAL(kind=KLOG),SAVE :: written=.false.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  The total number of forecast tasks.
!-----------------------------------------------------------------------
!
      INPES=int_state%INPES   !<-- I fcst tasks
      JNPES=int_state%JNPES   !<-- J fcst tasks
      NUM_PES=INPES*JNPES     !<-- # of fcst tasks
!
      MYPE=int_state%MYPE     !<-- the local PE
!
!-----------------------------------------------------------------------
!***  Extract the timestep count from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn_Run Gets Timestep from the ATM Clock"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &  !<-- The ESMF Clock
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- The number of times the clock has been advanced
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
!***  Update the current values inside the Dynamics internal state
!***  with the import state containing data from the Physics.
!***  In the first timestep we do not need to update since
!***  the Physics has not run yet.
!-----------------------------------------------------------------------
!
!d      btim=timef()
      CALL GET_VARS_FROM_STATE(int_state%VARS, int_state%NUM_VARS, IMP_STATE)
!
!d      IF(.NOT.int_state%FIRST)THEN
!d        CALL UPDATE_INTERNAL_STATE_DYN(IMP_STATE,INT_STATE)
!d      ENDIF 
!
!d      update_dyn_int_state_tim=update_dyn_int_state_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Do some work that only needs to be done once at the start of
!***  the Run step:  Dereference some variables and extract the
!***  horizontal diffusion flag.
!-----------------------------------------------------------------------
!
      firstpass: IF(FIRST_PASS)THEN
!
        FIRST_PASS=.FALSE.
!
        DDMPV=int_state%DDMPV
        DT=int_state%DT
        DYH=int_state%DYH
        DYV=int_state%DYV
        EF4T=int_state%EF4T
        GLOBAL=int_state%GLOBAL
        HYDRO=int_state%HYDRO
        IDTAD=int_state%IDTAD
        IDTADT=int_state%IDTADT
        IHRSTBC=int_state%IHRSTBC
        KSE=int_state%NUM_TRACERS_MET
        KSS=1
        LM=int_state%LM
        LNSAD=int_state%LNSAD
        LNSH=int_state%LNSH
        LNSV=int_state%LNSV
        LPT2=int_state%LPT2
        NBOCO=int_state%NBOCO
        PDTOP=int_state%PDTOP
        PT=int_state%PT
        RDYH=int_state%RDYH
        RDYV=int_state%RDYV
        RUNBC=int_state%RUNBC
        SECADV=int_state%SECADV
        TBOCO=int_state%TBOCO
!
        P_QV=int_state%P_QV
        P_QC=int_state%P_QC
        P_QR=int_state%P_QR
        P_QI=int_state%P_QI
        P_QS=int_state%P_QS
        P_QG=int_state%P_QG
!
        PARENT_CHILD_TIME_RATIO=int_state%PARENT_CHILD_TIME_RATIO
!
        IF(.NOT.ALLOCATED(DSG2))THEN
          ALLOCATE(DSG2(1:LM),STAT=ISTAT)
          ALLOCATE(PDSG1(1:LM),STAT=ISTAT)
          ALLOCATE(PSGML1(1:LM),STAT=ISTAT)
          ALLOCATE(SGML2(1:LM),STAT=ISTAT)
!
          ALLOCATE(SG1(1:LM+1),STAT=ISTAT)
          ALLOCATE(SG2(1:LM+1),STAT=ISTAT)
!
          ALLOCATE(CURV(JDS:JDE),STAT=ISTAT)
          ALLOCATE(DARE(JDS:JDE),STAT=ISTAT)
          ALLOCATE(DDMPU(JDS:JDE),STAT=ISTAT)
          ALLOCATE(DXV(JDS:JDE),STAT=ISTAT)
          ALLOCATE(FAD(JDS:JDE),STAT=ISTAT)
          ALLOCATE(FAH(JDS:JDE),STAT=ISTAT)
          ALLOCATE(FCP(JDS:JDE),STAT=ISTAT)
          ALLOCATE(FDIV(JDS:JDE),STAT=ISTAT)
          ALLOCATE(RARE(JDS:JDE),STAT=ISTAT)
          ALLOCATE(RDXH(JDS:JDE),STAT=ISTAT)
          ALLOCATE(RDXV(JDS:JDE),STAT=ISTAT)
          ALLOCATE(WPDAR(JDS:JDE),STAT=ISTAT)
!
          ALLOCATE(F(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(FIS(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(HDACX(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(HDACY(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(HDACVX(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(HDACVY(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(SICE(IMS:IME,JMS:JME),STAT=ISTAT)
          ALLOCATE(SM(IMS:IME,JMS:JME),STAT=ISTAT)
        ENDIF
!
        DO N=1,3
          IDATBC(N)=int_state%IDATBC(N)
        ENDDO
!
        DO L=1,LM
          DSG2(L)=int_state%DSG2(L)
          PDSG1(L)=int_state%PDSG1(L)
          PSGML1(L)=int_state%PSGML1(L)
          SGML2(L)=int_state%SGML2(L)
        ENDDO
!
        DO L=1,LM+1
          SG1(L)=int_state%SG1(L)
          SG2(L)=int_state%SG2(L)
        ENDDO
!
        DO J=JDS,JDE
          CURV(J)=int_state%CURV(J)
          DARE(J)=int_state%DARE(J)
          DDMPU(J)=int_state%DDMPU(J)
          DXV(J)=int_state%DXV(J)
          FAD(J)=int_state%FAD(J)
          FAH(J)=int_state%FAH(J)
          FCP(J)=int_state%FCP(J)
          FDIV(J)=int_state%FDIV(J)
          RARE(J)=int_state%RARE(J)
          RDXV(J)=int_state%RDXV(J)
          RDXH(J)=int_state%RDXH(J)
          WPDAR(J)=int_state%WPDAR(J)
        ENDDO
!
        DO J=JMS,JME
        DO I=IMS,IME
          F(I,J)=int_state%F(I,J)
          FIS(I,J)=int_state%FIS(I,J)
          HDACX(I,J)=int_state%HDACX(I,J)
          HDACY(I,J)=int_state%HDACY(I,J)
          HDACVX(I,J)=int_state%HDACVX(I,J)
          HDACVY(I,J)=int_state%HDACVY(I,J)
          SICE(I,J)=int_state%SICE(I,J)
          SM(I,J)=int_state%SM(I,J)
        ENDDO
        ENDDO
!
        N_PRINT_STATS=NINT(3600./DT)        !<-- Print layer statistics once per forecast hour
!
!-----------------------------------------------------------------------
!***  Extract the horizontal diffusion flag from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Dyn_Run Extracts Horizontal Diffusion Flag "
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &
                              ,name ='HDIFF'                            &
                              ,value=HDIFF_ON                           &
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF firstpass
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Begin the Dynamics calling sequence.
!***  Note that the first timestep begins differently
!***  than all subsequent timesteps.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      firststep: IF(int_state%FIRST.AND.                                &  !<--  The following block is used only for
                    .NOT.int_state%RESTART)THEN                            !     the first timestep and cold start
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%T,IMS,IME,JMS,JME,LM                              &
           ,INPES)
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%T                                                 &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          polehn_tim=polehn_tim+(timef()-btim)
!
        ENDIF
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  The pressure gradient routine.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (int_state%FIRST,int_state%RESTART,LM,DT,NTIMESTEP            &
          ,RDYV,DSG2,PDSG1,RDXV,WPDAR,FIS                               &
          ,int_state%PD                                                 &
          ,int_state%T,int_state%Q,int_state%CW                         &
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV)
!

        pgforce_tim=pgforce_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Divergence and horizontal pressure advection in thermo eqn
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL DHT                                                        &
          (LM,DYV,DSG2,PDSG1,DXV                                        &
          ,FCP,FDIV                                                     &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%U,int_state%V                                      &
          ,int_state%OMGALF                                             &
          ,int_state%PCNE,int_state%PCNW,int_state%PCX,int_state%PCY    &
          ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY    &
          ,int_state%DIV,int_state%TDIV)

!
        dht_tim=dht_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFHN                                                   &
           (LM                                                          &
           ,int_state%KHFILT                                            &
           ,int_state%HFILT                                             &
           ,int_state%DIV                                               &
#ifdef IBM
           ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3          &
           ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3          &
#else
           ,int_state%WFFTRH,int_state%NFFTRH                           &
#endif
           ,NUM_PES)
          fftfhn_tim=fftfhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
!
          CALL SWAPHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
!
          CALL POLEHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          polehn_tim=polehn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPWN                                                   &
            (int_state%U                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
!
          CALL SWAPWN                                                   &
            (int_state%V                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
          swapwn_tim=swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          polewn_tim=polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH                                                  &
         (int_state%T,LM                                                &
         ,int_state%U,LM                                                &
         ,int_state%V,LM                                                &
         ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF firststep
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      not_firststep: IF(.NOT.int_state%FIRST                            &  !<-- The following block is for all timesteps after
                        .OR.int_state%RESTART)THEN                         !    the first or all steps in restart case
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Horizontal diffusion (internal halo exchange for 4th order)
!-----------------------------------------------------------------------
!
        btim=timef()
!
        IF(HDIFF_ON>0)THEN
          CALL HDIFF                                                    &
            (GLOBAL,HYDRO                                               &
            ,INPES,JNPES,LM,LPT2                                        &
            ,DYH,RDYH                                                   &
            ,DXV,RARE,RDXH                                              &
            ,SICE,SM                                                    &
            ,HDACX,HDACY,HDACVX,HDACVY                                  &
            ,int_state%W,int_state%Z                                    &
            ,int_state%CW,int_state%Q,int_state%Q2                      &
            ,int_state%T,int_state%U,int_state%V)                          
        ENDIF
!
        hdiff_tim=hdiff_tim+(timef()-btim)

!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%T                                                &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q                                                &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%CW                                               &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q2                                               &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          poavhn_tim=poavhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polehn_tim=polehn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          swapwn_tim=swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polewn_tim=polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,1,1)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update the boundary mass points.
!
!***  For non-nested regional domains, read new boundary tendencies
!***  at the appropriate times.
!
!***  If this is a nested domain then unload the new boundary data
!***  from the Dynamics import state and compute the time tendencies.
!-----------------------------------------------------------------------
!
        bc_update: IF(.NOT.GLOBAL)THEN
!
          boundary_tendencies: IF(S_BDY.OR.N_BDY.OR.W_BDY.OR.E_BDY)THEN
!
!----------------------------------------------------------------
!***  Nests update boundary tendencies based on data from parent
!----------------------------------------------------------------
!
            IF(I_AM_A_NEST)THEN
              COMPUTE_BC=(NTIMESTEP==1.OR.                              &
                          MOD(NTIMESTEP,PARENT_CHILD_TIME_RATIO)==0)
!
              IF(COMPUTE_BC)THEN
!     call print_memory()
                CALL UPDATE_BC_TENDS(IMP_STATE                          &
                                    ,LM,LNSH,LNSV                       &
                                    ,PARENT_CHILD_TIME_RATIO,DT         &
                                    ,int_state%PDBS,int_state%PDBN      &
                                    ,int_state%PDBW,int_state%PDBE      &
                                    ,int_state%TBS,int_state%TBN        &
                                    ,int_state%TBW,int_state%TBE        &
                                    ,int_state%QBS,int_state%QBN        &
                                    ,int_state%QBW,int_state%QBE        &
                                    ,int_state%WBS,int_state%WBN        &
                                    ,int_state%WBW,int_state%WBE        &
                                    ,int_state%UBS,int_state%UBN        &
                                    ,int_state%UBW,int_state%UBE        &
                                    ,int_state%VBS,int_state%VBN        &
                                    ,int_state%VBW,int_state%VBE )

!d          if (mype==7) then
!d          DO J=JMS,JME
!d          write(301,*)J,int_state%PDBW(1,J,1),int_state%PDBW(1,J,2)
!d          END DO
!d          write(302,'(5F20.10)')int_state%TBW
!d          call flush_(301)
!d          call flush_(302)
!d          written=.true.
!d          end if

!     call print_memory()
              ENDIF
!
!---------------------------------------------------------------
!***  Single/uppermost domain reads its own boundary input data
!---------------------------------------------------------------
!
            ELSE
              READBC=(NTIMESTEP==1.OR.MOD(NTIMESTEP,NBOCO)==0)
!
              IF(READBC)THEN
                CALL READ_BC(LM,LNSH,LNSV,NTIMESTEP,DT                  &
                            ,RUNBC,IDATBC,IHRSTBC,TBOCO                 &
                            ,int_state%PDBS,int_state%PDBN              &
                            ,int_state%PDBW,int_state%PDBE              &
                            ,int_state%TBS,int_state%TBN                &
                            ,int_state%TBW,int_state%TBE                &
                            ,int_state%QBS,int_state%QBN                &
                            ,int_state%QBW,int_state%QBE                &
                            ,int_state%WBS,int_state%WBN                &
                            ,int_state%WBW,int_state%WBE                &
                            ,int_state%UBS,int_state%UBN                &
                            ,int_state%UBW,int_state%UBE                &
                            ,int_state%VBS,int_state%VBN                &
                            ,int_state%VBW,int_state%VBE                &
                            ,MY_DOMAIN_ID)
              ENDIF
!
            ENDIF
!
          ENDIF boundary_tendencies
!
          btim=timef()
!
          CALL BOCOH                                                    &
            (LM,LNSH,DT,PT,DSG2,PDSG1                                   &
             ,int_state%PD                                              &
             ,int_state%PDBE,int_state%PDBN                             &
             ,int_state%PDBS,int_state%PDBW                             &
             ,int_state%TBE,int_state%TBN                               &
             ,int_state%TBS,int_state%TBW                               &
             ,int_state%QBE,int_state%QBN                               &
             ,int_state%QBS,int_state%QBW                               &
             ,int_state%WBE,int_state%WBN                               &
             ,int_state%WBS,int_state%WBW                               &
             ,int_state%T,int_state%Q,int_state%CW                      &
             ,int_state%PINT)
!
          bocoh_tim=bocoh_tim+(timef()-btim)
!
        ENDIF bc_update
!
!-----------------------------------------------------------------------
!***  The pressure gradient routine.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (int_state%FIRST,int_state%RESTART,LM,DT,NTIMESTEP            &
          ,RDYV,DSG2,PDSG1,RDXV,WPDAR,FIS                               &
          ,int_state%PD                                                 &
          ,int_state%T,int_state%Q,int_state%CW                         &
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV)
!
        pgforce_tim=pgforce_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFUVN                                                  &
            (LM                                                         &
            ,int_state%KVFILT,int_state%VFILT                           &
            ,int_state%TCU,int_state%TCV                                &
#ifdef IBM
            ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3         &
            ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3         &
#else
            ,int_state%WFFTRW,int_state%NFFTRW                          &
#endif
            ,NUM_PES)
          fftfwn_tim=fftfwn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update the wind field.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL UPDATEUV                                                   &
         (LM                                                            &
         ,int_state%U,int_state%V                                       &
         ,int_state%TCU,int_state%TCV)
        updatet_tim=updatet_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          swapwn_tim=swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polewn_tim=polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update the boundary velocity points for the regional forecast.
!-----------------------------------------------------------------------
!
        IF(.NOT.GLOBAL)THEN
!
          btim=timef()
          CALL BOCOV                                                    &
            (LM,LNSV,DT                                                 &
            ,int_state%UBE,int_state%UBN,int_state%UBS,int_state%UBW    &
            ,int_state%VBE,int_state%VBN,int_state%VBS,int_state%VBW    &
            ,int_state%U,int_state%V)
          bocov_tim=bocov_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  The boundary winds have just been updated.  In order to replicate
!***  the integration of a restarted run compared to its free-forecast
!***  counterpart then we must save the wind data in the boundary
!***  arrays for the restart files at this place in the runstream.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP+1,int_state%NSTEPS_BC_RESTART)==0)THEN            !<-- Look ahead to the end of this timestep
          CALL SAVE_BC_DATA                                             &
            (LM,LNSV                                                    &
            ,int_state%UBS,int_state%UBN,int_state%UBW,int_state%UBE    &
            ,int_state%VBS,int_state%VBN,int_state%VBW,int_state%VBE    &
            ,int_state%NUM_WORDS_BC_SOUTH,int_state%RST_BC_DATA_SOUTH   &
            ,int_state%NUM_WORDS_BC_NORTH,int_state%RST_BC_DATA_NORTH   &
            ,int_state%NUM_WORDS_BC_WEST ,int_state%RST_BC_DATA_WEST    &
            ,int_state%NUM_WORDS_BC_EAST ,int_state%RST_BC_DATA_EAST    &
            ,EXP_STATE )
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Divergence and horizontal pressure advection in thermo eqn
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL DHT                                                        &
          (LM,DYV,DSG2,PDSG1,DXV                                        &
          ,FCP,FDIV                                                     &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%U,int_state%V                                      &
          ,int_state%OMGALF                                             &
          ,int_state%PCNE,int_state%PCNW,int_state%PCX,int_state%PCY    &
          ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY    &
          ,int_state%DIV,int_state%TDIV)
!
        dht_tim=dht_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFHN                                                   &
           (LM                                                          &
           ,int_state%KHFILT                                            &
           ,int_state%HFILT                                             &
           ,int_state%DIV                                               &
#ifdef IBM
           ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3          &
           ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3          &
#else
           ,int_state%WFFTRH,int_state%NFFTRH                           &
#endif
           ,NUM_PES)
          fftfhn_tim=fftfhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
!
          CALL SWAPHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
!
          CALL POLEHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          polehn_tim=polehn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,int_state%OMGALF,LM                              &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Divergence damping
!-----------------------------------------------------------------------
!
        btim=timef()
!
        IF(HDIFF_ON>0)THEN
          CALL DDAMP                                                    &
            (LM                                                         &
            ,DDMPV,PDTOP                                                &
            ,DSG2,PDSG1                                                 &
            ,SG1,SG2                                                    &
            ,DDMPU                                                      &
            ,int_state%FREERUN                                          &
            ,int_state%PD,int_state%PDO                                 &
            ,int_state%U,int_state%V                                    &
            ,int_state%DIV)
        ENDIF
!
        ddamp_tim=ddamp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPWN                                                   &
            (int_state%U                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
!
          CALL SWAPWN                                                   &
            (int_state%V                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
          swapwn_tim=swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          polewn_tim=polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%U,int_state%LM                         &
                      ,int_state%V,int_state%LM                         &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF not_firststep
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  The remainder of the Dynamics integration call sequence
!***  is the same for all timesteps.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      int_state%FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  Update the surface pressure.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL PDTSDT                                                       &
        (LM,DT,SG2                                                      &
        ,int_state%PD                                                   &
        ,int_state%PDO,int_state%PSDT                                   &
        ,int_state%PSGDT                                                &
!
!***  Temporary argument
!
       ,int_state%DIV,int_state%TDIV)
!
      pdtsdt_tim=pdtsdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
        btim=timef()
        CALL SWAPHN(int_state%PD,IMS,IME,JMS,JME,1,INPES)
        CALL SWAPHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%PD,IMS,IME,JMS,JME,1,INPES,JNPES)
        CALL POLEHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%PD,1                                     &
                    ,int_state%PSDT,1                                   &
                    ,int_state%PSGDT,LM-1                               &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Advection of T, U, and V
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL ADV1                                                         &
        (GLOBAL,SECADV                                                  &
        ,LM,LNSAD,INPES,JNPES                                           &
        ,DT,DYV,RDYH,RDYV                                               &
        ,DSG2,PDSG1                                                     &
        ,CURV,DXV,FAD,FAH,RDXH,RDXV,F                                   &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%OMGALF,int_state%PSGDT                               &
        ,int_state%T,int_state%U,int_state%V                            &
        ,int_state%TP,int_state%UP,int_state%VP                         &
!
!***  Temporary arguments
!
        ,int_state%PFNE,int_state%PFNW                                  &
        ,int_state%PFX,int_state%PFY                                    &
        ,int_state%TCT,int_state%TCU,int_state%TCV)

      adv1_tim=adv1_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Advection of tracers
!-----------------------------------------------------------------------
! 
      tracers: IF(ADVECT_TRACERS.AND.MOD(ABS(NTIMESTEP),IDTADT)==0)THEN
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
!wang 3-9-2010
        specadv=0
        IF(int_state%SPEC_ADV) specadv = 1   ! advect water 0-no, 1-yes
        IF(specadv == 1 .or. int_state%MICROPHYSICS=='wsm6') then
         kse1 = int_state%num_tracers_total
!!!         write(0,*)'spec adv is on, specadv=',specadv,kse1,kse
        ELSE
         kse1 = kse
        ENDIF

        CALL ADV2                                                       &
          (GLOBAL                                                       &
!          ,IDTADT,KSS,KSE,LM,LNSAD                                      &
          ,IDTADT,KSS,KSE1,LM,LNSAD                                     &
          ,DT,RDYH                                                      &
          ,DSG2,PDSG1                                                   &
          ,FAH,RDXH                                                     &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%PSGDT                                              &
          ,int_state%UP,int_state%VP                                    &
          ,int_state%Q2,int_state%INDX_Q2                               &
          ,int_state%TRACERS                                            &
          ,int_state%TRACERS_PREV                                       &
!
!***  Temporary arguments
!
          ,int_state%PFNE,int_state%PFNW                                &
          ,int_state%PFX,int_state%PFY                                  &
          ,int_state%TRACERS_SQRT                                       &
          ,int_state%TRACERS_TEND)
!
        adv2_tim=adv2_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
          IF(GLOBAL)THEN
!
            btim=timef()
!
 !             DO KS=KSS,KSE
              DO KS=KSS,KSE1
                CALL FFTFHN                                             &
                  (LM                                                   &
                  ,int_state%KHFILT                                     &
                  ,int_state%HFILT                                      &
                  ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,KS)      &
#ifdef IBM
                  ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3   &
                  ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3   &
#else
                  ,int_state%WFFTRH,int_state%NFFTRH                    &
#endif
                  ,NUM_PES)
              ENDDO
! 
            fftfhn_tim=fftfhn_tim+(timef()-btim)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Tracer monotonization
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL MONO                                                       &
!          (IDTADT,KSS,KSE,LM                                            &
          (IDTADT,KSS,KSE1,LM                                           &
          ,DSG2,PDSG1                                                   &
          ,DARE                                                         &
          ,int_state%PD                                                 &
          ,int_state%INDX_Q2                                            &
          ,int_state%TRACERS                                            &
          ,int_state%READ_GLOBAL_SUMS                                   &
          ,int_state%WRITE_GLOBAL_SUMS                                  &
!
!***  Temporary arguments
!
          ,int_state%TRACERS_SQRT                                       &
          ,int_state%TRACERS_TEND)
!
        mono_tim=mono_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update tracers
!-----------------------------------------------------------------------
!
        btim=timef()
!
!---------
!***  Q
!---------
!
        CALL UPDATES                                                    &
          (LM                                                           &
          ,int_state%Q                                                  &
!
!***  Temporary argument
!
          ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,int_state%INDX_Q))
!
!---------
!***  CW
!---------
!
        CALL UPDATES                                                    &
          (LM                                                           &
          ,int_state%CW                                                 &
!
!***  Temporary argument
!
          ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,int_state%INDX_CW))
!
!---------
!***  O3
!---------
!
        CALL UPDATES                                                    &
          (LM                                                           &
          ,int_state%O3                                                 &
!
!***  Temporary argument
!
          ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,int_state%INDX_O3))
!
!---------
!***  Q2
!---------
!
        CALL UPDATES                                                    &
          (LM                                                           &
          ,int_state%Q2                                                 &
!
!***  Temporary argument
!
          ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,int_state%INDX_Q2))
!
! wang 3-9-2010
        IF (specadv == 1) then
          DO KS=KSS,KSE1
           if(KS .ne. int_state%INDX_Q .and.                            &
              KS .ne. int_state%INDX_CW .and.                           &
              KS .ne. int_state%INDX_O3 .and.                           &
              KS .ne. int_state%INDX_Q2 ) then
            CALL UPDATES                                                &
              (LM                                                       &
              ,int_state%TRACERS(IMS:IME,JMS:JME,1:LM,KS)               &
              ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,KS))
           endif
          ENDDO
        ENDIF
! wang

        updatet_tim=updatet_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
! IF specadv == 1 , consider SWAPHN for other species. (in the future for GLOBAL)
!
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
! IF specadv == 1 , consider POLEHN for other species. (in the future for GLOBAL)
          polehn_tim=polehn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%O3,LM                                  &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
!
!wang 3-9-2010
        IF (specadv == 1) then
          DO KS=KSS,KSE1
           if(KS .ne. int_state%INDX_Q .and.                            &
              KS .ne. int_state%INDX_CW .and.                           &
              KS .ne. int_state%INDX_O3 .and.                           &
              KS .ne. int_state%INDX_Q2 ) then
            CALL HALO_EXCH(                                             &
               int_state%TRACERS(IMS:IME,JMS:JME,1:LM,KS),LM            &
              ,2,2)
           endif
          ENDDO
        ENDIF

!wang
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF tracers
!
!-----------------------------------------------------------------------
!***  Interface pressures and horizontal part of Omega-Alpha term
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL VTOA                                                         &
        (LM,DT,EF4T,PT,SG2                                              &
        ,int_state%PSDT                                                 &
        ,int_state%DWDT,int_state%RTOP                                  &
        ,int_state%OMGALF                                               &
        ,int_state%PINT                                                 &
!
!***  Temporary arguments
!
        ,int_state%TDIV,int_state%TCT)
!
      vtoa_tim=vtoa_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%TCT                                                &
#ifdef IBM
          ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3           &
          ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3           &
#else
          ,int_state%WFFTRH,int_state%NFFTRH                            &
#endif
          ,NUM_PES)
        fftfhn_tim=fftfhn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Update the temperature field.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL UPDATET                                                      &
        (LM                                                             &
        ,int_state%T                                                    &
!
!***  Temporary argument
!
        ,int_state%TCT)
!
      updatet_tim=updatet_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL SWAPHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%OMGALF,LM                                &
                    ,int_state%PINT,LM+1                                &
                    ,2,2)
      CALL HALO_EXCH(int_state%T,LM                                     &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Nonhydrostatic advection of height
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL CDZDT                                                        &
        (GLOBAL,HYDRO                                                   &
        ,LM,DT,DSG2,PDSG1,FAH,FIS                                       &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%PSGDT                                                &
        ,int_state%CW,int_state%Q,int_state%RTOP,int_state%T            & 
        ,int_state%PINT                                                 &
        ,int_state%DWDT,int_state%PDWDT,int_state%W                     &
        ,int_state%Z                                                    &
!
!***  temporary arguments
!
        ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY)
!
      cdzdt_tim=cdzdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%W                                                  &
#ifdef IBM
          ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3           &
          ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3           &
#else
          ,int_state%WFFTRH,int_state%NFFTRH                            &
#endif
          ,NUM_PES)
        fftfhn_tim=fftfhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%W,LM                                     &
                    ,3,3)
      exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Advection of W (with internal halo exchange)
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL CDWDT                                                        &
        (GLOBAL,HYDRO,int_state%RESTART                                 &
        ,INPES,JNPES,LM,ABS(NTIMESTEP)                                  &
        ,DT,G,DSG2,PDSG1,FAH                                            &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%PSGDT                                                &
        ,int_state%DWDT,int_state%PDWDT,int_state%W                     &
        ,int_state%PINT                                                 &
!
!***  External scratch areas
!
        ,int_state%PFX,int_state%PFY,int_state%PFNE,int_state%PFNW)
!
      cdwdt_tim=cdwdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%DWDT                                               &
#ifdef IBM
          ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3           &
          ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3           &
#else
          ,int_state%WFFTRH,int_state%NFFTRH                            &
#endif
          ,NUM_PES)
        fftfhn_tim=fftfhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%DWDT,LM                                  &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Vertically propagating fast waves
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL VSOUND                                                       &
        (GLOBAL,HYDRO,int_state%RESTART                                 &
        ,LM,ABS(NTIMESTEP)                                                   &
        ,CP,DT,PT,DSG2,PDSG1                                            &
        ,int_state%PD                                                   &
        ,int_state%CW,int_state%Q,int_state%RTOP                        &
        ,int_state%DWDT,int_state%T,int_state%W                         &
        ,int_state%PINT)
!
      vsound_tim=vsound_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%DWDT                                               &
          ,int_state%READ_GLOBAL_SUMS                                   &
          ,int_state%WRITE_GLOBAL_SUMS)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%W                                                  &
          ,int_state%READ_GLOBAL_SUMS                                   &
          ,int_state%WRITE_GLOBAL_SUMS)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%PINT                                               &
          ,int_state%READ_GLOBAL_SUMS                                   &
          ,int_state%WRITE_GLOBAL_SUMS)
        poavhn_tim=poavhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        swaphn_tim=swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        polehn_tim=polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%DWDT,LM                                  &
                    ,int_state%T,LM                                     &
                    ,2,2)
      CALL HALO_EXCH(int_state%W,LM                                     &
                    ,int_state%PINT,LM+1                                &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      passive_advec: IF(MOD(ABS(NTIMESTEP),IDTAD)==0.AND.OLD_PASSIVE)THEN
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Vertical advection of passive quantities 
!-----------------------------------------------------------------------
!
        btim=timef()
!
        vadv2_micro_check: IF(int_state%MICROPHYSICS=='fer')THEN
          CALL AVEQ2                                                    &
            (LM                                                         &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,int_state%PD                                               &
            ,int_state%Q2,int_state%E2                                  &
            ,1)
!
          CALL VADV2_SCAL                                               &
            (LM,IDTAD                                                   &
            ,DT,DSG2,PDSG1,PSGML1,SGML2                                 &
            ,int_state%PD,int_state%PSGDT                               &
            ,int_state%TRACERS                                          &
            ,int_state%NUM_TRACERS_MET,1,int_state%INDX_Q2)
!
        ELSE vadv2_micro_check
          CALL VADV2_SCAL                                               &
            (LM,IDTAD                                                   &
            ,DT,DSG2,PDSG1,PSGML1,SGML2                                 &
            ,int_state%PD,int_state%PSGDT                               &
            ,int_state%Q2                                               &
            ,1,1,int_state%INDX_Q2)
!
          CALL VADV2_SCAL                                               &
            (LM,IDTAD                                                   &
            ,DT,DSG2,PDSG1,PSGML1,SGML2                                 &
            ,int_state%PD,int_state%PSGDT                               &
            ,int_state%WATER                                            &
            ,int_state%NUM_WATER,2,int_state%INDX_Q2)
!
          DO K=1,LM
          DO J=JTS,JTE
          DO I=ITS,ITE
     !       int_state%Q=int_state%WATER(I,J,K,P_QV)                     &
     !                  /(1.+int_state%WATER(I,J,K,P_QV))
     !       int_state%Q(i,j,k)=int_state%WATER(I,J,K,P_QV)                     &
     !                  /(1.+int_state%WATER(I,J,K,P_QV))
          ENDDO
          ENDDO
          ENDDO
          int_state%Q(:,:,:)=int_state%WATER(:,:,:,P_QV)                     &
                      /(1.+int_state%WATER(:,:,:,P_QV))
           

!
        ENDIF vadv2_micro_check
        vadv2_tim=vadv2_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFHN                                                   &
            (LM                                                         &
            ,int_state%KHFILT                                           &
            ,int_state%HFILT                                            &
            ,int_state%CW                                               &
#ifdef IBM
            ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3         &
            ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3         &
#else
            ,int_state%WFFTRH,int_state%NFFTRH                          &
#endif
            ,NUM_PES)
!
          CALL FFTFHN                                                   &
            (LM                                                         &
            ,int_state%KHFILT                                           &
            ,int_state%HFILT                                            &
            ,int_state%Q                                                &
#ifdef IBM
            ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3         &
            ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3         &
#else
            ,int_state%WFFTRH,int_state%NFFTRH                          &
#endif
            ,NUM_PES)
!
          CALL FFTFHN                                                   &
            (LM                                                         &
            ,int_state%KHFILT                                           &
            ,int_state%HFILT                                            &
            ,int_state%E2                                               &
#ifdef IBM
            ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3         &
            ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3         &
#else
            ,int_state%WFFTRH,int_state%NFFTRH                          &
#endif
            ,NUM_PES)
!
          CALL FFTFHN                                                   &
            (LM                                                         &
            ,int_state%KHFILT                                           &
            ,int_state%HFILT                                            &
            ,int_state%O3                                               &
#ifdef IBM
            ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3         &
            ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3         &
#else
            ,int_state%WFFTRH,int_state%NFFTRH                          &
#endif
            ,NUM_PES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
!
            DO N=2,int_state%NUM_WATER
              CALL FFTFHN                                               &
                (LM                                                     &
                ,int_state%KHFILT                                       &
                ,int_state%HFILT                                        &
                ,int_state%WATER(:,:,:,N)                               &
#ifdef IBM
                ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3     &
                ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3     &
#else
                ,int_state%WFFTRH,int_state%NFFTRH                      &
#endif
                ,NUM_PES)
            ENDDO
!
          ENDIF
!
          fftfhn_tim=fftfhn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL SWAPHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES)
            ENDDO
          ENDIF
!
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL POLEHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES,JNPES)
            ENDDO
          ENDIF
!
          polehn_tim=polehn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%O3,LM                                  &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
!
        CALL HALO_EXCH(int_state%E2,LM                                  &
                      ,1,1)
!
        IF(int_state%MICROPHYSICS/='fer')THEN
          CALL HALO_EXCH(int_state%WATER,LM,int_state%NUM_WATER,2       &
                        ,2,2)
        ENDIF
!
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Horizontal advection of passive quantities (internal halo exchange)
!-----------------------------------------------------------------------
!
        btim=timef()
!
        hadv2_micro_check: IF(int_state%MICROPHYSICS=='fer')THEN
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,INPES,JNPES                                         &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%TRACERS                                          &
            ,int_state%NUM_TRACERS_MET,1,int_state%INDX_Q2              &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          CALL AVEQ2                                                    &
            (LM                                                         &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,int_state%PD                                               &
            ,int_state%Q2,int_state%E2                                  &
            ,2)
!
!***  Update the WATER array.
!***  Remember that WATER is used with the WRF physics and thus
!***  the P_QV slot (=2) is mixing ratio, not specific humidity.
!***  Although WATER is only used for physics in operations, it is
!***  updated here from Q every advection timestep for non-operational
!***  configurations where it may be used outside of the physics.
!
          IF(.NOT.OPERATIONAL_PHYSICS)THEN
              int_state%WATER(:,:,:,P_QV)=int_state%Q(:,:,:)/(1.-int_state%Q(:,:,:))
            DO K=1,LM
            KFLIP=LM+1-K
            DO J=JTS,JTE
            DO I=ITS,ITE
         !     int_state%WATER(I,J,K,P_QV)=int_state%Q(I,J,K)/(1.-int_state%Q(I,J,K))
         !     WC = int_state%CW(I,J,K)
              WC = int_state%CW(I-its+1,J-jts+1,K)
              QI = 0.
              QR = 0.
              QW = 0.
              FICE=int_state%F_ICE(I,J,KFLIP)
              FRAIN=int_state%F_RAIN(I,J,KFLIP)
!
              IF(FICE>=1.)THEN
                QI=WC
              ELSEIF(FICE<=0.)THEN
                QW=WC
              ELSE
                QI=FICE*WC
                QW=WC-QI
              ENDIF
!
              IF(QW>0..AND.FRAIN>0.)THEN
                IF(FRAIN>=1.)THEN
                  QR=QW
                  QW=0.
                ELSE
                  QR=FRAIN*QW
                  QW=QW-QR
                ENDIF
              ENDIF
                
!
  !            int_state%WATER(I,J,K,P_QC)=QW
  !            int_state%WATER(I,J,K,P_QR)=QR
  !            int_state%WATER(I,J,K,P_QI)=0.
  !            int_state%WATER(I,J,K,P_QS)=QI
  !            int_state%WATER(I,J,K,P_QG)=0.
              int_state%WATER(I-its+1,J-jts+1,K,P_QC)=QW
              int_state%WATER(I-its+1,J-jts+1,K,P_QR)=QR
              int_state%WATER(I-its+1,J-jts+1,K,P_QI)=0.
              int_state%WATER(I-its+1,J-jts+1,K,P_QS)=QI
              int_state%WATER(I-its+1,J-jts+1,K,P_QG)=0.
            ENDDO
            ENDDO
            ENDDO
          ENDIF
!
        ELSE hadv2_micro_check
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,INPES,JNPES                                         &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%Q2                                               &
            ,1,1,int_state%INDX_Q2                                      &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,INPES,JNPES                                         &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%WATER                                            &
            ,int_state%NUM_WATER,2,int_state%INDX_Q2                    &
            ,int_state%READ_GLOBAL_SUMS                                 &
            ,int_state%WRITE_GLOBAL_SUMS)
!
        ENDIF hadv2_micro_check
!
        hadv2_tim=hadv2_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL SWAPHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES)
            ENDDO
          ENDIF
!
          swaphn_tim=swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL POLEHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES,JNPES)
            ENDDO
          ENDIF
!
          polehn_tim=polehn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%O3,LM                                  &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
!
        IF(int_state%MICROPHYSICS/='fer')THEN
          CALL HALO_EXCH(int_state%WATER,LM,int_state%NUM_WATER,2       &
                        ,2,2)
        ENDIF
!
        exch_dyn_tim=exch_dyn_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF passive_advec
!
!-----------------------------------------------------------------------
!***  Close the file units used for Reads/Writes of global sums
!***  if the forecast is finished.
!-----------------------------------------------------------------------
!
      IF(ESMF_ClockIsStopTime(clock=CLOCK_ATM,rc=RC))THEN
        IF(int_state%WRITE_GLOBAL_SUMS.AND.MYPE==0)THEN
          CLOSE(IUNIT_ADVEC_SUMS)
          CLOSE(IUNIT_POLE_SUMS)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  NOTE:  The Dynamics export state is fully updated now
!***         because subroutine DYN_INITIALIZE inserted the 
!***         appropriate ESMF Fields into it.  Those Fields 
!***         contain pointers to the actual data and those
!***         pointers are never re-directed, i.e., no explicit
!***         action is needed to update the Dynamics export state.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!     if(mod(ntimestep,20)==0.or.ntimestep<=5)then
!       call twr(int_state%t,lm,'t_dyn',ntimestep,mype,num_pes,mpi_comm_comp &
!               ,ids,ide,jds,jde &
!               ,ims,ime,jms,jme &
!               ,its,ite,jts,jte)
!     endif
!
!-----------------------------------------------------------------------
!***  Write the layer statistics for temperature.
!-----------------------------------------------------------------------
!
      IF(MOD(ABS(NTIMESTEP)+1,N_PRINT_STATS)==0)THEN
!
        CALL FIELD_STATS(INT_STATE%T,MYPE,MPI_COMM_COMP,LM              &
                        ,ITS,ITE,JTS,JTE                                &
                        ,IMS,IME,JMS,JME                                &
                        ,IDS,IDE,JDS,JDE)
      ENDIF
!-----------------------------------------------------------------------
!
      RC=0
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYN RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      dyn_run_tim=dyn_run_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_FINALIZE(GRID_COMP                                 &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_ATM                                 &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE DYNAMICS COMPONENT.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM component's ESMF Clock.
!
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: RC,RC_FINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      MYPE=MYPE_SHARE
!
      IF(MYPE==0)THEN
        WRITE(0,*)' Dynamics Completed Normally.'
      ENDIF
!
!-----------------------------------------------------------------------
!***  DO NOT DEALLOCATE THE DYNAMICS INTERNAL STATE POINTER 
!***  WITHOUT DEALLOCATING ITS CONTENTS.
!-----------------------------------------------------------------------
!
!!!   DEALLOCATE(INT_STATE,stat=RC)
!
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'DYNAMICS FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYNAMICS FINALIZE STEP FAILED'
      ENDIF
!
!     IF(PRESENT(RC_FINALIZE))THEN
        RC_FINALIZE=RC_FINAL
!     ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_FINALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_BC_TENDS(IMP_STATE                              &
                                ,LM,LNSH,LNSV                           &
                                ,PARENT_CHILD_TIME_RATIO,DT             &
                                ,PDBS,PDBN,PDBW,PDBE                    &
                                ,TBS,TBN,TBW,TBE                        &
                                ,QBS,QBN,QBW,QBE                        &
                                ,WBS,WBN,WBW,WBE                        &
                                ,UBS,UBN,UBW,UBE                        &
                                ,VBS,VBN,VBW,VBE )
! 
!-----------------------------------------------------------------------
!***  THIS ROUTINE EXTRACTS BOUNDARY DATA FROM THE DYNAMICS IMPORT
!***  STATE OF NESTED DOMAINS THAT WAS RECEIVED FROM THEIR PARENTS.
!***  THIS DATA IS THEN USED TO UPDATE THE TIME TENDENCIES OF THE
!***  BOUNDARY VARIABLES.  THOSE TENDENCIES ARE VALID THROUGH EACH
!***  TIMESTEP OF THE NESTED DOMAIN'S PARENT.
!***  NOTE THAT THIS DATA WAS FIRST LOADED INTO THE EXPORT STATE OF
!***  THE PARENT-CHILD COUPLER IN SUBROUTINE EXPORT_CHILD_BOUNDARY.
!-----------------------------------------------------------------------
!
      USE MODULE_CONTROL,ONLY : E_BDY,N_BDY,S_BDY,W_BDY
!
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,INTENT(IN) :: LM                                          &  !<-- # of model layers
                           ,LNSH                                        &  !<-- # of boundary blending rows for H points
                           ,LNSV                                        &  !<-- # of boundary blending rows for V points
                           ,PARENT_CHILD_TIME_RATIO                        !<-- # of child timesteps per parent timestep
!
      REAL,INTENT(IN) :: DT                                                !<-- This domain's fundamental timestep
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                          !<-- Dynamics import state
!
      REAL,DIMENSION(IMS:IME,1:LNSH,     1:2),INTENT(INOUT) :: PDBS,PDBN   !<-- South/North PD values/tendencies
!
      REAL,DIMENSION(IMS:IME,1:LNSH,1:LM,1:2),INTENT(INOUT) :: TBS,TBN  &  !<-- South/North temperature values/tendencies
                                                              ,QBS,QBN  &  !<-- South/North specific humidity values/tendencies
                                                              ,WBS,WBN     !<-- South/North cloud condensate values/tendencies
!
      REAL,DIMENSION(IMS:IME,1:LNSV,1:LM,1:2),INTENT(INOUT) :: UBS,UBN  &  !<-- South/North U wind values/tendencies
                                                              ,VBS,VBN     !<-- South/North V wind values/tendencies
!
      REAL,DIMENSION(1:LNSH,JMS:JME,     1:2),INTENT(INOUT) :: PDBW,PDBE   !<-- West/East PD values/tendencies
!
      REAL,DIMENSION(1:LNSH,JMS:JME,1:LM,1:2),INTENT(INOUT) :: TBW,TBE  &  !<-- West/East temperature values/tendencies
                                                              ,QBW,QBE  &  !<-- West/East specific humidity values/tendencies
                                                              ,WBW,WBE     !<-- West/East cloud condensate values/tendencies
!
      REAL,DIMENSION(1:LNSV,JMS:JME,1:LM,1:2),INTENT(INOUT) :: UBW,UBE  &  !<-- West/East U wind values/tendencies
                                                              ,VBW,VBE     !<-- West/East V wind values/tendencies
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER,SAVE :: I1=-1,I2_H=-1,I2_V,J1=-1,J2_H=-1,J2_V
!
      INTEGER,SAVE :: KOUNT_S_H,KOUNT_S_V,KOUNT_N_H,KOUNT_N_V           &
                     ,KOUNT_W_H,KOUNT_W_V,KOUNT_E_H,KOUNT_E_V
!
      INTEGER      :: I,J,K,KOUNT,RC,RC_BCT
!
      REAL,SAVE :: RECIP
!
      REAL,DIMENSION(:),POINTER,SAVE :: BND_DATA_S_H                    &
                                       ,BND_DATA_S_V                    & 
                                       ,BND_DATA_N_H                    & 
                                       ,BND_DATA_N_V                    & 
                                       ,BND_DATA_W_H                    & 
                                       ,BND_DATA_W_V                    & 
                                       ,BND_DATA_E_H                    & 
                                       ,BND_DATA_E_V
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_BCT=ESMF_SUCCESS
!
      IF(I1<0)THEN                                                         !<-- Compute/allocate work variables first time only
!
!-----------------------------------------------------------------------
!***  GRIDPOINT INDEX LIMITS ALONG THE South/North AND West/East
!***  BOUNDARIES FOR MASS (H) AND VELOCITY (V) POINTS.  NOTE THAT
!***  THE BOUNDARY DATA GOES TWO POINTS INTO THE HALO.
!-----------------------------------------------------------------------
!
        I1  =MAX(ITS-2,IDS)
        I2_H=MIN(ITE+2,IDE)
        I2_V=MIN(ITE+2,IDE-1)
        J1  =MAX(JTS-2,JDS)
        J2_H=MIN(JTE+2,JDE)
        J2_V=MIN(JTE+2,JDE-1)
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING 'KOUNT' VARIABLES ARE THE NUMBER OF GRIDPOINTS
!***  ON THE GIVEN TASK SUBDOMAIN'S South/North/West/East BOUNDARIES
!***  FOR ALL QUANTITIES ON MASS AND VELOCITY POINTS.
!-----------------------------------------------------------------------
!
        KOUNT_S_H=(3*LM+1)*(I2_H-I1+1)*LNSH
        KOUNT_N_H=(3*LM+1)*(I2_H-I1+1)*LNSH
        KOUNT_S_V=2*LM*(I2_V-I1+1)*LNSV
        KOUNT_N_V=2*LM*(I2_V-I1+1)*LNSV
        KOUNT_W_H=(3*LM+1)*(J2_H-J1+1)*LNSH
        KOUNT_E_H=(3*LM+1)*(J2_H-J1+1)*LNSH
        KOUNT_W_V=2*LM*(J2_V-J1+1)*LNSV
        KOUNT_E_V=2*LM*(J2_V-J1+1)*LNSV
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE BOUNDARY POINTER ARRAYS INTO WHICH THE BOUBDARY
!***  DATA FROM THE DYNAMICS IMPORT STATE WILL BE UNLOADED.
!-----------------------------------------------------------------------
!
        ALLOCATE(BND_DATA_S_H(1:KOUNT_S_H))
        ALLOCATE(BND_DATA_S_V(1:KOUNT_S_V))
        ALLOCATE(BND_DATA_N_H(1:KOUNT_N_H))
        ALLOCATE(BND_DATA_N_V(1:KOUNT_N_V))
        ALLOCATE(BND_DATA_W_H(1:KOUNT_W_H))
        ALLOCATE(BND_DATA_W_V(1:KOUNT_W_V))
        ALLOCATE(BND_DATA_E_H(1:KOUNT_E_H))
        ALLOCATE(BND_DATA_E_V(1:KOUNT_E_V))
!
        NULLIFY(BND_DATA_S_H)
        NULLIFY(BND_DATA_S_V)
        NULLIFY(BND_DATA_N_H)
        NULLIFY(BND_DATA_N_V)
        NULLIFY(BND_DATA_W_H)
        NULLIFY(BND_DATA_W_V)
        NULLIFY(BND_DATA_E_H)
        NULLIFY(BND_DATA_E_V)
!
        RECIP=1./(DT*PARENT_CHILD_TIME_RATIO)
      ENDIF
!
!-----------------------------------------------------------------------
!***  UNLOAD THE BOUNDARY DATA FROM THE IMPORT STATE AND COMPUTE
!***  THE TIME TENDENCIES FOR THE TIME PERIOD SPANNING THE NUMBER
!***  OF THIS NEST'S TIMESTEPS NEEDED TO REACH THE END OF ITS 
!***  ITS PARENT'S TIMESTEP (FROM WHICH THE DATA WAS SENT).
!-----------------------------------------------------------------------
!
      south: IF(S_BDY)THEN
!
!-------------
!***  South H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South Boundary H Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='SOUTH_H'                      &  !<-- Name of desired boundary data
                              ,count    =KOUNT_S_H                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_S_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO J=1,LNSH
        DO I=I1,I2_H
          KOUNT=KOUNT+1
          PDBS(I,J,2)=(BND_DATA_S_H(KOUNT)-PDBS(I,J,1))*RECIP
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=1,LNSH
        DO I=I1,I2_H
          TBS(I,J,K,2)=(BND_DATA_S_H(KOUNT+1)-TBS(I,J,K,1))*RECIP
          QBS(I,J,K,2)=(BND_DATA_S_H(KOUNT+2)-QBS(I,J,K,1))*RECIP
          WBS(I,J,K,2)=(BND_DATA_S_H(KOUNT+3)-WBS(I,J,K,1))*RECIP
          KOUNT=KOUNT+3
        ENDDO
        ENDDO
        ENDDO
!
!-------------
!***  South V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South Boundary V Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='SOUTH_V'                      &  !<-- Name of desired boundary data
                              ,count    =KOUNT_S_V                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_S_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO K=1,LM
        DO J=1,LNSV
        DO I=I1,I2_V
          UBS(I,J,K,2)=(BND_DATA_S_V(KOUNT+1)-UBS(I,J,K,1))*RECIP
          VBS(I,J,K,2)=(BND_DATA_S_V(KOUNT+2)-VBS(I,J,K,1))*RECIP
          KOUNT=KOUNT+2
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF south
!
!-----------------------------------------------------------------------
!
      north: IF(N_BDY)THEN
!
!-------------
!***  North H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North Boundary H Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='NORTH_H'                      &  !<-- Name of desired boundary data
                              ,count    =KOUNT_N_H                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_N_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO J=1,LNSH
        DO I=I1,I2_H
          KOUNT=KOUNT+1
          PDBN(I,J,2)=(BND_DATA_N_H(KOUNT)-PDBN(I,J,1))*RECIP
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=1,LNSH
        DO I=I1,I2_H
          TBN(I,J,K,2)=(BND_DATA_N_H(KOUNT+1)-TBN(I,J,K,1))*RECIP
          QBN(I,J,K,2)=(BND_DATA_N_H(KOUNT+2)-QBN(I,J,K,1))*RECIP
          WBN(I,J,K,2)=(BND_DATA_N_H(KOUNT+3)-WBN(I,J,K,1))*RECIP
          KOUNT=KOUNT+3
        ENDDO
        ENDDO
        ENDDO
!
!-------------
!***  North V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North Boundary V Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='NORTH_V'                      &  !<-- Name of desired boundary data
                              ,count    =KOUNT_N_V                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_N_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO K=1,LM
        DO J=1,LNSV
        DO I=I1,I2_V
          UBN(I,J,K,2)=(BND_DATA_N_V(KOUNT+1)-UBN(I,J,K,1))*RECIP
          VBN(I,J,K,2)=(BND_DATA_N_V(KOUNT+2)-VBN(I,J,K,1))*RECIP
          KOUNT=KOUNT+2
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF north
!
!-----------------------------------------------------------------------
!
      west: IF(W_BDY)THEN
!
!------------
!***  West H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West Boundary H Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='WEST_H'                       &  !<-- Name of desired boundary data
                              ,count    =KOUNT_W_H                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_W_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO J=J1,J2_H
        DO I=1,LNSH
          KOUNT=KOUNT+1
          PDBW(I,J,2)=(BND_DATA_W_H(KOUNT)-PDBW(I,J,1))*RECIP
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J1,J2_H
        DO I=1,LNSH
          TBW(I,J,K,2)=(BND_DATA_W_H(KOUNT+1)-TBW(I,J,K,1))*RECIP
          QBW(I,J,K,2)=(BND_DATA_W_H(KOUNT+2)-QBW(I,J,K,1))*RECIP
          WBW(I,J,K,2)=(BND_DATA_W_H(KOUNT+3)-WBW(I,J,K,1))*RECIP
          KOUNT=KOUNT+3
        ENDDO
        ENDDO
        ENDDO
!
!------------
!***  West V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West Boundary V Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='WEST_V'                       &  !<-- Name of desired boundary data
                              ,count    =KOUNT_W_V                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_W_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO K=1,LM
        DO J=J1,J2_V
        DO I=1,LNSV
          UBW(I,J,K,2)=(BND_DATA_W_V(KOUNT+1)-UBW(I,J,K,1))*RECIP
          VBW(I,J,K,2)=(BND_DATA_W_V(KOUNT+2)-VBW(I,J,K,1))*RECIP
          KOUNT=KOUNT+2
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF west
!
!-----------------------------------------------------------------------
!
      east: IF(E_BDY)THEN
!
!------------
!***  East H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East Boundary H Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='EAST_H'                       &  !<-- Name of desired boundary data
                              ,count    =KOUNT_E_H                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_E_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO J=J1,J2_H
        DO I=1,LNSH
          KOUNT=KOUNT+1
          PDBE(I,J,2)=(BND_DATA_E_H(KOUNT)-PDBE(I,J,1))*RECIP
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J1,J2_H
        DO I=1,LNSH
          TBE(I,J,K,2)=(BND_DATA_E_H(KOUNT+1)-TBE(I,J,K,1))*RECIP
          QBE(I,J,K,2)=(BND_DATA_E_H(KOUNT+2)-QBE(I,J,K,1))*RECIP
          WBE(I,J,K,2)=(BND_DATA_E_H(KOUNT+3)-WBE(I,J,K,1))*RECIP
          KOUNT=KOUNT+3
        ENDDO
        ENDDO
        ENDDO
!
!------------
!***  East V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East Boundary V Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Dynamics import state
                              ,name     ='EAST_V'                       &  !<-- Name of desired boundary data
                              ,count    =KOUNT_E_V                      &  !<-- # of words in this boundary data
                              ,valueList=BND_DATA_E_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        DO K=1,LM
        DO J=J1,J2_V
        DO I=1,LNSV
          UBE(I,J,K,2)=(BND_DATA_E_V(KOUNT+1)-UBE(I,J,K,1))*RECIP
          VBE(I,J,K,2)=(BND_DATA_E_V(KOUNT+2)-VBE(I,J,K,1))*RECIP
          KOUNT=KOUNT+2
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF east
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_BC_TENDS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SAVE_BC_DATA(LM,LNSV                                   &
                             ,UBS,UBN,UBW,UBE                           &
                             ,VBS,VBN,VBW,VBE                           &
                             ,NUM_WORDS_BC_SOUTH,RST_BC_DATA_SOUTH      &
                             ,NUM_WORDS_BC_NORTH,RST_BC_DATA_NORTH      &
                             ,NUM_WORDS_BC_WEST ,RST_BC_DATA_WEST       &
                             ,NUM_WORDS_BC_EAST ,RST_BC_DATA_EAST       &
                             ,EXP_STATE_DYN                             &
                               )
! 
!-----------------------------------------------------------------------
!***  BOUNDARY ARRAY WINDS ARE NEEDED IN THE RESTART FILE IN ORDER TO
!***  ACHIEVE BIT IDENTICAL ANSWERS BETWEEN RESTARTED RUNS AND THEIR
!***  FREE-FORECAST ANALOGS.  THE BOUNDARY ARRAYS DO NOT SPAN THE
!***  INTEGRATION GRID THUS THEY CAN ONLY BE TRANSMITTED THROUGH
!***  ESMF STATES AS ATTRIBUTES.  NON-SCALAR ATTRIBUTES CAN ONLY
!***  CONTAIN ONE DIMENSION THEREFORE THE BOUNDARY DATA IS MOVED
!***  INTO 1-D ARRAYS IN THIS ROUTINE THEN INSERTED INTO THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
!---------------------
!***  Input Arguments
!---------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: LM                               &  !<-- # of model layers
                                      ,LNSV                             &  !<-- # of boundary blending rows for V points
                                      ,NUM_WORDS_BC_SOUTH               &  !<-- Total # of words in south bndry winds, this fcst task
                                      ,NUM_WORDS_BC_NORTH               &  !<-- Total # of words in north bndry winds, this fcst task
                                      ,NUM_WORDS_BC_WEST                &  !<-- Total # of words in west bndry winds, this fcst task
                                      ,NUM_WORDS_BC_EAST                   !<-- Total # of words in east bndry winds, this fcst task
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,1:LNSV,1:LM,1:2),INTENT(IN) ::  &
                                                               UBS,UBN  &  !<-- South/north boundary U
                                                              ,VBS,VBN     !<-- South/north boundary V
!
      REAL(kind=KFPT),DIMENSION(1:LNSV,JMS:JME,1:LM,1:2),INTENT(IN) ::  &
                                                               UBW,UBE  &  !<-- West/east boundary U
                                                              ,VBW,VBE     !<-- West/east boundary V
!
!---------------------
!***  Inout Arguments
!---------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_DYN                      !<-- The Dynamics export state
!
!----------------------
!***  Output Arguments
!----------------------
!
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_SOUTH),INTENT(OUT) ::    &
                                                     RST_BC_DATA_SOUTH     !<-- All south bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_NORTH),INTENT(OUT) ::    &
                                                     RST_BC_DATA_NORTH     !<-- All north bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_WEST ),INTENT(OUT) ::    &
                                                     RST_BC_DATA_WEST      !<-- All west bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_EAST ),INTENT(OUT) ::    &
                                                     RST_BC_DATA_EAST      !<-- All east bndry wind data on this fcst task
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IB,JB,KOUNT,L,RC,RC_SAVE
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  SOUTHERN BOUNDARY WINDS TO 1-D
!-----------------------------------------------------------------------
!
      IF(JTS==JDS)THEN                                                     !<-- Tasks on south boundary
        KOUNT=0
!
        DO L=1,LM
          DO JB=1,LNSV
            DO IB=ITS,MIN(ITE,IDE-1)
              RST_BC_DATA_SOUTH(KOUNT+1)=UBS(IB,JB,L,1)
              RST_BC_DATA_SOUTH(KOUNT+2)=UBS(IB,JB,L,2)
              RST_BC_DATA_SOUTH(KOUNT+3)=VBS(IB,JB,L,1)
              RST_BC_DATA_SOUTH(KOUNT+4)=VBS(IB,JB,L,2)
              KOUNT=KOUNT+4
            ENDDO
          ENDDO
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Dynamics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Dynamics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC South Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_SOUTH'            &  !<-- Name of 1-D string of south boundary values
                              ,count    =NUM_WORDS_BC_SOUTH             &  !<-- # of south boundary words on this fcst task
                              ,valueList=RST_BC_DATA_SOUTH              &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  NORTHERN BOUNDARY WINDS TO 1-D
!-----------------------------------------------------------------------
!
      IF(JTE==JDE)THEN                                                     !<-- Tasks on north boundary
        KOUNT=0
!
        DO L=1,LM
          DO JB=1,LNSV
            DO IB=ITS,MIN(ITE,IDE-1)
              RST_BC_DATA_NORTH(KOUNT+1)=UBN(IB,JB,L,1)
              RST_BC_DATA_NORTH(KOUNT+2)=UBN(IB,JB,L,2)
              RST_BC_DATA_NORTH(KOUNT+3)=VBN(IB,JB,L,1)
              RST_BC_DATA_NORTH(KOUNT+4)=VBN(IB,JB,L,2)
              KOUNT=KOUNT+4
            ENDDO
          ENDDO
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Dynamics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Dynamics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC North Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_NORTH'            &  !<-- Name of 1-D string of north boundary values
                              ,count    =NUM_WORDS_BC_NORTH             &  !<-- # of north boundary words on this fcst task
                              ,valueList=RST_BC_DATA_NORTH              &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  WESTERN BOUNDARY WINDS TO 1-D
!-----------------------------------------------------------------------
!
      IF(ITS==IDS)THEN                                                     !<-- Tasks on west boundary
        KOUNT=0
!
        DO L=1,LM
          DO JB=JTS,MIN(JTE,JDE-1)
            DO IB=1,LNSV
              RST_BC_DATA_WEST(KOUNT+1)=UBW(IB,JB,L,1)
              RST_BC_DATA_WEST(KOUNT+2)=UBW(IB,JB,L,2)
              RST_BC_DATA_WEST(KOUNT+3)=VBW(IB,JB,L,1)
              RST_BC_DATA_WEST(KOUNT+4)=VBW(IB,JB,L,2)
              KOUNT=KOUNT+4
            ENDDO
          ENDDO
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Dynamics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Dynamics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC West Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_WEST'             &  !<-- Name of 1-D string of west boundary values
                              ,count    =NUM_WORDS_BC_WEST              &  !<-- # of west boundary words on this fcst task
                              ,valueList=RST_BC_DATA_WEST               &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  EASTERN BOUNDARY WINDS TO 1-D
!-----------------------------------------------------------------------
!
      IF(ITE==IDE)THEN                                                     !<-- Tasks on east boundary
        KOUNT=0
!
        DO L=1,LM
          DO JB=JTS,MIN(JTE,JDE-1)
            DO IB=1,LNSV
              RST_BC_DATA_EAST(KOUNT+1)=UBE(IB,JB,L,1)
              RST_BC_DATA_EAST(KOUNT+2)=UBE(IB,JB,L,2)
              RST_BC_DATA_EAST(KOUNT+3)=VBE(IB,JB,L,1)
              RST_BC_DATA_EAST(KOUNT+4)=VBE(IB,JB,L,2)
              KOUNT=KOUNT+4
            ENDDO
          ENDDO
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Dynamics export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Dynamics export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC East Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_EAST'             &  !<-- Name of 1-D string of east boundary values
                              ,count    =NUM_WORDS_BC_EAST              &  !<-- # of east boundary words on this fcst task
                              ,valueList=RST_BC_DATA_EAST               &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SAVE_BC_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
