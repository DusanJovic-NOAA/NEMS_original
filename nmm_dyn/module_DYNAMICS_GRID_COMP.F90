!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_GRID_COMP
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
! HISTORY LOG:
!
!   2008-07-30  Janjic - Add CONVECTION='none' to OPERATIONAL_PHYSICS.
!               Janjic - Fix lower J limit in FFTFHN(WATER).
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_DYNAMICS_INTERNAL_STATE                                  !<-- Horizontal loop limits obtained here
!
      USE MODULE_DYNAMICS_FIELDS,ONLY : ARRAY_T                         &
                                       ,ARRAY_U,ARRAY_V                 &
                                       ,ARRAY_Q2,ARRAY_PD               &
                                       ,ARRAY_OMGALF                    &
                                       ,ARRAY_TRACERS                   &
                                       ,ALLOC_FIELDS_DYN
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &  
                                   ,MPI_COMM_COMP                       &
                                   ,MPI_COMM_INTER_ARRAY                &
                                   ,MYPE_SHARE

      USE MODULE_GET_CONFIG_DYN
!
      USE MODULE_CONTROL,ONLY : TIMEF
!
      USE MODULE_DIAGNOSE,ONLY : FIELD_STATS,TWR,VWR,EXIT
!
      USE MODULE_DYNAMICS_OUTPUT,ONLY: POINT_DYNAMICS_OUTPUT
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
      PUBLIC :: DYN_REGISTER
!
!
!-----------------------------------------------------------------------
      INCLUDE 'kind.inc'
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PUBLIC :: IM,JM,LM
      INTEGER(KIND=KINT)        :: MYPE,NUM_PES
!
      TYPE(INTERNAL_STATE),POINTER :: INT_STATE                           !<-- The Dynamics component internal state pointer.
!
      LOGICAL :: ADVECT_TRACERS                                         & !<-- Flag for advecting tracers
                ,OLD_PASSIVE                                            & !<-- Flag for old passive advection
                ,OPERATIONAL_PHYSICS                                      !<-- Flag to designate use of operational physics suite
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
      SUBROUTINE DYN_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Dynamics gridded component
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
      RC_REG=ESMF_SUCCESS                                                 !<-- Initialize error signal variable
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- The gridded component
                                     ,ESMF_SETINIT                      &  !<-- Predefined subroutine type
                                     ,DYN_INITIALIZE                    &  !<-- User's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
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
      SUBROUTINE DYN_INITIALIZE(GRID_COMP                               &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_ATM                               &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CARRY OUT ALL NECESSARY SETUPS FOR THE MODEL DYNAMICS.
!-----------------------------------------------------------------------
!
!     USE MODULE_ESMF_State
      USE MODULE_CONTROL,ONLY : DLMD,DPHD,DT,GLOBAL,HYDRO               &
                               ,ICYCLE                                  &
                               ,NHOURS_FCST,NPES                        &
                               ,READ_GLOBAL_SUMS,SBD,WBD                &
                               ,WRITE_GLOBAL_SUMS                       &
!
                               ,CONSTS,INIT
!
#ifdef IBM
      USE MODULE_FLTBNDS,ONLY : PREFFT
#else
      USE MODULE_FLTBNDS,ONLY : PREFFT, PRESMUD
#endif
!
!-----------------------------------------------------------------------
      INCLUDE 'mpif.h'
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                       !<-- The Dynamics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_ATM                       !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL,    INTENT(OUT)  :: RC_INIT
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!***  WRAP_INTERNAL_STATE IS DEFINED IN THE INTERNAL STATE MODULE.
!-----------------------------------------------------------------------
!
      TYPE(WRAP_INTERNAL_STATE) :: WRAP   ! This wrap is a derived type which contains
                                          ! only a pointer to the internal state.  It is needed
                                          ! for using different architectures or compilers.
!
      TYPE(ESMF_State)          :: IMP_STATE_WRITE                        !<-- The Dynamics import state  
      TYPE(ESMF_Grid)           :: GRID                                   !<-- The ESMF Grid
      TYPE(ESMF_VM)             :: VM                                     !<-- The ESMF Virtual Machine
      TYPE(ESMF_TimeInterval)   :: DT_ESMF                                !<-- The ESMF fundamental timestep (s)
      TYPE(ESMF_Time)           :: STARTTIME                              !<-- The ESMF start time  
!
      INTEGER :: IDENOMINATOR_DT,IERR,INTEGER_DT,KSE,KSS,L              &
                ,N,NUMERATOR_DT,RC
!
      LOGICAL(KIND=KLOG) :: RUN_LOCAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE DYNAMICS TIMERS.
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
!***  ALLOCATE THE INTERNAL STATE POINTER.
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE,STAT=RC)
!
!-----------------------------------------------------------------------
!***  ATTACH THE INTERNAL STATE TO THE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      WRAP%INT_STATE=>INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach Internal State to Dynamics Component"
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
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE IMPORT STATE OF THE WRITE GRIDDED COMPONENT
!***  FROM THE DYNAMICS EXPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Import State from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =EXP_STATE                          &  !<-- The Dynamics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Dynamics export state
                        ,nestedState=IMP_STATE_WRITE                    &  !<-- Extract write component import state from Dynamics export
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE LOCAL DOMAIN STARTING LIMITS AND THE HALO WIDTH INTO
!***  THE DYNAMICS INTERNAL STATE.
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
!***  USE ESMF UTILITIES TO GET INFORMATION FROM THE CONFIGURATION FILE.
!***  THE FUNCTION IS SIMILAR TO READING A NAMELIST.  THE GET_CONFIG
!***  ROUTINE IS THE USER'S.  IT EXTRACTS VALUES FROM THE CONFIG FILE
!***  AND PLACES THEM IN THE NAMELIST COMPONENTS OF THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
      CALL GET_CONFIG_DYN(GRID_COMP,INT_STATE,RC)                          !<-- User's routine to extract config file information
!
      IM=int_state%IM
      JM=int_state%JM
      LM=int_state%LM
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM TO OBTAIN THE TASK ID AND TOTAL NUMBER OF TASKS
!***  FOR THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM for Dynamics Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &  !<-- The Dynamics gridded component
                           ,vm      =VM                                 &
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs and Number from VM"
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
!***  int_state%NUM_PES TAKEN FROM VM IS THE TOTAL NUMBER OF TASKS 
!***  IN THE RUN INCLUDING I/O TASKS.  ACTUALLY WE WANT JUST THE
!***  NUMBER OF FORECAST TASKS.
!-----------------------------------------------------------------------
!
      int_state%NUM_PES=int_state%INPES*int_state%JNPES
      NUM_PES=int_state%NUM_PES
      MYPE=int_state%MYPE
!
      MYPE_SHARE=int_state%MYPE  ! This statement passes MYPE to
                                 ! module_DM_PARALLEL using
                                 ! MYPE_SHARE.
!
!-----------------------------------------------------------------------
!***  ONLY FORECAST TASKS ARE NEEDED FOR THE REMAINING
!***  INITIALIZATION PROCESS.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<NUM_PES)THEN
!
!-----------------------------------------------------------------------
!***  SET UP THE DYNAMICS INTERNAL STATE VARIABLES 
!***  AND ALLOCATE INTERNAL STATE ARRAYS.
!-----------------------------------------------------------------------
!
        CALL SET_INTERNAL_STATE_DYN(GRID_COMP,INT_STATE)
!
!-----------------------------------------------------------------------
!***  ASSIGN THE FUNDAMENTAL TIMESTEP RETRIEVED FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Fundamental Timestep from ATM's Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &
                          ,timeStep=DT_ESMF                             &
                          ,rc      =RC)
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
        int_state%DT=REAL(INTEGER_DT)+REAL(NUMERATOR_DT)                &
                                     /REAL(IDENOMINATOR_DT)
        DT=int_state%DT
!
!-----------------------------------------------------------------------
!***  PRIMARY INITIALIZATION OF SCALARS/ARRAYS.
!-----------------------------------------------------------------------
!
        KSS=1        
        KSE=int_state%NUM_TRACERS_MET
!
        btim=timef()
!
        CALL INIT(int_state%GLOBAL                                      &
                 ,KSS,KSE                                               &
                 ,int_state%PDTOP,int_state%PT,int_state%LPT2           &
                 ,int_state%SG1,int_state%DSG1                          &
                 ,int_state%PSG1,int_state%PDSG1                        &
                 ,int_state%SG2,int_state%DSG2,int_state%SGM            &
                 ,int_state%SGML1,int_state%PSGML1,int_state%SGML2      &
                 ,int_state%FIS,int_state%SM,int_state%SICE             &
                 ,int_state%PD,int_state%PDO,int_state%PINT             &
                 ,int_state%U,int_state%V,int_state%Q2,int_state%E2     &
                 ,int_state%T,int_state%Q,int_state%CW,int_state%PSGDT  &
                 ,int_state%TP,int_state%UP,int_state%VP                &
                 ,int_state%RRW,int_state%DWDT,int_state%W              &
                 ,int_state%OMGALF,int_state%DIV,int_state%Z            &
                 ,int_state%RTOP                                        &
                 ,int_state%TCU,int_state%TCV,int_state%TCT             &
                 ,int_state%TRACERS_PREV                                &
                 ,int_state%INDX_Q,int_state%INDX_CW                    &
                 ,int_state%INDX_RRW,int_state%INDX_Q2                  &
                 ,int_state%NTSTI,int_state%NTSTM                       &
                 ,int_state%IHR,int_state%IHRST,int_state%IDAT          &
                 ,RUN_LOCAL,int_state%RESTART                           &
                 ,int_state%NUM_WATER,int_state%WATER                   &
                 ,int_state%NUM_TRACERS_TOTAL,int_state%TRACERS         &
                 ,int_state%P_QV,int_state%P_QC,int_state%P_QR          &
                 ,int_state%P_QI,int_state%P_QS,int_state%P_QG)
!
!***  CHECK IF STARTING DATE/TIME IN INPUT DATA FILE AGREES WITH
!***  THE CONFIGURE FILE.
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
        init_tim=init_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  ASSIGN GRID-RELATED CONSTANTS AFTER DEREFERENCING NEEDED VARIABLES.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        DPHD=int_state%DPHD
        DLMD=int_state%DLMD
        SBD=int_state%SBD
        WBD=int_state%WBD
        NHOURS_FCST=int_state%NHOURS_FCST
        GLOBAL=int_state%GLOBAL
        READ_GLOBAL_SUMS=int_state%READ_GLOBAL_SUMS
        WRITE_GLOBAL_SUMS=int_state%WRITE_GLOBAL_SUMS
!
        CALL CONSTS(int_state%GLOBAL                                    &
                   ,int_state%IDTAD,int_state%IDTADT                    &
                   ,int_state%SECDIF                                    &
                   ,int_state%SMAG2,int_state%SMAG4                     &
                   ,int_state%CODAMP,int_state%WCOR                     &
                   ,int_state%PT                                        &
                   ,int_state%TPH0D,int_state%TLM0D                     &
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
                   ,int_state%NBOCO,int_state%TBOCO)
!
        consts_tim=consts_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE FILTERS.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
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
          prefft_tim=prefft_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
#ifndef IBM
          btim=timef()
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE POLAR FILTER FOR UNFILTERED VARIABLES.
!-----------------------------------------------------------------------
!
          CALL PRESMUD(int_state%DLMD,int_state%DPHD,int_state%SBD      &
                      ,int_state%NHSMUD)
!
          presmud_tim=presmud_tim+timef()-btim
#endif
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ESMF Grid THEN CREATE THE ESMF FIELDS ON THAT GRID
!***  FOR THE DYNAMICS IMPORT/EXPORT STATES.
!***  SEND SUBROUTINE ALLOC_FIELDS_DYN THE ENTIRE DYNAMICS INTERNAL STATE
!***  FROM WHICH THE DESIRED ARRAYS WILL BE EXTRACTED FOR
!***  INSERTION INTO THE IMPORT/EXPORT STATES.
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
        CALL ALLOC_FIELDS_DYN(GRID,INT_STATE)
!
!-----------------------------------------------------------------------
!***  ADD THE DESIRED ESMF Arrays TO THE DYNAMICS EXPORT STATE.
!***  THE POINTERS INSIDE THE Arrays ARE POINTING AT THE APPROPRIATE
!***  VARIABLES INSIDE THE INTERNAL STATE (see ALLOC_FIELDS_DYN
!***  in module_DYNAMICS_FIELDS.F).
!-----------------------------------------------------------------------
!
!-----------------------------------------------
!***  ADD THE 3D QUANTITIES TO THE EXPORT STATE.
!-----------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add 3-D Data to Dynamics Export State"
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
!!      CALL ESMF_StateAdd(state=EXP_STATE                              &
!!                        ,array=ARRAY_Q                                &  !<-- Specific humidity
!!                        ,rc   =RC)
!
!!      CALL ESMF_StateAdd(state=EXP_STATE                              &
!!                        ,array=ARRAY_CW                               &  !<-- Cloud condensate
!!                        ,rc   =RC)
!
        CALL ESMF_StateAdd(state=EXP_STATE                              & 
                          ,array=ARRAY_Q2                               &  !<-- TKE
                          ,rc   =RC)
!
        CALL ESMF_StateAdd(state=EXP_STATE                              & 
                          ,array=ARRAY_OMGALF                           &  !<-- Omega-alpha term
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------
!***  ADD THE 2D QUANTITIES TO THE EXPORT STATE.
!-----------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add 2-D Data to Dynamics Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
        MESSAGE_CHECK="Add 4-D Tracer Data to Dynamics Export State"
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
      dyn_init_tim=timef()-btim0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_RUN(GRID_COMP                                      &
                        ,IMP_STATE,EXP_STATE                            &
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
      USE MODULE_DYNAMICS_ROUTINES
!     USE MODULE_DYNAMICS_ROUTINES,ONLY: DHT                            &
!                                       ,IUNIT_ADVEC_SUMS               &
!                                       ,IUNIT_POLE_SUMS                &
!                                       ,PGFORCE,POLEHN_TIM             &
!                                       ,SWAPHN_TIM,WRITE_GLOBAL_SUMS
!
      USE MODULE_EXCHANGE,ONLY: HALO_EXCH
      USE MODULE_FLTBNDS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                     !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                     !<-- The Dynamics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                     !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM                     !<-- The ATM's ESMF Clock
!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_RUN
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: I,IER,INPES,IRTN,ISTAT,J,JNPES              &
                           ,K,KFLIP,KS,L,N,NSTEPS_HISTORY,NTIMESTEP,RC
!
      INTEGER(KIND=KINT),SAVE :: P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      LOGICAL(KIND=KLOG) :: READBC
!
!***  THE FOLLOWING SAVEs ARE FOR DEREFERENCED CONSTANT VARIABLES.
!
      INTEGER(KIND=KINT),SAVE :: IDTAD,IDTADT,IHRSTBC,KSE,KSS           &
                                ,LNSAD,LNSH,LNSV,LPT2,NBOCO
!
      INTEGER(KIND=KINT),DIMENSION(3),SAVE :: IDATBC
!
      REAL(KIND=KFPT) :: FICE,FRAIN,QI,QR,QW,SECONDS_TOTAL,WC
!
      REAL(KIND=KFPT),SAVE :: DDMPV,DT,DYH,DYV,EF4T,PDTOP,PT            &
                             ,RDYH,RDYV,TBOCO
!
      REAL(KIND=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: DSG2             &
                                                      ,PDSG1,PSGML1     &
                                                      ,SGML2
!
      REAL(KIND=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: SG2
!
      REAL(KIND=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: CURV             &
                                                      ,DARE,DDMPU,DXV   &
                                                      ,FAD,FAH          &
                                                      ,FCP,FDIV         &
                                                      ,RARE,RDXH,RDXV   &
                                                      ,WPDAR
!
      REAL(KIND=KFPT),DIMENSION(:,:),ALLOCATABLE,SAVE :: F,FIS          &
                                                        ,HDACX,HDACY    &
                                                        ,HDACVX,HDACVY  &
                                                        ,SICE,SM
!
      LOGICAL(KIND=KLOG),SAVE :: GLOBAL,HYDRO,RUNBC,SECADV,SECDIF
      LOGICAL(KIND=KLOG),SAVE :: FIRST_PASS=.TRUE.
!
      INTEGER(KIND=KINT),SAVE :: N_PRINT_STATS                            !<--- Timesteps between statistics prints
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
!***  THE TOTAL NUMBER OF FORECAST TASKS.
!-----------------------------------------------------------------------
!
      INPES=int_state%INPES   !<-- I fcst tasks
      JNPES=int_state%JNPES   !<-- J fcst tasks
      NUM_PES=INPES*JNPES     !<-- # of fcst tasks
!
      MYPE=int_state%MYPE     !<-- the local PE
!
!-----------------------------------------------------------------------
!***  EXTRACT THE TIMESTEP COUNT FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn_Run Gets ATM Clock"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &  !<-- The ESMF clock
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
!***  UPDATE THE CURRENT VALUES INSIDE THE DYNAMICS INTERNAL STATE
!***  WITH THE IMPORT STATE CONTAINING DATA FROM THE PHYSICS.
!***  IN THE FIRST TIMESTEP WE DO NOT NEED TO UPDATE SINCE
!***  THE PHYSICS HAS NOT RUN YET, I.E., THE IMPORT STATE
!***  IS EMPTY.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      IF(.NOT.int_state%FIRST)THEN
        CALL UPDATE_INTERNAL_STATE_DYN(IMP_STATE,INT_STATE)
      ENDIF 
!
      update_dyn_int_state_tim=update_dyn_int_state_tim+timef()-btim
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE MAIN DYNAMICS INTEGRATION LOOP.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  DEREFERENCE SOME VARIABLES FOR CONVENIENCE.
!-----------------------------------------------------------------------
!
      firstpass: IF(FIRST_PASS)THEN
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
        SECDIF=int_state%SECDIF
        TBOCO=int_state%TBOCO
!
        P_QV=int_state%P_QV
        P_QC=int_state%P_QC
        P_QR=int_state%P_QR
        P_QI=int_state%P_QI
        P_QS=int_state%P_QS
        P_QG=int_state%P_QG
!
        IF(.NOT.ALLOCATED(DSG2))THEN
          ALLOCATE(DSG2(1:LM),STAT=ISTAT)
          ALLOCATE(PDSG1(1:LM),STAT=ISTAT)
          ALLOCATE(PSGML1(1:LM),STAT=ISTAT)
          ALLOCATE(SGML2(1:LM),STAT=ISTAT)
!
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
      ENDIF firstpass
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  BEGIN THE DYNAMICS CALLING SEQUENCE.
!***  NOTE THAT THE FIRST TIMESTEP BEGINS DIFFERENTLY
!***  THAN ALL SUBSEQUENT TIMESTEPS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      firststep: IF(int_state%FIRST.AND.                                &  !<-- The following block is used only for
                    .NOT.int_state%RESTART)THEN                            !    the first timestep and cold start
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
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%T                                                 &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          polehn_tim=polehn_tim+timef()-btim
!
        ENDIF
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  THE PRESSURE GRADIENT ROUTINE.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (NTIMESTEP,int_state%FIRST,int_state%RESTART,LM,DT            &
          ,RDYV,DSG2,PDSG1,RDXV,WPDAR,FIS                               &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%T,int_state%Q,int_state%CW,int_state%DWDT          &
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV)
!
        pgforce_tim=pgforce_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL POAVHN                                                   &
           (IMS,IME,JMS,JME,LM                                          &
           ,int_state%RTOP,0)
          poavhn_tim=poavhn_tim+timef()-btim
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
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  DIVERGENCE AND HORIZONTAL PRESSURE ADVECTION IN THERMO EQN
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
        dht_tim=dht_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
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
          fftfhn_tim=fftfhn_tim+timef()-btim
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
          swaphn_tim=swaphn_tim+timef()-btim
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
          polehn_tim=polehn_tim+timef()-btim
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
          swapwn_tim=swapwn_tim+timef()-btim
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          polewn_tim=polewn_tim+timef()-btim
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
        exch_dyn_tim=exch_dyn_tim+timef()-btim
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
!***  HORIZONTAL DIFFUSION (Internal halo exchange for 4th order)
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL HDIFF                                                      &
          (GLOBAL,HYDRO,SECDIF                                          &
          ,INPES,JNPES,LM,LPT2                                          &
          ,DYH,RDYH                                                     &
          ,DXV,RARE,RDXH                                                &
          ,SICE,SM                                                      &
          ,HDACX,HDACY,HDACVX,HDACVY                                    &
          ,int_state%W,int_state%Z                                      &
          ,int_state%CW,int_state%Q,int_state%Q2                        &
          ,int_state%T,int_state%U,int_state%V)
!
        hdiff_tim=hdiff_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%T,NTIMESTEP)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q,NTIMESTEP)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%CW,NTIMESTEP)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q2,NTIMESTEP)
!
          poavhn_tim=poavhn_tim+timef()-btim
!
          btim=timef()
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polehn_tim=polehn_tim+timef()-btim
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          swapwn_tim=swapwn_tim+timef()-btim
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polewn_tim=polewn_tim+timef()-btim
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
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  READ BOUNDARY CONDITIONS FOR REGIONAL FORECAST
!***  AND UPDATE THE MASS POINTS.
!-----------------------------------------------------------------------
!
        bc_update: IF(.NOT.GLOBAL)THEN
!
          READBC=(NTIMESTEP==1.OR.MOD(NTIMESTEP,NBOCO)==0)
!
          IF(S_BDY.OR.N_BDY.OR.W_BDY.OR.E_BDY)THEN
            IF(READBC)THEN
              CALL READ_BC(LM,LNSH,LNSV,NTIMESTEP,DT                    &
                          ,RUNBC,IDATBC,IHRSTBC,TBOCO                   &
                          ,int_state%PDBS,int_state%PDBN                &
                          ,int_state%PDBW,int_state%PDBE                &
                          ,int_state%TBS,int_state%TBN                  &
                          ,int_state%TBW,int_state%TBE                  &
                          ,int_state%QBS,int_state%QBN                  &
                          ,int_state%QBW,int_state%QBE                  &
                          ,int_state%WBS,int_state%WBN                  &
                          ,int_state%WBW,int_state%WBE                  &
                          ,int_state%UBS,int_state%UBN                  &
                          ,int_state%UBW,int_state%UBE                  &
                          ,int_state%VBS,int_state%VBN                  &
                          ,int_state%VBW,int_state%VBE)
            ENDIF
          ENDIF
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
          bocoh_tim=bocoh_tim+timef()-btim
!
        ENDIF bc_update
!
!-----------------------------------------------------------------------
!***  THE PRESSURE GRADIENT ROUTINE.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (NTIMESTEP,int_state%FIRST,int_state%RESTART,LM,DT            &
          ,RDYV,DSG2,PDSG1,RDXV,WPDAR,FIS                               &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%T,int_state%Q,int_state%CW,int_state%DWDT          &
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV)
!
        pgforce_tim=pgforce_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
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
          fftfwn_tim=fftfwn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  UPDATE THE WIND FIELD.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL UPDATEUV                                                   &
         (LM                                                            &
         ,int_state%U,int_state%V                                       &
         ,int_state%TCU,int_state%TCV)
        updatet_tim=updatet_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%RTOP,NTIMESTEP)
          poavhn_tim=poavhn_tim+timef()-btim
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          swapwn_tim=swapwn_tim+timef()-btim
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          polewn_tim=polewn_tim+timef()-btim
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
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  UPDATE THE BOUNDARY VELOCITY POINTS FOR THE REGIONAL FORECAST.
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
          bocov_tim=bocov_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  DIVERGENCE AND HORIZONTAL PRESSURE ADVECTION IN THERMO EQN
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
        dht_tim=dht_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
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
          fftfhn_tim=fftfhn_tim+timef()-btim
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
          swaphn_tim=swaphn_tim+timef()-btim
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
          polehn_tim=polehn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,int_state%OMGALF,LM                              &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  DIVERGENCE DAMPING
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL DDAMP                                                      &
          (LM                                                           &
          ,DDMPV,PDTOP                                                  &
          ,DSG2,PDSG1                                                   &
          ,DDMPU                                                        &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%U,int_state%V                                      &
          ,int_state%DIV)
!
        ddamp_tim=ddamp_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR THE GLOBAL FORECAST.
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
          swapwn_tim=swapwn_tim+timef()-btim
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          polewn_tim=polewn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%U,int_state%LM                         &
                      ,int_state%V,int_state%LM                         &
                      ,2,2)
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF not_firststep
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE REMAINDER OF THE DYNAMICS INTEGRATION CALL SEQUENCE
!***  IS THE SAME FOR ALL TIMESTEPS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      int_state%FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  UPDATE THE SURFACE PRESSURE.
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
      pdtsdt_tim=pdtsdt_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
        btim=timef()
        CALL SWAPHN(int_state%PD,IMS,IME,JMS,JME,1,INPES)
        CALL SWAPHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%PD,IMS,IME,JMS,JME,1,INPES,JNPES)
        CALL POLEHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
!
        btim=timef()
        CALL SWAPHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%PD,1                                     &
                    ,int_state%PSDT,1                                   &
                    ,int_state%PSGDT,LM-1                               &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  ADVECTION OF T, U, AND V
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
!
      adv1_tim=adv1_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  ADVECTION OF TRACERS
!-----------------------------------------------------------------------
! 
      tracers: IF(ADVECT_TRACERS.AND.MOD(NTIMESTEP,IDTADT)==0)THEN
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL ADV2                                                       &
          (GLOBAL                                                       &
          ,IDTADT,KSS,KSE,LM,LNSAD                                      &
          ,DT,RDYH                                                      &
          ,DSG2,PDSG1                                                   &
          ,FAH,RDXH                                                     &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%PSGDT                                              &
          ,int_state%UP,int_state%VP                                    &
          ,int_state%Q2,int_state%INDX_Q2                               &
          ,int_state%TRACERS,int_state%TRACERS_PREV                     &
!
!***  Temporary arguments
!
          ,int_state%PFNE,int_state%PFNW                                &
          ,int_state%PFX,int_state%PFY                                  &
          ,int_state%TRACERS_SQRT,int_state%TRACERS_TEND)
!
        adv2_tim=adv2_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
!-----------------------------------------------------------------------
!
          IF(GLOBAL)THEN
!
            btim=timef()
!
              DO KS=KSS,KSE
                CALL FFTFHN                                             &
                  (LM                                                   &
                  ,int_state%KHFILT                                     &
                  ,int_state%HFILT                                      &
                  ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,KS)      &
#ifdef IBM
                  ,int_state%CRAUX1,int_state%CRAUX2,int_state%CRAUX3   &
                  ,int_state%RCAUX1,int_state%RCAUX2,int_state%RCAUX3   &
#else
          ,int_state%WFFTRH,int_state%NFFTRH                            &
#endif
                  ,NUM_PES)
              ENDDO
! 
            fftfhn_tim=fftfhn_tim+timef()-btim
!
          ENDIF
!-----------------------------------------------------------------------
!***  TRACER MONOTONIZATION
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL MONO                                                       &
          (IDTADT,KSS,KSE,LM                                            &
          ,DSG2,PDSG1                                                   &
          ,DARE                                                         &
          ,int_state%PD                                                 &
          ,int_state%INDX_Q2                                            &
          ,int_state%TRACERS                                            &
!
!***  Temporary arguments
!
          ,int_state%TRACERS_SQRT,int_state%TRACERS_TEND)
!
        mono_tim=mono_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  UPDATE TRACERS
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
!***  RRW
!---------
!
        CALL UPDATES                                                    &
          (LM                                                           &
          ,int_state%RRW                                                &
!
!***  Temporary argument
!
          ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,int_state%INDX_RRW))
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
        updatet_tim=updatet_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
          polehn_tim=polehn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%RRW,LM                                 &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
!
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF tracers
!
!-----------------------------------------------------------------------
!***  INTERFACE PRESSURES AND HORIZONTAL PART OF OMEGA-ALPHA TERM
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
      vtoa_tim=vtoa_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
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
        fftfhn_tim=fftfhn_tim+timef()-btim
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  UPDATE THE TEMPERATURE FIELD.
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
      updatet_tim=updatet_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL SWAPHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
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
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  NONHYDROSTATIC ADVECTION OF HEIGHT
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
      cdzdt_tim=cdzdt_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
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
        fftfhn_tim=fftfhn_tim+timef()-btim
!
        btim=timef()
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%W,LM                                     &
                    ,3,3)
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  ADVECTION OF W (WITH INTERNAL HALO EXCHANGE)
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL CDWDT                                                        &
        (GLOBAL,HYDRO,INPES,JNPES                                       &
        ,LM,NTIMESTEP,int_state%RESTART                                 &
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
      cdwdt_tim=cdwdt_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
!-----------------------------------------------------------------------
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL NHPOAV                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%DWDT,NTIMESTEP)
        poavhn_tim=poavhn_tim+timef()-btim
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
        fftfhn_tim=fftfhn_tim+timef()-btim
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%DWDT,LM                                  &
                    ,2,2)
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  VERTICALLY PROPAGATING FAST WAVES
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL VSOUND                                                       &
        (GLOBAL,HYDRO,int_state%RESTART                                 &
        ,LM,NTIMESTEP                                                   &
        ,CP,DT,PT,DSG2,PDSG1                                            &
        ,int_state%PD                                                   &
        ,int_state%CW,int_state%Q,int_state%RTOP                        &
        ,int_state%DWDT,int_state%T,int_state%W                         &
        ,int_state%PINT)
!
      vsound_tim=vsound_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%DWDT,NTIMESTEP)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%W,NTIMESTEP)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%PINT,NTIMESTEP)
        poavhn_tim=poavhn_tim+timef()-btim
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        polehn_tim=polehn_tim+timef()-btim
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
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      passive_advec: IF(MOD(NTIMESTEP,IDTAD)==0.AND.OLD_PASSIVE)THEN
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  VERTICAL ADVECTION OF PASSIVE QUANTITIES 
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
            int_state%Q=int_state%WATER(I,J,K,P_QV)                     &
                       /(1.+int_state%WATER(I,J,K,P_QV))
          ENDDO
          ENDDO
          ENDDO
!
        ENDIF vadv2_micro_check
        vadv2_tim=vadv2_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
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
            ,int_state%RRW                                              &
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
          fftfhn_tim=fftfhn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL SWAPHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES)
            ENDDO
          ENDIF
!
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL POLEHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES,JNPES)
            ENDDO
          ENDIF
!
          polehn_tim=polehn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%RRW,LM                                 &
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
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  HORIZONTAL ADVECTION OF PASSIVE QUANTITIES (INTERNAL HALO EXCHANGE)
!-----------------------------------------------------------------------
!
        btim=timef()
!
        hadv2_micro_check: IF(int_state%MICROPHYSICS=='fer')THEN
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,NTIMESTEP,INPES,JNPES                               &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%TRACERS                                          &
            ,int_state%NUM_TRACERS_MET,1,int_state%INDX_Q2)
!
          CALL AVEQ2                                                    &
            (LM                                                         &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,int_state%PD                                               &
            ,int_state%Q2,int_state%E2                                  &
            ,2)
!
!***  UPDATE WATER ARRAY.
!***  REMEMBER THAT WATER IS USED WITH THE WRF PHYSICS AND THUS
!***  THE P_QV SLOT (=2) IS MIXING RATIO, NOT SPECIFIC HUMIDITY.
!***  ALTHOUGH WATER IS ONLY USED FOR PHYSICS IN OPERATIONS, IT IS
!***  UPDATED HERE FROM Q EVERY ADVECTION TIMESTEP FOR NON-OPERATIONAL
!***  CONFIGURATIONS WHERE IT MAY BE USED OUTSIDE OF THE PHYSICS.
!
          IF(.NOT.OPERATIONAL_PHYSICS)THEN
            DO K=1,LM
            KFLIP=LM+1-K
            DO J=JTS,JTE
            DO I=ITS,ITE
              int_state%WATER(I,J,K,P_QV)=int_state%Q(I,J,K)/(1.-int_state%Q(I,J,K))
              WC = int_state%CW(I,J,K)
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
              int_state%WATER(I,J,K,P_QC)=QW
              int_state%WATER(I,J,K,P_QR)=QR
              int_state%WATER(I,J,K,P_QI)=0.
              int_state%WATER(I,J,K,P_QS)=QI
              int_state%WATER(I,J,K,P_QG)=0.
            ENDDO
            ENDDO
            ENDDO
          ENDIF
!
        ELSE hadv2_micro_check
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,NTIMESTEP,INPES,JNPES                               &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%Q2                                               &
            ,1,1,int_state%INDX_Q2)
!
          CALL HADV2_SCAL                                               &
            (GLOBAL,NTIMESTEP,INPES,JNPES                               &
            ,LM,IDTAD,DT,RDYH                                           &
            ,DSG2,PDSG1,PSGML1,SGML2                                    &
            ,DARE,RDXH                                                  &
            ,int_state%PD                                               &
            ,int_state%U,int_state%V                                    &
            ,int_state%WATER                                            &
            ,int_state%NUM_WATER,2,int_state%INDX_Q2)
!
        ENDIF hadv2_micro_check
!
        hadv2_tim=hadv2_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  FILTERING AND BOUNDARY CONDITIONS FOR GLOBAL FORECASTS
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL SWAPHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES)
            ENDDO
          ENDIF
!
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%RRW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!
          IF(int_state%MICROPHYSICS/='fer')THEN
            DO N=2,int_state%NUM_WATER
              CALL POLEHN(int_state%WATER(:,:,:,N)                      &
                         ,IMS,IME,JMS,JME,LM,INPES,JNPES)
            ENDDO
          ENDIF
!
          polehn_tim=polehn_tim+timef()-btim
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%RRW,LM                                 &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
!
        IF(int_state%MICROPHYSICS/='fer')THEN
          CALL HALO_EXCH(int_state%WATER,LM,int_state%NUM_WATER,2       &
                        ,2,2)
        ENDIF
!
        exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF passive_advec
!
!-----------------------------------------------------------------------
!***  CLOSE THE FILE UNITS USED FOR READS/WRITES OF GLOBAL SUMS
!***  IF THE FORECAST IS FINISHED.
!-----------------------------------------------------------------------
!
      IF(ESMF_ClockIsStopTime(clock=CLOCK_ATM,rc=RC))THEN
        IF(WRITE_GLOBAL_SUMS.AND.MYPE==0)THEN
          CLOSE(IUNIT_ADVEC_SUMS)
          CLOSE(IUNIT_POLE_SUMS)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  NOTE:  THE DYNAMICS EXPORT STATE IS FULLY UPDATED NOW
!***         BECAUSE SUBROUTINE DYN_INITIALIZE INSERTED THE 
!***         APPROPRIATE ESMF Fields INTO IT.  THOSE FIELDS 
!***         CONTAIN POINTERS TO THE ACTUAL DATA AND THOSE
!***         POINTERS ARE NEVER RE-DIRECTED, I.E., NO EXPLICIT
!***         ACTION IS NEEDED TO UPDATE THE DYNAMICS EXPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!     if(mod(ntimestep,20)==0.or.ntimestep<=5)then
!       call twr(int_state%t,lm,'t',ntimestep,mype,num_pes,mpi_comm_comp &
!               ,ids,ide,jds,jde &
!               ,ims,ime,jms,jme &
!               ,its,ite,jts,jte)
!     endif
!
!-----------------------------------------------------------------------
!***  WRITE THE LAYER STATISTICS
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP+1,N_PRINT_STATS)==0)THEN
!
        CALL FIELD_STATS(INT_STATE%T,MYPE,MPI_COMM_COMP,LM              &
                        ,ITS,ITE,JTS,JTE                                &
                        ,IMS,IME,JMS,JME                                &
                        ,IDS,IDE,JDS,JDE)
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        WRITE(0,25)NTIMESTEP,(NTIMESTEP+1)*DT/3600.
   25   FORMAT(' Finished Dyn Timestep ',i8,' ending at ',f10.3,' hours')
      ENDIF
!
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
      dyn_run_tim=dyn_run_tim+timef()-btim0
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
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK_ATM                                 &
                             ,RCFINAL)
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
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE                 !<-- The Dynamics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE_WRITE                 !<-- The Dynamics export state
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
      END SUBROUTINE DYN_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
