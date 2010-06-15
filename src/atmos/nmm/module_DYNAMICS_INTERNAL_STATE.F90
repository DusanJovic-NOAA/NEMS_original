!-----------------------------------------------------------------------

      MODULE MODULE_DYNAMICS_INTERNAL_STATE

!-----------------------------------------------------------------------
!***  This module declares the derived datatype called INTERNAL_STATE.
!***  For now the components of this datatype will be everything needed
!***  to advance the model integration, i.e. everything that would be
!***  part of a restart file.  Specifically this will include those
!***  quantities that evolve during the integration, the namelist
!***  variables, and the grid decomposition variables.
!-----------------------------------------------------------------------
!
      USE module_INCLUDE
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
      USE module_VARS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DYNAMICS_INTERNAL_STATE                                 &
               ,SET_INTERNAL_STATE_DYN_1                                &
               ,SET_INTERNAL_STATE_DYN_2                                &
               ,WRAP_DYN_INT_STATE 
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: MAX_VARS = 100
!
      TYPE DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Begin with the Configure File variables.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT), POINTER :: IM,JM,LM
        INTEGER(kind=KINT) :: INPES,JNPES                               &
                             ,MINUTES_HISTORY                           &
                             ,MINUTES_RESTART                           &
                             ,NHOURS_FCST                               &
                             ,NSTEPS_BC_RESTART                         &
                             ,NUM_TRACERS_MET                           &  !<-- Number of meteorological tracers (e.g. water)
                             ,NUM_TRACERS_CHEM                          &  !<-- Number of chem/aerosol tracers
                             ,START_YEAR                                &
                             ,START_MONTH                               &
                             ,START_DAY                                 &
                             ,START_HOUR                                &
                             ,START_MINUTE                              &
                             ,START_SECOND                              &
                             ,FILTER_METHOD                             
!
        REAL(kind=KFPT), POINTER :: DT,SBD,TSTART,TPH0D,TLM0D,WBD
        REAL(kind=KFPT) :: CODAMP                                       &
                          ,RUN_DURATION                                 &
                          ,SMAG2                                        &       
                          ,WCOR
!
        LOGICAL(kind=KLOG) :: ADIABATIC                                 &
                             ,ADVECT_TRACERS                            &
                             ,FREERUN                                   &
                             ,GLOBAL                                    &
                             ,HYDRO                                     &
                             ,NEMSIO_INPUT                              &
                             ,RESTART                                   &
                             ,SECADV                                    &
                             ,SPEC_ADV                                  &
                             ,READ_GLOBAL_SUMS                          &
                             ,WRITE_GLOBAL_SUMS
!
!-----------------------------------------------------------------------
!***  Distributed memory information.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: IHALO,JHALO,MYPE,NHALO,NUM_PES            &
                             ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP        
!
        INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE                           &
                             ,IMS,IME,JMS,JME                           &
                             ,IDS,IDE,JDS,JDE
!
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART         &
                                                  ,LOCAL_IEND           &
                                                  ,LOCAL_JSTART         &
                                                  ,LOCAL_JEND
!
!-----------------------------------------------------------------------
!***  Horizontal and vertical grid-related variables.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT), POINTER :: LPT2
!
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: KHFILT,KVFILT        &
                                                  ,NFFTRH,NFFTRW        &
                                                  ,NHSMUD
!
        REAL(kind=KFPT), POINTER :: DPHD,DLMD,DYH,PDTOP,PT
        REAL(kind=KFPT) :: DYV,RDYH,RDYV                                &
                          ,DDMPV                                        &
                          ,EF4T
!
        REAL(kind=KFPT),DIMENSION(:),POINTER :: SG1,PSG1                &
                                               ,DSG1,PDSG1              &
                                               ,SG2,DSG2                &
                                               ,SGML1,PSGML1            &
                                               ,SGML2,SGM               &
                                               ,DDMPU,WPDAR             &
                                               ,FCP,FDIV                &
                                               ,CURV,DDV                &
                                               ,DARE,RARE               &
                                               ,FAD,FAH                 &
                                               ,DXH,RDXH                &
                                               ,DXV,RDXV                &
                                               ,RDDV                    &
                                               ,WFFTRH,WFFTRW
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: F,FIS                 &
                                                 ,GLAT,GLON             &
                                                 ,HDACX,HDACY           &
                                                 ,HDACVX,HDACVY         &
                                                 ,HFILT,VFILT           &
                                                 ,VLAT,VLON
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: PD,PDO,PSDT,SICE,SM
!
!-----------------------------------------------------------------------
!***  Integration quantities.
!-----------------------------------------------------------------------
!
        LOGICAL(kind=KLOG) :: FIRST,READBC                              &
                             ,ADV_STANDARD,ADV_UPSTREAM
!
        INTEGER(kind=KINT), POINTER :: IHRST
        INTEGER(kind=KINT) :: NTSD,IDTAD,IDTADT,IHR,IHREND              &
                             ,LNSAD,NBOCO,NTSTI,NTSTM,NTSTM_MAX
!
        INTEGER(kind=KINT),DIMENSION(:), POINTER :: IDAT
!
        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: INSOIL,INVEG
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: PINT,RTOP           &
                                                   ,T,U,V               &
                                                   ,Q,CW,O3             &
                                                   ,Q2,E2               &
                                                   ,DIV                 &
                                                   ,TDIV                &
                                                   ,W,Z                 &
                                                   ,DWDT,PDWDT          &
                                                   ,OMGALF              &
                                                   ,PSGDT               &
                                                   ,TP,UP,VP            &
                                                   ,PCNE,PCNW           &
                                                   ,PFNE,PFNW           &
                                                   ,PCX,PCY             &
                                                   ,PFX,PFY             &
                                                   ,TCT,TCU,TCV
!
        LOGICAL(kind=KLOG) :: RUN
!
!-----------------------------------------------------------------------
!***  The general 4-D arrays for 3-D "tracers".
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: NUM_TRACERS_TOTAL                            !<-- Total number of "tracer" variables.
!
!-----------------------------------------------
!***  Declare indices of meteorological tracers
!-----------------------------------------------
!
        INTEGER(kind=KINT) :: INDX_Q                                    &  !<-- Location of Q in tracer arrays
                             ,INDX_CW                                   &  !<-- Location of CW in tracer arrays
                             ,INDX_O3                                   &  !<-- Location of O3 in tracer arrays
                             ,INDX_Q2                                      !<-- Location of Q2 in tracer arrays
!
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: TRACERS           &  !<-- Storage array for "tracer" variables.
                                                     ,TRACERS_SQRT      &  !<-- Sqrt of the tracer variables (for advection)
                                                     ,TRACERS_PREV      &  !<-- Values of tracer variables in prev timestep (for advection)
                                                     ,TRACERS_TEND         !<-- Tendencies of tracer variables (for advection)
!
!-----------------------------------------------------------------------
!***  Boundary conditions.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: IHRSTBC,LNSH,LNSV
!
        INTEGER(kind=KINT),DIMENSION(3) :: IDATBC
!
        REAL(kind=KFPT) :: TBOCO
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: PDBS,PDBN           &
                                                   ,PDBW,PDBE
!
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: TBS,TBN           & 
                                                     ,TBW,TBE           &
                                                     ,QBS,QBN           &
                                                     ,QBW,QBE           &
                                                     ,UBS,UBN           &
                                                     ,UBW,UBE           &
                                                     ,VBS,VBN           &
                                                     ,VBW,VBE           &
                                                     ,WBS,WBN           &
                                                     ,WBW,WBE
!
        LOGICAL(kind=KLOG) :: RUNBC
!
!----------------------------
!***  For 1-D restart output
!----------------------------
!
!
      INTEGER(kind=KINT) :: NUM_WORDS_BC_SOUTH                          &  !<-- Word counts of 1-D boundary data strings
                           ,NUM_WORDS_BC_NORTH                          &  !    for each side of the domain.
                           ,NUM_WORDS_BC_WEST                           &  !
                           ,NUM_WORDS_BC_EAST                              !<--
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RST_BC_DATA_SOUTH     &  !<-- 1-D strings of boundary data
                                                 ,RST_BC_DATA_NORTH     &  !    for each side of the domain.
                                                 ,RST_BC_DATA_WEST      &  !
                                                 ,RST_BC_DATA_EAST         !<--
!
!-----------------------------------------------------------------------
!***  FFT arrays.
!-----------------------------------------------------------------------
!
        REAL(kind=KDBL),DIMENSION(:),POINTER :: CRAUX1,CRAUX2,CRAUX3    &
                                               ,RCAUX1,RCAUX2,RCAUX3
!
!-----------------------------------------------------------------------
!***  Some physics variables are needed.
!-----------------------------------------------------------------------
!
        CHARACTER(99) :: SHORTWAVE                                      &
                        ,LONGWAVE                                       &   
                        ,SFC_LAYER                                      &
                        ,TURBULENCE                                     &
                        ,CONVECTION                                     &
                        ,MICROPHYSICS
!
        INTEGER(kind=KINT) :: NUM_WATER                                 &  !<-- 1 + types of water substance in microphysics
                             ,P_QV                                      &  !<-- Index for water vapor in WATER array
                             ,P_QC                                      &  !<-- Index for cloud water in WATER array
                             ,P_QR                                      &  !<-- Index for rain in WATER array
                             ,P_QI                                      &  !<-- Index for cloud ice in WATER array
                             ,P_QS                                      &  !<-- Index for snow in WATER array
                             ,P_QG                                         !<-- Index for graupel in WATER array
!
        INTEGER(kind=KINT) :: INDX_WATER_START                          &  !<-- Start index of the water in tracers array
                             ,INDX_WATER_END                               !<-- End index of the water in tracers array
!
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: WATER                !<-- Storage array for water substance
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: F_ICE,F_RAIN        &  !<-- Fractions of ice, rain, and rime
                                                   ,F_RIMEF
!
        LOGICAL :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
!-----------------------------------------------------------------------
!***  Nesting
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT), POINTER :: I_PAR_STA                        &  !<-- SW corner of nest domain on this parent I
                                      ,J_PAR_STA                           !<-- SW corner of nest domain on this parent J
        INTEGER(kind=KINT) :: PARENT_CHILD_TIME_RATIO                      !<-- # of child timesteps per parent timestep
!
!-----------------------------------------------------------------------
!***  Variable table
!-----------------------------------------------------------------------
!
        TYPE(VAR),DIMENSION(MAX_VARS) :: VARS
        INTEGER :: NUM_VARS = 0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END TYPE DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE INTERNAL_STATE TYPE IS SUPPORTED BY A C POINTER (NOT AN F90
!***  POINTER) AND THEREFORE THE FOLLOWING TYPE IS NEEDED.
!-----------------------------------------------------------------------
!
      TYPE WRAP_DYN_INT_STATE
        TYPE(DYNAMICS_INTERNAL_STATE),POINTER :: INT_STATE
      END TYPE WRAP_DYN_INT_STATE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_INTERNAL_STATE_DYN_1(INT_STATE,LM)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  Carry out all allocations of internal state quantities in this
!***  phase.  Put export variables into the export state.
!***  No initialization is done until phase 2.
!
!-----------------------------------------------------------------------
!
      USE module_DM_PARALLEL,ONLY: IDS,IDE,JDS,JDE                      &
                                  ,IMS,IME,JMS,JME                      &
                                  ,IHALO,JHALO
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(DYNAMICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- The Dynamics internal state
!
      INTEGER, INTENT(IN) :: LM                                            !<-- Number of model layers
!      
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_CYCLE,J                                 &
                           ,L,LNSH,LNSV,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The array called WATER is a special case needed to satisfy
!***  various WRF physics options.  The 4th dimension is set to
!***  1+Number_of_species including vapor.  The "1+" is needed
!***  because WRF never touches the first level.
!
!***  Set the P_ and F_ variables.
!***  V=>vapor; C=>cloudwater; R=>rain; I=>cloudice; S=>snow; G=>graupel
!***  Set the appropriate value for the logical F_ variables.
!***  For each species that is .TRUE., set the integer P_ variable
!***  to monotonically increasing values STARTING WITH 2.
!***  For the excluded species (F_*=.FALSE.), set the P_ variable to 1.
!-----------------------------------------------------------------------
!
      IF(TRIM(int_state%MICROPHYSICS)=='fer')THEN
        int_state%NUM_WATER=1+4
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=5
        int_state%P_QI=1
        int_state%P_QG=1
        int_state%F_QV=.TRUE.
        int_state%F_QC=.TRUE.
        int_state%F_QR=.TRUE.
        int_state%F_QS=.TRUE.
        int_state%F_QI=.FALSE.
        int_state%F_QG=.FALSE.
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='wsm3')THEN
        int_state%NUM_WATER=1+3
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=1
        int_state%P_QI=1
        int_state%P_QG=1
        int_state%F_QV=.TRUE.
        int_state%F_QC=.TRUE.
        int_state%F_QR=.TRUE.
        int_state%F_QS=.FALSE.
        int_state%F_QI=.FALSE.
        int_state%F_QG=.FALSE.
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='wsm6')THEN
        int_state%NUM_WATER=1+6
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=5
        int_state%P_QI=6
        int_state%P_QG=7
        int_state%F_QV=.TRUE.
        int_state%F_QC=.TRUE.
        int_state%F_QR=.TRUE.
        int_state%F_QS=.TRUE.
        int_state%F_QI=.TRUE.
        int_state%F_QG=.TRUE.
      ENDIF

      int_state%NUM_TRACERS_TOTAL=                                      &  !<-- # of 3-D arrays in 4-D TRACERS array
                                  int_state%NUM_TRACERS_MET             &  !<-- # of water, etc. tracers specified now (see below)
                                 +int_state%NUM_TRACERS_CHEM            &  !<-- # of specified scalars (chem, aerosol, etc.)
                                 +int_state%NUM_WATER                      !<-- # of water types

!------------------------------------------------
!***  Read and store the specifications for each
!***  internal state variable listed by the user
!***  in the Dynamics text file. 
!------------------------------------------------
!
      CALL READ_CONFIG('dyn_state.txt',int_state%VARS,int_state%NUM_VARS)
!
!-------------------------------------------------------------------
!***  Allocate appropriate memory within the Dynamics' composite 
!***  VARS array for all internal state variables that are 'Owned'
!***  by Dynamics and point those variables into that memory.
!***  In this step, all unowned variables will be pointed at NULL
!***  for the moment.
!-------------------------------------------------------------------
!
      CALL SET_DYN_VAR_PTR(INT_STATE, .TRUE., LM)
!
!-----------------------------------------------------------------------
!***  WE CAN RETRIEVE LM FROM THE INTERNAL STATE SINCE IT WAS
!***  PLACED THERE ALREADY FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
      I_CYCLE=IDE-3
!
      LNSH=int_state%LNSH
      LNSV=int_state%LNSV
!
!!!   int_state%NBOCO=NBOCO   !<-- Set later in call to CONSTS
!!!   int_state%TBOCO=TBOCO   !<-- Set later in call to CONSTS
!
!-----------------------------------------------------------------------
!***  Explicitly allocate the arrays of the Dynamics internal state
!***  that are never outside of the Dynamics component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Grid-related constants.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%PDSG1 (1:LM))                                    !<-- Thicknesses of pressure layers in press. range
      ALLOCATE(int_state%PSGML1(1:LM))                                    !<-- Pressure at midlayers in pressure range
!   
      ALLOCATE(int_state%PSG1(1:LM+1))                                    !<-- Pressure at interfaces in pressure range  (Pa)
!     
      ALLOCATE(int_state%KHFILT(JDS:JDE))                                 !<-- Polar filter, truncation wave #, h points
      ALLOCATE(int_state%KVFILT(JDS:JDE))                                 !<-- Polar filter, truncation wave #, v points
      ALLOCATE(int_state%NHSMUD(JMS:JME))                                 !<-- Polar smoother for unfiltered variables
!     
      ALLOCATE(int_state%HFILT(IDS:IDE,JDS:JDE))                          !<-- Polar filter, h points
      ALLOCATE(int_state%VFILT(IDS:IDE,JDS:JDE))                          !<-- Polar filter, v points
!     
      ALLOCATE(int_state%CURV (JDS:JDE))                                  !<-- Curvature term in coriolis force  (m-1)
      ALLOCATE(int_state%DARE (JDS:JDE))                                  !<-- Gridbox area  (m2)
      ALLOCATE(int_state%DDMPU(JDS:JDE))                                  !<-- Divergence damping coefficient, x direction  (m)
      ALLOCATE(int_state%DDV  (JDS:JDE))                                  !<-- Gridbox diagonal distance  (m)
      ALLOCATE(int_state%DXV  (JDS:JDE))                                  !<-- Delta x, v points  (m)
      ALLOCATE(int_state%FAD  (JDS:JDE))                                  !<-- Momentum advection factor
      ALLOCATE(int_state%FAH  (JDS:JDE))                                  !<-- z, w advection factor
      ALLOCATE(int_state%FCP  (JDS:JDE))                                  !<-- Temperature advection factor
      ALLOCATE(int_state%FDIV (JDS:JDE))                                  !<-- Divergence factor
      ALLOCATE(int_state%RARE (JDS:JDE))                                  !<-- 1 / gridbox area  (m-2)
      ALLOCATE(int_state%RDDV (JDS:JDE))                                  !<-- 1 / gridbox diagonal distance  (m-1)
      ALLOCATE(int_state%RDXH (JDS:JDE))                                  !<-- 1 / delta x, h points  (m-1)
      ALLOCATE(int_state%RDXV (JDS:JDE))                                  !<-- 1 / delta x, v points  (m-1)
      ALLOCATE(int_state%WPDAR(JDS:JDE))                                  !<-- Weight of grid separation filter
!     
      ALLOCATE(int_state%WFFTRH(1:2*I_CYCLE))                             !<-- FFT working field, h points
      ALLOCATE(int_state%WFFTRW(1:2*I_CYCLE))                             !<-- FFT working field, v points
      ALLOCATE(int_state%NFFTRH(1:15))                                    !<-- FFT working field, h points
      ALLOCATE(int_state%NFFTRW(1:15))                                    !<-- FFT working field, v points
!
      ALLOCATE(int_state%F     (IMS:IME,JMS:JME))                         !<-- Coriolis parameter  (s-1)
      ALLOCATE(int_state%HDACX (IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, h points  (s m-1)
      ALLOCATE(int_state%HDACY (IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, h points  (s m-1)
      ALLOCATE(int_state%HDACVX(IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, v points  (s m-1)
      ALLOCATE(int_state%HDACVY(IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, v points  (s m-1)
!
!-----------------------------------------------------------------------
!***  Local horizontal subdomain limits for all forecast tasks.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%LOCAL_ISTART(0:int_state%NUM_PES-1))
      ALLOCATE(int_state%LOCAL_IEND  (0:int_state%NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JSTART(0:int_state%NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JEND  (0:int_state%NUM_PES-1))
!
      int_state%IMS=IMS
      int_state%IME=IME
      int_state%JMS=JMS
      int_state%JME=JME
      int_state%IDS=IDS
      int_state%IDE=IDE
      int_state%JDS=JDS
      int_state%JDE=JDE
!
      int_state%IHALO=IHALO
      int_state%JHALO=JHALO
!
!-----------------------------------------------------------------------
!***  Fixed surface fields
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%SM(IMS:IME,JMS:JME))                            !<-- Sea mask
      ALLOCATE(int_state%SICE(IMS:IME,JMS:JME))                          !<-- Sea ice
!
!-----------------------------------------------------------------------
!***  Regional boundary conditions.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%UBN(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- U wind component at northern boundary  (m s-1)
      ALLOCATE(int_state%UBS(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- U wind component at southern boundary  (m s-1)
      ALLOCATE(int_state%VBN(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- V wind component at northern boundary  (m s-1)
      ALLOCATE(int_state%VBS(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- V wind component at southern boundary  (m s-1)
!
      ALLOCATE(int_state%UBE(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- U wind component at eastern boundary  (m s-1)
      ALLOCATE(int_state%UBW(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- U wind component at western boundary  (m s-1)
      ALLOCATE(int_state%VBE(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- V wind component at eastern boundary  (m s-1)
      ALLOCATE(int_state%VBW(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- V wind component at western boundary  (m s-1)
!
      IF(.NOT.int_state%GLOBAL)THEN
!
        ALLOCATE(int_state%PDBN(IMS:IME,1:LNSH,1:2))                       !<-- Pressure difference at northern boundary  (Pa)
        ALLOCATE(int_state%PDBS(IMS:IME,1:LNSH,1:2))                       !<-- Pressure difference at southern boundary  (Pa)
!
        ALLOCATE(int_state%PDBE(1:LNSH,JMS:JME,1:2))                       !<-- Pressure difference at eastern boundary  (Pa)
        ALLOCATE(int_state%PDBW(1:LNSH,JMS:JME,1:2))                       !<-- Pressure difference at western boundary  (Pa)
!
        ALLOCATE(int_state%QBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Specific humidity at northern boundary  (kg kg-1)
        ALLOCATE(int_state%QBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Specific humidity at southern boundary  (kg kg-1)
        ALLOCATE(int_state%TBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Temperature at northern boundary  (K)
        ALLOCATE(int_state%TBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Temperature at southern boundary  (K)
        ALLOCATE(int_state%WBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Vertical velocity at northern boundary  (m s-1)
        ALLOCATE(int_state%WBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Vertical velocity at southern boundary  (m s-1)
!
        ALLOCATE(int_state%QBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Specific humidity at eastern boundary  (kg kg-1)
        ALLOCATE(int_state%QBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Specific humidity at western boundary  (kg kg-1)
        ALLOCATE(int_state%TBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Temperature at eastern boundary  (K)
        ALLOCATE(int_state%TBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Temperature at western boundary  (K)
        ALLOCATE(int_state%WBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Vertical velocity at eastern boundary  (m s-1)
        ALLOCATE(int_state%WBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Vertical velocity at western boundary  (m s-1)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Atmospheric variables, hydrostatic (mostly)
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%PSDT(IMS:IME,JMS:JME))                            !<-- Hydrostatic surface pressure tendency  (Pa s-1)
!
      ALLOCATE(int_state%F_ICE(IMS:IME,JMS:JME,1:LM))                        !<-- Fraction of condensate as ice
      ALLOCATE(int_state%F_RAIN(IMS:IME,JMS:JME,1:LM))                       !<-- Fraction of condensate as rain
      ALLOCATE(int_state%F_RIMEF(IMS:IME,JMS:JME,1:LM))                      !<-- Fraction of condensate as rime
!
!-----------------------------------------------------------------------
!***  THE TRACERS ARRAY HOLDS ALL GENERAL "TRACER" VARIABLES INCLUDING
!***  WATER.  PLACE THE DESIRED VARIABLES AT THE TOP OF THE TRACERS
!***  ARRAY, LEVEL 1 THROUGH NUM_TRACERS_MET.  ALL OTHER SCALAR VARIABLES
!***  (e.g., chemistry and aerosols) WILL FOLLOW.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%TRACERS_SQRT(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- Sqrt of tracers (for advection)
      ALLOCATE(int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- Tendency of tracers (for advection)
!
!-----------------------------------------------------------------------
!***  Atmospheric variables, nonhydrostatic
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%PDWDT(IMS:IME,JMS:JME,1:LM))                     !<-- Correction factor, previous step  (m s-2)
!
!-----------------------------------------------------------------------
!***  Working arrays passed as arguments between subroutines.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%PCNE(IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, NE direction
      ALLOCATE(int_state%PCNW(IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, NW direction
      ALLOCATE(int_state%PCX (IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, X direction
      ALLOCATE(int_state%PCY (IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, Y direction
      ALLOCATE(int_state%PFNE(IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, NE direction
      ALLOCATE(int_state%PFNW(IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, NW direction
      ALLOCATE(int_state%PFX (IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, X direction
      ALLOCATE(int_state%PFY (IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, Y direction
      ALLOCATE(int_state%TDIV(IMS:IME,JMS:JME,1:LM))                      !<-- Integrated horizontal mass divergence
!
!-----------------------------------------------------------------------
!***  FFT arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE%CRAUX1(1:25000))                                 !<-- FFT working field
      ALLOCATE(INT_STATE%CRAUX2(1:20000))                                 !<-- FFT working field
      ALLOCATE(INT_STATE%CRAUX3(1:1    ))                                 !<-- FFT working field
      ALLOCATE(INT_STATE%RCAUX1(1:25000))                                 !<-- FFT working field
      ALLOCATE(INT_STATE%RCAUX2(1:20000))                                 !<-- FFT working field
      ALLOCATE(INT_STATE%RCAUX3(1:1    ))                                 !<-- FFT working field
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_INTERNAL_STATE_DYN_1
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_INTERNAL_STATE_DYN_2(INT_STATE,LM)
!
!-----------------------------------------------------------------------
!***  Point all unallocated internal state variables into allocated
!***  memory from Physics.  Do fundamental initialization.
!-----------------------------------------------------------------------
!
      USE module_DM_PARALLEL,ONLY: IDS,IDE,JDS,JDE                      &
                                  ,IMS,IME,JMS,JME                      &
                                  ,IHALO,JHALO
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(DYNAMICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- Dynamics internal state
!
      INTEGER(kind=KINT),INTENT(IN) :: LM                                  !<-- Number of model layers
!      
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_CYCLE,INDX,J,L,LL,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the Dynamics' internal state variables listed in
!***  the user's text file and re-point them into the composite VARS
!***  array.  In phase 1 of this step, all unowned variables were
!***  pointed at NULL (unallocated pointers).  Now in phase 2 the
!***  locations in VARS associated with unowned variables are pointing
!***  into allocated memory of the owned variables in Physics.
!-----------------------------------------------------------------------
!
      CALL SET_DYN_VAR_PTR(INT_STATE, .FALSE., LM)
!
!-----------------------------------------------------------------------
!***  Point Q at level 1 of the TRACERS array.
!-----------------------------------------------------------------------
!
      int_state%INDX_Q=1
      CALL FIND_VAR_INDX('Q',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_Q)
      int_state%Q=>int_state%VARS(I)%R3D
!
!-----------------------------------------------------------------------
!***  Additional tracers.
!-----------------------------------------------------------------------
!
!--------------------------------
!***  Combined cloud water array
!--------------------------------
!
      int_state%INDX_CW=2
      CALL FIND_VAR_INDX('CW',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_CW)
      int_state%CW=>int_state%VARS(I)%R3D
!
!--------------------------------
!***  Turbulence kinetic energy
!--------------------------------
!
      int_state%INDX_Q2=3
      CALL FIND_VAR_INDX('E2',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_Q2)
      int_state%E2=>int_state%VARS(I)%R3D
!
!--------------------------------
!***  General tracer for testing
!--------------------------------
!
      int_state%INDX_O3=4
      CALL FIND_VAR_INDX('O3',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_O3)
      int_state%O3=>int_state%VARS(I)%R3D
!
!--------------------------------
!***  Water tracers
!--------------------------------
!
      int_state%INDX_WATER_START = int_state%NUM_TRACERS_MET + int_state%NUM_TRACERS_CHEM + 1
      int_state%INDX_WATER_END = int_state%INDX_WATER_START + int_state%NUM_WATER - 1
      int_state%WATER=>int_state%TRACERS(:,:,:,int_state%INDX_WATER_START:int_state%INDX_WATER_END)
!
!
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%F(I,J)=-999.
	int_state%GLAT(I,J)=-999.
	int_state%GLON(I,J)=-999.
	int_state%VLAT(I,J)=-999.
	int_state%VLON(I,J)=-999.
        int_state%HDACX(I,J)=-999.
        int_state%HDACY(I,J)=-999.
        int_state%HDACVX(I,J)=-999.
        int_state%HDACVY(I,J)=-999.
      ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=-1.E6
      ENDDO
      ENDDO
!
      DO N=1,2
      DO L=1,LM
      DO LL=1,int_state%LNSV
      DO I=IMS,IME
        int_state%UBN(I,LL,L,N)=-1.E6
        int_state%UBS(I,LL,L,N)=-1.E6
        int_state%VBN(I,LL,L,N)=-1.E6
        int_state%VBS(I,LL,L,N)=-1.E6
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
      DO N=1,2
      DO L=1,LM
      DO J=JMS,JME
      DO LL=1,int_state%LNSV
        int_state%UBE(LL,J,L,N)=-1.E6
        int_state%UBW(LL,J,L,N)=-1.E6
        int_state%VBE(LL,J,L,N)=-1.E6
        int_state%VBW(LL,J,L,N)=-1.E6
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
      IF(.NOT.int_state%GLOBAL)THEN
!
        DO N=1,2
        DO LL=1,int_state%LNSH
        DO I=IMS,IME
          int_state%PDBN(I,LL,N)=0.
          int_state%PDBS(I,LL,N)=0.
        ENDDO
        ENDDO
        ENDDO
!
        DO N=1,2
        DO J=JMS,JME
        DO LL=1,int_state%LNSH
          int_state%PDBE(LL,J,N)=0.
          int_state%PDBW(LL,J,N)=0.
        ENDDO
        ENDDO
        ENDDO
!
        DO N=1,2
        DO L=1,LM
        DO LL=1,int_state%LNSH
        DO I=IMS,IME
          int_state%QBN(I,LL,L,N)=-999.
          int_state%QBS(I,LL,L,N)=-999.
          int_state%TBN(I,LL,L,N)=-999.
          int_state%TBS(I,LL,L,N)=-999.
          int_state%WBN(I,LL,L,N)=-999.
          int_state%WBS(I,LL,L,N)=-999.
        ENDDO
        ENDDO
        ENDDO
        ENDDO
!
        DO N=1,2
        DO L=1,LM
        DO J=JMS,JME
        DO LL=1,int_state%LNSH
          int_state%QBE(LL,J,L,N)=-999.
          int_state%QBW(LL,J,L,N)=-999.
          int_state%TBE(LL,J,L,N)=-999.
          int_state%TBW(LL,J,L,N)=-999.
          int_state%WBE(LL,J,L,N)=-999.
          int_state%WBW(LL,J,L,N)=-999.
        ENDDO
        ENDDO
        ENDDO
        ENDDO
!
        int_state%NUM_WORDS_BC_SOUTH=-1                                    !<-- Word counts of 1-D boundary data strings
        int_state%NUM_WORDS_BC_NORTH=-1                                    !
        int_state%NUM_WORDS_BC_WEST =-1                                    !
        int_state%NUM_WORDS_BC_EAST =-1                                    !<--
!
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
        int_state%PDO(I,J)=0.
        int_state%PSDT(I,J)=-1.E6
      ENDDO
      ENDDO
!
      DO L=1,LM-1
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PSGDT(I,J,L)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%F_ICE(I,J,L)=0.
        int_state%F_RAIN(I,J,L)=0.
        int_state%F_RIMEF(I,J,L)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO N=1,int_state%NUM_TRACERS_MET
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%TRACERS     (I,J,L,N)=1.E-20
        int_state%TRACERS_SQRT(I,J,L,N)=1.E-20
        int_state%TRACERS_PREV(I,J,L,N)=1.E-20
        int_state%TRACERS_TEND(I,J,L,N)=1.E-20
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%OMGALF(I,J,L)=-1.E6
        int_state%Q2(I,J,L)    = 0.02
        int_state%T(I,J,L)     =-1.E6
        int_state%TP(I,J,L)    =-1.E6
        int_state%U(I,J,L)     =-1.E6
        int_state%UP(I,J,L)    =-1.E6
        int_state%V(I,J,L)     =-1.E6
        int_state%VP(I,J,L)    =-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%DWDT(I,J,L)=-1.E6
        int_state%PDWDT(I,J,L)=-1.E6
        int_state%RTOP(I,J,L)=-1.E6
        int_state%W(I,J,L)=-1.E6
        int_state%Z(I,J,L)=-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%DIV(I,J,L) =-1.E6
        int_state%PCNE(I,J,L)=-1.E6
        int_state%PCNW(I,J,L)=-1.E6
        int_state%PCX(I,J,L) =-1.E6
        int_state%PCY(I,J,L) =-1.E6
        int_state%PFNE(I,J,L)=-1.E6
        int_state%PFNW(I,J,L)=-1.E6
        int_state%PFX(I,J,L) =-1.E6
        int_state%PFY(I,J,L) =-1.E6
        int_state%TDIV(I,J,L)=-1.E6
        int_state%TCT(I,J,L) =-1.E6
        int_state%TCU(I,J,L) =-1.E6
        int_state%TCV(I,J,L) =-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      int_state%I_PAR_STA=0
      int_state%J_PAR_STA=0
!
      int_state%PARENT_CHILD_TIME_RATIO=-999
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_INTERNAL_STATE_DYN_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_DYN_VAR_PTR(INT_STATE,ALLOC_FLAG,LM)

!---------------------------------------------------------------
!***  We must allocate memory within the composite VARS array
!***  for each Dynamics internal state variable that is 'Owned'
!***  when ALLOC_FLAG is true in phase 1 of Dynamics Init.
!***  In phase 2 ALLOC_FLAG is false, no allocation takes
!***  place, and the designated internal state variables are
!***  re-pointed into VARS whose locations associated with
!***  unowned variables will be pointing at allocated memory
!***  for the owned variable in Physics.
!---------------------------------------------------------------
!
      USE module_DM_PARALLEL,ONLY: IMS,IME,JMS,JME,JDS,JDE
!
      IMPLICIT NONE
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(DYNAMICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- Dynamics internal state
      LOGICAL, INTENT(IN) :: ALLOC_FLAG                                    !<-- Do we want to allocate the variables?
      INTEGER, INTENT(IN) :: LM                                            !<-- Number of model layers
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'IM'        ,int_state%IM        )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'JM'        ,int_state%JM        )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'LM'        ,int_state%LM        )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'IHRST'     ,int_state%IHRST     )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'I_PAR_STA' ,int_state%I_PAR_STA )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'J_PAR_STA' ,int_state%J_PAR_STA )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'LPT2'      ,int_state%LPT2      )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DT'        ,int_state%DT        )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DYH'       ,int_state%DYH       )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PDTOP'     ,int_state%PDTOP     )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PT'        ,int_state%PT        )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TLM0D'     ,int_state%TLM0D     )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TPH0D'     ,int_state%TPH0D     )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TSTART'    ,int_state%TSTART    )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DPHD'      ,int_state%DPHD      )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DLMD'      ,int_state%DLMD      )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SBD'       ,int_state%SBD       )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'WBD'       ,int_state%WBD       )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'IDAT'      ,int_state%IDAT    ,1 ,3 )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DXH'       ,int_state%DXH     ,JDS, JDE )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SG1'       ,int_state%SG1     ,1, LM+1  )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SG2'       ,int_state%SG2     ,1, LM+1  )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DSG1'      ,int_state%DSG1    ,1, LM    )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DSG2'      ,int_state%DSG2    ,1, LM    )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SGML1'     ,int_state%SGML1   ,1, LM    )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SGML2'     ,int_state%SGML2   ,1, LM    )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SGM'       ,int_state%SGM     ,1, LM+1  )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'FIS'       ,int_state%FIS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'GLAT'      ,int_state%GLAT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'GLON'      ,int_state%GLON    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PD'        ,int_state%PD      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'VLAT'      ,int_state%VLAT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'VLON'      ,int_state%VLON    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PDO'       ,int_state%PDO     ,(/ IMS,JMS /),(/ IME,JME /) )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'T'         ,int_state%T       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Q'         ,int_state%Q       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'U'         ,int_state%U       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'V'         ,int_state%V       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Q2'        ,int_state%Q2      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CW'        ,int_state%CW      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'W'         ,int_state%W       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DWDT'      ,int_state%DWDT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PINT'      ,int_state%PINT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM+1 /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'OMGALF'    ,int_state%OMGALF  ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'O3'        ,int_state%O3      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'DIV'       ,int_state%DIV     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RTOP'      ,int_state%RTOP    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TCU'       ,int_state%TCU     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TCV'       ,int_state%TCV     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TCT'       ,int_state%TCT     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TP'        ,int_state%TP      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'UP'        ,int_state%UP      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'VP'        ,int_state%VP      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'E2'        ,int_state%E2      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PSGDT'     ,int_state%PSGDT   ,(/ IMS,JMS,1 /),(/ IME,JME,LM-1 /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Z'         ,int_state%Z       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /) )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TRACERS'     ,int_state%TRACERS      ,(/ IMS,JMS,1,1 /),(/ IME,JME,LM,int_state%NUM_TRACERS_TOTAL /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TRACERS_PREV',int_state%TRACERS_PREV ,(/ IMS,JMS,1,1 /),(/ IME,JME,LM,int_state%NUM_TRACERS_TOTAL /) )

      DO N=1,int_state%NUM_VARS
        IF (int_state%VARS(N)%TKR==0) THEN
          write(0,*)' Error in SET_DYN_VAR_PTR. '
          write(0,*)' Variable ',TRIM(int_state%VARS(N)%VBL_NAME),' is not associated to an internal state fortran pointer'
          STOP 
        END IF
      END DO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_DYN_VAR_PTR
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE module_DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
