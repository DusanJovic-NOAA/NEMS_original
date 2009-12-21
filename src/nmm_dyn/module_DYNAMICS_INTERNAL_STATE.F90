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
      USE ESMF_Mod
      USE MODULE_INCLUDE
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
      PUBLIC :: DYNAMICS_INTERNAL_STATE                                 &
               ,SET_INTERNAL_STATE_DYN                                  &
               ,UPDATE_INTERNAL_STATE_DYN                               &
               ,WRAP_DYN_INT_STATE 
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      TYPE DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Begin with the Configure File variables.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: IM,JM,LM                                  &
                             ,INPES,JNPES                               &
                             ,NHOURS_FCST                               &
                             ,NHOURS_HISTORY                            &
                             ,NHOURS_RESTART                            &
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
        REAL(KIND=KFPT) :: CODAMP,DT                                    &
                          ,RUN_DURATION                                 &
                          ,SBD                                          &
                          ,SMAG2                                        &       
                          ,TPH0D,TLM0D                                  &
                          ,TSTART                                       &
                          ,WBD,WCOR
!
        LOGICAL(KIND=KLOG) :: ADIABATIC                                 &
                             ,ADVECT_TRACERS                            &
                             ,HYDRO                                     &
                             ,GLOBAL                                    &
                             ,RESTART                                   &
                             ,SECADV                                    &
                             ,READ_GLOBAL_SUMS                          &
                             ,WRITE_GLOBAL_SUMS                         &
                             ,NEMSIO_INPUT
!
        TYPE(ESMF_Logical) :: GLOBAL_E
!
!-----------------------------------------------------------------------
!***  Distributed memory information.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: IHALO,JHALO,MYPE,NHALO,NUM_PES            &
                             ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP        
!
        INTEGER(KIND=KINT) :: ITS,ITE,JTS,JTE                           &
                             ,IMS,IME,JMS,JME                           &
                             ,IDS,IDE,JDS,JDE
!
        INTEGER(KIND=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART         &
                                                  ,LOCAL_IEND           &
                                                  ,LOCAL_JSTART         &
                                                  ,LOCAL_JEND
!
!-----------------------------------------------------------------------
!***  Horizontal and vertical grid-related variables.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: LPT2
!
        INTEGER(KIND=KINT),DIMENSION(:),POINTER :: KHFILT,KVFILT        &
                                                  ,NFFTRH,NFFTRW        &
                                                  ,NHSMUD
!
        REAL(KIND=KFPT) :: DPHD,DLMD                                    &
                          ,DYH,DYV,RDYH,RDYV                            &
                          ,DDMPV                                        &
                          ,EF4T                                         &
                          ,PDTOP,PT
!
        REAL(KIND=KFPT),DIMENSION(:),POINTER :: SG1,PSG1                &
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
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: F,FIS                 &
                                                 ,GLAT,GLON             &
                                                 ,HDACX,HDACY           &
                                                 ,HDACVX,HDACVY         &
                                                 ,HFILT,VFILT           &
                                                 ,VLAT,VLON
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: PD,PDO,PSDT,SICE,SM
!
!-----------------------------------------------------------------------
!***  Integration quantities.
!-----------------------------------------------------------------------
!
        LOGICAL(KIND=KLOG) :: FIRST,READBC                              &
                             ,ADV_STANDARD,ADV_UPSTREAM
!
        INTEGER(KIND=KINT) :: NTSD,IDTAD,IDTADT,IHR,IHRST,IHREND        &
                             ,LNSAD,NBOCO,NTSTI,NTSTM,NTSTM_MAX
!
        INTEGER(KIND=KINT),DIMENSION(3) :: IDAT
!
        INTEGER(KIND=KINT),DIMENSION(:,:),POINTER :: INSOIL,INVEG
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: PINT,RTOP           &
                                                   ,T,U,V               &
                                                   ,Q,CW,RRW            &
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
        LOGICAL(KIND=KLOG) :: RUN
!!!     TYPE(ESMF_Logical) :: RUN
!
!-----------------------------------------------------------------------
!***  The general 4-D arrays for 3-D "tracers".
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: NUM_TRACERS_TOTAL                            !<-- Total number of "tracer" variables.
!
!-----------------------------------------------
!***  Declare indices of meteorological tracers
!-----------------------------------------------
!
        INTEGER(KIND=KINT) :: INDX_Q                                    &  !<-- Location of Q in tracer arrays
                             ,INDX_CW                                   &  !<-- Location of CW in tracer arrays
                             ,INDX_RRW                                  &  !<-- Location of RRW in tracer arrays
                             ,INDX_Q2                                      !<-- Location of Q2 in tracer arrays
!
        REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: TRACERS           &  !<-- Storage array for "tracer" variables.
                                                     ,TRACERS_SQRT      &  !<-- Sqrt of the tracer variables (for advection)
                                                     ,TRACERS_PREV      &  !<-- Values of tracer variables in prev timestep (for advection)
                                                     ,TRACERS_TEND         !<-- Tendencies of tracer variables (for advection)
!
!-----------------------------------------------------------------------
!***  Boundary conditions.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: IHRSTBC,LNSH,LNSV
!
        INTEGER(KIND=KINT),DIMENSION(3) :: IDATBC
!
        REAL(KIND=KFPT) :: TBOCO
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: PDBS,PDBN           &
                                                   ,PDBW,PDBE
!
        REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: TBS,TBN           & 
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
        LOGICAL(KIND=KLOG) :: RUNBC
!
!----------------------------
!***  For 1-D restart output
!----------------------------
!
!
      INTEGER :: NUM_WORDS_BC_SOUTH                                     &  !<-- Word counts of 1-D boundary data strings
                ,NUM_WORDS_BC_NORTH                                     &  !    for each side of the domain.
                ,NUM_WORDS_BC_WEST                                      &  !
                ,NUM_WORDS_BC_EAST                                         !<--
!
      REAL,DIMENSION(:),ALLOCATABLE :: RST_BC_DATA_SOUTH                &  !<-- 1-D strings of boundary data
                                      ,RST_BC_DATA_NORTH                &  !    for each side of the domain.
                                      ,RST_BC_DATA_WEST                 &  !
                                      ,RST_BC_DATA_EAST                    !<--
!
!-----------------------------------------------------------------------
!***  FFT arrays.
!-----------------------------------------------------------------------
!
        REAL(KIND=KDBL),DIMENSION(:),POINTER :: CRAUX1,CRAUX2,CRAUX3    &
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
        INTEGER(KIND=KINT) :: NUM_WATER                                 &  !<-- 1 + types of water substance in microphysics
                             ,P_QV                                      &  !<-- Index for water vapor in WATER array
                             ,P_QC                                      &  !<-- Index for cloud water in WATER array
                             ,P_QR                                      &  !<-- Index for rain in WATER array
                             ,P_QI                                      &  !<-- Index for cloud ice in WATER array
                             ,P_QS                                      &  !<-- Index for snow in WATER array
                             ,P_QG                                         !<-- Index for graupel in WATER array
!
        INTEGER(KIND=KINT) :: INDX_WATER_START                          &  !<-- Start index of the water in tracers array
                             ,INDX_WATER_END                               !<-- End index of the water in tracers array
!
        REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: WATER                !<-- Storage array for water substance
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: F_ICE,F_RAIN        &  !<-- Fractions of ice, rain, and rime
                                                   ,F_RIMEF
!
        LOGICAL :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
!-----------------------------------------------------------------------
!***  Nesting
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: I_PAR_STA                                 &  !<-- SW corner of nest domain on this parent I
                             ,J_PAR_STA                                 &  !<-- SW corner of nest domain on this parent J
                             ,PARENT_CHILD_TIME_RATIO                      !<-- # of child timesteps per parent timestep
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END TYPE DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
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
      SUBROUTINE SET_INTERNAL_STATE_DYN(GRID_COMP,INT_STATE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL,ONLY: IDS,IDE,JDS,JDE                      &
                                  ,IMS,IME,JMS,JME                      &
                                  ,IHALO,JHALO
      USE MODULE_CONSTANTS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)          ,INTENT(IN)    :: GRID_COMP             !<-- The Dynamics gridded component
      TYPE(DYNAMICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- The Dynamics internal state
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: LNSH_BC=1
!
      INTEGER(KIND=KINT) :: I,I_CYCLE,J,KS                              &
                           ,L,LL,LM,LNSH,LNSV,N,NUM_PES
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  WE CAN RETRIEVE LM FROM THE INTERNAL STATE SINCE IT WAS
!***  PLACED THERE ALREADY FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
      LM=int_state%LM
!
      I_CYCLE=IDE-3
!
      LNSH=LNSH_BC
      int_state%LNSH=LNSH
      int_state%LNSV=LNSH
      LNSV=int_state%LNSV
!
!!!   int_state%NBOCO=NBOCO   !<-- Set later in call to CONSTS
!!!   int_state%TBOCO=TBOCO   !<-- Set later in call to CONSTS
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS OF THE DYNAMICS INTERNAL STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  GRID-RELATED CONSTANTS.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DSG1  (1:LM))                                    !<-- Thicknesses of sigma layers in pressure range
      ALLOCATE(int_state%DSG2  (1:LM))                                    !<-- Thicknesses of sigma layers in sigma range
      ALLOCATE(int_state%PDSG1 (1:LM))                                    !<-- Thicknesses of pressure layers in press. range
      ALLOCATE(int_state%PSGML1(1:LM))                                    !<-- Pressure at midlayers in pressure range
      ALLOCATE(int_state%SGML1 (1:LM))                                    !<-- Sigma at midlayers in pressure range
      ALLOCATE(int_state%SGML2 (1:LM))                                    !<-- Sigma at midlayers in sigma range
!   
      ALLOCATE(int_state%PSG1(1:LM+1))                                    !<-- Pressure at interfaces in pressure range  (Pa)
      ALLOCATE(int_state%SG1 (1:LM+1))                                    !<-- Sigma at interfaces in pressure range
      ALLOCATE(int_state%SG2 (1:LM+1))                                    !<-- Sigma at interfaces in sigma range
      ALLOCATE(int_state%SGM (1:LM+1))                                    !<-- Sigma at interfaces
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
      ALLOCATE(int_state%DXH  (JDS:JDE))                                  !<-- Delta x, h points  (m)
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
      ALLOCATE(int_state%GLAT  (IMS:IME,JMS:JME))                         !<-- Latitudes of h points  (radians)
      ALLOCATE(int_state%GLON  (IMS:IME,JMS:JME))                         !<-- Longitudes of h points (radians)
      ALLOCATE(int_state%VLAT  (IMS:IME,JMS:JME))                         !<-- Latitudes of V points  (radians)
      ALLOCATE(int_state%VLON  (IMS:IME,JMS:JME))                         !<-- Longitudes of V points (radians)
      ALLOCATE(int_state%HDACX (IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, h points  (s m-1)
      ALLOCATE(int_state%HDACY (IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, h points  (s m-1)
      ALLOCATE(int_state%HDACVX(IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, v points  (s m-1)
      ALLOCATE(int_state%HDACVY(IMS:IME,JMS:JME))                         !<-- Lateral diffusion coefficient, v points  (s m-1)
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
!-----------------------------------------------------------------------
!***  LOCAL HORIZONTAL SUBDOMAIN LIMITS FOR ALL FORECAST TASKS.
!-----------------------------------------------------------------------
!
      NUM_PES=int_state%NUM_PES
      ALLOCATE(int_state%LOCAL_ISTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_IEND  (0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JSTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JEND  (0:NUM_PES-1))
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
!***  FIXED SURFACE FIELDS
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%FIS(IMS:IME,JMS:JME))                           !<-- Surface geopotential (m2 s-2)
      ALLOCATE(int_state%SM(IMS:IME,JMS:JME))                            !<-- Sea mask
      ALLOCATE(int_state%SICE(IMS:IME,JMS:JME))                          !<-- Sea ice
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=-1.E6
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  REGIONAL BOUNDARY CONDITIONS.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%UBN(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- U wind component at northern boundary  (m s-1)
      ALLOCATE(int_state%UBS(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- U wind component at southern boundary  (m s-1)
      ALLOCATE(int_state%VBN(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- V wind component at northern boundary  (m s-1)
      ALLOCATE(int_state%VBS(IMS:IME,1:LNSV,1:LM,1:2))                     !<-- V wind component at southern boundary  (m s-1)
!
      DO N=1,2
      DO L=1,LM
      DO LL=1,LNSV
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
      ALLOCATE(int_state%UBE(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- U wind component at eastern boundary  (m s-1)
      ALLOCATE(int_state%UBW(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- U wind component at western boundary  (m s-1)
      ALLOCATE(int_state%VBE(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- V wind component at eastern boundary  (m s-1)
      ALLOCATE(int_state%VBW(1:LNSV,JMS:JME,1:LM,1:2))                     !<-- V wind component at western boundary  (m s-1)
!
      DO N=1,2
      DO L=1,LM
      DO J=JMS,JME
      DO LL=1,LNSV
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
        ALLOCATE(int_state%PDBN(IMS:IME,1:LNSH,1:2))                       !<-- Pressure difference at northern boundary  (Pa)
        ALLOCATE(int_state%PDBS(IMS:IME,1:LNSH,1:2))                       !<-- Pressure difference at southern boundary  (Pa)
!
        DO N=1,2
        DO LL=1,LNSH
        DO I=IMS,IME
          int_state%PDBN(I,LL,N)=0.
          int_state%PDBS(I,LL,N)=0.
        ENDDO
        ENDDO
        ENDDO
!
        ALLOCATE(int_state%PDBE(1:LNSH,JMS:JME,1:2))                       !<-- Pressure difference at eastern boundary  (Pa)
        ALLOCATE(int_state%PDBW(1:LNSH,JMS:JME,1:2))                       !<-- Pressure difference at western boundary  (Pa)
!
        DO N=1,2
        DO J=JMS,JME
        DO LL=1,LNSH
          int_state%PDBE(LL,J,N)=0.
          int_state%PDBW(LL,J,N)=0.
        ENDDO
        ENDDO
        ENDDO
!
        ALLOCATE(int_state%QBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Specific humidity at northern boundary  (kg kg-1)
        ALLOCATE(int_state%QBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Specific humidity at southern boundary  (kg kg-1)
        ALLOCATE(int_state%TBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Temperature at northern boundary  (K)
        ALLOCATE(int_state%TBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Temperature at southern boundary  (K)
        ALLOCATE(int_state%WBN(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Vertical velocity at northern boundary  (m s-1)
        ALLOCATE(int_state%WBS(IMS:IME,1:LNSH,1:LM,1:2))                   !<-- Vertical velocity at southern boundary  (m s-1)
!
        DO N=1,2
        DO L=1,LM
        DO LL=1,LNSH
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
        ALLOCATE(int_state%QBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Specific humidity at eastern boundary  (kg kg-1)
        ALLOCATE(int_state%QBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Specific humidity at western boundary  (kg kg-1)
        ALLOCATE(int_state%TBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Temperature at eastern boundary  (K)
        ALLOCATE(int_state%TBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Temperature at western boundary  (K)
        ALLOCATE(int_state%WBE(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Vertical velocity at eastern boundary  (m s-1)
        ALLOCATE(int_state%WBW(1:LNSH,JMS:JME,1:LM,1:2))                   !<-- Vertical velocity at western boundary  (m s-1)
!
        DO N=1,2
        DO L=1,LM
        DO J=JMS,JME
        DO LL=1,LNSH
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
!-----------------------------------------------------------------------
!***  ATMOSPHERIC VARIABLES, HYDROSTATIC (mostly)
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%PD  (IMS:IME,JMS:JME))                            !<-- Pressure difference, sigma range  (Pa)
      ALLOCATE(int_state%PDO (IMS:IME,JMS:JME))                            !<-- Previous pressure difference, sigma range  (Pa)
      ALLOCATE(int_state%PSDT(IMS:IME,JMS:JME))                            !<-- Hydrostatic surface pressure tendency  (Pa s-1)
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
        int_state%PDO(I,J)=0.
        int_state%PSDT(I,J)=-1.E6
      ENDDO
      ENDDO
!
      ALLOCATE(int_state%PINT(IMS:IME,JMS:JME,1:LM+1))                     !<-- Nonhydrostatic interface pressure  (Pa)
      ALLOCATE(int_state%PSGDT(IMS:IME,JMS:JME,1:LM-1))                    !<-- Specific volume  (m3 kg-1)
!
      DO L=1,LM-1
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PSGDT(I,J,L)=0.
      ENDDO
      ENDDO
      ENDDO
!
      ALLOCATE(int_state%OMGALF(IMS:IME,JMS:JME,1:LM))                     !<-- Omega-alpha  (K)
      ALLOCATE(int_state%Q2    (IMS:IME,JMS:JME,1:LM))                     !<-- 2*tke  (m2 s-2)
      ALLOCATE(int_state%T     (IMS:IME,JMS:JME,1:LM))                     !<-- Sensible temperature  (K)
      ALLOCATE(int_state%TP    (IMS:IME,JMS:JME,1:LM))                     !<-- Sensible temperature, previous step  (K)
      ALLOCATE(int_state%U     (IMS:IME,JMS:JME,1:LM))                     !<-- U wind component  (m s-1)
      ALLOCATE(int_state%UP    (IMS:IME,JMS:JME,1:LM))                     !<-- U wind component, previous step  (m s-1)
      ALLOCATE(int_state%V     (IMS:IME,JMS:JME,1:LM))                     !<-- V wind component  (m s-1)
      ALLOCATE(int_state%VP    (IMS:IME,JMS:JME,1:LM))                     !<-- V wind component, previous step  (m s-1)
!
!-----------------------------------------------------------------------
!***  THE ARRAY CALLED WATER IS A SPECIAL CASE NEEDED TO SATISFY
!***  VARIOUS WRF PHYSICS OPTIONS.  THE IS SET TO 1+Number_of_species
!***  INCLUDING VAPOR.  THE "1+" IS NEEDED BECAUSE WRF NEVER TOUCHES
!***  THE FIRST LEVEL.
!
!***  SET THE P_ and F_ VARIABLES.
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
      ENDIF
!
      ALLOCATE(int_state%F_ICE(IMS:IME,JMS:JME,1:LM))                        !<-- Fraction of condensate as ice
      ALLOCATE(int_state%F_RAIN(IMS:IME,JMS:JME,1:LM))                       !<-- Fraction of condensate as rain
      ALLOCATE(int_state%F_RIMEF(IMS:IME,JMS:JME,1:LM))                      !<-- Fraction of condensate as rime
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
!-----------------------------------------------------------------------
!***  THE TRACERS ARRAY HOLDS ALL GENERAL "TRACER" VARIABLES INCLUDING
!***  WATER.  PLACE THE DESIRED VARIABLES AT THE TOP OF THE TRACERS
!***  ARRAY, LEVEL 1 THROUGH NUM_TRACERS_MET.  ALL OTHER SCALAR VARIABLES
!***  (e.g., chemistry and aerosols) WILL FOLLOW.
!-----------------------------------------------------------------------
!
      int_state%NUM_TRACERS_TOTAL=                                      &  !<-- # of 3-D arrays in 4-D TRACERS array
                                  int_state%NUM_TRACERS_MET             &  !<-- # of water, etc. tracers specified now (see below)
                                 +int_state%NUM_TRACERS_CHEM            &  !<-- # of specified scalars (chem, aerosol, etc.)
                                 +int_state%NUM_WATER                      !<-- # of water types
!
      ALLOCATE(int_state%TRACERS     (IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- All tracer variables
      ALLOCATE(int_state%TRACERS_SQRT(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- Sqrt of tracers (for advection)
      ALLOCATE(int_state%TRACERS_PREV(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- Tracers in previous timestep (for advection)
      ALLOCATE(int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- Tendency of tracers (for advection)
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
!-----------------------------------------------------------------------
!***  POINT Q AT LEVEL 1 OF THE TRACERS ARRAY.
!-----------------------------------------------------------------------
!
      int_state%INDX_Q=1                                                    !<-- Water vapor is always in location 1 of TRACERS arrays
      int_state%Q=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_Q) !<-- Water vapor is always in location 1 of TRACERS arrays
!
!-----------------------------------------------------------------------
!***  ADDITIONAL TRACERS.
!-----------------------------------------------------------------------
!
!--------------------------------
!***  Combined cloud water array
!--------------------------------
!
      int_state%INDX_CW=2
      int_state%CW=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_CW)
!
!--------------------------------
!***  Turbulence kinetic energy
!--------------------------------
!
      int_state%INDX_Q2=3                                                 !<-- Used for the SQRT, PREV, and TEND arrays
      int_state%E2=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_Q2)
!
!--------------------------------
!***  General tracer for testing
!--------------------------------
!
      int_state%INDX_RRW=4
      int_state%RRW=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_RRW)
!
!--------------------------------
!***  Water tracers
!--------------------------------
!
      int_state%INDX_WATER_START = int_state%NUM_TRACERS_MET + int_state%NUM_TRACERS_CHEM + 1
      int_state%INDX_WATER_END = int_state%INDX_WATER_START + int_state%NUM_WATER - 1
      int_state%WATER=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_WATER_START:int_state%INDX_WATER_END)
!
!-----------------------------------------------------------------------
!***  NOTE:  Remaining scalar variables in the Tracer arrays will 
!***         begin at location NUM_TRACERS_MET+1.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!***  Atmospheric variables, nonhydrostatic
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DWDT (IMS:IME,JMS:JME,1:LM))                     !<-- Vertical acceleration, correction factor  (m s-2)
      ALLOCATE(int_state%PDWDT(IMS:IME,JMS:JME,1:LM))                     !<-- Correction factor, previous step  (m s-2)
      ALLOCATE(int_state%RTOP (IMS:IME,JMS:JME,1:LM))                     !<-- RT/P, specific volume  (m3 kg-1)
      ALLOCATE(int_state%W    (IMS:IME,JMS:JME,1:LM))                     !<-- Vertical velocity at midlayers  (m s-1)
      ALLOCATE(int_state%Z    (IMS:IME,JMS:JME,1:LM))                     !<-- Height at midlayers  (m)
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
!-----------------------------------------------------------------------
!***  Working arrays passed as arguments between subroutines.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DIV (IMS:IME,JMS:JME,1:LM))                      !<-- Horizontal mass divergence
      ALLOCATE(int_state%PCNE(IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, NE direction
      ALLOCATE(int_state%PCNW(IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, NW direction
      ALLOCATE(int_state%PCX (IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, X direction
      ALLOCATE(int_state%PCY (IMS:IME,JMS:JME,1:LM))                      !<-- 2nd term of pgf, Y direction
      ALLOCATE(int_state%PFNE(IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, NE direction
      ALLOCATE(int_state%PFNW(IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, NW direction
      ALLOCATE(int_state%PFX (IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, X direction
      ALLOCATE(int_state%PFY (IMS:IME,JMS:JME,1:LM))                      !<-- Mass flux, Y direction
      ALLOCATE(int_state%TDIV(IMS:IME,JMS:JME,1:LM))                      !<-- Integrated horizontal mass divergence
      ALLOCATE(int_state%TCT (IMS:IME,JMS:JME,1:LM))                      !<-- Time change of T (K s-1)
      ALLOCATE(int_state%TCU (IMS:IME,JMS:JME,1:LM))                      !<-- Time change of U (m s-2)
      ALLOCATE(int_state%TCV (IMS:IME,JMS:JME,1:LM))                      !<-- Time change of V (m s-2)
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
!***  Initialize nesting quantities
!-----------------------------------------------------------------------
!
      int_state%I_PAR_STA=0
      int_state%J_PAR_STA=0
!
      int_state%PARENT_CHILD_TIME_RATIO=-999
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_INTERNAL_STATE_DYN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_INTERNAL_STATE_DYN(IMP_STATE,INT_STATE)
!
!-----------------------------------------------------------------------
!***  UPDATE THE DYNAMICS INTERNAL STATE WITH THE DATA
!***  IN THE IMPORT STATE SENT FROM THE PHYSICS.
!-----------------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL,ONLY: IHALO,JHALO
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State)             ,INTENT(INOUT) :: IMP_STATE             !<-- The Dynamics import state
      TYPE(DYNAMICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- The Dynamics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,RC,RC_UPD
      INTEGER(KIND=KINT) :: IMS,IME,JMS,JME,LM
!
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
      REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D
!
      CHARACTER(20) :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  FOR EACH VARIABLE, EXTRACT ITS ESMF Field FROM THE IMPORT STATE.
!***  EACH Field CONTAINS A POINTER THAT POINTS AT THAT VARIABLE IN THE
!***  INTERNAL STATE OF THE DYNAMICS COMPONENT.
!***  EXTRACT THE POINTER FROM THE ESMF Field AND USE IT TO UPDATE
!***  THE DYNAMICS COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS
!
      IMS=int_state%IMS
      IME=int_state%IME
      JMS=int_state%JMS
      JME=int_state%JME
      LM =int_state%LM
!
!-----------------------------------------------------------------------
!***  FIRST UPDATE THE 3-D Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   T   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      FIELD_NAME='T'
      NULLIFY(HOLD_3D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Temperature Field from Dynamics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Field
                        ,itemName=FIELD_NAME                            &  !<-- Name of the Field we want
                        ,field   =HOLD_FIELD                            &  !<-- Put extracted Field here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farray   =HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO L=1,LM
        JJ=-JHALO
        DO J=JMS,JME
          II=-IHALO
          JJ=JJ+1
          DO I=IMS,IME
            II=II+1
            int_state%T(I,J,L)=HOLD_3D(II,JJ,L)                            !<-- Update Temperature
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   U   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      FIELD_NAME='U'
      NULLIFY(HOLD_3D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract U Wind Field from Dynamics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Field
                        ,itemName=FIELD_NAME                            &  !<-- Name of the Field we want
                        ,field   =HOLD_FIELD                            &  !<-- Put extracted Field here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn Update: Extract U Wind Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farray   =HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO L=1,LM
        JJ=-JHALO
        DO J=JMS,JME
          II=-IHALO
          JJ=JJ+1
          DO I=IMS,IME
            II=II+1
            int_state%U(I,J,L)=HOLD_3D(II,JJ,L)                            !<-- Update U wind component
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   V   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      FIELD_NAME='V'
      NULLIFY(HOLD_3D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract V Wind Field from Dynamics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Field
                        ,itemName=FIELD_NAME                            &  !<-- Name of the Field we want
                        ,field   =HOLD_FIELD                            &  !<-- Put extracted Field here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn Update: Extract V Wind Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farray   =HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO L=1,LM
        JJ=-JHALO
        DO J=JMS,JME
          II=-IHALO
          JJ=JJ+1
          DO I=IMS,IME
            II=II+1
            int_state%V(I,J,L)=HOLD_3D(II,JJ,L)                            !<-- Update V wind component
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - -   Q2   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      FIELD_NAME='Q2'
      NULLIFY(HOLD_3D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract TKE Field from Dynamics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Field
                        ,itemName=FIELD_NAME                            &  !<-- Name of the Field we want
                        ,field   =HOLD_FIELD                            &  !<-- Put extracted Field here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn Update: Extract TKE Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farray   =HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO L=1,LM
        JJ=-JHALO
        DO J=JMS,JME
          II=-IHALO
          JJ=JJ+1
          DO I=IMS,IME
            II=II+1
            int_state%Q2(I,J,L)=HOLD_3D(II,JJ,L)                           !<-- Update TKE
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOW UPDATE THE 4-D TRACERS FIELD.  
!***  THE NUMBER OF 3-D QUANTITIES WITHIN IT IS KNOWN THROUGH THE
!***  NUM_TRACERS_TOTAL VARIABLE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - TRACERS - - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      FIELD_NAME='TRACERS'
      NULLIFY(HOLD_4D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract 4-D Tracers Field from Dynamics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Field
                        ,itemName=FIELD_NAME                            &  !<-- Name of the Field we want
                        ,field   =HOLD_FIELD                            &  !<-- Put extracted Field here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dyn Update: Extract Tracer Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farray   =HOLD_4D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,int_state%NUM_TRACERS_TOTAL
        DO L=1,LM
          JJ=-JHALO
          DO J=JMS,JME
            II=-IHALO
            JJ=JJ+1
            DO I=IMS,IME
              II=II+1
              int_state%TRACERS(I,J,L,N)=HOLD_4D(II,JJ,L,N)                !<-- Update tracer variables
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      IF(RC_UPD==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYNAMICS UPDATE SUCCEEDED'
      ELSE
        WRITE(0,*)'DYNAMICS UPDATE FAILED RC_UPD=',RC_UPD
      ENDIF
!
!-----------------------------------------------------------------------
!
!
      END SUBROUTINE UPDATE_INTERNAL_STATE_DYN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
