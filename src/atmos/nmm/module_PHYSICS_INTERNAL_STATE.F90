!-----------------------------------------------------------------------

      MODULE MODULE_PHYSICS_INTERNAL_STATE

!-----------------------------------------------------------------------
!***  This module declares the derived datatype called INTERNAL_STATE.
!***  For now the components of this datatype will be everything needed
!***  to advance the model integration, i.e. everything that would be
!***  part of a restart file.  Specifically this will include those
!***  quantities that evolve during the integration, the namelist
!***  variables, and the grid decomposition variables.
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!                   
!   2008-07-28  Vasic - Precompute accumulation counters.
!   2009-08-11  Ferrier - Changed accumulation counters to 2D arrays (ESMF issue)
!
!-----------------------------------------------------------------------
!

      USE MODULE_INCLUDE
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE
!
      USE MODULE_LS_NOAHLSM  ,ONLY: NUM_SOIL_LAYERS,SLDPTH
      USE MODULE_MICROPHYSICS_NMM,ONLY: MICRO_RESTART
      USE MODULE_ERR_MSG     ,ONLY: ERR_MSG,MESSAGE_CHECK
      USE MODULE_VARS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: PHYSICS_INTERNAL_STATE                                  &
               ,SET_INTERNAL_STATE_PHY_1                                &
               ,SET_INTERNAL_STATE_PHY_2                                &
               ,WRAP_PHY_INT_STATE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: MAX_VARS = 200
!
      TYPE PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Begin with the namelist variables.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT), POINTER :: NSOIL
        INTEGER(kind=KINT) :: IM,JM,LM                                  &
                             ,MINUTES_HISTORY                           &
                             ,MINUTES_RESTART                           &
                             ,NUM_TRACERS_MET                           &  !<-- Number of meteorological tracers (e.g. water)
                             ,NUM_TRACERS_CHEM                          &  !<-- Number of chem/aerosol tracers
                             ,START_YEAR,START_MONTH,START_DAY          &
                             ,START_HOUR,START_MINUTE,START_SECOND
!
        INTEGER(kind=KINT), POINTER :: NPHS
        INTEGER(kind=KINT) :: DT_INT,NPRECIP,NRADL,NRADS                &
                             ,PCPHR,UCMCALL,IGBP
!
        REAL(kind=KFPT) :: DT,SBD,WBD,TPH0D,TLM0D
!
        LOGICAL :: GLOBAL,GWDFLG,HYDRO,NEMSIO_INPUT,NESTED,NHRS_UDEF    &
                  ,PCPFLG,RESTART,SPECIFIED,WRITE_PREC_ADJ              &
                  ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP,RST_OUT_00     &
                  ,SPEC_ADV                                                ! Cloud water species advection option
!
!-----------------------------------------------------------------------
!***  Distributed memory information.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: INPES                                     &  !<-- Forecast PE's in the E-W direction
                             ,JNPES                                     &  !<-- Forecast PE's in the N-S direction
                             ,MYPE                                      &  !<-- This task's ID
                             ,NUM_PES                                      !<-- Total number of forecast tasks
!
        INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE                           &
                             ,IMS,IME,JMS,JME
!
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART         &  !<-- Each task's local starting I
                                                  ,LOCAL_IEND           &  !<-- Each task's local ending I
                                                  ,LOCAL_JSTART         &  !<-- Each task's local starting J
                                                  ,LOCAL_JEND              !<-- Each task's local ending J
!
!-----------------------------------------------------------------------
!***  Horizontal/Vertical grid
!-----------------------------------------------------------------------
!
        REAL(kind=KFPT) :: DLMD,DPHD                                    &
                          ,DYH,DYV                                      &
                          ,FRES,FR,FSL,FSS                              &
                          ,PDTOP,PT                                     &
                          ,RDYH,RDYV
!
        REAL(kind=KFPT),DIMENSION(:),POINTER :: DSG2,DXH,DXV            &
                                               ,PDSG1,PSG1,PSGML1       &
                                               ,RDXH,RDXV               &
                                               ,SG1,SG2,SGM,SGML2
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: GLAT,GLON
!
!-----------------------------------------------------------------------
!***  Integration quantities.
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: NTSD                                         !<-- Internal timestep counter
!
        INTEGER(kind=KINT) :: NHRS_CLOD                                 &  !<-- Fcst hours cloud is accumulated
                             ,NHRS_HEAT                                 &  !<-- Fcst hours heating is accumulated
                             ,NHRS_PREC                                 &  !<-- Fcst hours precip is accumulated
                             ,NHRS_RDLW                                 &  !<-- Fcst hours LW radiation is accumulated
                             ,NHRS_RDSW                                 &  !<-- Fcst hours SW radiation is accumulated
                             ,NHRS_SRFC                                    !<-- Fcst hours sfc evap/flux is accumulated
!
        INTEGER(kind=KINT), POINTER :: NCLOD                            &  !<-- # of fundamental timesteps cloud is accumulated
                                      ,NHEAT                            &  !<-- # of fundamental timesteps latent heating is accumulated
                                      ,NPREC                            &  !<-- # of fundamental timesteps precip is accumulated
                                      ,NRDLW                            &  !<-- # of fundamental timesteps LW radiation is accumulated
                                      ,NRDSW                            &  !<-- # of fundamental timesteps SW radiation is accumulated
                                      ,NSRFC                            &  !<-- # of fundamental timesteps sfc evap/flux is accumulated
                                      ,AVGMAXLEN                           !<-- Fcst sec over which avg/max diags are accumulated
!
        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: ISLTYP,IVGTYP      &
                                                    ,LPBL               &
                                                    ,NCFRCV,NCFRST
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: T,U,V               &
                                                   ,DUDT,DVDT,DWDT      &
                                                   ,Q,CW                &
                                                   ,Q2,O3               &
                                                   ,OMGALF              &
                                                   ,PINT                &
                                                   ,W,Z
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: RLWTT,RSWTT
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: EXCH_H,PPTDAT
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: RQVBLTEN,RTHBLTEN   &
                                                   ,TCUCN,W0AVG
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: CLDFRA              &
                                                   ,F_ICE,F_RAIN        &
                                                   ,F_RIMEF             &
                                                   ,TRAIN,XLEN_MIX
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: SH2O,SMC,STC
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: CROT,SROT             &
                                                 ,HANGL,HANIS,HASYS     &
                                                 ,HASYSW,HASYNW,HASYW   &
                                                 ,HCNVX                 &
                                                 ,HLENNW,HLENSW         &
                                                 ,HLENW,HLENS           &
                                                 ,HSLOP,HSTDV,HZMAX
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ACFRCV,ACFRST         &
                                                 ,AKHS,AKHS_OUT         &
                                                 ,AKMS,AKMS_OUT         &
                                                 ,CFRACH,CFRACL,CFRACM  &
                                                 ,CMC,CNVBOT,CNVTOP     &
                                                 ,CPRATE,CUPPT          &
                                                 ,CZMEAN,CZEN,LSPA      &
                                                 ,DDATA,EPSR,FIS        &
                                                 ,HBOT,HBOTD,HBOTS      &
                                                 ,HTOP,HTOPD,HTOPS      &
                                                 ,MIXHT,PBLH,PD,QSH,QZ0 &
                                                 ,RLWIN,RSWIN,RSWINC    &
                                                 ,RSWOUT,RLWTOA,RSWTOA  &
                                                 ,SICE,SIGT4,SM         &
                                                 ,SST,STDH              & !zj
                                                 ,THS,THZ0,USTAR        & !zj
                                                 ,UZ0,VZ0               &
                                                 ,Z0,Z0BASE
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ALBASE,ALBEDO         &
                                                 ,ALWIN,ALWOUT,ALWTOA   &
                                                 ,ASWIN,ASWOUT,ASWTOA   &
                                                 ,BGROFF,GRNFLX         &
                                                 ,MAVAIL,MXSNAL         &
                                                 ,POTEVP,POTFLX         &
                                                 ,QCG,QSG,QVG,QWBS      &
                                                 ,RADOT,RMOL            &
                                                 ,SFCEVP,SFCEXC         &
                                                 ,SFCLHX,SFCSHX         &
                                                 ,SHDMAX,SHDMIN         &
                                                 ,SI,SMSTAV,SMSTOT      &
                                                 ,SNO,SNOAVG,SNOPCX     &
                                                 ,SOILT1,SOILTB         &
                                                 ,SR,SSROFF,SUBSHX      &
                                                 ,TG,TSKIN,TSNAV,TWBS   &
                                                 ,VEGFRC
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ACPREC,ACSNOM,ACSNOW  &
                                                 ,CUPREC,CLDEFI         &
                                                 ,PREC,PSHLTR,P10,Q02   &
                                                 ,Q10,QSHLTR,PSFC       &
                                                 ,PSFCAVG               &
                                                 ,T2,TH02,TH10,TSHLTR   &
                                                 ,T02MAX,T02MIN         &
                                                 ,RH02MAX,RH02MIN       &
                                                 ,U10,V10,T10,T10AVG    &
                                                 ,U10MAX,V10MAX,SPD10MAX &
                                                 ,TLMIN,TLMAX           &
                                                 ,UPVVELMAX,DNVVELMAX   &
                                                 ,UPHLMAX,REFDMAX       &
                                                 ,AKHSAVG,AKMSAVG
!
!***  The following are 2-D arrays needed only to hold scalars.
!***  This is done since ESMF does not permit scalar Attributes
!***  to evolve.
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ACUTIM                &  !<-- Counter for cloud processes, used by post0 code
                                                 ,APHTIM                &  !<-- Counter for other processes, used by post0 code
                                                 ,ARDLW                 &  !<-- Counter in summing LW radiation flux
                                                 ,ARDSW                 &  !<-- Counter in summing SW radiation flux
                                                 ,ASRFC                 &  !<-- Counter in summing sfc flux
                                                 ,AVRAIN                &  !<-- Counter in summing latent heating from grid microphysics
                                                 ,AVCNVC                   !<-- Counter in summing latent heating from convection
!
        REAL(kind=KFPT),DIMENSION(:),POINTER :: MP_RESTART_STATE        &
                                               ,SLDPTH                  &
                                               ,TBPVS_STATE,TBPVS0_STATE
!
!-----------------------------------------------------------------------
!***  GFS physics additional arrays
!-----------------------------------------------------------------------
!
        REAL(kind=KDBL)                              :: CDEC,SDEC       &
                                                       ,SLAG,SOLCON
        INTEGER        ,DIMENSION(:)      ,POINTER   :: JINDX1,JINDX2
        REAL(kind=KDBL),DIMENSION(:)      ,POINTER   :: DDY
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: TMPMIN,TMPMAX
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: DUGWD,DVGWD
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: SEMIS,SFALB     &
                                                       ,SFCDLW,SFCDSW   &
                                                       ,SFCNSW,TSFLW 
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: SICFCS,SIHFCS   &
                                                       ,SLPFCS,SOTFCS   &
                                                       ,TG3FCS          &
                                                       ,VEGFCS,VETFCS   &
                                                       ,ZORFCS
        REAL(kind=KDBL),DIMENSION(:,:,:)  ,POINTER   :: ALBFC1,ALFFC1
        REAL(kind=KDBL),DIMENSION(:,:,:)  ,POINTER   :: PHY_F2DV   ! save last time step 2Ds
        REAL(kind=KDBL),DIMENSION(:,:,:,:),POINTER   :: PHY_F3DV   ! save last time step 3Ds
        REAL(kind=KDBL),DIMENSION(:,:,:,:),POINTER   :: OZPLIN
!_-----------------------------------------------------------------------------
!***  gfs microphysics additional arrays saving surface pressure, Temperature,water vapor
!     at previous time steps, Weiguo Wang 11-22-2010
!-------------------------------------------------------------------------------
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: TP1,QP1
        REAL(kind=KFPT),DIMENSION(:,:),  POINTER :: PSP1
!
!-----------------------------------------------------------------------
!***  Physics Options
!-----------------------------------------------------------------------
!
        CHARACTER(99) :: CONVECTION,MICROPHYSICS                        &
                        ,LONGWAVE,SHORTWAVE                             &
                        ,SFC_LAYER,TURBULENCE,LAND_SURFACE
!
        LOGICAL :: GFS
!
        INTEGER :: CO2TF
!
!-----------------------------------------------------------------------
!***  Microphysics Indices and Water Substance Storage Array
!-----------------------------------------------------------------------
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
!
        LOGICAL :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
!-----------------------------------------------------------------------
!***  The general 4-D array of 3-D "tracers".
!-----------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: NUM_TRACERS_TOTAL                            !<-- Total number of "tracer" variables.
!
        INTEGER(kind=KINT) :: INDX_Q                                    &  !<-- Location of Q in tracer arrays
                             ,INDX_CW                                   &  !<-- Location of CW in tracer arrays
                             ,INDX_O3                                   &  !<-- Location of O3 in tracer arrays
                             ,INDX_Q2                                      !<-- Location of Q2 in tracer arrays
!
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: TRACERS              !<-- Storage array for "tracer" variables.
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
      END TYPE PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  The internal_state type is supported by a C pointer (not an F90
!***  pointer) and therefore the following type is needed.
!-----------------------------------------------------------------------
!
      TYPE WRAP_PHY_INT_STATE
        TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE
      END TYPE WRAP_PHY_INT_STATE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_INTERNAL_STATE_PHY_1(INT_STATE,LM)
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
      USE MODULE_CONSTANTS
!
!-----------------------------------------------------------------------
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
      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- The Physics internal state
      INTEGER, INTENT(IN) :: LM
!      
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: I,J,L,N,NSTEPS_PER_HOUR,NUM_PES
      INTEGER :: LATSOZP,TIMEOZ,LEVOZP,PL_COEFF,KOZPL=28
!
!-----------------------------------------------------------------------
!***********************************************************************
!

!-----------------------------------------------------------------------
!***  The array called WATER is a special case needed to satisfy
!***  various WRF physics options.  The 4th index limit is set to
!***  1+Number_of_species INCLUDING VAPOR.  The "1+" is needed because
!***  WRF never touches the first level.
!
!***  Set the P_ AND F_ variables.
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
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='fer_hires')THEN
        int_state%NUM_WATER=1+4
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=5
        int_state%P_QI=1
        int_state%P_QG=1
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='wsm3')THEN
        int_state%NUM_WATER=1+3
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=1
        int_state%P_QI=1
        int_state%P_QG=1
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='wsm6')THEN
        int_state%NUM_WATER=1+6
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=4
        int_state%P_QS=5
        int_state%P_QI=6
        int_state%P_QG=7
      ELSEIF(TRIM(int_state%MICROPHYSICS)=='gfs')THEN
        int_state%NUM_WATER=1+3
        int_state%P_QV=2
        int_state%P_QC=3
        int_state%P_QR=1
        int_state%P_QS=1
        int_state%P_QI=4
        int_state%P_QG=1
      ENDIF
!
!-----------------------------------------------------------------------
!***  Moved the initiation of the int_state%F_Q logicals to 
!***  be defined by a series of IF tests below.
!-----------------------------------------------------------------------
!
      IF(int_state%P_QV <= 1) THEN
         int_state%F_QV=.FALSE.
      ELSE
         int_state%F_QV=.TRUE.
      ENDIF

      IF(int_state%P_QC <= 1) THEN
         int_state%F_QC=.FALSE.
      ELSE
         int_state%F_QC=.TRUE.
      ENDIF

      IF(int_state%P_QR <= 1) THEN
         int_state%F_QR=.FALSE.
      ELSE
         int_state%F_QR=.TRUE.
      ENDIF

      IF(int_state%P_QS <= 1) THEN
         int_state%F_QS=.FALSE.
      ELSE
         int_state%F_QS=.TRUE.
      ENDIF

      IF(int_state%P_QI <= 1) THEN
         int_state%F_QI=.FALSE.
      ELSE
         int_state%F_QI=.TRUE.
      ENDIF

      IF(int_state%P_QG <= 1) THEN
         int_state%F_QG=.FALSE.
      ELSE
         int_state%F_QG=.TRUE.
      ENDIF

!-----
      int_state%NUM_TRACERS_TOTAL=                                      &  !<-- # of 3-D arrays in 4-D TRACERS array
                                  int_state%NUM_TRACERS_MET             &  !<-- # of meteorological tracers such as water (see below)
                                 +int_state%NUM_TRACERS_CHEM            &  !<-- # of specified scalars (chem, aerosol, etc.)
                                 +int_state%NUM_WATER                      !<-- # of water types

!------------------------------------------------
!***  Read and store the specifications for each
!***  internal state variable listed by the user
!***  in the Physics text file.
!------------------------------------------------

      CALL READ_CONFIG('phy_state.txt',int_state%VARS,int_state%NUM_VARS)

!------------------------------------------------------------------
!***  Allocate appropriate memory within the Physics' composite
!***  VARS array for all internal state variables that are 'Owned'
!***  by Physics and point those variables into that memory.
!***  In this step, all unowned variables will be pointed at NULL
!***  for the moment.
!------------------------------------------------------------------

      CALL SET_PHY_VAR_PTR(INT_STATE, .TRUE., LM)
!
!-----------------------------------------------------------------------
!***  The array memory limits need to be set in the internal state.
!-----------------------------------------------------------------------
!
      int_state%IMS=IMS
      int_state%IME=IME
      int_state%JMS=JMS
      int_state%JME=JME
!
!-----------------------------------------------------------------------
!***  Initialize
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the arrays of the internal state.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DSG2(1:LM))            ;int_state%DSG2   = R4_IN ! Delta sigma (bottom domain)
      ALLOCATE(int_state%PDSG1(1:LM))           ;int_state%PDSG1  = R4_IN ! Delta pressure (top domain)
      ALLOCATE(int_state%PSGML1(1:LM))          ;int_state%PSGML1 = R4_IN ! Midlayer pressure (top domain)
      ALLOCATE(int_state%SGML2(1:LM))           ;int_state%SGML2  = R4_IN ! Midlayer sigma (bottom domain)
!
      ALLOCATE(INT_STATE%SG1(1:LM+1))           ;int_state%SG1    = R4_IN ! First hybrid component
      ALLOCATE(INT_STATE%PSG1(1:LM+1))          ;int_state%PSG1   = R4_IN ! First hybrid component press.
      ALLOCATE(INT_STATE%SG2(1:LM+1))           ;int_state%SG2    = R4_IN ! Second hybrid component
      ALLOCATE(INT_STATE%SGM(1:LM+1))           ;int_state%SGM    = R4_IN ! Reference sigma
!
      ALLOCATE(int_state%RDXH(JDS:JDE))         ;int_state%RDXH   = R4_IN ! 1./DX for H point rows !zj
      ALLOCATE(int_state%RDXV(JDS:JDE))         ;int_state%RDXV   = R4_IN ! 1./DX for V point rows !zj
!
      ALLOCATE(int_state%DXH(JDS:JDE))          ;int_state%DXH    = R4_IN ! DX for H point rows
      ALLOCATE(int_state%DXV(JDS:JDE))          ;int_state%DXV    = R4_IN ! DX for V point rows
!
!-----------------------------------------------------------------------
!***  Local horizontal subdomain limits for all forecast tasks.
!-----------------------------------------------------------------------
!
      NUM_PES=int_state%NUM_PES
      ALLOCATE(int_state%LOCAL_ISTART(0:NUM_PES-1)) ;int_state%LOCAL_ISTART = I4_IN 
      ALLOCATE(int_state%LOCAL_IEND  (0:NUM_PES-1)) ;int_state%LOCAL_IEND   = I4_IN 
      ALLOCATE(int_state%LOCAL_JSTART(0:NUM_PES-1)) ;int_state%LOCAL_JSTART = I4_IN 
      ALLOCATE(int_state%LOCAL_JEND  (0:NUM_PES-1)) ;int_state%LOCAL_JEND   = I4_IN 
!
!-----------------------------------------------------------------------
!***  Prognostic arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DUDT(IMS:IME,JMS:JME,1:LM))       ;int_state%DUDT     = R4_IN ! U wind component tendency  (m s-2)
      ALLOCATE(int_state%DVDT(IMS:IME,JMS:JME,1:LM))       ;int_state%DVDT     = R4_IN ! V wind component tendency  (m s-2)
!
      ALLOCATE(int_state%RQVBLTEN(IMS:IME,JMS:JME,1:LM+1)) ;int_state%RQVBLTEN = R4_IN ! Specific humidity tendency from turbulence  (kg kg-1 s-1)
      ALLOCATE(int_state%RTHBLTEN(IMS:IME,JMS:JME,1:LM+1)) ;int_state%RTHBLTEN = R4_IN ! Theta tendency from turbulence  (K s-1)
!
      ALLOCATE(int_state%W0AVG(IMS:IME,1:LM+1,JMS:JME))    ;int_state%W0AVG    = R4_IN ! Time-averaged vertical velocity (for K-F)  (m s-1)
!
      ALLOCATE(int_state%PPTDAT(IMS:IME,JMS:JME,1:int_state%PCPHR)) ;int_state%PPTDAT = R4_IN 

!-----------------------------------------------------------------------
!*** gfs microphysics, wang, 11-22-2010
!-----------------------------------------------------------------------
        ALLOCATE(int_state%TP1(IMS:IME,JMS:JME,1:LM))
        ALLOCATE(int_state%QP1(IMS:IME,JMS:JME,1:LM))
        ALLOCATE(int_state%PSP1(IMS:IME,JMS:JME))
        DO I=IMS,IME
        DO J=JMS,JME
          int_state%PSP1(I,J) = -999.0
          DO L=1,LM
           int_state%TP1(I,J,L) = -999.0
           int_state%QP1(I,J,L) = -999.0
          ENDDO
        ENDDO
        ENDDO
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Only for GFS physics
!-----------------------------------------------------------------------
!
      IF ( int_state%GFS ) THEN
        REWIND (KOZPL)
        READ (KOZPL) PL_COEFF, LATSOZP, LEVOZP, TIMEOZ
        ALLOCATE(int_state%OZPLIN(LATSOZP,LEVOZP,PL_COEFF,TIMEOZ)) ;int_state%OZPLIN = R8_IN 
      ENDIF
!
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%LPBL(IMS:IME,JMS:JME))   ;int_state%LPBL   = I4_IN ! Model layer containing top of the PBL
      ALLOCATE(int_state%DDATA(IMS:IME,JMS:JME))  ;int_state%DDATA  = R4_IN ! Observed precip to each physics timestep (kg m-2)
      ALLOCATE(int_state%MAVAIL(IMS:IME,JMS:JME)) ;int_state%MAVAIL = R4_IN ! Moisture availability
      ALLOCATE(int_state%QCG(IMS:IME,JMS:JME))    ;int_state%QCG    = R4_IN ! Cloud water mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%QSG(IMS:IME,JMS:JME))    ;int_state%QSG    = R4_IN ! Surface saturation water vapor mixing ratio  (kg kg-1)
      ALLOCATE(int_state%QVG(IMS:IME,JMS:JME))    ;int_state%QVG    = R4_IN ! Water vapor mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%SHDMAX(IMS:IME,JMS:JME)) ;int_state%SHDMAX = R4_IN ! Maximum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SHDMIN(IMS:IME,JMS:JME)) ;int_state%SHDMIN = R4_IN ! Minimum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SOILT1(IMS:IME,JMS:JME)) ;int_state%SOILT1 = R4_IN ! Snow temperature  (K)
      ALLOCATE(int_state%STDH(IMS:IME,JMS:JME))   ;int_state%STDH   = R4_IN ! Standard deviation of topography height (m) !zj
      ALLOCATE(int_state%TSNAV(IMS:IME,JMS:JME))  ;int_state%TSNAV  = R4_IN ! Average snow temperature  (K)
      ALLOCATE(int_state%CROT(IMS:IME,JMS:JME))   ;int_state%CROT   = R4_IN ! Cosine of the angle between Earth and model coordinates
      ALLOCATE(int_state%SROT(IMS:IME,JMS:JME))   ;int_state%SROT   = R4_IN ! Sine of the angle between Earth and model coordinates
      ALLOCATE(int_state%HSTDV(IMS:IME,JMS:JME))  ;int_state%HSTDV  = R4_IN ! Standard deviation of the height (m)
      ALLOCATE(int_state%HCNVX(IMS:IME,JMS:JME))  ;int_state%HCNVX  = R4_IN ! Orographic convexity
      ALLOCATE(int_state%HASYW(IMS:IME,JMS:JME))  ;int_state%HASYW  = R4_IN ! Orographic asymmetry, west wind direction
      ALLOCATE(int_state%HASYS(IMS:IME,JMS:JME))  ;int_state%HASYS  = R4_IN ! Orographic asymmetry, south wind direction
      ALLOCATE(int_state%HASYSW(IMS:IME,JMS:JME)) ;int_state%HASYSW = R4_IN ! Orographic asymmetry, southwest wind direction
      ALLOCATE(int_state%HASYNW(IMS:IME,JMS:JME)) ;int_state%HASYNW = R4_IN ! Orographic asymmetry, northwest wind direction
      ALLOCATE(int_state%HLENW(IMS:IME,JMS:JME))  ;int_state%HLENW  = R4_IN ! Orographic length scale, west wind direction
      ALLOCATE(int_state%HLENS(IMS:IME,JMS:JME))  ;int_state%HLENS  = R4_IN ! Orographic length scale, south wind direction
      ALLOCATE(int_state%HLENSW(IMS:IME,JMS:JME)) ;int_state%HLENSW = R4_IN ! Orographic length scale, southwest wind direction
      ALLOCATE(int_state%HLENNW(IMS:IME,JMS:JME)) ;int_state%HLENNW = R4_IN ! Orographic length scale, northwest wind direction
      ALLOCATE(int_state%HANGL(IMS:IME,JMS:JME))  ;int_state%HANGL  = R4_IN ! Angle of mountain range with respect to east
      ALLOCATE(int_state%HANIS(IMS:IME,JMS:JME))  ;int_state%HANIS  = R4_IN ! Anisotropy/aspect ratio
      ALLOCATE(int_state%HSLOP(IMS:IME,JMS:JME))  ;int_state%HSLOP  = R4_IN ! Slope of orography
      ALLOCATE(int_state%HZMAX(IMS:IME,JMS:JME))  ;int_state%HZMAX  = R4_IN ! Maximum height about mean terrain
      ALLOCATE(int_state%Q02(IMS:IME,JMS:JME))    ;int_state%Q02    = R4_IN ! Specific humidity at 2-m  (kg k-1)
      ALLOCATE(int_state%TH02(IMS:IME,JMS:JME))   ;int_state%TH02   = R4_IN ! Theta at 2-m  (K)
!
!-----------------------------------------------------------------------
!***  GFS physics
!-----------------------------------------------------------------------
!
      gfs_physics: IF(int_state%GFS)THEN
!
        ALLOCATE(int_state%DDY              (JTS:JTE))    ;int_state%DDY    = R8_IN     !
        ALLOCATE(int_state%JINDX1           (JTS:JTE))    ;int_state%JINDX1 = I4_IN     !
        ALLOCATE(int_state%JINDX2           (JTS:JTE))    ;int_state%JINDX2 = I4_IN     !
!
        ALLOCATE(int_state%DUGWD    (IMS:IME,JMS:JME))    ;int_state%DUGWD  = R8_IN     ! U comp. GWD tend (m s-1)
        ALLOCATE(int_state%DVGWD    (IMS:IME,JMS:JME))    ;int_state%DVGWD  = R8_IN     ! V comp. GWD tend (m s-1)
!
        ALLOCATE(int_state%TMPMIN   (IMS:IME,JMS:JME))    ;int_state%TMPMIN = R8_IN     ! Max temp (K)
        ALLOCATE(int_state%TMPMAX   (IMS:IME,JMS:JME))    ;int_state%TMPMAX = R8_IN     ! Min temp (K)
!
        ALLOCATE(int_state%SFALB    (IMS:IME,JMS:JME))    ;int_state%SFALB  = R8_IN     !
        ALLOCATE(int_state%TSFLW    (IMS:IME,JMS:JME))    ;int_state%TSFLW  = R8_IN     !
        ALLOCATE(int_state%SEMIS    (IMS:IME,JMS:JME))    ;int_state%SEMIS  = R8_IN     !
        ALLOCATE(int_state%SFCDLW   (IMS:IME,JMS:JME))    ;int_state%SFCDLW = R8_IN     !
        ALLOCATE(int_state%SFCDSW   (IMS:IME,JMS:JME))    ;int_state%SFCDSW = R8_IN     !
        ALLOCATE(int_state%SFCNSW   (IMS:IME,JMS:JME))    ;int_state%SFCNSW = R8_IN     !
!
        ALLOCATE(int_state%ZORFCS   (IMS:IME,JMS:JME))    ;int_state%ZORFCS = R8_IN     !
        ALLOCATE(int_state%SIHFCS   (IMS:IME,JMS:JME))    ;int_state%SIHFCS = R8_IN     !
        ALLOCATE(int_state%SICFCS   (IMS:IME,JMS:JME))    ;int_state%SICFCS = R8_IN     !
        ALLOCATE(int_state%SLPFCS   (IMS:IME,JMS:JME))    ;int_state%SLPFCS = R8_IN     !
        ALLOCATE(int_state%TG3FCS   (IMS:IME,JMS:JME))    ;int_state%TG3FCS = R8_IN     !
        ALLOCATE(int_state%VEGFCS   (IMS:IME,JMS:JME))    ;int_state%VEGFCS = R8_IN     !
        ALLOCATE(int_state%VETFCS   (IMS:IME,JMS:JME))    ;int_state%VETFCS = R8_IN     !
        ALLOCATE(int_state%SOTFCS   (IMS:IME,JMS:JME))    ;int_state%SOTFCS = R8_IN     !
!
        ALLOCATE(int_state%ALBFC1   (IMS:IME,JMS:JME,4))  ;int_state%ALBFC1 = R8_IN     !
        ALLOCATE(int_state%ALFFC1   (IMS:IME,JMS:JME,2))  ;int_state%ALFFC1 = R8_IN     !
!
        ALLOCATE(int_state%PHY_F2DV (IMS:IME,JMS:JME,3))    ;int_state%PHY_F2DV = R8_IN ! for Zhao =3, Ferr=1
        ALLOCATE(int_state%PHY_F3DV (IMS:IME,JMS:JME,LM,4)) ;int_state%PHY_F3DV = R8_IN ! for Zhao =4, Ferr=3
!
      ENDIF  gfs_physics
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_INTERNAL_STATE_PHY_1
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_INTERNAL_STATE_PHY_2(INT_STATE,LM)

!-----------------------------------------------------------------------
!***  Point all unallocated internal state variables into allocated
!***  memory from Dynamics.  Do fundamental initialization.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE
      INTEGER, INTENT(IN) :: LM
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: INDX, I, J, L
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  Loop through the Physics' internal state variables listed in
!***  the user's text file and re-point them into the composite VARS
!***  array.  In phase 1 of this step, all unowned variables were
!***  pointed at NULL (unallocated pointers).  Now in phase 2 the
!***  locations in VARS associated with unowned variables are pointing
!***  into allocated memory of the owned variables in Dynamics.
!-----------------------------------------------------------------------
!
      CALL SET_PHY_VAR_PTR(INT_STATE, .FALSE., LM)
!
!-----------------------------------------------------------------------
!***  Point Q at level 1 of the Tracers array.
!-----------------------------------------------------------------------
!
      int_state%INDX_Q=1
      CALL FIND_VAR_INDX('Q',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_Q)
      int_state%Q=>int_state%VARS(I)%R3D
!
!-----------------------------------------------------------------------
!***  Additional tracers:
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
      END SUBROUTINE SET_INTERNAL_STATE_PHY_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_PHY_VAR_PTR(INT_STATE,AF,LM)
!
!--------------------------------------------------------------
!***  We must allocate memory within the composite VARS array
!***  for each Physics internal state variable that is 'Owned'
!***  when ALLOC_FLAG is true in phase 1 of Physics Init.
!***  In phase 2 ALLOC_FLAG is false, no allocation takes
!***  place, and the designated internal state variables are
!***  re-pointed into VARS whose locations associated with 
!***  unowned variables will be pointing at allocated memory
!***  for the owned variable in Dynamics.
!--------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL,ONLY: IMS,IME,JMS,JME
!
!--------------------------------------------------------------
!
      IMPLICIT NONE
!
!--------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE
      LOGICAL, INTENT(IN) :: AF                                  ! ALLOC_FLAG
      INTEGER, INTENT(IN) :: LM
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: N,NV,RC
!
!--------------------------------------------------------------
!**************************************************************
!--------------------------------------------------------------

      NV=int_state%NUM_VARS

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NSOIL'      ,int_state%NSOIL ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NPHS'       ,int_state%NPHS  )  
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NCLOD'      ,int_state%NCLOD ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NHEAT'      ,int_state%NHEAT ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NPREC'      ,int_state%NPREC ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NRDLW'      ,int_state%NRDLW ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NRDSW'      ,int_state%NRDSW ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NSRFC'      ,int_state%NSRFC ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AVGMAXLEN'  ,int_state%AVGMAXLEN ) 

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ISLTYP'     ,int_state%ISLTYP   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'IVGTYP'     ,int_state%IVGTYP   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NCFRCV'     ,int_state%NCFRCV   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'NCFRST'     ,int_state%NCFRST   ,(/ IMS,JMS /),(/ IME,JME /) )  

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SLDPTH'     ,int_state%SLDPTH	      ,1 ,NUM_SOIL_LAYERS  ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'MP_RESTART' ,int_state%MP_RESTART_STATE  ,1 ,MICRO_RESTART    ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TBPVS_STAT' ,int_state%TBPVS_STATE       ,1 ,MICRO_RESTART    ) 
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TBPVS0_STA' ,int_state%TBPVS0_STATE      ,1 ,MICRO_RESTART    ) 

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACFRCV'     ,int_state%ACFRCV   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACFRST'     ,int_state%ACFRST   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACPREC'     ,int_state%ACPREC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACSNOM'     ,int_state%ACSNOM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACSNOW'     ,int_state%ACSNOW   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKHS_OUT'   ,int_state%AKHS_OUT ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKHSAVG'    ,int_state%AKHSAVG  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKMS_OUT'   ,int_state%AKMS_OUT ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKMSAVG'    ,int_state%AKMSAVG  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ALBASE'     ,int_state%ALBASE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ALBEDO'     ,int_state%ALBEDO   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ALWIN'      ,int_state%ALWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ALWOUT'     ,int_state%ALWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ALWTOA'     ,int_state%ALWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ASWIN'      ,int_state%ASWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ASWOUT'     ,int_state%ASWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ASWTOA'     ,int_state%ASWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'BGROFF'     ,int_state%BGROFF   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CFRACH'     ,int_state%CFRACH   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CFRACL'     ,int_state%CFRACL   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CFRACM'     ,int_state%CFRACM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CLDEFI'     ,int_state%CLDEFI   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CMC'        ,int_state%CMC      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CNVBOT'     ,int_state%CNVBOT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CNVTOP'     ,int_state%CNVTOP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CPRATE'     ,int_state%CPRATE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CUPPT'      ,int_state%CUPPT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CUPREC'     ,int_state%CUPREC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CZEN'       ,int_state%CZEN     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CZMEAN'     ,int_state%CZMEAN   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'DNVVELMAX'  ,int_state%DNVVELMAX,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'EPSR'       ,int_state%EPSR     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'FIS'        ,int_state%FIS      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'GLAT'       ,int_state%GLAT     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'GLON'       ,int_state%GLON     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'GRNFLX'     ,int_state%GRNFLX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HBOTD'      ,int_state%HBOTD    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HBOTS'      ,int_state%HBOTS    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HTOPD'      ,int_state%HTOPD    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HTOPS'      ,int_state%HTOPS    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'MIXHT'      ,int_state%MIXHT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'MXSNAL'     ,int_state%MXSNAL   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PBLH'       ,int_state%PBLH     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PD'         ,int_state%PD       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'POTEVP'     ,int_state%POTEVP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PREC'       ,int_state%PREC     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PSFCAVG'    ,int_state%PSFCAVG  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PSHLTR'     ,int_state%PSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'P10'        ,int_state%P10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Q10'        ,int_state%Q10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'QSH'        ,int_state%QSH      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'QSHLTR'     ,int_state%QSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'QWBS'       ,int_state%QWBS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'QZ0'        ,int_state%QZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RADOT'      ,int_state%RADOT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'REFDMAX'    ,int_state%REFDMAX  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RH02MAX'    ,int_state%RH02MAX  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RH02MIN'    ,int_state%RH02MIN  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RLWIN'      ,int_state%RLWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RLWTOA'     ,int_state%RLWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RSWIN'      ,int_state%RSWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RSWINC'     ,int_state%RSWINC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RSWOUT'     ,int_state%RSWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SFCEVP'     ,int_state%SFCEVP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SFCEXC'     ,int_state%SFCEXC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SFCLHX'     ,int_state%SFCLHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SFCSHX'     ,int_state%SFCSHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SI'         ,int_state%SI       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SICE'       ,int_state%SICE     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SIGT4'      ,int_state%SIGT4    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SM'         ,int_state%SM       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SMSTAV'     ,int_state%SMSTAV   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SMSTOT'     ,int_state%SMSTOT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SNO'        ,int_state%SNO      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SNOAVG'     ,int_state%SNOAVG   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SNOPCX'     ,int_state%SNOPCX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SOILTB'     ,int_state%SOILTB   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SPD10MAX'   ,int_state%SPD10MAX ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SR'         ,int_state%SR       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SSROFF'     ,int_state%SSROFF   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SST'        ,int_state%SST      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SUBSHX'     ,int_state%SUBSHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T02MAX'     ,int_state%T02MAX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T02MIN'     ,int_state%T02MIN   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T10'        ,int_state%T10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T10AVG'     ,int_state%T10AVG   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TG'         ,int_state%TG       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TH10'       ,int_state%TH10     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'THS'        ,int_state%THS      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'THZ0'       ,int_state%THZ0     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TSHLTR'     ,int_state%TSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TWBS'       ,int_state%TWBS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'UPHLMAX'    ,int_state%UPHLMAX  ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'UPVVELMAX'  ,int_state%UPVVELMAX,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'U10'        ,int_state%U10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'U10MAX'     ,int_state%U10MAX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'USTAR'      ,int_state%USTAR    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'UZ0'        ,int_state%UZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'V10'        ,int_state%V10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'V10MAX'     ,int_state%V10MAX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'VEGFRC'     ,int_state%VEGFRC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'VZ0'        ,int_state%VZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Z0'         ,int_state%Z0       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TSKIN'      ,int_state%TSKIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKHS'       ,int_state%AKHS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AKMS'       ,int_state%AKMS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HBOT'       ,int_state%HBOT     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'HTOP'       ,int_state%HTOP     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RSWTOA'     ,int_state%RSWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'POTFLX'     ,int_state%POTFLX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RMOL'       ,int_state%RMOL     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T2'         ,int_state%T2       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Z0BASE'     ,int_state%Z0BASE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PSFC'       ,int_state%PSFC     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TLMIN'      ,int_state%TLMIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TLMAX'      ,int_state%TLMAX    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'LSPA'       ,int_state%LSPA     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ACUTIM'     ,int_state%ACUTIM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'APHTIM'     ,int_state%APHTIM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ARDLW'      ,int_state%ARDLW    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ARDSW'      ,int_state%ARDSW    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'ASRFC'      ,int_state%ASRFC    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AVRAIN'     ,int_state%AVRAIN   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'AVCNVC'     ,int_state%AVCNVC   ,(/ IMS,JMS /),(/ IME,JME /) )

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Q2'         ,int_state%Q2       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'OMGALF'     ,int_state%OMGALF   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'T'          ,int_state%T	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'PINT'       ,int_state%PINT     ,(/ IMS,JMS,1 /),(/ IME,JME,LM+1 /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'DWDT'       ,int_state%DWDT     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Q'          ,int_state%Q	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'U'          ,int_state%U	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'V'          ,int_state%V	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'W'          ,int_state%W	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'Z'          ,int_state%Z	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RLWTT'      ,int_state%RLWTT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'RSWTT'      ,int_state%RSWTT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'EXCH_H'     ,int_state%EXCH_H   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CLDFRA'     ,int_state%CLDFRA   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'CW'         ,int_state%CW       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'F_ICE'      ,int_state%F_ICE    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'F_RAIN'     ,int_state%F_RAIN   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'F_RIMEF'    ,int_state%F_RIMEF  ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TRAIN'      ,int_state%TRAIN    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'XLEN_MIX'   ,int_state%XLEN_MIX ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TCUCN'      ,int_state%TCUCN    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SMC'        ,int_state%SMC      ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'STC'        ,int_state%STC      ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'SH2O'       ,int_state%SH2O     ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'O3'         ,int_state%O3       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))

      CALL SET_VAR_PTR(int_state%VARS,NV,AF,'TRACERS'    ,int_state%TRACERS,(/ IMS,JMS,1,1 /),(/ IME,JME,LM,int_state%NUM_TRACERS_TOTAL /) )

      DO N=1,NV
        IF (int_state%VARS(N)%TKR==0) THEN
          write(0,*)' Error in SET_PHY_VAR_PTR. '
          write(0,*)' Variable ',TRIM(int_state%VARS(N)%VBL_NAME),' is not associated to an internal state fortran pointer'
          STOP
        END IF
      END DO

      END SUBROUTINE SET_PHY_VAR_PTR
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
