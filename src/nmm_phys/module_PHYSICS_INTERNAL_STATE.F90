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
                             ,NHOURS_HISTORY                            &
                             ,NHOURS_RESTART                            &
                             ,NUM_TRACERS_MET                           &  !<-- Number of meteorological tracers (e.g. water)
                             ,NUM_TRACERS_CHEM                          &  !<-- Number of chem/aerosol tracers
                             ,START_YEAR,START_MONTH,START_DAY          &
                             ,START_HOUR,START_MINUTE,START_SECOND
!
        INTEGER(kind=KINT), POINTER :: NPHS
        INTEGER(kind=KINT) :: DT_INT,NPRECIP,NRADL,NRADS                &
                             ,PCPHR,UCMCALL
!
        REAL(kind=KFPT) :: DT,SBD,WBD,TPH0D,TLM0D
!
        LOGICAL :: GLOBAL,GWDFLG,HYDRO,NEMSIO_INPUT,NESTED,NHRS_UDEF    &
                  ,PCPFLG,RESTART,SPECIFIED,WRITE_PREC_ADJ              &
                  ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP
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
                                      ,NSRFC                               !<-- # of fundamental timesteps sfc evap/flux is accumulated
!
        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: ISLTYP,IVGTYP      &
                                                    ,LPBL               &
                                                    ,NCFRCV,NCFRST
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: T,U,V               &
                                                   ,DUDT,DVDT           &
                                                   ,Q,CW                &
                                                   ,Q2,OMGALF           &
                                                   ,RRW
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: RLWTT,RSWTT
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: EXCH_H,PPTDAT
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: RQVBLTEN,RTHBLTEN   &
                                                   ,TCUCN,W0AVG,WINT
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
                                                 ,SNO,SNOPCX            &
                                                 ,SOILT1,SOILTB         &
                                                 ,SR,SSROFF,SUBSHX      &
                                                 ,TG,TSKIN,TSNAV,TWBS   &
                                                 ,VEGFRC
!
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ACPREC,ACSNOM,ACSNOW  &
                                                 ,CUPREC,CLDEFI         &
                                                 ,PREC,PSHLTR,Q02,Q10   &
                                                 ,QSHLTR,PSFC           &
                                                 ,T2,TH02,TH10,TSHLTR   &
                                                 ,U10,V10,TLMIN,TLMAX
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
        REAL(kind=KDBL)                              :: SOLCON, SLAG, SDEC, CDEC
        INTEGER        ,DIMENSION(:)      ,POINTER   :: JINDX1,JINDX2
        REAL(kind=KDBL),DIMENSION(:)      ,POINTER   :: DDY
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: TMPMIN,TMPMAX
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: DUGWD,DVGWD
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: SFCNSW,SFCDSW,SFALB,SFCDLW   &
                                                       ,TSFLW
        REAL(kind=KDBL),DIMENSION(:,:)    ,POINTER   :: ZORFCS,SIHFCS,SICFCS,SLPFCS  &
                                                       ,TG3FCS,VEGFCS,VETFCS,SOTFCS
        REAL(kind=KDBL),DIMENSION(:,:,:)  ,POINTER   :: ALBFC1,ALFFC1
        REAL(kind=KDBL),DIMENSION(:,:,:)  ,POINTER   :: SWH,HLW
        REAL(kind=KDBL),DIMENSION(:,:,:)  ,POINTER   :: PHY_F2DV   ! save last time step 2Ds
        REAL(kind=KDBL),DIMENSION(:,:,:,:),POINTER   :: PHY_F3DV   ! save last time step 3Ds
        REAL(kind=KDBL),DIMENSION(:,:,:,:),POINTER   :: OZPLIN
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
                             ,INDX_RRW                                  &  !<-- Location of RRW in tracer arrays
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
!***  THE INTERNAL_STATE TYPE IS SUPPORTED BY A C POINTER (NOT AN F90
!***  POINTER) AND THEREFORE THE FOLLOWING TYPE IS NEEDED.
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
      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE             !<-- The Physics internal state
      INTEGER, INTENT(IN) :: LM
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,L,N,NSTEPS_PER_HOUR,NUM_PES
      INTEGER :: LATSOZP,TIMEOZ,LEVOZP,PL_COEFF,KOZPL=28
!
!-----------------------------------------------------------------------
!***********************************************************************
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
      ALLOCATE(int_state%DSG2(1:LM))     ! Delta sigma (bottom domain)
      ALLOCATE(int_state%PDSG1(1:LM))    ! Delta pressure (top domain)
      ALLOCATE(int_state%PSGML1(1:LM))   ! Midlayer pressure (top domain)
      ALLOCATE(int_state%SGML2(1:LM))    ! Midlayer sigma (bottom domain)
!
      ALLOCATE(INT_STATE%SG1(1:LM+1))    ! First hybrid component
      ALLOCATE(INT_STATE%PSG1(1:LM+1))   ! First hybrid component press.
      ALLOCATE(INT_STATE%SG2(1:LM+1))    ! Second hybrid component
      ALLOCATE(INT_STATE%SGM(1:LM+1))    ! Reference sigma
!
      ALLOCATE(int_state%RDXH(JDS:JDE))  ! 1./DX for H point rows !zj
      ALLOCATE(int_state%RDXV(JDS:JDE))  ! 1./DX for V point rows !zj
!
      ALLOCATE(int_state%DXH(JDS:JDE))   ! DX for H point rows
      ALLOCATE(int_state%DXV(JDS:JDE))   ! DX for V point rows
!
      ALLOCATE(int_state%GLAT(IMS:IME,JMS:JME))    ! Geographic latitude (radians)
      ALLOCATE(int_state%GLON(IMS:IME,JMS:JME))    ! Geographic longitude (radians, positive east)
!
!-----------------------------------------------------------------------
!***  Local horizontal subdomain limits for all forecast tasks.
!-----------------------------------------------------------------------
!
      NUM_PES=int_state%NUM_PES
      ALLOCATE(int_state%LOCAL_ISTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_IEND  (0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JSTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JEND  (0:NUM_PES-1))
!
!-----------------------------------------------------------------------
!***  Prognostic arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DUDT(IMS:IME,JMS:JME,1:LM))      ! U wind component tendency  (m s-2)
      ALLOCATE(int_state%DVDT(IMS:IME,JMS:JME,1:LM))      ! V wind component tendency  (m s-2)
!
      ALLOCATE(int_state%RQVBLTEN(IMS:IME,JMS:JME,1:LM+1)) ! Specific humidity tendency from turbulence  (kg kg-1 s-1)
      ALLOCATE(int_state%RTHBLTEN(IMS:IME,JMS:JME,1:LM+1)) ! Theta tendency from turbulence  (K s-1)
!
      ALLOCATE(int_state%W0AVG(IMS:IME,1:LM+1,JMS:JME))   ! Time-averaged vertical velocity (for K-F)  (m s-1)
      ALLOCATE(int_state%WINT(IMS:IME,1:LM+1,JMS:JME))    ! Interface vertical velocity (m s-1)
!
      ALLOCATE(int_state%PPTDAT(IMS:IME,JMS:JME,1:int_state%PCPHR))

!-----------------------------------------------------------------------
!***  Only for GFS physics
!-----------------------------------------------------------------------
!
      IF ( int_state%GFS ) THEN
        REWIND (KOZPL)
        READ (KOZPL) PL_COEFF, LATSOZP, LEVOZP, TIMEOZ
        ALLOCATE(int_state%OZPLIN(LATSOZP,LEVOZP,PL_COEFF,TIMEOZ))
      ENDIF
!
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%LPBL(IMS:IME,JMS:JME))     ! Model layer containing top of the PBL
      ALLOCATE(int_state%DDATA(IMS:IME,JMS:JME))    ! Observed precip to each physics timestep (kg m-2)
      ALLOCATE(int_state%MAVAIL(IMS:IME,JMS:JME))   ! Moisture availability
      ALLOCATE(int_state%MIXHT(IMS:IME,JMS:JME))    ! Mixed layer height  (m)
      ALLOCATE(int_state%QCG(IMS:IME,JMS:JME))      ! Cloud water mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%QSG(IMS:IME,JMS:JME))      ! Surface saturation water vapor mixing ratio  (kg kg-1)
      ALLOCATE(int_state%QVG(IMS:IME,JMS:JME))      ! Water vapor mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%SHDMAX(IMS:IME,JMS:JME))   ! Maximum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SHDMIN(IMS:IME,JMS:JME))   ! Minimum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SOILT1(IMS:IME,JMS:JME))   ! Snow temperature  (K)
      ALLOCATE(int_state%STDH(IMS:IME,JMS:JME))     ! Standard deviation of topography height (m) !zj
      ALLOCATE(int_state%TSNAV(IMS:IME,JMS:JME))    ! Average snow temperature  (K)
      ALLOCATE(int_state%CROT(IMS:IME,JMS:JME))     ! Cosine of the angle between Earth and model coordinates
      ALLOCATE(int_state%SROT(IMS:IME,JMS:JME))     ! Sine of the angle between Earth and model coordinates
      ALLOCATE(int_state%HSTDV(IMS:IME,JMS:JME))    ! Standard deviation of the height (m)
      ALLOCATE(int_state%HCNVX(IMS:IME,JMS:JME))    ! Orographic convexity
      ALLOCATE(int_state%HASYW(IMS:IME,JMS:JME))    ! Orographic asymmetry, west wind direction
      ALLOCATE(int_state%HASYS(IMS:IME,JMS:JME))    ! Orographic asymmetry, south wind direction
      ALLOCATE(int_state%HASYSW(IMS:IME,JMS:JME))   ! Orographic asymmetry, southwest wind direction
      ALLOCATE(int_state%HASYNW(IMS:IME,JMS:JME))   ! Orographic asymmetry, northwest wind direction
      ALLOCATE(int_state%HLENW(IMS:IME,JMS:JME))    ! Orographic length scale, west wind direction
      ALLOCATE(int_state%HLENS(IMS:IME,JMS:JME))    ! Orographic length scale, south wind direction
      ALLOCATE(int_state%HLENSW(IMS:IME,JMS:JME))   ! Orographic length scale, southwest wind direction
      ALLOCATE(int_state%HLENNW(IMS:IME,JMS:JME))   ! Orographic length scale, northwest wind direction
      ALLOCATE(int_state%HANGL(IMS:IME,JMS:JME))    ! Angle of mountain range with respect to east
      ALLOCATE(int_state%HANIS(IMS:IME,JMS:JME))    ! Anisotropy/aspect ratio
      ALLOCATE(int_state%HSLOP(IMS:IME,JMS:JME))    ! Slope of orography
      ALLOCATE(int_state%HZMAX(IMS:IME,JMS:JME))    ! Maximum height about mean terrain
      ALLOCATE(int_state%Q02(IMS:IME,JMS:JME))     ! Specific humidity at 2-m  (kg k-1)
      ALLOCATE(int_state%TH02(IMS:IME,JMS:JME))    ! Theta at 2-m  (K)
!
!-----------------------------------------------------------------------
!***  GFS physics
!-----------------------------------------------------------------------
!
      gfs_physics: IF(int_state%GFS)THEN
!
        ALLOCATE(int_state%DDY              (JTS:JTE))         !
        ALLOCATE(int_state%JINDX1           (JTS:JTE))         !
        ALLOCATE(int_state%JINDX2           (JTS:JTE))         !
!
        ALLOCATE(int_state%DUGWD    (IMS:IME,JMS:JME))         ! U comp. GWD tend (m s-1)
        ALLOCATE(int_state%DVGWD    (IMS:IME,JMS:JME))         ! V comp. GWD tend (m s-1)
!
        ALLOCATE(int_state%TMPMIN   (IMS:IME,JMS:JME))         ! Max temp (K)
        ALLOCATE(int_state%TMPMAX   (IMS:IME,JMS:JME))         ! Min temp (K)
!
        ALLOCATE(int_state%SFCNSW   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SFCDSW   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SFALB    (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SFCDLW   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%TSFLW    (IMS:IME,JMS:JME))         !
!
        ALLOCATE(int_state%ZORFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SIHFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SICFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SLPFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%TG3FCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%VEGFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%VETFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%SOTFCS   (IMS:IME,JMS:JME))         !
        ALLOCATE(int_state%ALBFC1   (IMS:IME,JMS:JME,4))       !
        ALLOCATE(int_state%ALFFC1   (IMS:IME,JMS:JME,2))       !
!
        ALLOCATE(int_state%SWH      (IMS:IME,JMS:JME,LM))      !
        ALLOCATE(int_state%HLW      (IMS:IME,JMS:JME,LM))      !
!
        ALLOCATE(int_state%PHY_F2DV (IMS:IME,JMS:JME,3))       ! for Zhao =3, Ferr=1
        ALLOCATE(int_state%PHY_F3DV (IMS:IME,JMS:JME,LM,4))    ! for Zhao =4, Ferr=3
!
        int_state%SOLCON=0.0D0
        int_state%SLAG  =0.0D0
        int_state%SDEC  =0.0D0
        int_state%CDEC  =0.0D0
!
        DO J=JTS,JTE
          int_state%DDY   (J)=0.0D0
          int_state%JINDX1(J)=0
          int_state%JINDX2(J)=0
        ENDDO
!
        DO J=JMS,JME
        DO I=IMS,IME
!
          int_state%DUGWD  (I,J)=0.0D0
          int_state%DVGWD  (I,J)=0.0D0
!
          int_state%TMPMIN (I,J)=373.0D0
          int_state%TMPMAX (I,J)=173.0D0
!
          int_state%SFCNSW (I,J)=0.0D0
          int_state%SFCDSW (I,J)=0.0D0
          int_state%SFALB  (I,J)=0.0D0
          int_state%SFCDLW (I,J)=0.0D0
          int_state%TSFLW  (I,J)=0.0D0
          int_state%ZORFCS (I,J)=-1.D6
          int_state%SIHFCS (I,J)=-1.D6
          int_state%SICFCS (I,J)=-1.D6
          int_state%SLPFCS (I,J)=-1.D6
          int_state%TG3FCS (I,J)=-1.D6
          int_state%VEGFCS (I,J)=-1.D6
          int_state%VETFCS (I,J)=-1.D6
          int_state%SOTFCS (I,J)=-1.D6
!
          DO N=1,4
            int_state%ALBFC1(I,J,N)=-1.D6
          ENDDO
!
          DO N=1,2
            int_state%ALFFC1(I,J,N)=-1.D6
          ENDDO
!
          DO L=1,LM
            int_state%SWH  (I,J,L)=-1.D6
            int_state%HLW  (I,J,L)=-1.D6
          ENDDO
!
          DO N=1,3                                 ! for Zhao =3, Ferr=1
            int_state%PHY_F2DV (I,J,N)=0.0D0
          ENDDO
!
          DO N=1,4                                 ! for Zhao =4, Ferr=3
          DO L=1,LM
            int_state%PHY_F3DV (I,J,L,N)=0.0D0
          ENDDO
          ENDDO
!
        ENDDO
        ENDDO

        DO N=1,TIMEOZ
        DO L=1,PL_COEFF
        DO J=1,LEVOZP
        DO I=1,LATSOZP
          int_state%OZPLIN(I,J,L,N)=-1.D6
        ENDDO
        ENDDO
        ENDDO
        ENDDO
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

      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE
      INTEGER, INTENT(IN) :: LM

      INTEGER :: INDX, I, J, L

!-----------------------------------------------------------------------
!
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
!***  ADDITIONAL TRACERS.
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
      int_state%INDX_RRW=4
      CALL FIND_VAR_INDX('RRW',int_state%VARS,int_state%NUM_VARS,I)
      int_state%VARS(I)%R3D=>int_state%TRACERS(:,:,:,int_state%INDX_RRW)
      int_state%RRW=>int_state%VARS(I)%R3D
!
!--------------------------------
!***  Water tracers
!--------------------------------
!
      int_state%INDX_WATER_START = int_state%NUM_TRACERS_MET + int_state%NUM_TRACERS_CHEM + 1
      int_state%INDX_WATER_END = int_state%INDX_WATER_START + int_state%NUM_WATER - 1
      int_state%WATER=>int_state%TRACERS(:,:,:,int_state%INDX_WATER_START:int_state%INDX_WATER_END)
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%Q2(I,J,L)=0.02
        int_state%T(I,J,L)=-1.E6
        int_state%U(I,J,L)=-1.E6
        int_state%V(I,J,L)=-1.E6
        int_state%DUDT(I,J,L)=-1.E6
        int_state%DVDT(I,J,L)=-1.E6

        int_state%RLWTT(I,J,L)=0.
        int_state%RSWTT(I,J,L)=0.

        int_state%EXCH_H(I,J,L)=0.
        int_state%XLEN_MIX(I,J,L)=-1.E6

        int_state%CLDFRA(I,J,L)=0.
        int_state%TRAIN(I,J,L) =0.
        int_state%TCUCN(I,J,L) =0.
      ENDDO
      ENDDO
      ENDDO

      DO J=JMS,JME
      DO L=1,LM+1
      DO I=IMS,IME
        int_state%RQVBLTEN(I,J,L)=-1.E6
        int_state%RTHBLTEN(I,J,L)=-1.E6
        int_state%W0AVG(I,L,J)=-1.E6
        int_state%WINT(I,L,J)=-1.E6
      ENDDO
      ENDDO
      ENDDO

      DO L=1,NUM_SOIL_LAYERS
        int_state%SLDPTH(L)=SLDPTH(L)
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SMC(I,J,L)=-1.E6
        int_state%STC(I,J,L)=-1.E6
        int_state%SH2O(I,J,L)=-1.E6
      ENDDO
      ENDDO
      ENDDO
  
      DO L=1,int_state%PCPHR
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PPTDAT(I,J,L)=-1.E6
      ENDDO
      ENDDO
      ENDDO

      int_state%NSOIL=NUM_SOIL_LAYERS

      DO J=JMS,JME
      DO I=IMS,IME
        int_state%LPBL(I,J)    =-999
        int_state%NCFRCV(I,J)  =-999
        int_state%NCFRST(I,J)  =-999
        int_state%ACFRCV(I,J)  =-1.E6
        int_state%ACFRST(I,J)  =-1.E6
        int_state%AKHS(I,J)    = 0.
        int_state%AKHS_OUT(I,J)= 0.
        int_state%AKMS(I,J)    = 0.
        int_state%AKMS_OUT(I,J)= 0.
        int_state%ALBASE(I,J)  =-1.E6
        int_state%ALBEDO(I,J)  =-1.E6
        int_state%ALWIN(I,J)   =-1.E6
        int_state%ALWOUT(I,J)  =-1.E6
        int_state%ALWTOA(I,J)  =-1.E6
        int_state%ASWIN(I,J)   =-1.E6
        int_state%ASWOUT(I,J)  =-1.E6
        int_state%ASWTOA(I,J)  =-1.E6
        int_state%BGROFF(I,J)  =-1.E6
        int_state%CFRACL(I,J)  =-1.E6
        int_state%CFRACM(I,J)  =-1.E6
        int_state%CFRACL(I,J)  =-1.E6
        int_state%CNVBOT(I,J)  =0.0
        int_state%CNVTOP(I,J)  =0.0
        int_state%CMC(I,J)     =-1.E6
        int_state%CPRATE(I,J)  =-1.E6
        int_state%CUPPT(I,J)   =-1.E6
        int_state%CZMEAN(I,J)  =-1.E6
        int_state%CZEN(I,J)    =-1.E6
        int_state%DDATA(I,J)   =-1.E6
        int_state%LSPA(I,J)    =-1.E6
        int_state%EPSR(I,J)    =-1.E6
        int_state%FIS(I,J)     =-1.E6
        int_state%HBOT(I,J)    =-1.E6
        int_state%HBOTD(I,J)   =-1.E6
        int_state%HBOTS(I,J)   =-1.E6
        int_state%HTOP(I,J)    =-1.E6
        int_state%HTOPD(I,J)   =-1.E6
        int_state%HTOPS(I,J)   =-1.E6
        int_state%GLAT(I,J)    =-1.E6
        int_state%GLON(I,J)    =-1.E6
        int_state%GRNFLX(I,J)  = 0.
        int_state%MAVAIL(I,J)  = 1.
        int_state%MXSNAL(I,J)  =-1.E6
        int_state%PBLH(I,J)    =-1.E6
        int_state%MIXHT(I,J)   =-1.E6
        int_state%PD(I,J)      =-1.E6
        int_state%POTEVP(I,J)  = 0.
        int_state%POTFLX(I,J)  =-1.E6
        int_state%QCG(I,J)     =-1.E6
        int_state%QSG(I,J)     =-1.E6
        int_state%QSH(I,J)     = 0.
        int_state%QVG(I,J)     =-1.E6
        int_state%QWBS(I,J)    =-1.E6
        int_state%QZ0(I,J)     = 0.
        int_state%RADOT(I,J)   = 0.
        int_state%RLWIN(I,J)   = 0.
        int_state%RMOL(I,J)    =-1.E6
        int_state%RSWIN(I,J)   = 0.
        int_state%RSWINC(I,J)  = 0.
        int_state%RSWOUT(I,J)  = 0.
        int_state%RLWTOA(I,J)  = 0.
        int_state%RSWTOA(I,J)  = 0.
        int_state%SFCEVP(I,J)  = 0.
        int_state%SFCEXC(I,J)  = 0.
        int_state%SFCLHX(I,J)  =-1.E6
        int_state%SFCSHX(I,J)  =-1.E6
        int_state%SHDMAX(I,J)  =-1.E6
        int_state%SHDMIN(I,J)  =-1.E6
        int_state%SICE(I,J)    =-1.E6
        int_state%SIGT4(I,J)   =-1.E6
        int_state%SM(I,J)      =-1.E6
        int_state%SMSTAV(I,J)  = 0.
        int_state%SMSTOT(I,J)  = 0.
        int_state%SNO(I,J)     = 0.
        int_state%SNOPCX(I,J)  =-1.E6
        int_state%SOILT1(I,J)  =-1.E6
        int_state%SOILTB(I,J)  =-1.E6
        int_state%SR(I,J)      =-1.E6
        int_state%SSROFF(I,J)  = 0.
        int_state%SST(I,J)     =-1.E6
        int_state%STDH(I,J)    =-1.E6 !zj
        int_state%SUBSHX(I,J)  =-1.E6
        int_state%THS(I,J)     =-1.E6
        int_state%THZ0(I,J)    = 0.
        int_state%TSKIN(I,J)   =-1.E6
        int_state%TSNAV(I,J)   =-1.E6
        int_state%TWBS(I,J)    =-1.E6
        int_state%USTAR(I,J)   = 0.1
        int_state%UZ0(I,J)     = 0.
        int_state%VEGFRC(I,J)  =-1.E6
        int_state%VZ0(I,J)     = 0.
        int_state%Z0(I,J)      =-1.E6
        int_state%Z0BASE(I,J)  =-1.E6
        int_state%CROT(I,J)    = 0.
        int_state%SROT(I,J)    = 0.
        int_state%HSTDV(I,J)   = 0.
        int_state%HCNVX(I,J)   = 0.
        int_state%HASYW(I,J)   = 0.
        int_state%HASYS(I,J)   = 0.
        int_state%HASYSW(I,J)  = 0.
        int_state%HASYNW(I,J)  = 0.
        int_state%HLENW(I,J)   = 0.
        int_state%HLENS(I,J)   = 0.
        int_state%HLENSW(I,J)  = 0.
        int_state%HLENNW(I,J)  = 0.
        int_state%HANGL(I,J)   = 0.
        int_state%HANIS(I,J)   = 0.
        int_state%HSLOP(I,J)   = 0.
        int_state%HZMAX(I,J)   = 0.
      ENDDO
      ENDDO
      
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)= 0.
        int_state%CUPREC(I,J)= 0.
        int_state%PREC(I,J)  = 0.
        int_state%CLDEFI(I,J)=-1.E6
        int_state%PSHLTR(I,J)=-1.E6
        int_state%PSFC(I,J)  =-1.E6
        int_state%Q02(I,J)   = 0.
        int_state%Q10(I,J)   = 0.
        int_state%QSHLTR(I,J)= 0.
        int_state%T2(I,J)    = 0.
        int_state%TH02(I,J)  = 0.
        int_state%TH10(I,J)  = 0.
        int_state%TSHLTR(I,J)= 0.
        int_state%U10(I,J)   = 0.
        int_state%V10(I,J)   = 0.
        int_state%TLMIN(I,J) = 0.
        int_state%TLMAX(I,J) = 0.

        int_state%ACUTIM(I,J) = 0.
        int_state%APHTIM(I,J) = 0.
        int_state%ARDLW(I,J)  = 0.
        int_state%ARDSW(I,J)  = 0.
        int_state%ASRFC(I,J)  = 0.
        int_state%AVRAIN(I,J) = 0.
        int_state%AVCNVC(I,J) = 0.
      ENDDO
      ENDDO

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
      END SUBROUTINE SET_INTERNAL_STATE_PHY_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_PHY_VAR_PTR(INT_STATE,ALLOC_FLAG,LM)

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

      USE MODULE_DM_PARALLEL,ONLY: IMS,IME,JMS,JME

      IMPLICIT NONE

      TYPE(PHYSICS_INTERNAL_STATE),INTENT(INOUT) :: INT_STATE
      LOGICAL, INTENT(IN) :: ALLOC_FLAG
      INTEGER, INTENT(IN) :: LM

      INTEGER :: N

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NSOIL'      ,int_state%NSOIL ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NPHS'       ,int_state%NPHS  )  
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NCLOD'      ,int_state%NCLOD ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NHEAT'      ,int_state%NHEAT ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NPREC'      ,int_state%NPREC ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NRDLW'      ,int_state%NRDLW ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NRDSW'      ,int_state%NRDSW ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NSRFC'      ,int_state%NSRFC ) 

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ISLTYP'     ,int_state%ISLTYP   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'IVGTYP'     ,int_state%IVGTYP   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NCFRCV'     ,int_state%NCFRCV   ,(/ IMS,JMS /),(/ IME,JME /) )  
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'NCFRST'     ,int_state%NCFRST   ,(/ IMS,JMS /),(/ IME,JME /) )  

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SLDPTH'     ,int_state%SLDPTH	      ,1 ,NUM_SOIL_LAYERS  ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'MP_RESTART' ,int_state%MP_RESTART_STATE  ,1 ,MICRO_RESTART    ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TBPVS_STAT' ,int_state%TBPVS_STATE       ,1 ,MICRO_RESTART    ) 
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TBPVS0_STA' ,int_state%TBPVS0_STATE      ,1 ,MICRO_RESTART    ) 

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACFRCV'     ,int_state%ACFRCV   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACFRST'     ,int_state%ACFRST   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACPREC'     ,int_state%ACPREC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACSNOM'     ,int_state%ACSNOM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACSNOW'     ,int_state%ACSNOW   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AKHS_OUT'   ,int_state%AKHS_OUT ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AKMS_OUT'   ,int_state%AKMS_OUT ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ALBASE'     ,int_state%ALBASE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ALBEDO'     ,int_state%ALBEDO   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ALWIN'      ,int_state%ALWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ALWOUT'     ,int_state%ALWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ALWTOA'     ,int_state%ALWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ASWIN'      ,int_state%ASWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ASWOUT'     ,int_state%ASWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ASWTOA'     ,int_state%ASWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'BGROFF'     ,int_state%BGROFF   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CFRACH'     ,int_state%CFRACH   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CFRACL'     ,int_state%CFRACL   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CFRACM'     ,int_state%CFRACM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CLDEFI'     ,int_state%CLDEFI   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CMC'        ,int_state%CMC      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CNVBOT'     ,int_state%CNVBOT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CNVTOP'     ,int_state%CNVTOP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CPRATE'     ,int_state%CPRATE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CUPPT'      ,int_state%CUPPT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CUPREC'     ,int_state%CUPREC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CZEN'       ,int_state%CZEN     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CZMEAN'     ,int_state%CZMEAN   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'EPSR'       ,int_state%EPSR     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'FIS'        ,int_state%FIS      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'GRNFLX'     ,int_state%GRNFLX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HBOTD'      ,int_state%HBOTD    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HBOTS'      ,int_state%HBOTS    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HTOPD'      ,int_state%HTOPD    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HTOPS'      ,int_state%HTOPS    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'MXSNAL'     ,int_state%MXSNAL   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PBLH'       ,int_state%PBLH     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PD'         ,int_state%PD       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'POTEVP'     ,int_state%POTEVP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PREC'       ,int_state%PREC     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PSHLTR'     ,int_state%PSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Q10'        ,int_state%Q10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'QSH'        ,int_state%QSH      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'QSHLTR'     ,int_state%QSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'QWBS'       ,int_state%QWBS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'QZ0'        ,int_state%QZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RADOT'      ,int_state%RADOT    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RLWIN'      ,int_state%RLWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RLWTOA'     ,int_state%RLWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RSWIN'      ,int_state%RSWIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RSWINC'     ,int_state%RSWINC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RSWOUT'     ,int_state%RSWOUT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SFCEVP'     ,int_state%SFCEVP   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SFCEXC'     ,int_state%SFCEXC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SFCLHX'     ,int_state%SFCLHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SFCSHX'     ,int_state%SFCSHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SI'         ,int_state%SI       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SICE'       ,int_state%SICE     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SIGT4'      ,int_state%SIGT4    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SM'         ,int_state%SM       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SMSTAV'     ,int_state%SMSTAV   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SMSTOT'     ,int_state%SMSTOT   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SNO'        ,int_state%SNO      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SNOPCX'     ,int_state%SNOPCX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SOILTB'     ,int_state%SOILTB   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SR'         ,int_state%SR       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SSROFF'     ,int_state%SSROFF   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SST'        ,int_state%SST      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SUBSHX'     ,int_state%SUBSHX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TG'         ,int_state%TG       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TH10'       ,int_state%TH10     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'THS'        ,int_state%THS      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'THZ0'       ,int_state%THZ0     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TSHLTR'     ,int_state%TSHLTR   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TWBS'       ,int_state%TWBS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'U10'        ,int_state%U10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'USTAR'      ,int_state%USTAR    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'UZ0'        ,int_state%UZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'V10'        ,int_state%V10      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'VEGFRC'     ,int_state%VEGFRC   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'VZ0'        ,int_state%VZ0      ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Z0'         ,int_state%Z0       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TSKIN'      ,int_state%TSKIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AKHS'       ,int_state%AKHS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AKMS'       ,int_state%AKMS     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HBOT'       ,int_state%HBOT     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'HTOP'       ,int_state%HTOP     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RSWTOA'     ,int_state%RSWTOA   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'POTFLX'     ,int_state%POTFLX   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RMOL'       ,int_state%RMOL     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'T2'         ,int_state%T2       ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Z0BASE'     ,int_state%Z0BASE   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'PSFC'       ,int_state%PSFC     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TLMIN'      ,int_state%TLMIN    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TLMAX'      ,int_state%TLMAX    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'LSPA'       ,int_state%LSPA     ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ACUTIM'     ,int_state%ACUTIM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'APHTIM'     ,int_state%APHTIM   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ARDLW'      ,int_state%ARDLW    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ARDSW'      ,int_state%ARDSW    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'ASRFC'      ,int_state%ASRFC    ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AVRAIN'     ,int_state%AVRAIN   ,(/ IMS,JMS /),(/ IME,JME /) )
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'AVCNVC'     ,int_state%AVCNVC   ,(/ IMS,JMS /),(/ IME,JME /) )

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Q2'         ,int_state%Q2       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'OMGALF'     ,int_state%OMGALF   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'T'          ,int_state%T	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'Q'          ,int_state%Q	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'U'          ,int_state%U	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'V'          ,int_state%V	     ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RLWTT'      ,int_state%RLWTT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RSWTT'      ,int_state%RSWTT    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'EXCH_H'     ,int_state%EXCH_H   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CLDFRA'     ,int_state%CLDFRA   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'CW'         ,int_state%CW       ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'F_ICE'      ,int_state%F_ICE    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'F_RAIN'     ,int_state%F_RAIN   ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'F_RIMEF'    ,int_state%F_RIMEF  ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TRAIN'      ,int_state%TRAIN    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'XLEN_MIX'   ,int_state%XLEN_MIX ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TCUCN'      ,int_state%TCUCN    ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SMC'        ,int_state%SMC      ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'STC'        ,int_state%STC      ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'SH2O'       ,int_state%SH2O     ,(/ IMS,JMS,1 /),(/ IME,JME,NUM_SOIL_LAYERS /))
      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'RRW'        ,int_state%RRW      ,(/ IMS,JMS,1 /),(/ IME,JME,LM /))

      CALL SET_VAR_PTR(int_state%VARS,int_state%NUM_VARS,ALLOC_FLAG,'TRACERS'    ,int_state%TRACERS,(/ IMS,JMS,1,1 /),(/ IME,JME,LM,int_state%NUM_TRACERS_TOTAL /) )

      DO N=1,int_state%NUM_VARS
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
