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
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
      USE MODULE_INCLUDE
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                    
!
      USE MODULE_LANDSURFACE ,ONLY: NUM_SOIL_LAYERS
      USE MODULE_MICROPHYSICS,ONLY: MICRO_RESTART
      USE MODULE_ERR_MSG     ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: INTERNAL_STATE,SET_INTERNAL_STATE_PHY                   &
               ,UPDATE_INTERNAL_STATE_PHY,WRAP_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      TYPE INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Begin with the namelist variables.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: IM,JM,LM                                  &
                             ,NHOURS_HISTORY                            &
                             ,NHOURS_RESTART                            &
                             ,NSOIL                                     &  !<-- # of soil layers
                             ,NUM_TRACERS_MET                           &  !<-- Number of meteorological tracers (e.g. water)
                             ,NUM_TRACERS_CHEM                          &  !<-- Number of chem/aerosol tracers
                             ,START_YEAR,START_MONTH,START_DAY          &
                             ,START_HOUR,START_MINUTE,START_SECOND
!
        INTEGER(KIND=KINT) :: NPHS,NPRECIP,NRADL,NRADS,UCMCALL,DT_INT
!
        REAL(KIND=KFPT) :: DT,SBD,WBD,TPH0D,TLM0D
!
        LOGICAL :: GLOBAL,GWDFLG,HYDRO,NESTED,PCPFLG,RESTART,SPECIFIED  &
                  ,NHRS_UDEF
!
!-----------------------------------------------------------------------
!***  Distributed memory information.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: INPES                                     &  !<-- Forecast PE's in the E-W direction
                             ,JNPES                                     &  !<-- Forecast PE's in the N-S direction
                             ,MYPE                                      &  !<-- This task's ID
                             ,NUM_PES                                      !<-- Total number of forecast tasks
!
        INTEGER(KIND=KINT) :: ITS,ITE,JTS,JTE                           &
                             ,IMS,IME,JMS,JME
!
        INTEGER(KIND=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART         &  !<-- Each task's local starting I
                                                  ,LOCAL_IEND           &  !<-- Each task's local ending I
                                                  ,LOCAL_JSTART         &  !<-- Each task's local starting J
                                                  ,LOCAL_JEND              !<-- Each task's local ending J
!
!-----------------------------------------------------------------------
!***  Horizontal/Vertical grid
!-----------------------------------------------------------------------
!
        REAL(KIND=KFPT) :: DYH,DYV,RDYH,RDYV,PDTOP,PT
!
        REAL(KIND=KFPT),DIMENSION(:),POINTER :: DSG2,DXH,DXV,RDXH,RDXV  &
                                               ,PDSG1,PSGML1,SGML2
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: GLAT,GLON
!
!-----------------------------------------------------------------------
!***  Integration quantities.
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: NTSD                                         !<-- Internal timestep counter
!
        INTEGER(KIND=KINT) :: NHRS_CLOD                                 &  !<-- Fcst hours cloud is accumulated
                             ,NHRS_HEAT                                 &  !<-- Fcst hours heating is accumulated
                             ,NHRS_PREC                                 &  !<-- Fcst hours precip is accumulated
                             ,NHRS_RDLW                                 &  !<-- Fcst hours LW radiation is accumulated
                             ,NHRS_RDSW                                 &  !<-- Fcst hours SW radiation is accumulated
                             ,NHRS_SRFC                                    !<-- Fcst hours sfc evap/flux is accumulated
!
        INTEGER(KIND=KINT) :: NCLOD                                     &  !<-- # of fundamental timesteps cloud is accumulated
                             ,NHEAT                                     &  !<-- # of fundamental timesteps latent heating is accumulated
                             ,NPREC                                     &  !<-- # of fundamental timesteps precip is accumulated
                             ,NRDLW                                     &  !<-- # of fundamental timesteps LW radiation is accumulated
                             ,NRDSW                                     &  !<-- # of fundamental timesteps SW radiation is accumulated
                             ,NSRFC                                        !<-- # of fundamental timesteps sfc evap/flux is accumulated
!
        INTEGER(KIND=KINT),DIMENSION(:,:),POINTER :: ISLTYP,IVGTYP      &
                                                    ,LPBL               &
                                                    ,NCFRCV,NCFRST
!
        REAL(KIND=KFPT) :: ACUTIM                                       &  
                          ,APHTIM                                       &
                          ,ARDLW                                        &  !<-- Counter in summing LW radiation flux
                          ,ARDSW                                        &  !<-- Counter in summing SW radiation flux
                          ,ASRFC                                        &  !<-- Counter in summing sfc flux
                          ,AVRAIN                                       &  !<-- Counter in summing latent heating from all precip
                          ,AVCNVC                                          !<-- Counter in summing latent heating from cnvc precip
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: T,U,V               &
                                                   ,DUDT,DVDT           &
                                                   ,Q,CW                &
                                                   ,Q2,OMGALF           &
                                                   ,RRW
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: RLWTT,RSWTT
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: EXCH_H
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: TCUCN,W0AVG,WINT
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: CLDFRA              &
                                                   ,F_ICE,F_RAIN        &
                                                   ,F_RIMEF             &
                                                   ,RQVBLTEN,RTHBLTEN   &
                                                   ,TRAIN,XLEN_MIX
!
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: SMC,STC,SH2O
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: CROT,SROT             &
                         ,HANGL,HANIS,HASYS,HASYSW,HASYNW,HASYW,HCNVX   &
                         ,HLENNW,HLENSW,HLENW,HLENS,HSLOP,HSTDV,HZMAX
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ACFRCV,ACFRST         &
                                                 ,AKHS,AKHS_OUT         &
                                                 ,AKMS,AKMS_OUT         &
                                                 ,CFRACH,CFRACL,CFRACM  &
                                                 ,CMC,CNVBOT,CNVTOP     &
                                                 ,CPRATE,CUPPT          &
                                                 ,CZMEAN,CZEN           &
                                                 ,DDATA,EPSR,FIS        &
                                                 ,HBOT,HBOTD,HBOTS      &
                                                 ,HTOP,HTOPD,HTOPS      &
                                                 ,MIXHT,PBLH,PD,QSH,QZ0 &
                                                 ,RLWIN,RSWIN,RSWINC    &
                                                 ,RSWOUT,RLWTOA,RSWTOA  &
                                                 ,SICE,SIGT4,SM         &
!zj                                                 ,SST,THS,THZ0,USTAR    &
                                                 ,SST,STDH              & !zj
                                                 ,THS,THZ0,USTAR        & !zj
                                                 ,UZ0,UZ0H,VZ0,VZ0H     &
                                                 ,Z0,Z0BASE
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ALBASE,ALBEDO         &
                                                 ,ALWIN,ALWOUT,ALWTOA   &
                                                 ,ASWIN,ASWOUT,ASWTOA   &
                                                 ,BGROFF,GRNFLX         &
                                                 ,MAVAIL,MXSNAL         &
                                                 ,POTEVP,POTFLX         &
                                                 ,QCG,QSG,QVG,QWBS      &
                                                 ,RADOT,RMOL            &
                                                 ,SFCEVP,SFCEXC         &
                                                 ,SFCLHX,SFCSHX         &
!                                                 ,SH2O,SHDMAX,SHDMIN    &
                                                 ,SHDMAX,SHDMIN         &
                                                 ,SI,SMSTAV,SMSTOT      &
                                                 ,SNO,SNOPCX            &
                                                 ,SOILT1,SOILTB         &
                                                 ,SR,SSROFF,SUBSHX      &
                                                 ,TG,TSKIN,TSNAV,TWBS   &
                                                 ,VEGFRC
!
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ACPREC,ACSNOM,ACSNOW  &
                                                 ,CUPREC,CLDEFI         &
                                                 ,PREC,PSHLTR,Q02,Q10   &
                                                 ,QSHLTR                &
                                                 ,T2,TH02,TH10,TSHLTR   &
                                                 ,U10,V10
!
        REAL(KIND=KFPT),DIMENSION(:),POINTER :: MP_RESTART_STATE        &
                                               ,SLDPTH                  &
                                               ,TBPVS_STATE,TBPVS0_STATE
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
!
        LOGICAL :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
!-----------------------------------------------------------------------
!***  The general 4-D array of 3-D "tracers".
!-----------------------------------------------------------------------
!
        INTEGER(KIND=KINT) :: NUM_TRACERS_TOTAL                            !<-- Total number of "tracer" variables.
!
        INTEGER(KIND=KINT) :: INDX_Q                                    &  !<-- Location of Q in tracer arrays
                             ,INDX_CW                                   &  !<-- Location of CW in tracer arrays
                             ,INDX_RRW                                  &  !<-- Location of RRW in tracer arrays
                             ,INDX_Q2                                      !<-- Location of Q2 in tracer arrays
!
        REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: TRACERS              !<-- Storage array for "tracer" variables.
!
!-----------------------------------------------------------------------
!
      END TYPE INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  THE INTERNAL_STATE TYPE IS SUPPORTED BY A C POINTER (NOT AN F90
!***  POINTER) AND THEREFORE THE FOLLOWING TYPE IS NEEDED.
!-----------------------------------------------------------------------
!
      TYPE WRAP_INTERNAL_STATE
        TYPE(INTERNAL_STATE),POINTER :: INT_STATE
      END TYPE WRAP_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_INTERNAL_STATE_PHY(GRID_COMP,INT_STATE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_CONSTANTS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp), INTENT(IN)    :: GRID_COMP                     !<-- The Physics gridded component
      TYPE(INTERNAL_STATE),INTENT(INOUT) :: INT_STATE                     !<-- The Physics internal state
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,L,LM,N,NSTEPS_PER_HOUR,NUM_PES
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  WE CAN RETRIEVE LM FROM THE INTERNAL STATE SINCE IT WAS 
!***  PLACED THERE ALREADY FROM THE CONFIG FILE.
!***  THE ARRAY MEMORY LIMITS NEED TO BE SET IN THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
      LM=int_state%LM
!      
      int_state%IMS=IMS
      int_state%IME=IME
      int_state%JMS=JMS
      int_state%JME=JME
!
      int_state%NSOIL=NUM_SOIL_LAYERS
!
!-----------------------------------------------------------------------
!***  INITIALIZE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS OF THE INTERNAL STATE.
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%DSG2(1:LM))     ! Delta sigma (bottom domain)
      ALLOCATE(int_state%PDSG1(1:LM))    ! Delta pressure (top domain)
      ALLOCATE(int_state%PSGML1(1:LM))   ! Midlayer pressure (top domain)
      ALLOCATE(int_state%SGML2(1:LM))    ! Midlayer sigma (bottom domain)
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
!***  LOCAL HORIZONTAL SUBDOMAIN LIMITS FOR ALL FORECAST TASKS.
!-----------------------------------------------------------------------
!
      NUM_PES=int_state%NUM_PES
      ALLOCATE(int_state%LOCAL_ISTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_IEND  (0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JSTART(0:NUM_PES-1))
      ALLOCATE(int_state%LOCAL_JEND  (0:NUM_PES-1))
!
!-----------------------------------------------------------------------
!***  PROGNOSTIC ARRAYS
!-----------------------------------------------------------------------
!
      ALLOCATE(int_state%Q2  (IMS:IME,JMS:JME,1:LM))      ! 2*tke  (m2 s-2)
      ALLOCATE(int_state%OMGALF(IMS:IME,JMS:JME,1:LM))    ! Omega alpha from the dynamics  (K)
      ALLOCATE(int_state%T   (IMS:IME,JMS:JME,1:LM))      ! Sensible temperature  (K)
      ALLOCATE(int_state%U   (IMS:IME,JMS:JME,1:LM))      ! U wind component  (m s-1)
      ALLOCATE(int_state%V   (IMS:IME,JMS:JME,1:LM))      ! V wind component  (m s-1)
      ALLOCATE(int_state%DUDT(IMS:IME,JMS:JME,1:LM))      ! U wind component tendency  (m s-2)
      ALLOCATE(int_state%DVDT(IMS:IME,JMS:JME,1:LM))      ! V wind component tendency  (m s-2)
!
      ALLOCATE(int_state%RLWTT(IMS:IME,JMS:JME,1:LM))     ! Radiative T tendency, longwave  (K s-1)
      ALLOCATE(int_state%RSWTT(IMS:IME,JMS:JME,1:LM))     ! Radiative T tendency, shortwave (K s-1)
!
      ALLOCATE(int_state%RQVBLTEN(IMS:IME,JMS:JME,1:LM+1)) ! Specific humidity tendency from turbulence  (kg kg-1 s-1)
      ALLOCATE(int_state%RTHBLTEN(IMS:IME,JMS:JME,1:LM+1)) ! Theta tendency from turbulence  (K s-1)
!
      ALLOCATE(int_state%EXCH_H(IMS:IME,JMS:JME,1:LM))    ! Turbulent exchange coefficient for heat  (m2 s-1)
!
      ALLOCATE(int_state%CLDFRA(IMS:IME,JMS:JME,1:LM))    ! Cloud fraction
      ALLOCATE(int_state%F_ICE(IMS:IME,JMS:JME,1:LM))     ! Fraction of condensate that is ice
      ALLOCATE(int_state%F_RAIN(IMS:IME,JMS:JME,1:LM))    ! Rain fraction of liquid part of CWM
      ALLOCATE(int_state%F_RIMEF(IMS:IME,JMS:JME,1:LM))   ! Rime factor
      ALLOCATE(int_state%TRAIN(IMS:IME,JMS:JME,1:LM))     ! Accumulated stratiform T tendency  (K s-1)
      ALLOCATE(int_state%XLEN_MIX(IMS:IME,JMS:JME,1:LM))  ! Mixing length  (m)
!
      ALLOCATE(int_state%TCUCN(IMS:IME,JMS:JME,1:LM))     ! Accumulated convective T tendency  (K s-1)
      ALLOCATE(int_state%W0AVG(IMS:IME,1:LM+1,JMS:JME))   ! Time-averaged vertical velocity (for K-F)  (m s-1)
      ALLOCATE(int_state%WINT(IMS:IME,1:LM+1,JMS:JME))    ! Interface vertical velocity (m s-1)
!
      ALLOCATE(int_state%SLDPTH(1:NUM_SOIL_LAYERS))                       ! Thickness of soil layers (m) from top
      ALLOCATE(int_state%SMC(IMS:IME,JMS:JME,1:NUM_SOIL_LAYERS))          ! Soil moisture volume fraction
      ALLOCATE(int_state%STC(IMS:IME,JMS:JME,1:NUM_SOIL_LAYERS))          ! Soil temperature  (K)
      ALLOCATE(int_state%SH2O(IMS:IME,JMS:JME,1:NUM_SOIL_LAYERS))         ! Soil non-frozen moisture
!
!-----------------------------------------------------------------------
!***  THE ARRAY CALLED WATER IS A SPECIAL CASE NEEDED TO SATISFY
!***  VARIOUS WRF PHYSICS OPTIONS.  THE IS SET TO 1+Number_of_species
!***  INCLUDING VAPOR.  THE "1+" IS NEEDED BECAUSE WRF NEVER TOUCHES
!***  THE FIRST LEVEL.
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
                                  int_state%NUM_TRACERS_MET             &  !<-- # of meteorological tracers such as water (see below)
                                 +int_state%NUM_TRACERS_CHEM            &  !<-- # of specified scalars (chem, aerosol, etc.)
                                 +int_state%NUM_WATER                      !<-- # of water types
!
      ALLOCATE(int_state%TRACERS(IMS:IME,JMS:JME,1:LM,1:int_state%NUM_TRACERS_TOTAL))  !<-- All tracer variables
!
      DO N=1,int_state%NUM_TRACERS_MET
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%TRACERS(I,J,L,N)=1.E-20
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  POINT Q AT LEVEL 1 OF THE TRACERS ARRAY.
!-----------------------------------------------------------------------
!
      int_state%INDX_Q=1                                                    !<-- Water vapor is always in location 1 of TRACERS array
      int_state%Q=>int_state%TRACERS(IMS:IME,JMS:JME,1:LM,int_state%INDX_Q) !<-- Water vapor is always in location 1 of TRACERS array
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
      int_state%INDX_Q2=3
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
!***  NOTE:  ADDITIONAL SCALAR VARIABLES IN THE TRACER ARRAYS WILL
!***         BEGIN AT LOCATION NUM_TRACERS_MET+1.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      DO L=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%Q2(I,J,L)=-1.E6
        int_state%T(I,J,L)=-1.E6
        int_state%U(I,J,L)=-1.E6
        int_state%V(I,J,L)=-1.E6
        int_state%DUDT(I,J,L)=-1.E6
        int_state%DVDT(I,J,L)=-1.E6
!
        int_state%RLWTT(I,J,L)=0.
        int_state%RSWTT(I,J,L)=0.
!
        int_state%EXCH_H(I,J,L)=0.
        int_state%XLEN_MIX(I,J,L)=-1.E6
!
        int_state%CLDFRA(I,J,L)=0.
        int_state%TRAIN(I,J,L) =0.
        int_state%TCUCN(I,J,L) =0.
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,LM+1
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%RQVBLTEN(I,J,L)=-1.E6
        int_state%RTHBLTEN(I,J,L)=-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO L=1,LM+1
      DO I=IMS,IME
        int_state%W0AVG(I,L,J)=-1.E6
        int_state%WINT(I,L,J)=-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      DO L=1,NUM_SOIL_LAYERS
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SMC(I,J,L)=-1.E6
        int_state%STC(I,J,L)=-1.E6
        int_state%SH2O(I,J,L)=-1.E6
      ENDDO
      ENDDO
      ENDDO
!
      ALLOCATE(int_state%ISLTYP(IMS:IME,JMS:JME))   ! Soil type
      ALLOCATE(int_state%IVGTYP(IMS:IME,JMS:JME))   ! Vegetation type
      ALLOCATE(int_state%LPBL(IMS:IME,JMS:JME))     ! Model layer containing top of the PBL
      ALLOCATE(int_state%NCFRCV(IMS:IME,JMS:JME))   ! Number of times convective cloud fraction >0 between outputs
      ALLOCATE(int_state%NCFRST(IMS:IME,JMS:JME))   ! Number of times stratiform cloud fraction >0 between outputs
      ALLOCATE(int_state%ACFRCV(IMS:IME,JMS:JME))   ! Accumulated convective cloud fractions between outputs
      ALLOCATE(int_state%ACFRST(IMS:IME,JMS:JME))   ! Accumulated statiform cloud fractions between outputs
      ALLOCATE(int_state%AKHS(IMS:IME,JMS:JME))     ! Sfc exch coeff T and q divided by delta Z  (m s-1)
      ALLOCATE(int_state%AKHS_OUT(IMS:IME,JMS:JME)) ! Sfc exch coeff T and q  (m2 s-1)
      ALLOCATE(int_state%AKMS(IMS:IME,JMS:JME))     ! Sfc exch coeff momentum divided by delta Z  (m s-1)
      ALLOCATE(int_state%AKMS_OUT(IMS:IME,JMS:JME)) ! Sfc exch coeff momentum  (m2 s-1)
      ALLOCATE(int_state%ALBASE(IMS:IME,JMS:JME))   ! Base albedo
      ALLOCATE(int_state%ALBEDO(IMS:IME,JMS:JME))   ! Dynamic albedo
      ALLOCATE(int_state%ALWIN(IMS:IME,JMS:JME))    ! Accumulated LW down at surface  (W m-2)
      ALLOCATE(int_state%ALWOUT(IMS:IME,JMS:JME))   ! Accumulated LW from the ground  (W m-2)
      ALLOCATE(int_state%ALWTOA(IMS:IME,JMS:JME))   ! Accumulated LW at top of atmosphere  (W m-2)
      ALLOCATE(int_state%ASWIN(IMS:IME,JMS:JME))    ! Accumulated SW down at surface  (W m-2)
      ALLOCATE(int_state%ASWOUT(IMS:IME,JMS:JME))   ! Accumulated SW up at ground  (W m-2)
      ALLOCATE(int_state%ASWTOA(IMS:IME,JMS:JME))   ! Accumulated SW at top of atmosphere  (W m-2)
      ALLOCATE(int_state%BGROFF(IMS:IME,JMS:JME))   ! Subsurface runoff  (mm)
      ALLOCATE(int_state%CFRACL(IMS:IME,JMS:JME))   ! Low atmosphere cloud fraction
      ALLOCATE(int_state%CFRACM(IMS:IME,JMS:JME))   ! Mid atmosphere cloud fraction
      ALLOCATE(int_state%CFRACH(IMS:IME,JMS:JME))   ! High atmosphere cloud fraction
      ALLOCATE(int_state%CMC(IMS:IME,JMS:JME))      ! Canopy moisture  (m)
      ALLOCATE(int_state%CNVBOT(IMS:IME,JMS:JME))   ! Lowest convec cloud bottom model lyr between outputs
      ALLOCATE(int_state%CNVTOP(IMS:IME,JMS:JME))   ! Highest convec cloud top model lyr between outputs
      ALLOCATE(int_state%CPRATE(IMS:IME,JMS:JME))   ! Instantaneous convective precipitation
      ALLOCATE(int_state%CUPPT(IMS:IME,JMS:JME))    ! Convective precip between radiation calls  (m)
      ALLOCATE(int_state%CZMEAN(IMS:IME,JMS:JME))   ! Mean cosine of solar zenith angle between radiation calls
      ALLOCATE(int_state%CZEN(IMS:IME,JMS:JME))     ! Current cosine of solar zenith angle
      ALLOCATE(int_state%DDATA(IMS:IME,JMS:JME))    ! Observed precip to each physics timestep (kg m-2)
      ALLOCATE(int_state%EPSR(IMS:IME,JMS:JME))     !
      ALLOCATE(int_state%FIS(IMS:IME,JMS:JME))      ! Surface geopotential  (m2 s-2)
      ALLOCATE(int_state%GRNFLX(IMS:IME,JMS:JME))   ! Deep sloil heat flux  (W m-2)
      ALLOCATE(int_state%HBOT(IMS:IME,JMS:JME))     ! Bottom model layer of convective cloud for radiation
      ALLOCATE(int_state%HBOTD(IMS:IME,JMS:JME))    ! Bottom model layer of deep convective cloud for output
      ALLOCATE(int_state%HBOTS(IMS:IME,JMS:JME))    ! Bottom model layer of shallow convective cloud for output
      ALLOCATE(int_state%HTOP(IMS:IME,JMS:JME))     ! Top model layer of convective cloud for radiation
      ALLOCATE(int_state%HTOPD(IMS:IME,JMS:JME))    ! Top model layer of deep convective cloud for output
      ALLOCATE(int_state%HTOPS(IMS:IME,JMS:JME))    ! Top model layer of shallow convective cloud for output
      ALLOCATE(int_state%MAVAIL(IMS:IME,JMS:JME))   ! Moisture availability
      ALLOCATE(int_state%MXSNAL(IMS:IME,JMS:JME))   ! Maximum deep snow albedo
      ALLOCATE(int_state%PBLH(IMS:IME,JMS:JME))     ! PBL height  (m)
      ALLOCATE(int_state%MIXHT(IMS:IME,JMS:JME))    ! Mixed layer height  (m)
      ALLOCATE(int_state%PD(IMS:IME,JMS:JME))       ! Mass in the sigma domain  (Pa)
      ALLOCATE(int_state%POTEVP(IMS:IME,JMS:JME))   ! Accumulated potential evaporation (m)
      ALLOCATE(int_state%POTFLX(IMS:IME,JMS:JME))   ! Energy equivalent of POTEVP (W m-2)
      ALLOCATE(int_state%QCG(IMS:IME,JMS:JME))      ! Cloud water mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%QSG(IMS:IME,JMS:JME))      ! Surface saturation water vapor mixing ratio  (kg kg-1)
      ALLOCATE(int_state%QVG(IMS:IME,JMS:JME))      ! Water vapor mixing ratio at the surface  (kg kg-1)
      ALLOCATE(int_state%QSH(IMS:IME,JMS:JME))      ! Surface specific humidity  (kg kg-1)
      ALLOCATE(int_state%QWBS(IMS:IME,JMS:JME))     ! Instantaneous latent heat flux (W m-2)
      ALLOCATE(int_state%QZ0(IMS:IME,JMS:JME))      ! Specific humidity at top of viscous sublayer  (kg kg-1)
      ALLOCATE(int_state%RADOT(IMS:IME,JMS:JME))    ! Longwave from the ground  (W m-2)
      ALLOCATE(int_state%RLWIN(IMS:IME,JMS:JME))    ! Longwave down at ground  (W m-2)
      ALLOCATE(int_state%RMOL(IMS:IME,JMS:JME))     ! Reciprocal of Monin-Obukhov length  (m-1)
      ALLOCATE(int_state%RSWIN(IMS:IME,JMS:JME))    ! Shortwave down at ground  (W m-2)
      ALLOCATE(int_state%RSWINC(IMS:IME,JMS:JME))   ! Clear-sky shortwave down at ground  (W m-2)
      ALLOCATE(int_state%RSWOUT(IMS:IME,JMS:JME))   ! Shortwave up at ground  (W m-2)
      ALLOCATE(int_state%RLWTOA(IMS:IME,JMS:JME))   ! Longwave at top of atmosphere  (W m-2)
      ALLOCATE(int_state%RSWTOA(IMS:IME,JMS:JME))   ! Shortwave at top of atmosphere  (W m-2)
      ALLOCATE(int_state%SFCEVP(IMS:IME,JMS:JME))   ! Surface evaporation  (?)
      ALLOCATE(int_state%SFCEXC(IMS:IME,JMS:JME))   ! Another surface exchange coefficient for T and q (?) (m2 s-1) (see AKHS_OUT)
      ALLOCATE(int_state%SFCLHX(IMS:IME,JMS:JME))   ! Accumulated sfc latent heat flux (W m-2)
      ALLOCATE(int_state%SFCSHX(IMS:IME,JMS:JME))   ! Accumulated sfc sensible heat flux (W m-2)
!      ALLOCATE(int_state%SH2O(IMS:IME,JMS:JME))     ! Unfrozen soil moisture volume fraction
      ALLOCATE(int_state%SHDMAX(IMS:IME,JMS:JME))   ! Maximum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SHDMIN(IMS:IME,JMS:JME))   ! Minimum areal fractional coverage of annual green vegetation
      ALLOCATE(int_state%SI(IMS:IME,JMS:JME))       ! Snow depth  (m)
      ALLOCATE(int_state%SICE(IMS:IME,JMS:JME))     ! Sea ice fraction
      ALLOCATE(int_state%SIGT4(IMS:IME,JMS:JME))    ! Sigma*T**4  (W m-2)
      ALLOCATE(int_state%SM(IMS:IME,JMS:JME))       ! Sea mask (1=>sea ; 0=>land)
      ALLOCATE(int_state%SMSTAV(IMS:IME,JMS:JME))   ! Soil moisture availability for evapotranspiration
      ALLOCATE(int_state%SMSTOT(IMS:IME,JMS:JME))   ! Total soil moisture  (?)
      ALLOCATE(int_state%SNO(IMS:IME,JMS:JME))      ! Liquid water snow amount  (m)
      ALLOCATE(int_state%SNOPCX(IMS:IME,JMS:JME))   ! Snow phase change heat flux  (W m-2)
      ALLOCATE(int_state%SOILT1(IMS:IME,JMS:JME))   ! Snow temperature  (K)
      ALLOCATE(int_state%SOILTB(IMS:IME,JMS:JME))   ! Deep ground soil temperature  (K)
      ALLOCATE(int_state%SR(IMS:IME,JMS:JME))       ! Timestep mass ratio of snow:precip
      ALLOCATE(int_state%SSROFF(IMS:IME,JMS:JME))   ! Surface runoff  (mm)
      ALLOCATE(int_state%SUBSHX(IMS:IME,JMS:JME))   ! Accumulated deep soil heat flux (W m-2)
      ALLOCATE(int_state%SST(IMS:IME,JMS:JME))      ! Sea sfc temperature  (K)
      ALLOCATE(int_state%STDH(IMS:IME,JMS:JME))     ! Standard deviation of topography height (m) !zj
      ALLOCATE(int_state%TG(IMS:IME,JMS:JME))       ! Deep ground soil temperature  (K)
      ALLOCATE(int_state%THS(IMS:IME,JMS:JME))      ! Surface theta  (K)
      ALLOCATE(int_state%THZ0(IMS:IME,JMS:JME))     ! Theta at top of viscous sublayer  (K)
      ALLOCATE(int_state%TSKIN(IMS:IME,JMS:JME))    ! Skin temperature  (K)
      ALLOCATE(int_state%TSNAV(IMS:IME,JMS:JME))    ! Average snow temperature  (K)
      ALLOCATE(int_state%TWBS(IMS:IME,JMS:JME))     ! Instantaneous sensible heat flux (W m-2)
      ALLOCATE(int_state%USTAR(IMS:IME,JMS:JME))    ! Friction velocity  (m s-1)
      ALLOCATE(int_state%UZ0(IMS:IME,JMS:JME))      ! U component at top of viscous sublayer  (m s-1)
      ALLOCATE(int_state%UZ0H(IMS:IME,JMS:JME))     ! UZ0 on mass points  (m s-1)
      ALLOCATE(int_state%VEGFRC(IMS:IME,JMS:JME))   ! Vegetation fraction
      ALLOCATE(int_state%VZ0(IMS:IME,JMS:JME))      ! V component at top of viscous sublayer  (m s-1)
      ALLOCATE(int_state%VZ0H(IMS:IME,JMS:JME))     ! VZ0 on mass points  (m s-1)
      ALLOCATE(int_state%Z0(IMS:IME,JMS:JME))       ! Roughness length  (m)
      ALLOCATE(int_state%Z0BASE(IMS:IME,JMS:JME))   ! Background roughness length  (m)
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
!
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
        int_state%CMC(I,J)     =-1.E6
        int_state%CPRATE(I,J)  =-1.E6
        int_state%CUPPT(I,J)   =-1.E6
        int_state%CZMEAN(I,J)  =-1.E6
        int_state%CZEN(I,J)    =-1.E6
        int_state%DDATA(I,J)   =-1.E6
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
        int_state%MAVAIL(I,J)  =-1.E6
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
        int_state%RADOT(I,J)   =-1.E6
        int_state%RLWIN(I,J)   = 0.
        int_state%RMOL(I,J)    =-1.E6
        int_state%RSWIN(I,J)   = 0.
        int_state%RSWINC(I,J)  = 0.
        int_state%RSWOUT(I,J)  = 0.
        int_state%RLWTOA(I,J)  = 0.
        int_state%RSWTOA(I,J)  = 0.
        int_state%SFCEVP(I,J)  =-1.E6
        int_state%SFCEXC(I,J)  =-1.E6
        int_state%SFCLHX(I,J)  =-1.E6
        int_state%SFCSHX(I,J)  =-1.E6
!        int_state%SH2O(I,J)   =-1.E6
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
        int_state%USTAR(I,J)   =-1.E6
        int_state%UZ0(I,J)     =-1.E6
        int_state%UZ0H(I,J)    =-1.E6
        int_state%VEGFRC(I,J)  =-1.E6
        int_state%VZ0(I,J)     =-1.E6
        int_state%VZ0H(I,J)    =-1.E6
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
!
      ALLOCATE(int_state%ACPREC(IMS:IME,JMS:JME))  ! Accumulated precip  (m)
      ALLOCATE(int_state%ACSNOM(IMS:IME,JMS:JME))  ! Accumulated melted snow  (?)
      ALLOCATE(int_state%ACSNOW(IMS:IME,JMS:JME))  ! Accumulated snow  (?)
      ALLOCATE(int_state%CUPREC(IMS:IME,JMS:JME))  ! Conv precip  (m)
      ALLOCATE(int_state%CLDEFI(IMS:IME,JMS:JME))  ! Convective cloud efficiency
      ALLOCATE(int_state%PREC(IMS:IME,JMS:JME))    ! Precip within physics timestep  (m)
      ALLOCATE(int_state%pSHLTR(IMS:IME,JMS:JME))  ! Pressure at 2-m  (Pa)
      ALLOCATE(int_state%Q02(IMS:IME,JMS:JME))     ! Specific humidity at 2-m  (kg k-1)
      ALLOCATE(int_state%Q10(IMS:IME,JMS:JME))     ! Specific humidity at 10-m  (kg k-1)
      ALLOCATE(int_state%QSHLTR(IMS:IME,JMS:JME))  ! Specific humidity at 2-m  (kg kg-1)
      ALLOCATE(int_state%T2(IMS:IME,JMS:JME))      ! Temperature 2-m  (K)
      ALLOCATE(int_state%TH02(IMS:IME,JMS:JME))    ! Theta at 2-m  (K)
      ALLOCATE(int_state%TH10(IMS:IME,JMS:JME))    ! Theta at 10-m  (K)
      ALLOCATE(int_state%TSHLTR(IMS:IME,JMS:JME))  ! Theta at 2-m again  (K)
      ALLOCATE(int_state%U10(IMS:IME,JMS:JME))     ! U at 10-m  (m s-1)
      ALLOCATE(int_state%V10(IMS:IME,JMS:JME))     ! V at 10-m  (m s-1)
!
      ALLOCATE(int_state%MP_RESTART_STATE(MICRO_RESTART))  ! For Ferrier restart
      ALLOCATE(int_state%TBPVS_STATE(MICRO_RESTART))       ! For Ferrier restart
      ALLOCATE(int_state%TBPVS0_STATE(MICRO_RESTART))      ! For Ferrier restart
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)= 0.
        int_state%CUPREC(I,J)= 0.
        int_state%PREC(I,J)  = 0.
        int_state%CLDEFI(I,J)=-1.E6
        int_state%PSHLTR(I,J)=-1.E6
        int_state%Q02(I,J)   = 0.
        int_state%Q10(I,J)   = 0.
        int_state%QSHLTR(I,J)= 0.
        int_state%T2(I,J)    = 0.
        int_state%TH02(I,J)  = 0.
        int_state%TH10(I,J)  = 0.
        int_state%TSHLTR(I,J)= 0.
        int_state%U10(I,J)   = 0.
        int_state%V10(I,J)   = 0.
      ENDDO
      ENDDO
!
      int_state%AVRAIN=3600*int_state%NHRS_PREC/(int_state%NPRECIP*int_state%DT_INT)
      int_state%AVCNVC=3600*int_state%NHRS_CLOD/(int_state%NPRECIP*int_state%DT_INT)
      int_state%ASRFC =3600*int_state%NHRS_SRFC/(int_state%NPHS   *int_state%DT_INT)
      int_state%ARDSW =3600*int_state%NHRS_RDSW/(int_state%NPHS   *int_state%DT_INT)
      int_state%ARDLW =3600*int_state%NHRS_RDLW/(int_state%NPHS   *int_state%DT_INT)
      int_state%APHTIM=3600*int_state%NHRS_HEAT/(int_state%NPHS   *int_state%DT_INT)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_INTERNAL_STATE_PHY
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_INTERNAL_STATE_PHY(IMP_STATE,INT_STATE)
!
!-----------------------------------------------------------------------
!***  UPDATE THE PHYSICS INTERNAL STATE WITH THE DATA
!***  IN THE IMPORT STATE SENT FROM THE DYNAMICS.
!-----------------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL,ONLY: IHALO,JHALO
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State)    ,INTENT(INOUT) :: IMP_STATE                     ! The Physics import state
      TYPE(INTERNAL_STATE),INTENT(INOUT) :: INT_STATE                     ! The Physics internal state
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
      CHARACTER(20) :: ARRAY_NAME
!
      TYPE(ESMF_Array) :: HOLD_ARRAY
!
!-----------------------------------------------------------------------
!***********************************************************************
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
!***  FOR EACH VARIABLE, EXTRACT ITS ESMF Array FROM THE IMPORT STATE.
!***  EACH Array CONTAINS A POINTER THAT POINTS AT THAT VARIABLE IN THE
!***  INTERNAL STATE OF THE COMPONENT FROM WHICH THE IMPORT STATE CAME.
!***  EXTRACT THE POINTER FROM THE ESMF Array AND USE IT TO UPDATE
!***  THIS COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST UPDATE THE 3-D Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   T   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='T'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Temperature Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract Temperature Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
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
            int_state%T(I,J,L)=HOLD_3D(II,JJ,L)                            !<-- Update Temperature
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   U   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='U'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract U Wind Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract U Wind Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
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
      ARRAY_NAME='V'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract V Wind Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract V Wind Pointer from Field"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
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
!- - - - - - - - - - - - - - - -   Q2  - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='Q2'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract TKE Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract TKE Pointer from Array"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
                        ,rc       =RC)
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
!- - - - - - - - - - - - - - -  OMGALF - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='OMGALF'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Omega-alpha Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract Omega-alpha Pointer from Array"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_3D                              &  !<-- Put the pointer here
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
            int_state%OMGALF(I,J,L)=HOLD_3D(II,JJ,L)                           !<-- Update OMGALF
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  UPDATE THE 2-D Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -  PD - - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='PD'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Sigma Domain Pressure Depth Array from Physics import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract Sigma Domain Pressure Depth Pointer from Array"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_2D                              &  !<-- Put the pointer here
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      JJ=-JHALO
      DO J=JMS,JME
        II=-IHALO
        JJ=JJ+1
        DO I=IMS,IME
          II=II+1
          int_state%PD(I,J)=HOLD_2D(II,JJ)                                 !<-- Update PD
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOW UPDATE THE 4-D TRACERS ARRAY.
!***  THE NUMBER OF 3-D QUANTITIES WITHIN IT IS KNOWN THROUGH THE
!***  NUM_TRACERS_TOTAL VARIABLE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - TRACERS - - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
      ARRAY_NAME='TRACERS'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract 4-D Tracers Array from Physics Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- State that holds the Array
                        ,itemName=ARRAY_NAME                            &  !<-- Name of the Array we want
                        ,array   =HOLD_ARRAY                            &  !<-- Put extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phy Update: Extract Tracers Pointer from Array"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ArrayGet(array    =HOLD_ARRAY                           &  !<-- Array that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=HOLD_4D                              &  !<-- Put the pointer here
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
!       WRITE(0,*)'PHYSICS UPDATE SUCCEEDED'
      ELSE
        WRITE(0,*)'PHYSICS UPDATE FAILED RC_UPD=',RC_UPD
      ENDIF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_INTERNAL_STATE_PHY
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
