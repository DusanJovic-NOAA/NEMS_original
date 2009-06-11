!-----------------------------------------------------------------------
!
      MODULE MODULE_MICROPHYSICS_NMM
!
!-----------------------------------------------------------------------
!
!***  THE MICROPHYSICS DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
!     USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
!                                  ,IMS,IME,JMS,JME                     &
!                                  ,ITS,ITE,JTS,JTE                     &
      USE MODULE_DM_PARALLEL,ONLY : ITS_B1,ITE_B1,ITE_B2                &
                                   ,ITS_B1_H1,ITE_B1_H1,ITE_B1_H2       &
                                   ,ITS_B1_H2,ITE_H2                    &
                                   ,JTS_B1,JTE_B1,JTE_B2                &
                                   ,JTS_B1_H1,JTE_B1_H1,JTE_B1_H2       &
                                   ,JTS_B1_H2,JTE_H2                    &
                                   ,MPI_COMM_COMP                       &
                                   ,MYPE_SHARE,NUM_TILES
!
      USE MODULE_CONSTANTS,ONLY : CICE,CLIQ,CPV,EP_1,EP_2,EPSILON,G     &
                                 ,P608,PSAT,R_D,R_V,RHOAIR0,RHOWATER    &
                                 ,SVPT0,XLF,XLV
!
      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: FERRIER_INIT,FPVS,GSMDRIVE,WSM3INIT
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE MICROPHYSICS OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: KESSLERSCHEME=1                   &
                                     ,LINSCHEME=2                       &
                                     ,WSM3SCHEME=3                      &
                                     ,ETAMPNEW=5                        & !<--- (Ferrier, WRF)
                                     ,THOMPSON=8
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR FERRIER MICROPHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     REAL,PRIVATE,SAVE ::  ABFR, CBFR, CIACW, CIACR, C_N0r0,            &
     &  CN0r0, CN0r_DMRmin, CN0r_DMRmax, CRACW, CRAUT, ESW0,            &
     &  RFmax, RQR_DR1, RQR_DR2, RQR_DR3, RQR_DRmin,                    &
     &  RQR_DRmax, RR_DRmin, RR_DR1, RR_DR2, RR_DR3, RR_DRmax
!
      INTEGER, PRIVATE,PARAMETER :: MY_T1=1, MY_T2=35
      REAL,PRIVATE,DIMENSION(MY_T1:MY_T2),SAVE :: MY_GROWTH_NMM
!
      REAL, PRIVATE,PARAMETER :: DMImin=.05e-3, DMImax=1.e-3,           &
     &      DelDMI=1.e-6,XMImin=1.e6*DMImin
      INTEGER, PUBLIC,PARAMETER :: XMImax=1.e6*DMImax, XMIexp=.0536,    &
     &                             MDImin=XMImin, MDImax=XMImax
      REAL, PRIVATE,DIMENSION(MDImin:MDImax) ::                         &
     &      ACCRI,SDENS,VSNOWI,VENTI1,VENTI2
!
      REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3, DMRmax=.45e-3,          &
     &      DelDMR=1.e-6,XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax
      INTEGER, PRIVATE,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
      REAL, PRIVATE,DIMENSION(MDRmin:MDRmax)::                          &
     &      ACCRR,MASSR,RRATE,VRAIN,VENTR1,VENTR2
!
      INTEGER, PRIVATE,PARAMETER :: Nrime=40
      REAL, DIMENSION(2:9,0:Nrime),PRIVATE,SAVE :: VEL_RF
!
      INTEGER,PARAMETER :: NX=7501
      REAL, PARAMETER :: XMIN=180.0,XMAX=330.0
      REAL, DIMENSION(NX),PRIVATE,SAVE :: TBPVS,TBPVS0
      REAL, PRIVATE,SAVE :: C1XPVS0,C2XPVS0,C1XPVS,C2XPVS
!
      REAL, PRIVATE,PARAMETER ::                                        &
!--- Physical constants follow:
     &   CP=1004.6, EPSQ=1.E-12, GRAV=9.806, RHOL=1000., RD=287.04      &
     &  ,RV=461.5, T0C=273.15, XLS=2.834E6                              &
!--- Derived physical constants follow:
     &  ,EPS=RD/RV, EPS1=RV/RD-1., EPSQ1=1.001*EPSQ                     &
     &  ,RCP=1./CP, RCPRV=RCP/RV, RGRAV=1./GRAV, RRHOL=1./RHOL          &
     &  ,XLS1=XLS*RCP, XLS2=XLS*XLS*RCPRV, XLS3=XLS*XLS/RV              &
!--- Constants specific to the parameterization follow:
!--- CLIMIT/CLIMIT1 are lower limits for treating accumulated precipitation
     &  ,CLIMIT=10.*EPSQ, CLIMIT1=-CLIMIT                               &
     &  ,C1=1./3.                                                       &
     &  ,DMR1=.1E-3, DMR2=.2E-3, DMR3=.32E-3                            &
     &  ,XMR1=1.e6*DMR1, XMR2=1.e6*DMR2, XMR3=1.e6*DMR3
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3
!
! ======================================================================
!--- Important tunable parameters that are exported to other modules
!  * RHgrd - threshold relative humidity for onset of condensation
!  * T_ICE - temperature (C) threshold at which all remaining liquid water
!            is glaciated to ice
!  * T_ICE_init - maximum temperature (C) at which ice nucleation occurs
!  * NLImax - maximum number concentrations (m**-3) of large ice (snow/graupel/sleet)
!  * NLImin - minimum number concentrations (m**-3) of large ice (snow/graupel/sleet)
!  * N0r0 - assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 0.45 mm
!  * N0rmin - minimum intercept (m**-4) for rain drops
!  * NCW - number concentrations of cloud droplets (m**-3)
!  * FLARGE1, FLARGE2 - number fraction of large ice to total (large+snow) ice
!          at T>0C and in presence of sublimation (FLARGE1), otherwise in
!          presence of ice saturated/supersaturated conditions
! ======================================================================
      REAL, PUBLIC,PARAMETER ::                                         &
     &  RHgrd=1.                                                        &
     & ,T_ICE=-40.                                                      &
     & ,T_ICEK=T0C+T_ICE                                                &
     & ,T_ICE_init=-15.                                                 &
     & ,NLImax=5.E3                                                     &
     & ,NLImin=1.E3                                                     &
     & ,N0r0=8.E6                                                       &
     & ,N0rmin=1.E4                                                     &
     & ,NCW=100.E6                                                      &
     & ,FLARGE1=1.                                                      &
     & ,FLARGE2=.03
!--- Other public variables passed to other routines:
      REAL,PUBLIC,SAVE ::  QAUT0
      REAL, PUBLIC,DIMENSION(MDImin:MDImax) :: MASSI
      INTEGER,PUBLIC,PARAMETER :: MICRO_RESTART=7501
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR WSM3 MICROPHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   REAL, PARAMETER, PRIVATE :: dtcldcr     = 120.
   REAL, PARAMETER, PRIVATE :: n0r = 8.e6
   REAL, PARAMETER, PRIVATE :: avtr = 841.9
   REAL, PARAMETER, PRIVATE :: bvtr = 0.8
   REAL, PARAMETER, PRIVATE :: r0 = .8e-5 ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER, PRIVATE :: peaut = .55   ! collection efficiency
   REAL, PARAMETER, PRIVATE :: xncr = 3.e8   ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER, PRIVATE :: xmyu = 1.718e-5 ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER, PRIVATE :: avts = 11.72
   REAL, PARAMETER, PRIVATE :: bvts = .41
   REAL, PARAMETER, PRIVATE :: n0smax =  1.e11 ! t=-90C unlimited
   REAL, PARAMETER, PRIVATE :: lamdarmax = 8.e4
   REAL, PARAMETER, PRIVATE :: lamdasmax = 1.e5
   REAL, PARAMETER, PRIVATE :: lamdagmax = 6.e4
   REAL, PARAMETER, PRIVATE :: betai = .6
   REAL, PARAMETER, PRIVATE :: xn0 = 1.e-2
   REAL, PARAMETER, PRIVATE :: dicon = 11.9
   REAL, PARAMETER, PRIVATE :: di0 = 12.9e-6
   REAL, PARAMETER, PRIVATE :: dimax = 500.e-6
   REAL, PARAMETER, PRIVATE :: n0s = 2.e6             ! temperature dependent n0s
   REAL, PARAMETER, PRIVATE :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER, PRIVATE :: qcrmin = 1.e-9
   REAL, SAVE ::                                     &
             qc0, qck1,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,&
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,   &
             precr1,precr2,xm0,xmmax,roqimax,bvts1,  &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,    &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r,&
             pidn0s,xlv1,                            &
             rslopermax,rslopesmax,rslopegmax,       &
             rsloperbmax,rslopesbmax,rslopegbmax,    &
             rsloper2max,rslopes2max,rslopeg2max,    &
             rsloper3max,rslopes3max,rslopeg3max
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE GSMDRIVE(ITIMESTEP,DT,NPHS,NUM_WATER                   &
                         ,DX,DY,SM,FIS                                  &
                         ,DSG2,SGML2,PDSG1,PSGML1,PT,PD                 &
                         ,T,Q,CWM,OMGALF,WATER                          &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG                 &
                         ,F_QV,F_QC,F_QR,F_QI,F_QS,F_QG                 &
                         ,PREC,ACPREC,AVRAIN                            &
                         ,MP_RESTART_STATE                              &
                         ,TBPVS_STATE,TBPVS0_STATE                      &
                         ,SPECIFIED,NESTED                              &
                         ,MICROPHYSICS                                  &
                         ,IDS,IDE,JDS,JDE,LM                            &
                         ,IMS,IME,JMS,JME                               &
                         ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    GSMDRIVE    MICROPHYSICS OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-26
!
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NMMB PHYSICS COMPONENT
!     AND THE WRF MICROPHYSICS DRIVER.
!
! PROGRAM HISTORY LOG:
!   02-03-26  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-07-31  BLACK      - BUILT INTO NMMB PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!
! USAGE: CALL GSMDRIVE FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!$$$
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: ITIMESTEP,NPHS,NUM_WATER                    &
                           ,IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE 
!
      INTEGER,INTENT(IN) :: P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
!
      REAL,INTENT(IN) :: DT,DX,DY,PT
!
      REAL,INTENT(INOUT) :: AVRAIN
!
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS,PD,SM
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: OMGALF
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACPREC,PREC
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM,Q,T     &
     &                                                     ,TRAIN
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,NUM_WATER),INTENT(INOUT) :: WATER
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: F_ICE     &
     &                                                     ,F_RAIN    &
     &                                                     ,F_RIMEF
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: SR
!
      CHARACTER(99),INTENT(IN) :: MICROPHYSICS
!
      LOGICAL,INTENT(IN) :: NESTED,SPECIFIED
!
      LOGICAL,INTENT(IN) :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
!***  State Variables for ETAMPNEW Microphysics 
!
      REAL,DIMENSION(:),INTENT(INOUT) :: MP_RESTART_STATE               &
     &                                  ,TBPVS_STATE,TBPVS0_STATE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,IJ,J,K,KFLIP,MP_PHYSICS,N,NTSD
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: LOWLYR
!
      REAL :: CAPA,DPL,DTPHS,FICE,FRAIN,PCPCOL,PDSL,PLYR,QI,QR,QW       &
             ,RDTPHS,RG,TNEW,WC
!
      REAL,DIMENSION(1:LM) :: QL,TL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUBOT,CUTOP,RAINNC,RAINNCV     &
                                        ,SNOWNC,SNOWNCV,XLAND 
!
      REAL,DIMENSION(IMS:IME,1:LM+1,JMS:JME) :: DZ,CWM_PHY              &
                                               ,F_ICE_PHY               &
                                               ,F_RAIN_PHY              &
                                               ,F_RIMEF_PHY             &
                                               ,P8W,P_PHY,PI_PHY        &
                                               ,RR,T_PHY,TH_PHY,WMID
!
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: WATER_TRANS
!
      LOGICAL :: WARM_RAIN,F_QT
!
!-----------------------------------------------------------------------
!
      INTEGER,DIMENSION(NUM_TILES) :: I_START,I_END,J_START,J_END
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(WATER_TRANS))THEN
        ALLOCATE(WATER_TRANS(IMS:IME,1:LM+1,JMS:JME,NUM_WATER))
      ENDIF
!
!-----------------------------------------------------------------------
!
!!!   NTSD=ITIMESTEP-1
      NTSD=ITIMESTEP
      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      CAPA=RD/CP
      RG=1./G
      AVRAIN=AVRAIN+1.
!
!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!!!   DO J=JTS_B1,JTE_B1
!!!   DO I=ITS_B1,ITE_B1
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,pdsl,kflip,dpl,plyr,ql,tl)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)
        P8W(I,LM+1,J)=PT
        LOWLYR(I,J)=1
        XLAND(I,J)=SM(I,J)+1.
!
!-----------------------------------------------------------------------
!***   FILL RAINNC WITH ZERO (NORMALLY CONTAINS THE NONCONVECTIVE
!***                          ACCUMULATED RAIN BUT NOT YET USED BY NMM)
!***   COULD BE OBTAINED FROM ACPREC AND CUPREC (ACPREC-CUPREC)
!-----------------------------------------------------------------------
!
        RAINNC(I,J)=0.
        SNOWNC(I,J)=0.
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1   ! We are moving down from the top in the flipped arrays
          KFLIP=LM+1-K ! The flipped index will thus be for the NMMB arrays
!
          DPL=PDSG1(KFLIP)+DSG2(KFLIP)*PDSL
          PLYR=SGML2(KFLIP)*PDSL+PSGML1(KFLIP)
          TL(K)=T(I,J,KFLIP)
          QL(K)=AMAX1(Q(I,J,KFLIP),EPSQ)
!
          RR(I,K,J)=PLYR/(R_D*TL(K)*(P608*QL(K)+1.))
          T_PHY(I,K,J)=TL(K)
          CWM_PHY(I,K,J)=CWM(I,J,KFLIP)
!zj flipped twice          WATER_TRANS(I,K,J,P_QV)=WATER(I,J,KFLIP,P_QV)
          F_ICE_PHY(I,K,J)=F_ICE(I,J,KFLIP)
          F_RAIN_PHY(I,K,J)=F_RAIN(I,J,KFLIP)
          F_RIMEF_PHY(I,K,J)=F_RIMEF(I,J,KFLIP)
!
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA
          TH_PHY(I,K,J)=TL(K)/PI_PHY(I,K,J)
          P8W(I,K,J)=P8W(I,K+1,J)+PDSG1(KFLIP)+DSG2(KFLIP)*PDSL
          P_PHY(I,K,J)=PLYR
          DZ(I,K,J)=DPL*RG/RR(I,K,J)
!
          WMID(I,K,J)=-OMGALF(I,J,KFLIP)*CP/(G*DT)
        ENDDO
!
        WMID(I,LM,J)=0.   !<---  W in the top model layer must equal zero for WSM3.
        WMID(I,LM+1,J)=0.
!
        F_ICE_PHY(I,LM+1,J)=0.
        F_RAIN_PHY(I,LM+1,J)=0.
        F_RIMEF_PHY(I,LM+1,J)=0.
!
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  SYNCHRONIZE MIXING RATIO IN WATER ARRAY WITH SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
      DO K=1,LM
        DO J=JMS,JME
          DO I=IMS,IME
            WATER(I,J,K,P_QV)=Q(I,J,K)/(1.-Q(I,J,K))
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  TRANSPOSE THE MOIST ARRAY (IJK) FOR THE PHYSICS (IKJ).
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(n,k,j,i,kflip)
!.......................................................................
      DO N=1,NUM_WATER
        DO K=1,LM
        KFLIP=LM+1-K
        DO J=JMS,JME
        DO I=IMS,IME
          WATER_TRANS(I,K,J,N)=WATER(I,J,KFLIP,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
!***  CALL MICROPHYSICS
!
!-----------------------------------------------------------------------
!
!!!   CALL SET_TILES(GRID,IDS+1,IDE-1,JDS+2,JDE-2,ITS,ITE,JTS,JTE)
      DO K=1,NUM_TILES
!!!     I_START(K)=ITS
!!!     I_END(K)=ITE
!!!     J_START(K)=JTS
!!!     J_END(K)=JTE
        I_START(K)=ITS_B1
        I_END(K)=ITE_B1
        J_START(K)=JTS_B1
        J_END(K)=JTE_B1
      ENDDO
!
!-----------------------------------------------------------------------
!***  TRANSLATE THE MICROPHYSICS OPTIONS IN THE CONFIG FILE TO THEIR
!***  ANALOGS IN THE WRF REGISTRY SO THAT THE WRF MICROPHYSICS DRIVER
!***  REMAINS UNTOUCHED.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(MICROPHYSICS))
        CASE ('fer')
          MP_PHYSICS=5
        CASE ('kes')
          MP_PHYSICS=1
        CASE ('lin')
          MP_PHYSICS=2
        CASE ('wsm3')
          MP_PHYSICS=3
        CASE ('tho')
          MP_PHYSICS=8
        CASE DEFAULT
          WRITE(0,*)' User selected MICROPHYSICS=',MICROPHYSICS
          WRITE(0,*)' Improper selection of Microphysics scheme in GSMDRIVE'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
      END SELECT
!
      CALL MICROPHYSICS_DRIVER(                                         &
     &                  TH=TH_PHY                                       &
     &                 ,RHO=RR,PI_PHY=PI_PHY,P=P_PHY                    &
     &                 ,RAINNC=RAINNC                                   &
     &                 ,RAINNCV=RAINNCV                                 &
     &                 ,SNOWNC=SNOWNC                                   &
     &                 ,SNOWNCV=SNOWNCV                                 &
     &                 ,DZ8W=DZ,P8W=P8W,W=WMID                          &
     &                 ,DT=DTPHS,DX=DX,DY=DY                            &
     &                 ,MP_PHYSICS=MP_PHYSICS                           &
     &                 ,SPECIFIED=SPECIFIED.OR.NESTED                   &
     &                 ,SPEC_ZONE=0,WARM_RAIN=WARM_RAIN                 &
     &                 ,XLAND=XLAND,ITIMESTEP=NTSD                      &
     &                 ,F_ICE_PHY=F_ICE_PHY                             &
     &                 ,F_RAIN_PHY=F_RAIN_PHY                           &
     &                 ,F_RIMEF_PHY=F_RIMEF_PHY                         &
     &                 ,LOWLYR=LOWLYR,SR=SR                             &
     &                 ,QV_CURR=WATER_TRANS(IMS,1,JMS,P_QV),F_QV=F_QV   &
     &                 ,QC_CURR=WATER_TRANS(IMS,1,JMS,P_QC),F_QC=F_QC   &
     &                 ,QR_CURR=WATER_TRANS(IMS,1,JMS,P_QR),F_QR=F_QR   &
     &                 ,QI_CURR=WATER_TRANS(IMS,1,JMS,P_QI),F_QI=F_QI   &
     &                 ,QS_CURR=WATER_TRANS(IMS,1,JMS,P_QS),F_QS=F_QS   &
     &                 ,QG_CURR=WATER_TRANS(IMS,1,JMS,P_QG),F_QG=F_QG   &
     &                 ,QT_CURR=CWM_PHY,F_QT=F_QT                       &
     &                 ,MP_RESTART_STATE=MP_RESTART_STATE               &
     &                 ,TBPVS_STATE=TBPVS_STATE                         &
     &                 ,TBPVS0_STATE=TBPVS0_STATE                       &
     &                 ,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=1,KDE=LM+1  &
     &                 ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=1,KME=LM+1  &
     &                 ,I_START=I_START,I_END=I_END                     &
     &                 ,J_START=J_START,J_END=J_END                     &
     &                 ,KTS=1,KTE=LM,NUM_TILES=NUM_TILES                &
     &                                                        )
!            
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE FOLLOWING MUST BE RECONCILED WHEN THREADING IS TURNED ON.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(ij)
!     DO IJ=1,NUM_TILES
!       CALL MICROPHYSICS_ZERO_OUT(                                     &
!                    WATER,N_MOIST,CONFIG_FLAGS                         &
!                   ,IDS,IDE,JDS,JDE,KDS,KDE                            &
!                   ,IMS,IME,JMS,JME,KMS,KME                            &
!                   ,GRID%I_START(IJ),GRID%I_END(IJ)                    &
!                   ,GRID%J_START(IJ),GRID%J_END(IJ)                    &
!                   ,KTS,KTE                                       )
!     ENDDO
!
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,kflip,tnew)
!.......................................................................
      DO J=JTS_B1,JTE_B1
        DO K=1,LM
        KFLIP=LM+1-K
        DO I=ITS_B1,ITE_B1
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          TNEW=TH_PHY(I,K,J)*PI_PHY(I,K,J)
          TRAIN(I,J,KFLIP)=TRAIN(I,J,KFLIP)+(TNEW-T(I,J,KFLIP))*RDTPHS
          T(I,J,KFLIP)=TNEW
          Q(I,J,KFLIP)=WATER_TRANS(I,K,J,P_QV)/(1.+WATER_TRANS(I,K,J,P_QV)) !To s.h.
          CWM(I,J,KFLIP)=CWM_PHY(I,K,J)
!
!      THE FOLLOWING LINE IS NO LONGER NEEDED BECAUSE CWM IS NOT IN WATER
!          CWM(I,J,KFLIP)=WATER_TRANS(I,K,J,P_QC)
!
        F_ICE(I,J,KFLIP)=F_ICE_PHY(I,K,J)
        F_RAIN(I,J,KFLIP)=F_RAIN_PHY(I,K,J)
        F_RIMEF(I,J,KFLIP)=F_RIMEF_PHY(I,K,J)
!
        ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  UPDATE PRECIPITATION
!-----------------------------------------------------------------------
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,pcpcol)
      DO J=JTS_B1,JTE_B1
      DO I=ITS_B1,ITE_B1
        PCPCOL=RAINNCV(I,J)*1.E-3
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  REFILL THE WATER ARRAY.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(n,j,k,i,kflip)
!.......................................................................
      DO N=2,NUM_WATER
        DO J=JMS,JME
        DO K=1,LM
        KFLIP=LM+1-K
        DO I=IMS,IME
          WATER(I,J,KFLIP,N)=WATER_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(WATER_TRANS)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GSMDRIVE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
SUBROUTINE microphysics_driver(                                          &
                       th, rho, pi_phy, p                                &
                      ,ht, dz8w, p8w, dt,dx,dy                           &
                      ,mp_physics, spec_zone                             &
                      ,specified, channel_switch                         &
                      ,warm_rain                                         &
                      ,xland,itimestep                                   &
                      ,f_ice_phy,f_rain_phy,f_rimef_phy                  &
                      ,lowlyr,sr                                         &
                      ,ids,ide, jds,jde, kds,kde                         &
                      ,ims,ime, jms,jme, kms,kme                         &
                      ,i_start,i_end,j_start,j_end,kts,kte,num_tiles     &
                      ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr   &
                      ,qni_curr                                          &
                      ,f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qni               &
                      ,qt_curr,f_qt                                      &
                      ,mp_restart_state,tbpvs_state,tbpvs0_state          & ! for etampnew
                      ,w ,z                                              &
                      ,rainnc, rainncv                                   &
                      ,snownc, snowncv                                   &
                      ,graupelnc, graupelncv                             &
                                                                         )
! Framework
!!!USE module_state_description, ONLY :                                  &
!!!                  KESSLERSCHEME, LINSCHEME, WSM3SCHEME, WSM5SCHEME    &
!!!                 ,WSM6SCHEME, ETAMPNEW, NCEPCLOUD3, NCEPCLOUD5, THOMPSON


! Model Layer
!!!USE module_model_constants
!!!USE module_wrf_error

! *** add new modules of schemes here

!!!USE module_mp_kessler
!!!USE module_mp_lin
!!!USE module_mp_ncloud3
!!!USE module_mp_ncloud5
!!!USE module_mp_wsm3
!!!USE module_mp_wsm5
!!!USE module_mp_wsm6
!!!USE module_mp_etanew
!!!USE module_mp_thompson
    
!----------------------------------------------------------------------
   ! This driver calls subroutines for the microphys.
   !
   ! Schemes
   !
   ! Kessler scheme
   ! Lin et al. (1983), Rutledge and Hobbs (1984)
   ! WRF Single-Moment 3-class, Hong, Dudhia and Chen (2004)
   ! WRF Single-Moment 5-class, Hong, Dudhia and Chen (2004)
   ! WRF Single-Moment 6-class, Lim and Hong (2003 WRF workshop)
   ! Eta Grid-scale Cloud and Precipitation scheme (EGCP01, Ferrier)
   ! NCEP cloud3, Hong et al. (1998) with some mod, Dudhia (1989)
   ! NCEP cloud5, Hong et al. (1998) with some mod, Rutledge and Hobbs (1984)
   ! 
!----------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------  
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------  
! In WRF, kms (smallest number) is the bottom level and kme (largest 
! number) is the top level.  In your scheme, if 1 is at the top level, 
! then you have to reverse the order in the k direction.
!                 
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
! Qni        cloud ice number concentration (#/kg)
!
!----------------------------------------------------------------------
!-- th        potential temperature    (K)
!-- moist_new     updated moisture array   (kg/kg)
!-- moist_old     Old moisture array       (kg/kg)
!-- rho           density of air           (kg/m^3)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- rhowater      water density                      (kg/m^3)
!-- rhosnow       snow density               (kg/m^3)
!-- F_ICE_PHY     Fraction of ice.
!-- F_RAIN_PHY    Fraction of rain.
!-- F_RIMEF_PHY   Mass ratio of rimed ice (rime factor)
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QS          species index for snow
!-- P_QG          species index for graupel
!-- P_QNI         species index for cloud ice number concentration
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!======================================================================

   INTEGER,    INTENT(IN   )    :: mp_physics
   LOGICAL,    INTENT(IN   )    :: specified
!
   INTEGER,      INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
   INTEGER,      INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
   INTEGER,      INTENT(IN   )    ::                         kts,kte
   INTEGER,      INTENT(IN   )    ::     itimestep,num_tiles,spec_zone
   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   LOGICAL,      INTENT(IN   )    ::   warm_rain
!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                    &
         INTENT(INOUT) ::                                     th
!

!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                    &
         INTENT(IN   ) ::                                             &
                                                                 rho, &
                                                                dz8w, &
                                                                 p8w, &
                                                              pi_phy, &
                                                               p
!

   REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::     &
                                     F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY

!

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)   :: XLAND

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(OUT)   :: SR

   REAL, INTENT(IN   ) :: dt,dx,dy

   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: LOWLYR

!
! Optional
!
   LOGICAL,  OPTIONAL,   INTENT(IN   )    :: channel_switch
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         OPTIONAL,                                                &
         INTENT(INOUT ) ::                                        &
                  w, z                                            & 
                 ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
                 ,qt_curr,qni_curr

!
   REAL, DIMENSION( ims:ime , jms:jme ),                          &
         INTENT(INOUT),                                           &
         OPTIONAL   ::                                            &
                                                           RAINNC &
                                                         ,RAINNCV &
                                                          ,SNOWNC &
                                                         ,SNOWNCV &
                                                       ,GRAUPELNC &
                                                      ,GRAUPELNCV

   REAL , DIMENSION( ims:ime , jms:jme ) , OPTIONAL ,             &
         INTENT(IN)   ::                                       ht

   REAL, DIMENSION (:), OPTIONAL, INTENT(INOUT) :: mp_restart_state &
                                         ,tbpvs_state,tbpvs0_state
!

   LOGICAL, OPTIONAL ::      f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qni,f_qt


! LOCAL  VAR

   INTEGER :: i,j,k,its,ite,jts,jte,ij,sz,n
   LOGICAL :: channel

!---------------------------------------------------------------------
!  check for microphysics type.  We need a clean way to 
!  specify these things!
!---------------------------------------------------------------------

   channel = .FALSE.
   IF ( PRESENT ( channel_switch ) ) channel = channel_switch

   if (mp_physics .eq. 0) return
   IF( specified ) THEN
     sz = spec_zone
   ELSE
     sz = 0
   ENDIF

!jaa   !$OMP PARALLEL DO   &
!jaa   !$OMP PRIVATE ( ij, its, ite, jts, jte, i,j,k,n )

   DO ij = 1 , num_tiles

       IF (channel) THEN
         its = max(i_start(ij),ids)
         ite = min(i_end(ij),ide-1)
       ELSE
         its = max(i_start(ij),ids+sz)
         ite = min(i_end(ij),ide-1-sz)
       ENDIF
       jts = max(j_start(ij),jds+sz)
       jte = min(j_end(ij),jde-1-sz)

!-----------

     micro_select: SELECT CASE(mp_physics)

!!!     CASE (KESSLERSCHEME)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling kessler' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT( QC_CURR ) .AND.  &
!!!                                        PRESENT( QR_CURR ) .AND.  &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!                                        PRESENT( Z       ))  THEN
!!!            CALL kessler(                                        &
!!!               T=th                                              &
!!!              ,QV=qv_curr                                        &
!!!              ,QC=qc_curr                                        &
!!!              ,QR=qr_curr                                        &
!!!              ,RHO=rho, PII=pi_phy,DT_IN=dt, Z=z, XLV=xlv, CP=cp &
!!!              ,EP2=ep_2,SVP1=svp1,SVP2=svp2                      &
!!!              ,SVP3=svp3,SVPT0=svpt0,RHOWATER=rhowater           &
!!!              ,DZ8W=dz8w                                         &
!!!              ,RAINNC=rainnc,RAINNCV=rainncv                     &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )
!!!          ELSE 
!!!             CALL wrf_error_fatal ( 'arguments not present for calling kessler' )
!!!          ENDIF

!
!!!     CASE (THOMPSON)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling thompson et al' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!               PRESENT( QS_CURR ) .AND. PRESENT ( QG_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!                                        PRESENT ( QNI_CURR ).AND.  &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( Z       ) .AND. PRESENT ( W       )  ) THEN
!!!          CALL mp_thomp(                              &
!!!                  ITIMESTEP=itimestep,                &
!!!                  TH=th,                              &
!!!                  QV=qv_curr,                         &
!!!                  QL=qc_curr,                         &
!!!                  QR=qr_curr,                         &
!!!                  QI=qi_curr,                         &
!!!                  QS=qs_curr,                         &
!!!                  QG=qg_curr,                         &
!!!                  QNI=qni_curr,                       &
!!!                  RHO=rho,                            &
!!!                  PII=pi_phy,                         &
!!!                  P=p,                                &
!!!                  DT_IN=dt,                           &
!!!                  Z=z,                                &
!!!                  HT=ht,                              &
!!!                  DZ8W=dz8w,                          &
!!!                  W=w                                 &
!!!                 ,RAINNC=RAINNC                       &
!!!                 ,RAINNCV=RAINNCV                     &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )
!!!          ELSE 
!!!             CALL wrf_error_fatal ( 'arguments not present for calling thompson_et_al' )
!!!          ENDIF
!
!!!     CASE (LINSCHEME)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling lin_et_al' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!               PRESENT( QS_CURR )                           .AND.  &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( Z       ) .AND. PRESENT ( W       )  ) THEN
!!!            CALL lin_et_al(                                      &
!!!               TH=th                                             &
!!!              ,QV=qv_curr                                        &
!!!              ,QL=qc_curr                                        &
!!!              ,QR=qr_curr                                        &
!!!              ,QI=qi_curr                                        &
!!!              ,QS=qs_curr                                        &
!!!              ,RHO=rho, PII=pi_phy, P=p, DT_IN=dt, Z=z           &
!!!              ,HT=ht, DZ8W=dz8w, GRAV=G,  CP=cp                  &
!!!              ,RAIR=r_d, RVAPOR=R_v                              &
!!!              ,XLS=xls, XLV=xlv, XLF=xlf                         &
!!!              ,RHOWATER=rhowater, RHOSNOW=rhosnow                &
!!!              ,EP2=ep_2,SVP1=svp1,SVP2=svp2                      &
!!!              ,SVP3=svp3,SVPT0=svpt0                             &
!!!              ,RAINNC=rainnc, RAINNCV=rainncv                    &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!              ,F_QG=f_qg                                         &
!!!              ,QG=qg_curr                                        &
!!!                                                                 )
!!!          ELSE 
!!!             CALL wrf_error_fatal ( 'arguments not present for calling lin_et_al' )
!!!          ENDIF

        CASE (WSM3SCHEME)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling wsm3' )
             IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
                  PRESENT( QR_CURR ) .AND.                            &
                  PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
                  PRESENT( SNOWNC  ) .AND. PRESENT ( SNOWNCV ) .AND.  &
                  PRESENT( W       )                            ) THEN
             CALL wsm3(                                             &
                  TH=th                                             &
                 ,Q=qv_curr                                         &
                 ,QCI=qc_curr                                       &
                 ,QRS=qr_curr                                       &
                 ,W=w,DEN=rho,PII=pi_phy,P=p,DELZ=dz8w              &
                 ,DELT=dt,G=g,CPD=cp,CPV=cpv                        &
                 ,RD=r_d,RV=r_v,T0C=svpt0                           &
                 ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
                 ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
                 ,DEN0=rhoair0, DENR=rhowater                       &
                 ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
                 ,RAIN=rainnc ,RAINNCV=rainncv                      &
                 ,SNOW=snownc ,SNOWNCV=snowncv                      &
                 ,SR=sr                                             &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                 ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                                                                    )
             ELSE 
!!!             CALL wrf_error_fatal ( 'arguments not present for calling wsm3' )
                WRITE(0,*)'arguments not present for calling wsm3'
             ENDIF

!!!     CASE (WSM5SCHEME)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling wsm5' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!               PRESENT( QS_CURR ) .AND.                            &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( SNOWNC  ) .AND. PRESENT ( SNOWNCV ) .AND.  &
!!!               PRESENT( W       )                            ) THEN
!!!          CALL wsm5(                                             &
!!!               TH=th                                             &
!!!              ,Q=qv_curr                                         &
!!!              ,QC=qc_curr                                        &
!!!              ,QR=qr_curr                                        &
!!!              ,QI=qi_curr                                        &
!!!              ,QS=qs_curr                                        &
!!!              ,W=w,DEN=rho,PII=pi_phy,P=p,DELZ=dz8w              &
!!!              ,DELT=dt,G=g,CPD=cp,CPV=cpv                        &
!!!              ,RD=r_d,RV=r_v,T0C=svpt0                           &
!!!              ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
!!!              ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
!!!              ,DEN0=rhoair0, DENR=rhowater                       &
!!!              ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
!!!              ,RAIN=rainnc ,RAINNCV=rainncv                      &
!!!              ,SNOW=snownc ,SNOWNCV=snowncv                      &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )
!!!          ELSE
!!!             CALL wrf_error_fatal ( 'arguments not present for calling wsm5' )
!!!          ENDIF

!!!     CASE (WSM6SCHEME)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling wsm6' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!               PRESENT( QS_CURR ) .AND. PRESENT ( QG_CURR ) .AND.  &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( SNOWNC  ) .AND. PRESENT ( SNOWNCV ) .AND.  &
!!!               PRESENT( GRAUPELNC  ) .AND. PRESENT ( GRAUPELNCV ) .AND.  &
!!!               PRESENT( W       )                            ) THEN
!!!          CALL wsm6(                                             &
!!!               TH=th                                             &
!!!              ,Q=qv_curr                                         &
!!!              ,QC=qc_curr                                        &
!!!              ,QR=qr_curr                                        &
!!!              ,QI=qi_curr                                        &
!!!              ,QS=qs_curr                                        &
!!!              ,QG=qg_curr                                        &
!!!              ,W=w,DEN=rho,PII=pi_phy,P=p,DELZ=dz8w              &
!!!              ,DELT=dt,G=g,CPD=cp,CPV=cpv                        &
!!!              ,RD=r_d,RV=r_v,T0C=svpt0                           &
!!!              ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
!!!              ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
!!!              ,DEN0=rhoair0, DENR=rhowater                       &
!!!              ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
!!!              ,RAIN=rainnc ,RAINNCV=rainncv                      &
!!!              ,SNOW=snownc ,SNOWNCV=snowncv                      &
!!!              ,GRAUPEL=graupelnc ,GRAUPELNCV=graupelncv          &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )
!!!          ELSE
!!!             CALL wrf_error_fatal ( 'arguments not present for calling wsm6' )
!!!          ENDIF

        CASE (ETAMPNEW)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling etampnew')

             IF ( PRESENT( qv_curr ) .AND. PRESENT( qt_curr ) .AND. &
                  PRESENT( qc_curr ) .AND. PRESENT( qr_curr ) .AND. &
                  PRESENT( qs_curr ) .AND.                          &
                  PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
                  PRESENT( mp_restart_state )                  .AND. &
                  PRESENT( tbpvs_state )                      .AND. &
                  PRESENT( tbpvs0_state )                       ) THEN
               CALL ETAMP_NEW(                                      &
                  ITIMESTEP=itimestep,DT=dt,DX=dx,DY=dy             &
                 ,DZ8W=dz8w,RHO_PHY=rho,P_PHY=p,PI_PHY=pi_phy,TH_PHY=th &
                 ,QV=qv_curr                                        &
                 ,QC=qc_curr                                        &
                 ,QS=qs_curr                                        &
                 ,QR=qr_curr                                        &
                 ,QT=qt_curr                                        &
                 ,LOWLYR=LOWLYR,SR=SR                               &
                 ,F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY         &
                 ,F_RIMEF_PHY=F_RIMEF_PHY                           &
                 ,RAINNC=rainnc,RAINNCV=rainncv                     &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                 ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                 ,MP_RESTART_STATE=mp_restart_state                 &
                 ,TBPVS_STATE=tbpvs_state,TBPVS0_STATE=tbpvs0_state &
                                                                    )
             ELSE
!!!             CALL wrf_error_fatal ( 'arguments not present for calling etampnew' )
                WRITE(0,*)' arguments not present for calling etampnew'
!!!             CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
                CALL NMMB_FINALIZE
             ENDIF

!!!     CASE (NCEPCLOUD3)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling ncloud3' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND.                            &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( W       )                            ) THEN
!!!          CALL ncloud3(                                          &
!!!               TH=th                                             &
!!!              ,Q=qv_curr                                         &
!!!              ,QCI=qc_curr                                       &
!!!              ,QRS=qr_curr                                       &
!!!              ,W=w, DEN=rho, PII=pi_phy, P=p, DELZ=dz8w          &
!!!              ,DELT=dt,G=g,CPD=cp,CPV=cpv                        &
!!!              ,RD=r_d,RV=r_v,T0C=SVPT0                           &
!!!              ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
!!!              ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
!!!              ,DEN0=rhoair0, DENR=rhowater                       &
!!!              ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
!!!              ,RAIN=RAINNC,RAINNCV=RAINNCV                       &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )

!!!          ELSE
!!!             CALL wrf_error_fatal ( 'arguments not present for calling ncepcloud3' )
!!!          ENDIF

!!!     CASE (NCEPCLOUD5)
!!!          CALL wrf_debug ( 100 , 'microphysics_driver: calling ncloud5' )
!!!          IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!!!               PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!!!               PRESENT( QS_CURR ) .AND. PRESENT ( QG_CURR ) .AND.  &
!!!               PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!!!               PRESENT( W       )                            ) THEN
!!!          CALL ncloud5(                                          &
!!!               TH=th                                             &
!!!              ,Q=qv_curr                                         &
!!!              ,QC=qc_curr                                        &
!!!              ,QR=qr_curr                                        &
!!!              ,QI=qi_curr                                        &
!!!              ,QS=qs_curr                                        &
!!!              ,W=w, DEN=rho, PII=pi_phy, P=p, DELZ=dz8w          &
!!!              ,DELT=dt,G=g,CPD=cp,CPV=cpv                        &
!!!              ,RD=r_d,RV=r_v,T0C=SVPT0                           &
!!!              ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
!!!              ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
!!!              ,DEN0=rhoair0, DENR=rhowater                       &
!!!              ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
!!!              ,RAIN=RAINNC,RAINNCV=RAINNCV                       &
!!!              ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
!!!              ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
!!!              ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
!!!                                                                 )
!!!          ELSE
!!!             CALL wrf_error_fatal ( 'arguments not present for calling ncepcloud5' )
!!!          ENDIF


      CASE DEFAULT 

!!!      WRITE( wrf_err_message , * ) 'The microphysics option does not exist: mp_physics = ', mp_physics
!!!      CALL wrf_error_fatal ( wrf_err_message )
         WRITE(0,*)' The microphysics option does not exist: mp_physics = ', mp_physics
!!!      CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
         CALL NMMB_FINALIZE

      END SELECT micro_select 

   ENDDO
!jaa   !$OMP END PARALLEL DO

!!!CALL wrf_debug ( 200 , 'microphysics_driver: returning from' )

   RETURN

   END SUBROUTINE microphysics_driver

!!!END MODULE module_microphysics_driver

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE ETAMP_NEW (itimestep,DT,DX,DY,                         &
     &                      dz8w,rho_phy,p_phy,pi_phy,th_phy,qv,qt,     &
     &                      LOWLYR,SR,                                  &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,           &
     &                      QC,QR,QS,                                   &
     &                      mp_restart_state,tbpvs_state,tbpvs0_state,  &
     &                      RAINNC,RAINNCV,                             &
     &                      ids,ide, jds,jde, kds,kde,		        &
     &                      ims,ime, jms,jme, kms,kme,		        &
     &                      its,ite, jts,jte, kts,kte )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40

      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,ITIMESTEP

      REAL, INTENT(IN) 	    :: DT,DX,DY
      REAL, INTENT(IN),     DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      th_phy,qv,qt
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      qc,qr,qs
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
     &                                                   RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR
!
      REAL,DIMENSION(*),INTENT(INOUT) :: MP_RESTART_STATE
!
      REAL,DIMENSION(nx),INTENT(INOUT) :: TBPVS_STATE,TBPVS0_STATE
!
      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR

!-----------------------------------------------------------------------
!     LOCAL VARS
!-----------------------------------------------------------------------

!     NSTATS,QMAX,QTOT are diagnostic vars

      INTEGER,DIMENSION(ITLO:ITHI,4) :: NSTATS
      REAL,   DIMENSION(ITLO:ITHI,5) :: QMAX
      REAL,   DIMENSION(ITLO:ITHI,22):: QTOT

!     SOME VARS WILL BE USED FOR DATA ASSIMILATION (DON'T NEED THEM NOW). 
!     THEY ARE TREATED AS LOCAL VARS, BUT WILL BECOME STATE VARS IN THE 
!     FUTURE. SO, WE DECLARED THEM AS MEMORY SIZES FOR THE FUTURE USE

!     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
!     the microphysics scheme. Instead, they will be used by Eta precip 
!     assimilation.

      REAL,  DIMENSION( ims:ime, kms:kme, jms:jme ) ::                  &
     &       TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      REAL,  DIMENSION(its:ite, kts:kte, jts:jte):: t_phy

      INTEGER :: I,J,K,KFLIP
      REAL :: WC
!
!-----------------------------------------------------------------------
!**********************************************************************
!-----------------------------------------------------------------------
!
      MY_GROWTH_NMM(MY_T1:MY_T2)=MP_RESTART_STATE(MY_T1:MY_T2)
!
      C1XPVS0=MP_RESTART_STATE(MY_T2+1)
      C2XPVS0=MP_RESTART_STATE(MY_T2+2)
      C1XPVS =MP_RESTART_STATE(MY_T2+3)
      C2XPVS =MP_RESTART_STATE(MY_T2+4)
      CIACW  =MP_RESTART_STATE(MY_T2+5)
      CIACR  =MP_RESTART_STATE(MY_T2+6)
      CRACW  =MP_RESTART_STATE(MY_T2+7)
      CRAUT  =MP_RESTART_STATE(MY_T2+8)
!
      TBPVS(1:NX) =TBPVS_STATE(1:NX)
      TBPVS0(1:NX)=TBPVS0_STATE(1:NX)
!
!     write(0,*)' in ETAMP_NEW its=',its,' ite=',ite,' jts=',jts,' jte=',jte,' kts=',kts,' kte=',kte
!     write(0,*)' ims=',ims,' ime=',ime,' jms=',jms,' jme=',jme,' kms=',kms,' kme=',kme
!.......................................................................
!$omp parallel do private(j,k,i)
!.......................................................................
      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
        t_phy(i,k,j) = th_phy(i,k,j)*pi_phy(i,k,j)
        qv(i,k,j)=qv(i,k,j)/(1.+qv(i,k,j)) !Convert to specific humidity
!     if(i==101.and.j==163)then
!       write(0,*)' enter etamp_new k=',k,' t=',t_phy(i,k,j),' th=',th_phy(i,k,j),' pii=',pi_phy(i,k,j)
!     endif
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

!     initial diagnostic variables and data assimilation vars
!     (will need to delete this part in the future)

      DO k = 1,4
      DO i = ITLO,ITHI
         NSTATS(i,k)=0. 
      ENDDO
      ENDDO

      DO k = 1,5
      DO i = ITLO,ITHI
         QMAX(i,k)=0.
      ENDDO
      ENDDO

      DO k = 1,22
      DO i = ITLO,ITHI
         QTOT(i,k)=0.
      ENDDO
      ENDDO

! initial data assimilation vars (will need to delete this part in the future)

!.......................................................................
!$omp parallel do private(j,k,i)
!.......................................................................
      DO j = jts,jte
       DO i = its,ite
         ACPREC(i,j)=0.
         APREC (i,j)=0.
         PREC  (i,j)=0.
         SR    (i,j)=0.
       ENDDO
       DO k = kts,kte
       DO i = its,ite
	 TLATGS_PHY (i,k,j)=0.
	 TRAIN_PHY  (i,k,j)=0.
       ENDDO
       ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

!-----------------------------------------------------------------------

      CALL EGCP01DRV(DT,LOWLYR,                                         &
     &               APREC,PREC,ACPREC,SR,NSTATS,QMAX,QTOT,	        &
     &               dz8w,rho_phy,qt,t_phy,qv,F_ICE_PHY,P_PHY,          &
     &               F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY,TRAIN_PHY,       &
     &               ids,ide, jds,jde, kds,kde,		                &
     &               ims,ime, jms,jme, kms,kme,		                &
     &               its,ite, jts,jte, kts,kte		          )
!-----------------------------------------------------------------------

!.......................................................................
!$omp parallel do private(j,k,i,wc)
!.......................................................................
     DO j = jts,jte
        DO k = kts,kte
	DO i = its,ite
  	   th_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)
           qv(i,k,j)=qv(i,k,j)/(1.-qv(i,k,j))  !Convert to mixing ratio
           WC=qt(I,K,J)
           QS(I,K,J)=0.
           QR(I,K,J)=0.
           QC(I,K,J)=0.
           IF(F_ICE_PHY(I,K,J)>=1.)THEN
             QS(I,K,J)=WC
           ELSEIF(F_ICE_PHY(I,K,J)<=0.)THEN
             QC(I,K,J)=WC
           ELSE
             QS(I,K,J)=F_ICE_PHY(I,K,J)*WC
             QC(I,K,J)=WC-QS(I,K,J)
           ENDIF
!
           IF(QC(I,K,J)>0..AND.F_RAIN_PHY(I,K,J)>0.)THEN
             IF(F_RAIN_PHY(I,K,J).GE.1.)THEN
               QR(I,K,J)=QC(I,K,J)
               QC(I,K,J)=0.
             ELSE
               QR(I,K,J)=F_RAIN_PHY(I,K,J)*QC(I,K,J)
               QC(I,K,J)=QC(I,K,J)-QR(I,K,J)
             ENDIF
           ENDIF
	ENDDO
        ENDDO
     ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

! 
! update rain (from m to mm)

       DO j=jts,jte
       DO i=its,ite
          RAINNC(i,j)=APREC(i,j)*1000.+RAINNC(i,j)
          RAINNCV(i,j)=APREC(i,j)*1000.
       ENDDO
       ENDDO
!
     MP_RESTART_STATE(MY_T1:MY_T2)=MY_GROWTH_NMM(MY_T1:MY_T2)
     MP_RESTART_STATE(MY_T2+1)=C1XPVS0
     MP_RESTART_STATE(MY_T2+2)=C2XPVS0
     MP_RESTART_STATE(MY_T2+3)=C1XPVS
     MP_RESTART_STATE(MY_T2+4)=C2XPVS
     MP_RESTART_STATE(MY_T2+5)=CIACW
     MP_RESTART_STATE(MY_T2+6)=CIACR
     MP_RESTART_STATE(MY_T2+7)=CRACW
     MP_RESTART_STATE(MY_T2+8)=CRAUT
!
     TBPVS_STATE(1:NX) =TBPVS(1:NX)
     TBPVS0_STATE(1:NX)=TBPVS0(1:NX)

!-----------------------------------------------------------------------

  END SUBROUTINE ETAMP_NEW

!-----------------------------------------------------------------------
      SUBROUTINE EGCP01DRV(                                             &
     &  DTPH,LOWLYR,APREC,PREC,ACPREC,SR,                               &
     &  NSTATS,QMAX,QTOT,                                               &
     &  dz8w,RHO_PHY,CWM_PHY,T_PHY,Q_PHY,F_ICE_PHY,P_PHY,               &
     &  F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY,TRAIN_PHY,                    &
     &  ids,ide, jds,jde, kds,kde,                                      &
     &  ims,ime, jms,jme, kms,kme,                                      &
     &  its,ite, jts,jte, kts,kte)
!-----------------------------------------------------------------------
! DTPH           Physics time step (s)
! CWM_PHY (qt)   Mixing ratio of the total condensate. kg/kg
! Q_PHY          Mixing ratio of water vapor. kg/kg
! F_RAIN_PHY     Fraction of rain. 
! F_ICE_PHY      Fraction of ice.
! F_RIMEF_PHY    Mass ratio of rimed ice (rime factor).
!
!TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related the
!micrphysics sechme. Instead, they will be used by Eta precip assimilation.
!
!NSTATS,QMAX,QTOT are used for diagnosis purposes.
!
!-----------------------------------------------------------------------
!--- Variables APREC,PREC,ACPREC,SR are calculated for precip assimilation
!    and/or ZHAO's scheme in Eta and are not required by this microphysics 
!    scheme itself.  
!--- NSTATS,QMAX,QTOT are used for diagnosis purposes only.  They will be 
!    printed out when PRINT_diag is true.
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40
      LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
!     VARIABLES PASSED IN/OUT
      INTEGER,INTENT(IN ) :: ids,ide, jds,jde, kds,kde                  &
     &                      ,ims,ime, jms,jme, kms,kme                  &
     &                      ,its,ite, jts,jte, kts,kte
      REAL,INTENT(IN) :: DTPH
      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR
      INTEGER,DIMENSION(ITLO:ITHI,4),INTENT(INOUT) :: NSTATS
      REAL,DIMENSION(ITLO:ITHI,5),INTENT(INOUT) :: QMAX
      REAL,DIMENSION(ITLO:ITHI,22),INTENT(INOUT) :: QTOT
      REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) ::                  &
     &                         APREC,PREC,ACPREC,SR
      REAL,DIMENSION( its:ite, kts:kte, jts:jte ),INTENT(INOUT) :: t_phy
      REAL,DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN) ::         &
     &                                             dz8w,P_PHY,RHO_PHY
      REAL,DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(INOUT) ::      &
     &   CWM_PHY, F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,TLATGS_PHY           &
     &   ,Q_PHY,TRAIN_PHY
!
!-----------------------------------------------------------------------
!LOCAL VARIABLES
!-----------------------------------------------------------------------
!
#define CACHE_FRIENDLY_MP_ETANEW
#ifdef CACHE_FRIENDLY_MP_ETANEW
#  define TEMP_DIMS  kts:kte,its:ite,jts:jte
#  define TEMP_DEX   L,I,J
#else
#  define TEMP_DIMS  its:ite,jts:jte,kts:kte
#  define TEMP_DEX   I,J,L
#endif
!
      INTEGER :: LSFC,I,J,I_index,J_index,L,K,KFLIP
      REAL,DIMENSION(TEMP_DIMS) :: CWM,T,Q,TRAIN,TLATGS,P
      REAL,DIMENSION(kts:kte,its:ite,jts:jte) :: F_ice,F_rain,F_RimeF       
      INTEGER,DIMENSION(its:ite,jts:jte) :: LMH
      REAL :: TC,WC,QI,QR,QW,Fice,Frain,DUM,ASNOW,ARAIN
      REAL,DIMENSION(kts:kte) :: P_col,Q_col,T_col,QV_col,WC_col,       &
         RimeF_col,QI_col,QR_col,QW_col, THICK_col,DPCOL 
      REAL,DIMENSION(2) :: PRECtot,PRECmax
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,i,l,kflip)
!.......................................................................
        DO J=JTS,JTE    
         DO I=ITS,ITE  
          LMH(I,J) = KTE-LOWLYR(I,J)+1
           DO L=KTS,KTE
             KFLIP=KTE+1-L
             CWM(TEMP_DEX)=CWM_PHY(I,KFLIP,J)
             T(TEMP_DEX)=T_PHY(I,KFLIP,J)
             Q(TEMP_DEX)=Q_PHY(I,KFLIP,J)
             P(TEMP_DEX)=P_PHY(I,KFLIP,J)
             TLATGS(TEMP_DEX)=TLATGS_PHY(I,KFLIP,J)
             TRAIN(TEMP_DEX)=TRAIN_PHY(I,KFLIP,J)
             F_ice(L,I,J)=F_ice_PHY(I,KFLIP,J)
             F_rain(L,I,J)=F_rain_PHY(I,KFLIP,J)
             F_RimeF(L,I,J)=F_RimeF_PHY(I,KFLIP,J)
           ENDDO
         ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do                                                       &
!$omp private (j,i,k,lsfc,kflip,dpcol,l,p_col,thick_col,t_col,tc,qv_col,&
!$omp          wc_col,wc,qi,qr,qw,fice,frain,rimef_col,qi_col,qr_col,   &
!$omp          qw_col,i_index,j_index,arain,asnow,dum,prectot,precmax,  &
!$omp          qmax,qtot,nstats),SCHEDULE(dynamic)  
!.......................................................................
       DO J=JTS,JTE    
        DO I=ITS,ITE  
          LSFC=LMH(I,J)                      ! "L" of surface
!
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            DPCOL(K)=RHO_PHY(I,KFLIP,J)*GRAV*dz8w(I,KFLIP,J)
          ENDDO
!   
   !
   !--- Initialize column data (1D arrays)
   !
        L=1
        IF (CWM(TEMP_DEX) .LE. EPSQ) CWM(TEMP_DEX)=EPSQ
          F_ice(1,I,J)=1.
          F_rain(1,I,J)=0.
          F_RimeF(1,I,J)=1.
          DO L=1,LSFC
      !
      !--- Pressure (Pa) = (Psfc-Ptop)*(ETA/ETA_sfc)+Ptop
      !
            P_col(L)=P(TEMP_DEX)
      !
      !--- Layer thickness = RHO*DZ = -DP/G = (Psfc-Ptop)*D_ETA/(G*ETA_sfc)
      !
            THICK_col(L)=DPCOL(L)*RGRAV
            T_col(L)=T(TEMP_DEX)
            TC=T_col(L)-T0C
            QV_col(L)=max(EPSQ, Q(TEMP_DEX))
            IF (CWM(TEMP_DEX) .LE. EPSQ1) THEN
              WC_col(L)=0.
              IF (TC .LT. T_ICE) THEN
                F_ice(L,I,J)=1.
              ELSE
                F_ice(L,I,J)=0.
              ENDIF
              F_rain(L,I,J)=0.
              F_RimeF(L,I,J)=1.
            ELSE
              WC_col(L)=CWM(TEMP_DEX)
            ENDIF
      !
      !--- Determine composition of condensate in terms of 
      !      cloud water, ice, & rain
      !
            WC=WC_col(L)
            QI=0.
            QR=0.
            QW=0.
            Fice=F_ice(L,I,J)
            Frain=F_rain(L,I,J)
            IF (Fice .GE. 1.) THEN
              QI=WC
            ELSE IF (Fice .LE. 0.) THEN
              QW=WC
            ELSE
              QI=Fice*WC
              QW=WC-QI
            ENDIF
            IF (QW.GT.0. .AND. Frain.GT.0.) THEN
              IF (Frain .GE. 1.) THEN
                QR=QW
                QW=0.
              ELSE
                QR=Frain*QW
                QW=QW-QR
              ENDIF
            ENDIF
!ratko-wrf
            IF (QI .LE. 0.) F_RimeF(L,I,J)=1.
!ratko-wrf
            RimeF_col(L)=F_RimeF(L,I,J)               ! (real)
            QI_col(L)=QI
            QR_col(L)=QR
            QW_col(L)=QW
          ENDDO
!
!#######################################################################
   !
   !--- Perform the microphysical calculations in this column
   !
          I_index=I
          J_index=J
       CALL EGCP01COLUMN ( ARAIN, ASNOW, DTPH, I_index, J_index, LSFC,  &
     & P_col, QI_col, QR_col, QV_col, QW_col, RimeF_col, T_col,         &
     & THICK_col, WC_col,KTS,KTE,NSTATS,QMAX,QTOT )


   !
!#######################################################################
!
   !
   !--- Update storage arrays
   !
          DO L=1,LSFC
            TRAIN(TEMP_DEX)=(T_col(L)-T(TEMP_DEX))/DTPH
            TLATGS(TEMP_DEX)=T_col(L)-T(TEMP_DEX)
            T(TEMP_DEX)=T_col(L)
            Q(TEMP_DEX)=QV_col(L)
            CWM(TEMP_DEX)=WC_col(L)
      !
      !--- REAL*4 array storage
      !
!ratko-wrf            F_RimeF(L,I,J)=MAX(1., RimeF_col(L))
            IF (QI_col(L) .LE. EPSQ) THEN
              F_ice(L,I,J)=0.
              IF (T_col(L) .LT. T_ICEK) F_ice(L,I,J)=1.
!ratko-wrf
              F_RimeF(L,I,J)=1.
!ratko-wrf
            ELSE
              F_ice(L,I,J)=MAX( 0., MIN(1., QI_col(L)/WC_col(L)) )
!ratko-wrf
              F_RimeF(L,I,J)=MAX(1., RimeF_col(L))
!ratko-wrf
            ENDIF
            IF (QR_col(L) .LE. EPSQ) THEN
              DUM=0
            ELSE
              DUM=QR_col(L)/(QR_col(L)+QW_col(L))
            ENDIF
            F_rain(L,I,J)=DUM
      !
          ENDDO
   !
   !--- Update accumulated precipitation statistics
   !
   !--- Surface precipitation statistics; SR is fraction of surface 
   !    precipitation (if >0) associated with snow
   !
        APREC(I,J)=(ARAIN+ASNOW)*RRHOL       ! Accumulated surface precip (depth in m)  !<--- Ying
        PREC(I,J)=PREC(I,J)+APREC(I,J)
        ACPREC(I,J)=ACPREC(I,J)+APREC(I,J)
        IF(APREC(I,J) .LT. 1.E-8) THEN
          SR(I,J)=0.
        ELSE
          SR(I,J)=RRHOL*ASNOW/APREC(I,J)
        ENDIF
   !
   !--- Debug statistics 
   !
        IF (PRINT_diag) THEN
          PRECtot(1)=PRECtot(1)+ARAIN
          PRECtot(2)=PRECtot(2)+ASNOW
          PRECmax(1)=MAX(PRECmax(1), ARAIN)
          PRECmax(2)=MAX(PRECmax(2), ASNOW)
        ENDIF


!#######################################################################
!#######################################################################
!
    enddo                          ! End "I" loop
    enddo                          ! End "J" loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do private(j,i,l,kflip)
!.......................................................................
        DO J=JTS,JTE    
         DO I=ITS,ITE  
           DO L=KTS,KTE
              KFLIP=KTE+1-L
             CWM_PHY(I,KFLIP,J)=CWM(TEMP_DEX)
             T_PHY(I,KFLIP,J)=T(TEMP_DEX)
             Q_PHY(I,KFLIP,J)=Q(TEMP_DEX)
             TLATGS_PHY(I,KFLIP,J)=TLATGS(TEMP_DEX)
             TRAIN_PHY(I,KFLIP,J)=TRAIN(TEMP_DEX)
             F_ice_PHY(I,KFLIP,J)=F_ice(L,I,J)
             F_rain_PHY(I,KFLIP,J)=F_rain(L,I,J)
             F_RimeF_PHY(I,KFLIP,J)=F_RimeF(L,I,J)
           ENDDO
        enddo
       enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      END SUBROUTINE EGCP01DRV
!
!###############################################################################
! ***** VERSION OF MICROPHYSICS DESIGNED FOR HIGHER RESOLUTION MESO ETA MODEL
!       (1) Represents sedimentation by preserving a portion of the precipitation
!           through top-down integration from cloud-top.  Modified procedure to
!           Zhao and Carr (1997).
!       (2) Microphysical equations are modified to be less sensitive to time
!           steps by use of Clausius-Clapeyron equation to account for changes in
!           saturation mixing ratios in response to latent heating/cooling.  
!       (3) Prevent spurious temperature oscillations across 0C due to 
!           microphysics.
!       (4) Uses lookup tables for: calculating two different ventilation 
!           coefficients in condensation and deposition processes; accretion of
!           cloud water by precipitation; precipitation mass; precipitation rate
!           (and mass-weighted precipitation fall speeds).
!       (5) Assumes temperature-dependent variation in mean diameter of large ice
!           (Houze et al., 1979; Ryan et al., 1996).
!        -> 8/22/01: This relationship has been extended to colder temperatures
!           to parameterize smaller large-ice particles down to mean sizes of MDImin,
!           which is 50 microns reached at -55.9C.
!       (6) Attempts to differentiate growth of large and small ice, mainly for
!           improved transition from thin cirrus to thick, precipitating ice
!           anvils.
!        -> 8/22/01: This feature has been diminished by effectively adjusting to
!           ice saturation during depositional growth at temperatures colder than
!           -10C.  Ice sublimation is calculated more explicitly.  The logic is
!           that sources of are either poorly understood (e.g., nucleation for NWP) 
!           or are not represented in the Eta model (e.g., detrainment of ice from 
!           convection).  Otherwise the model is too wet compared to the radiosonde
!           observations based on 1 Feb - 18 March 2001 retrospective runs.  
!       (7) Top-down integration also attempts to treat mixed-phase processes,
!           allowing a mixture of ice and water.  Based on numerous observational
!           studies, ice growth is based on nucleation at cloud top &
!           subsequent growth by vapor deposition and riming as the ice particles 
!           fall through the cloud.  Effective nucleation rates are a function
!           of ice supersaturation following Meyers et al. (JAM, 1992).  
!        -> 8/22/01: The simulated relative humidities were far too moist compared 
!           to the rawinsonde observations.  This feature has been substantially 
!           diminished, limited to a much narrower temperature range of 0 to -10C.  
!       (8) Depositional growth of newly nucleated ice is calculated for large time
!           steps using Fig. 8 of Miller and Young (JAS, 1979), at 1 deg intervals
!           using their ice crystal masses calculated after 600 s of growth in water
!           saturated conditions.  The growth rates are normalized by time step
!           assuming 3D growth with time**1.5 following eq. (6.3) in Young (1993).
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!       (9) Ice precipitation rates can increase due to increase in response to
!           cloud water riming due to (a) increased density & mass of the rimed
!           ice, and (b) increased fall speeds of rimed ice.
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!###############################################################################
!###############################################################################
!
      SUBROUTINE EGCP01COLUMN ( ARAIN, ASNOW, DTPH, I_index, J_index,   &
     & LSFC, P_col, QI_col, QR_col, QV_col, QW_col, RimeF_col, T_col,   &
     & THICK_col, WC_col ,KTS,KTE,NSTATS,QMAX,QTOT)                          
!
!###############################################################################
!###############################################################################
!
!-------------------------------------------------------------------------------
!----- NOTE:  Code is currently set up w/o threading!  
!-------------------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:  Grid-scale microphysical processes - condensation & precipitation
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: 08-2001
!   PRGRMMR: Jin  (Modification for WRF structure)
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Merges original GSCOND & PRECPD subroutines.   
!   * Code has been substantially streamlined and restructured.
!   * Exchange between water vapor & small cloud condensate is calculated using
!     the original Asai (1965, J. Japan) algorithm.  See also references to
!     Yau and Austin (1979, JAS), Rutledge and Hobbs (1983, JAS), and Tao et al.
!     (1989, MWR).  This algorithm replaces the Sundqvist et al. (1989, MWR)
!     parameterization.  
!-------------------------------------------------------------------------------
!     
! USAGE: 
!   * CALL EGCP01COLUMN FROM SUBROUTINE EGCP01DRV
!
! INPUT ARGUMENT LIST:
!   DTPH       - physics time step (s)
!   I_index    - I index
!   J_index    - J index
!   LSFC       - Eta level of level above surface, ground
!   P_col      - vertical column of model pressure (Pa)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!   THICK_col  - vertical column of model mass thickness (density*height increment)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   
!
! OUTPUT ARGUMENT LIST: 
!   ARAIN      - accumulated rainfall at the surface (kg)
!   ASNOW      - accumulated snowfall at the surface (kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!     
! OUTPUT FILES:
!     NONE
!     
! Subprograms & Functions called:
!   * Real Function CONDENSE  - cloud water condensation
!   * Real Function DEPOSIT   - ice deposition (not sublimation)
!
! UNIQUE: NONE
!  
! LIBRARY: NONE
!  
! COMMON BLOCKS:  
!     CMICRO_CONS  - key constants initialized in GSMCONST
!     CMICRO_STATS - accumulated and maximum statistics
!     CMY_GROWTH   - lookup table for growth of ice crystals in 
!                    water saturated conditions (Miller & Young, 1979)
!     IVENT_TABLES - lookup tables for ventilation effects of ice
!     IACCR_TABLES - lookup tables for accretion rates of ice
!     IMASS_TABLES - lookup tables for mass content of ice
!     IRATE_TABLES - lookup tables for precipitation rates of ice
!     IRIME_TABLES - lookup tables for increase in fall speed of rimed ice
!     RVENT_TABLES - lookup tables for ventilation effects of rain
!     RACCR_TABLES - lookup tables for accretion rates of rain
!     RMASS_TABLES - lookup tables for mass content of rain
!     RVELR_TABLES - lookup tables for fall speeds of rain
!     RRATE_TABLES - lookup tables for precipitation rates of rain
!   
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
!
!------------------------------------------------------------------------- 
!--------------- Arrays & constants in argument list --------------------- 
!------------------------------------------------------------------------- 
!
      IMPLICIT NONE
!    
      INTEGER,INTENT(IN) :: KTS,KTE,I_index, J_index, LSFC
      REAL,INTENT(INOUT) ::  ARAIN, ASNOW
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) ::  P_col, QI_col,QR_col    &
     & ,QV_col ,QW_col, RimeF_col, T_col, THICK_col,WC_col
!
!------------------------------------------------------------------------- 
!-------------- Common blocks for microphysical statistics ---------------
!------------------------------------------------------------------------- 
!
!------------------------------------------------------------------------- 
!--------- Common blocks for constants initialized in GSMCONST ----------
!
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40
      INTEGER,INTENT(INOUT) :: NSTATS(ITLO:ITHI,4)
      REAL,INTENT(INOUT) :: QMAX(ITLO:ITHI,5),QTOT(ITLO:ITHI,22) 
!
!------------------------------------------------------------------------- 
!--------------- Common blocks for various lookup tables -----------------
!
!--- Discretized growth rates of small ice crystals after their nucleation 
!     at 1 C intervals from -1 C to -35 C, based on calculations by Miller 
!     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
!     are multiplied by physics time step in GSMCONST.
!
!------------------------------------------------------------------------- 
!
!--- Mean ice-particle diameters varying from 50 microns to 1000 microns
!      (1 mm), assuming an exponential size distribution.  
!
!---- Meaning of the following arrays: 
!        - mdiam - mean diameter (m)
!        - VENTI1 - integrated quantity associated w/ ventilation effects 
!                   (capacitance only) for calculating vapor deposition onto ice
!        - VENTI2 - integrated quantity associated w/ ventilation effects 
!                   (with fall speed) for calculating vapor deposition onto ice
!        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
!        - MASSI  - integrated quantity associated w/ ice mass 
!        - VSNOWI - mass-weighted fall speed of snow (large ice), used to calculate 
!                   precipitation rates
!
!
!------------------------------------------------------------------------- 
!
!--- VEL_RF - velocity increase of rimed particles as functions of crude
!      particle size categories (at 0.1 mm intervals of mean ice particle
!      sizes) and rime factor (different values of Rime Factor of 1.1**N, 
!      where N=0 to Nrime).
!
!------------------------------------------------------------------------- 
!
!--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns 
!      (0.45 mm), assuming an exponential size distribution.  
!
!------------------------------------------------------------------------- 
!------- Key parameters, local variables, & important comments ---------
!-----------------------------------------------------------------------
!
!--- TOLER => Tolerance or precision for accumulated precipitation 
!
      REAL, PARAMETER :: TOLER=5.E-7, C2=1./6., RHO0=1.194, Xratio=.025                                           
!
!--- If BLEND=1:
!      precipitation (large) ice amounts are estimated at each level as a 
!      blend of ice falling from the grid point above and the precip ice
!      present at the start of the time step (see TOT_ICE below).
!--- If BLEND=0:
!      precipitation (large) ice amounts are estimated to be the precip
!      ice present at the start of the time step.
!
!--- Extended to include sedimentation of rain on 2/5/01 
!
      REAL, PARAMETER :: BLEND=1.
!
!--- This variable is for debugging purposes (if .true.)
!
      LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
!
!-----------------------------------------------------------------------
!--- Local variables
!-----------------------------------------------------------------------
!
      REAL EMAIRI, N0r, NLICE, NSmICE
      LOGICAL CLEAR, ICE_logical, DBG_logical, RAIN_logical
      INTEGER :: IDR,INDEX_MY,INDEXR,INDEXR1,INDEXS,IPASS,ITDX,IXRF,    &
     &           IXS,LBEF,L
!
      REAL :: ABI,ABW,AIEVP,ARAINnew,ASNOWnew,BLDTRH,BUDGET,            &
     &        CREVP,DELI,DELR,DELT,DELV,DELW,DENOMF,                    &
     &        DENOMI,DENOMW,DENOMWI,DIDEP,                              &
     &        DIEVP,DIFFUS,DLI,DTPH,DTRHO,DUM,DUM1,                     &
     &        DUM2,DWV0,DWVI,DWVR,DYNVIS,ESI,ESW,FIR,FLARGE,FLIMASS,    &
     &        FSMALL,FWR,FWS,GAMMAR,GAMMAS,                             &
     &        PCOND,PIACR,PIACW,PIACWI,PIACWR,PICND,PIDEP,PIDEP_max,    &
     &        PIEVP,PILOSS,PIMLT,PP,PRACW,PRAUT,PREVP,PRLOSS,           &
     &        QI,QInew,QLICE,QR,QRnew,QSI,QSIgrd,QSInew,QSW,QSW0,       &
     &        QSWgrd,QSWnew,QT,QTICE,QTnew,QTRAIN,QV,QW,QW0,QWnew,      &
     &        RFACTOR,RHO,RIMEF,RIMEF1,RQR,RR,RRHO,SFACTOR,             &
     &        TC,TCC,TFACTOR,THERM_COND,THICK,TK,TK2,TNEW,              &
     &        TOT_ICE,TOT_ICEnew,TOT_RAIN,TOT_RAINnew,                  &
     &        VEL_INC,VENTR,VENTIL,VENTIS,VRAIN1,VRAIN2,VRIMEF,VSNOW,   &
     &        WC,WCnew,WSgrd,WS,WSnew,WV,WVnew,WVQW,                    &
     &        XLF,XLF1,XLI,XLV,XLV1,XLV2,XLIMASS,XRF,XSIMASS          
!
!#######################################################################
!########################## Begin Execution ############################
!#######################################################################
!
!
      ARAIN=0.                ! Accumulated rainfall into grid box from above (kg/m**2)
      ASNOW=0.                ! Accumulated snowfall into grid box from above (kg/m**2)
!
!-----------------------------------------------------------------------
!------------ Loop from top (L=1) to surface (L=LSFC) ------------------
!-----------------------------------------------------------------------
!

      DO 10 L=1,LSFC

!--- Skip this level and go to the next lower level if no condensate 
!      and very low specific humidities
!
        IF (QV_col(L).LE.EPSQ .AND. WC_col(L).LE.EPSQ) GO TO 10
!
!-----------------------------------------------------------------------
!------------ Proceed with cloud microphysics calculations -------------
!-----------------------------------------------------------------------
!
          TK=T_col(L)         ! Temperature (deg K)
          TC=TK-T0C           ! Temperature (deg C)
          PP=P_col(L)         ! Pressure (Pa)
          QV=QV_col(L)        ! Specific humidity of water vapor (kg/kg)
          WV=QV/(1.-QV)       ! Water vapor mixing ratio (kg/kg)
          WC=WC_col(L)        ! Grid-scale mixing ratio of total condensate (water or ice; kg/kg)
!
!-----------------------------------------------------------------------
!--- Moisture variables below are mixing ratios & not specifc humidities
!-----------------------------------------------------------------------
!
          CLEAR=.TRUE.
!    
!--- This check is to determine grid-scale saturation when no condensate is present
!    
          ESW=1000.*FPVS0(TK)              ! Saturation vapor pressure w/r/t water
          QSW=EPS*ESW/(PP-ESW)             ! Saturation mixing ratio w/r/t water
          WS=QSW                           ! General saturation mixing ratio (water/ice)
          IF (TC .LT. 0.) THEN
            ESI=1000.*FPVS(TK)             ! Saturation vapor pressure w/r/t ice
            QSI=EPS*ESI/(PP-ESI)           ! Saturation mixing ratio w/r/t water
            WS=QSI                         ! General saturation mixing ratio (water/ice)
          ENDIF
!
!--- Effective grid-scale Saturation mixing ratios
!
          QSWgrd=RHgrd*QSW
          QSIgrd=RHgrd*QSI
          WSgrd=RHgrd*WS
!
!--- Check if air is subsaturated and w/o condensate
!
          IF (WV.GT.WSgrd .OR. WC.GT.EPSQ) CLEAR=.FALSE.
!
!--- Check if any rain is falling into layer from above
!
          IF (ARAIN .GT. CLIMIT) THEN
            CLEAR=.FALSE.
          ELSE
            ARAIN=0.
          ENDIF
!
!--- Check if any ice is falling into layer from above
!
!--- NOTE that "SNOW" in variable names is synonomous with 
!    large, precipitation ice particles
!
          IF (ASNOW .GT. CLIMIT) THEN
            CLEAR=.FALSE.
          ELSE
            ASNOW=0.
          ENDIF
!
!-----------------------------------------------------------------------
!-- Loop to the end if in clear, subsaturated air free of condensate ---
!-----------------------------------------------------------------------
!
          IF (CLEAR) GO TO 10
!
!-----------------------------------------------------------------------
!--------- Initialize RHO, THICK & microphysical processes -------------
!-----------------------------------------------------------------------
!
!
!--- Virtual temperature, TV=T*(1./EPS-1)*Q, Q is specific humidity;
!    (see pp. 63-65 in Fleagle & Businger, 1963)
!
          RHO=PP/(RD*TK*(1.+EPS1*QV))   ! Air density (kg/m**3)
          RRHO=1./RHO                ! Reciprocal of air density
          DTRHO=DTPH*RHO             ! Time step * air density
          BLDTRH=BLEND*DTRHO         ! Blend parameter * time step * air density
          THICK=THICK_col(L)         ! Layer thickness = RHO*DZ = -DP/G = (Psfc-Ptop)*D_ETA/(G*ETA_sfc)
!
          ARAINnew=0.                ! Updated accumulated rainfall
          ASNOWnew=0.                ! Updated accumulated snowfall
          QI=QI_col(L)               ! Ice mixing ratio
          QInew=0.                   ! Updated ice mixing ratio
          QR=QR_col(L)               ! Rain mixing ratio
          QRnew=0.                   ! Updated rain ratio
          QW=QW_col(L)               ! Cloud water mixing ratio
          QWnew=0.                   ! Updated cloud water ratio
!
          PCOND=0.                   ! Condensation (>0) or evaporation (<0) of cloud water (kg/kg)
          PIDEP=0.                   ! Deposition (>0) or sublimation (<0) of ice crystals (kg/kg)
          PIACW=0.                   ! Cloud water collection (riming) by precipitation ice (kg/kg; >0)
          PIACWI=0.                  ! Growth of precip ice by riming (kg/kg; >0)
          PIACWR=0.                  ! Shedding of accreted cloud water to form rain (kg/kg; >0)
          PIACR=0.                   ! Freezing of rain onto large ice at supercooled temps (kg/kg; >0)
          PICND=0.                   ! Condensation (>0) onto wet, melting ice (kg/kg)
          PIEVP=0.                   ! Evaporation (<0) from wet, melting ice (kg/kg)
          PIMLT=0.                   ! Melting ice (kg/kg; >0)
          PRAUT=0.                   ! Cloud water autoconversion to rain (kg/kg; >0)
          PRACW=0.                   ! Cloud water collection (accretion) by rain (kg/kg; >0)
          PREVP=0.                   ! Rain evaporation (kg/kg; <0)
!
!--- Double check input hydrometeor mixing ratios
!
!          DUM=WC-(QI+QW+QR)
!          DUM1=ABS(DUM)
!          DUM2=TOLER*MIN(WC, QI+QW+QR)
!          IF (DUM1 .GT. DUM2) THEN
!            WRITE(6,"(/2(a,i4),a,i2)") '{@ i=',I_index,' j=',J_index,
!     &                                 ' L=',L
!            WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{@ TCold=',TC,'P=',.01*PP,'DIFF=',DUM,'WCold=',WC,
!     & '{@ QIold=',QI,'QWold=',QW,'QRold=',QR
!          ENDIF
!
!***********************************************************************
!*********** MAIN MICROPHYSICS CALCULATIONS NOW FOLLOW! ****************
!***********************************************************************
!
!--- Calculate a few variables, which are used more than once below
!
!--- Latent heat of vaporization as a function of temperature from
!      Bolton (1980, JAS)
!
          XLV=3.148E6-2370*TK        ! Latent heat of vaporization (Lv)
          XLF=XLS-XLV                ! Latent heat of fusion (Lf)
          XLV1=XLV*RCP               ! Lv/Cp
          XLF1=XLF*RCP               ! Lf/Cp
          TK2=1./(TK*TK)             ! 1./TK**2
          XLV2=XLV*XLV*QSW*TK2/RV    ! Lv**2*Qsw/(Rv*TK**2)
          DENOMW=1.+XLV2*RCP         ! Denominator term, Clausius-Clapeyron correction
!
!--- Basic thermodynamic quantities
!      * DYNVIS - dynamic viscosity  [ kg/(m*s) ]
!      * THERM_COND - thermal conductivity  [ J/(m*s*K) ]
!      * DIFFUS - diffusivity of water vapor  [ m**2/s ]
!
          TFACTOR=TK**1.5/(TK+120.)
          DYNVIS=1.496E-6*TFACTOR
          THERM_COND=2.116E-3*TFACTOR
          DIFFUS=8.794E-5*TK**1.81/PP
!
!--- Air resistance term for the fall speed of ice following the
!      basic research by Heymsfield, Kajikawa, others 
!
          GAMMAS=(1.E5/PP)**C1
!
!--- Air resistance for rain fall speed (Beard, 1985, JAS, p.470)
!
          GAMMAR=(RHO0/RHO)**.4
!
!----------------------------------------------------------------------
!-------------  IMPORTANT MICROPHYSICS DECISION TREE  -----------------
!----------------------------------------------------------------------
!
!--- Determine if conditions supporting ice are present
!
          IF (TC.LT.0. .OR. QI.GT. EPSQ .OR. ASNOW.GT.CLIMIT) THEN
            ICE_logical=.TRUE.
          ELSE
            ICE_logical=.FALSE.
            QLICE=0.
            QTICE=0.
          ENDIF
!
!--- Determine if rain is present
!
          RAIN_logical=.FALSE.
          IF (ARAIN.GT.CLIMIT .OR. QR.GT.EPSQ) RAIN_logical=.TRUE.
!
          IF (ICE_logical) THEN
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!---
!  * FLARGE  - ratio of number of large ice to total (large & small) ice
!  * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!---
!  * TOT_ICE - total mass (small & large) ice before microphysics,
!              which is the sum of the total mass of large ice in the 
!              current layer and the input flux of ice from above
!  * PILOSS  - greatest loss (<0) of total (small & large) ice by
!              sublimation, removing all of the ice falling from above
!              and the ice within the layer
!  * RimeF1  - Rime Factor, which is the mass ratio of total (unrimed & rimed) 
!              ice mass to the unrimed ice mass (>=1)
!  * VrimeF  - the velocity increase due to rime factor or melting (ratio, >=1)
!  * VSNOW   - Fall speed of rimed snow w/ air resistance correction
!  * EMAIRI  - equivalent mass of air associated layer and with fall of snow into layer
!  * XLIMASS - used for calculating large ice mixing ratio
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE   - time-averaged number concentration of large ice
!  * NSmICE  - number concentration of small ice crystals at current level
!---
!--- Assumed number fraction of large ice particles to total (large & small) 
!    ice particles, which is based on a general impression of the literature.
!
            WVQW=WV+QW                ! Water vapor & cloud water
!


            IF (TC.GE.0. .OR. WVQW.LT.QSIgrd) THEN
   !
   !--- Eliminate small ice particle contributions for melting & sublimation
   !
              FLARGE=FLARGE1
            ELSE
   !
   !--- Enhanced number of small ice particles during depositional growth
   !    (effective only when 0C > T >= T_ice [-10C] )
   !
              FLARGE=FLARGE2
   !
   !--- Larger number of small ice particles due to rime splintering
   !
              IF (TC.GE.-8. .AND. TC.LE.-3.) FLARGE=.5*FLARGE
!
            ENDIF            ! End IF (TC.GE.0. .OR. WVQW.LT.QSIgrd)
            FSMALL=(1.-FLARGE)/FLARGE
            XSIMASS=RRHO*MASSI(MDImin)*FSMALL
            IF (QI.LE.EPSQ .AND. ASNOW.LE.CLIMIT) THEN
              INDEXS=MDImin
              TOT_ICE=0.
              PILOSS=0.
              RimeF1=1.
              VrimeF=1.
              VEL_INC=GAMMAS
              VSNOW=0.
              EMAIRI=THICK
              XLIMASS=RRHO*RimeF1*MASSI(INDEXS)
              FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
              QLICE=0.
              QTICE=0.
              NLICE=0.
              NSmICE=0.
            ELSE
   !
   !--- For T<0C mean particle size follows Houze et al. (JAS, 1979, p. 160), 
   !    converted from Fig. 5 plot of LAMDAs.  Similar set of relationships 
   !    also shown in Fig. 8 of Ryan (BAMS, 1996, p. 66).
   !
              DUM=XMImax*EXP(.0536*TC)
              INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
              TOT_ICE=THICK*QI+BLEND*ASNOW
              PILOSS=-TOT_ICE/THICK
              LBEF=MAX(1,L-1)
              DUM1=RimeF_col(LBEF)
              DUM2=RimeF_col(L)
              RimeF1=(DUM2*THICK*QI+DUM1*BLEND*ASNOW)/TOT_ICE
              RimeF1=MIN(RimeF1, RFmax)
              DO IPASS=0,1
                IF (RimeF1 .LE. 1.) THEN
                  RimeF1=1.
                  VrimeF=1.
                ELSE
                  IXS=MAX(2, MIN(INDEXS/100, 9))
                  XRF=10.492*ALOG(RimeF1)
                  IXRF=MAX(0, MIN(INT(XRF), Nrime))
                  IF (IXRF .GE. Nrime) THEN
                    VrimeF=VEL_RF(IXS,Nrime)
                  ELSE
                    VrimeF=VEL_RF(IXS,IXRF)+(XRF-FLOAT(IXRF))*          &
     &                    (VEL_RF(IXS,IXRF+1)-VEL_RF(IXS,IXRF))
                  ENDIF
                ENDIF            ! End IF (RimeF1 .LE. 1.)
                VEL_INC=GAMMAS*VrimeF
                VSNOW=VEL_INC*VSNOWI(INDEXS)
                EMAIRI=THICK+BLDTRH*VSNOW
                XLIMASS=RRHO*RimeF1*MASSI(INDEXS)
                FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
                QTICE=TOT_ICE/EMAIRI
                QLICE=FLIMASS*QTICE
                NLICE=QLICE/XLIMASS
                NSmICE=Fsmall*NLICE
   !
                IF ( (NLICE.GE.NLImin .AND. NLICE.LE.NLImax)            &
     &                .OR. IPASS.EQ.1) THEN
                  EXIT
                ELSE
                  IF (TC < 0) THEN
                    XLI=RHO*(QTICE/DUM-XSIMASS)/RimeF1
                    IF (XLI .LE. MASSI(MDImin) ) THEN
                      INDEXS=MDImin
                    ELSE IF (XLI .LE. MASSI(450) ) THEN
                      DLI=9.5885E5*XLI**.42066         ! DLI in microns
                      INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                    ELSE IF (XLI .LE. MASSI(MDImax) ) THEN
                      DLI=3.9751E6*XLI**.49870         ! DLI in microns
                      INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                    ELSE
                      INDEXS=MDImax
                    ENDIF             ! End IF (XLI .LE. MASSI(MDImin) )
                  ENDIF               ! End IF (TC < 0)
        !
        !--- Reduce excessive accumulation of ice at upper levels
        !    associated with strong grid-resolved ascent
        !
        !--- Force NLICE to be between NLImin and NLImax
        !
        !
        !--- 8/22/01: Increase density of large ice if maximum limits 
        !    are reached for number concentration (NLImax) and mean size 
        !    (MDImax).  Done to increase fall out of ice.
        !
                  DUM=MAX(NLImin, MIN(NLImax, NLICE) )
                  IF (DUM.GE.NLImax .AND. INDEXS.GE.MDImax)             &
     &                RimeF1=RHO*(QTICE/NLImax-XSIMASS)/MASSI(INDEXS)
!            WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{$ TC=',TC,'P=',.01*PP,'NLICE=',NLICE,'DUM=',DUM,
!     & '{$ XLI=',XLI,'INDEXS=',FLOAT(INDEXS),'RHO=',RHO,'QTICE=',QTICE,
!     & '{$ XSIMASS=',XSIMASS,'RimeF1=',RimeF1
                ENDIF                  ! End IF ( (NLICE.GE.NLImin .AND. NLICE.LE.NLImax) ...
              ENDDO                    ! End DO IPASS=0,1
            ENDIF                      ! End IF (QI.LE.EPSQ .AND. ASNOW.LE.CLIMIT)
          ENDIF                        ! End IF (ICE_logical)
!
!----------------------------------------------------------------------
!--------------- Calculate individual processes -----------------------
!----------------------------------------------------------------------
!
!--- Cloud water autoconversion to rain and collection by rain
!
          IF (QW.GT.EPSQ .AND. TC.GE.T_ICE) THEN
   !
   !--- QW0 could be modified based on land/sea properties, 
   !      presence of convection, etc.  This is why QAUT0 and CRAUT
   !      are passed into the subroutine as externally determined
   !      parameters.  Can be changed in the future if desired.
   !
            QW0=QAUT0*RRHO
            PRAUT=MAX(0., QW-QW0)*CRAUT
            IF (QLICE .GT. EPSQ) THEN
      !
      !--- Collection of cloud water by large ice particles ("snow")
      !    PIACWI=PIACW for riming, PIACWI=0 for shedding
      !
              FWS=MIN(1., CIACW*VEL_INC*NLICE*ACCRI(INDEXS)/PP**C1)
              PIACW=FWS*QW
              IF (TC .LT. 0.) PIACWI=PIACW    ! Large ice riming
            ENDIF           ! End IF (QLICE .GT. EPSQ)
          ENDIF             ! End IF (QW.GT.EPSQ .AND. TC.GE.T_ICE)
!
!----------------------------------------------------------------------
!--- Loop around some of the ice-phase processes if no ice should be present
!----------------------------------------------------------------------
!
          IF (ICE_logical .EQV. .FALSE.) GO TO 20
!
!--- Now the pretzel logic of calculating ice deposition
!
          IF (TC.LT.T_ICE .AND. (WV.GT.QSIgrd .OR. QW.GT.EPSQ)) THEN
   !
   !--- Adjust to ice saturation at T<T_ICE (-10C) if supersaturated.
   !    Sources of ice due to nucleation and convective detrainment are
   !    either poorly understood, poorly resolved at typical NWP 
   !    resolutions, or are not represented (e.g., no detrained 
   !    condensate in BMJ Cu scheme).
   !    
            PCOND=-QW
            DUM1=TK+XLV1*PCOND                 ! Updated (dummy) temperature (deg K)
            DUM2=WV+QW                         ! Updated (dummy) water vapor mixing ratio
            DUM=1000.*FPVS(DUM1)               ! Updated (dummy) saturation vapor pressure w/r/t ice
            DUM=RHgrd*EPS*DUM/(PP-DUM)         ! Updated (dummy) saturation mixing ratio w/r/t ice
            IF (DUM2 .GT. DUM) PIDEP=DEPOSIT (PP, DUM1, DUM2)
            DWVi=0.    ! Used only for debugging
   !
          ELSE IF (TC .LT. 0.) THEN
   !
   !--- These quantities are handy for ice deposition/sublimation
   !    PIDEP_max - max deposition or minimum sublimation to ice saturation
   !
            DENOMI=1.+XLS2*QSI*TK2
            DWVi=MIN(WVQW,QSW)-QSI
            PIDEP_max=MAX(PILOSS, DWVi/DENOMI)
            IF (QTICE .GT. 0.) THEN
      !
      !--- Calculate ice deposition/sublimation
      !      * SFACTOR - [VEL_INC**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
      !        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
      !      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
      !               VENTIL, VENTIS - m**-2 ;  VENTI1 - m ;  
      !               VENTI2 - m**2/s**.5 ; DIDEP - unitless
      !
              SFACTOR=VEL_INC**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
              ABI=1./(RHO*XLS3*QSI*TK2/THERM_COND+1./DIFFUS)
      !
      !--- VENTIL - Number concentration * ventilation factors for large ice
      !--- VENTIS - Number concentration * ventilation factors for small ice
      !
      !--- Variation in the number concentration of ice with time is not
      !      accounted for in these calculations (could be in the future).
      !
              VENTIL=(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))*NLICE
              VENTIS=(VENTI1(MDImin)+SFACTOR*VENTI2(MDImin))*NSmICE
              DIDEP=ABI*(VENTIL+VENTIS)*DTPH
      !
      !--- Account for change in water vapor supply w/ time
      !
              IF (DIDEP .GE. Xratio)then
                DIDEP=(1.-EXP(-DIDEP*DENOMI))/DENOMI
              endif
              IF (DWVi .GT. 0.) THEN
                PIDEP=MIN(DWVi*DIDEP, PIDEP_max)
              ELSE IF (DWVi .LT. 0.) THEN
                PIDEP=MAX(DWVi*DIDEP, PIDEP_max)
              ENDIF
      !
            ELSE IF (WVQW.GT.QSI .AND. TC.LE.T_ICE_init) THEN
      !
      !--- Ice nucleation in near water-saturated conditions.  Ice crystal
      !    growth during time step calculated using Miller & Young (1979, JAS).
      !--- These deposition rates could drive conditions below water saturation,
      !    which is the basis of these calculations.  Intended to approximate
      !    more complex & computationally intensive calculations.
      !
              INDEX_MY=MAX(MY_T1, MIN( INT(.5-TC), MY_T2 ) )
      !
      !--- DUM1 is the supersaturation w/r/t ice at water-saturated conditions
      !
      !--- DUM2 is the number of ice crystals nucleated at water-saturated 
      !    conditions based on Meyers et al. (JAM, 1992).
      !
      !--- Prevent unrealistically large ice initiation (limited by PIDEP_max)
      !      if DUM2 values are increased in future experiments
      !
              DUM1=QSW/QSI-1.      
              DUM2=1.E3*EXP(12.96*DUM1-.639)
              PIDEP=MIN(PIDEP_max, DUM2*MY_GROWTH_NMM(INDEX_MY)*RRHO)
      !
            ENDIF       ! End IF (QTICE .GT. 0.)
   !
          ENDIF         ! End IF (TC.LT.T_ICE .AND. (WV.GT.QSIgrd .OR. QW.GT.EPSQ))
!
!------------------------------------------------------------------------
!
20      CONTINUE     ! Jump here if conditions for ice are not present


!
!------------------------------------------------------------------------
!
!--- Cloud water condensation
!
          IF (TC.GE.T_ICE .AND. (QW.GT.EPSQ .OR. WV.GT.QSWgrd)) THEN
            IF (PIACWI.EQ.0. .AND. PIDEP.EQ.0.) THEN
              PCOND=CONDENSE (PP, QW, TK, WV)
            ELSE
   !
   !--- Modify cloud condensation in response to ice processes
   !
              DUM=XLV*QSWgrd*RCPRV*TK2
              DENOMWI=1.+XLS*DUM
              DENOMF=XLF*DUM
              DUM=MAX(0., PIDEP)
              PCOND=(WV-QSWgrd-DENOMWI*DUM-DENOMF*PIACWI)/DENOMW
              DUM1=-QW
              DUM2=PCOND-PIACW
              IF (DUM2 .LT. DUM1) THEN
      !
      !--- Limit cloud water sinks
      !
                DUM=DUM1/DUM2
                PCOND=DUM*PCOND
                PIACW=DUM*PIACW
                PIACWI=DUM*PIACWI
              ENDIF        ! End IF (DUM2 .LT. DUM1)
            ENDIF          ! End IF (PIACWI.EQ.0. .AND. PIDEP.EQ.0.)
          ENDIF            ! End IF (TC.GE.T_ICE .AND. (QW.GT.EPSQ .OR. WV.GT.QSWgrd))
!
!--- Limit freezing of accreted rime to prevent temperature oscillations,
!    a crude Schumann-Ludlam limit (p. 209 of Young, 1993). 
!
          TCC=TC+XLV1*PCOND+XLS1*PIDEP+XLF1*PIACWI
          IF (TCC .GT. 0.) THEN
            PIACWI=0.
            TCC=TC+XLV1*PCOND+XLS1*PIDEP
          ENDIF
          IF (TC.GT.0. .AND. TCC.GT.0. .AND. ICE_logical) THEN
   !
   !--- Calculate melting and evaporation/condensation
   !      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
   !               VENTIL - m**-2 ;  VENTI1 - m ;  
   !               VENTI2 - m**2/s**.5 ; CIEVP - /s
   !
            SFACTOR=VEL_INC**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
            VENTIL=NLICE*(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))
            AIEVP=VENTIL*DIFFUS*DTPH
            IF (AIEVP .LT. Xratio) THEN
              DIEVP=AIEVP
            ELSE
              DIEVP=1.-EXP(-AIEVP)
            ENDIF
            QSW0=EPS*ESW0/(PP-ESW0)
            DWV0=MIN(WV,QSW)-QSW0
            DUM=QW+PCOND
            IF (WV.LT.QSW .AND. DUM.LE.EPSQ) THEN
   !
   !--- Evaporation from melting snow (sink of snow) or shedding
   !    of water condensed onto melting snow (source of rain)
   !
              DUM=DWV0*DIEVP
              PIEVP=MAX( MIN(0., DUM), PILOSS)
              PICND=MAX(0., DUM)
            ENDIF            ! End IF (WV.LT.QSW .AND. DUM.LE.EPSQ)
            PIMLT=THERM_COND*TCC*VENTIL*RRHO*DTPH/XLF
   !
   !--- Limit melting to prevent temperature oscillations across 0C
   !
            DUM1=MAX( 0., (TCC+XLV1*PIEVP)/XLF1 )
            PIMLT=MIN(PIMLT, DUM1)
   !
   !--- Limit loss of snow by melting (>0) and evaporation
   !
            DUM=PIEVP-PIMLT
            IF (DUM .LT. PILOSS) THEN
              DUM1=PILOSS/DUM
              PIMLT=PIMLT*DUM1
              PIEVP=PIEVP*DUM1
            ENDIF           ! End IF (DUM .GT. QTICE)
          ENDIF             ! End IF (TC.GT.0. .AND. TCC.GT.0. .AND. ICE_logical) 
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!  * TOT_RAIN - total mass of rain before microphysics, which is the sum of
!               the total mass of rain in the current layer and the input 
!               flux of rain from above
!  * VRAIN1   - fall speed of rain into grid from above (with air resistance correction)
!  * QTRAIN   - time-averaged mixing ratio of rain (kg/kg)
!  * PRLOSS   - greatest loss (<0) of rain, removing all rain falling from
!               above and the rain within the layer
!  * RQR      - rain content (kg/m**3)
!  * INDEXR   - mean size of rain drops to the nearest 1 micron in size
!  * N0r      - intercept of rain size distribution (typically 10**6 m**-4)
!
          TOT_RAIN=0.
          VRAIN1=0.
          QTRAIN=0.
          PRLOSS=0.
          RQR=0.
          N0r=0.
          INDEXR=MDRmin
          INDEXR1=INDEXR    !-- For debugging only
          IF (RAIN_logical) THEN
            IF (ARAIN .LE. 0.) THEN
              INDEXR=MDRmin
              VRAIN1=0.
            ELSE
   !
   !--- INDEXR (related to mean diameter) & N0r could be modified 
   !      by land/sea properties, presence of convection, etc.
   !
   !--- Rain rate normalized to a density of 1.194 kg/m**3
   !
              RR=ARAIN/(DTPH*GAMMAR)
   !
              IF (RR .LE. RR_DRmin) THEN
        !
        !--- Assume fixed mean diameter of rain (0.2 mm) for low rain rates, 
        !      instead vary N0r with rain rate
        !
                INDEXR=MDRmin
              ELSE IF (RR .LE. RR_DR1) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.05 and 0.10 mm:
        !      V(Dr)=5.6023e4*Dr**1.136, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*5.6023e4*Dr**(4+1.136) = 1.408e15*Dr**5.136,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.123e-3*RR**.1947 -> Dr (microns) = 1.123e3*RR**.1947
        !
                INDEXR=INT( 1.123E3*RR**.1947 + .5 )
                INDEXR=MAX( MDRmin, MIN(INDEXR, MDR1) )
              ELSE IF (RR .LE. RR_DR2) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.10 and 0.20 mm:
        !      V(Dr)=1.0867e4*Dr**.958, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*1.0867e4*Dr**(4+.958) = 2.731e14*Dr**4.958,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.225e-3*RR**.2017 -> Dr (microns) = 1.225e3*RR**.2017
        !
                INDEXR=INT( 1.225E3*RR**.2017 + .5 )
                INDEXR=MAX( MDR1, MIN(INDEXR, MDR2) )
              ELSE IF (RR .LE. RR_DR3) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.20 and 0.32 mm:
        !      V(Dr)=2831.*Dr**.80, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*2831.*Dr**(4+.80) = 7.115e13*Dr**4.80, 
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.3006e-3*RR**.2083 -> Dr (microns) = 1.3006e3*RR**.2083
        !
                INDEXR=INT( 1.3006E3*RR**.2083 + .5 )
                INDEXR=MAX( MDR2, MIN(INDEXR, MDR3) )
              ELSE IF (RR .LE. RR_DRmax) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.32 and 0.45 mm:
        !      V(Dr)=944.8*Dr**.6636, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*944.8*Dr**(4+.6636) = 2.3745e13*Dr**4.6636,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.355e-3*RR**.2144 -> Dr (microns) = 1.355e3*RR**.2144
        !
                INDEXR=INT( 1.355E3*RR**.2144 + .5 )
                INDEXR=MAX( MDR3, MIN(INDEXR, MDRmax) )
              ELSE 
        !
        !--- Assume fixed mean diameter of rain (0.45 mm) for high rain rates, 
        !      instead vary N0r with rain rate
        !
                INDEXR=MDRmax
              ENDIF              ! End IF (RR .LE. RR_DRmin) etc. 
              VRAIN1=GAMMAR*VRAIN(INDEXR)
            ENDIF              ! End IF (ARAIN .LE. 0.)
            INDEXR1=INDEXR     ! For debugging only
            TOT_RAIN=THICK*QR+BLEND*ARAIN
            QTRAIN=TOT_RAIN/(THICK+BLDTRH*VRAIN1)
            PRLOSS=-TOT_RAIN/THICK
            RQR=RHO*QTRAIN
   !
   !--- RQR - time-averaged rain content (kg/m**3)
   !
            IF (RQR .LE. RQR_DRmin) THEN
              N0r=MAX(N0rmin, CN0r_DMRmin*RQR)
              INDEXR=MDRmin
            ELSE IF (RQR .GE. RQR_DRmax) THEN
              N0r=CN0r_DMRmax*RQR
              INDEXR=MDRmax
            ELSE
              N0r=N0r0
              INDEXR=MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
            ENDIF
   !
            IF (TC .LT. T_ICE) THEN
              PIACR=-PRLOSS
            ELSE
              DWVr=WV-PCOND-QSW
              DUM=QW+PCOND
              IF (DWVr.LT.0. .AND. DUM.LE.EPSQ) THEN
      !
      !--- Rain evaporation
      !
      !    * RFACTOR - [GAMMAR**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
      !        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
      !
      !    * Units: RFACTOR - s**.5/m ;  ABW - m**2/s ;  VENTR - m**-2 ;  
      !             N0r - m**-4 ;  VENTR1 - m**2 ;  VENTR2 - m**3/s**.5 ;
      !             CREVP - unitless
      !
                RFACTOR=GAMMAR**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                ABW=1./(RHO*XLV2/THERM_COND+1./DIFFUS)
      !
      !--- Note that VENTR1, VENTR2 lookup tables do not include the 
      !      1/Davg multiplier as in the ice tables
      !
                VENTR=N0r*(VENTR1(INDEXR)+RFACTOR*VENTR2(INDEXR))
                CREVP=ABW*VENTR*DTPH
                IF (CREVP .LT. Xratio) THEN
                  DUM=DWVr*CREVP
                ELSE
                  DUM=DWVr*(1.-EXP(-CREVP*DENOMW))/DENOMW
                ENDIF
                PREVP=MAX(DUM, PRLOSS)
              ELSE IF (QW .GT. EPSQ) THEN
                FWR=CRACW*GAMMAR*N0r*ACCRR(INDEXR)
                PRACW=MIN(1.,FWR)*QW
              ENDIF           ! End IF (DWVr.LT.0. .AND. DUM.LE.EPSQ)
      !
              IF (TC.LT.0. .AND. TCC.LT.0.) THEN
         !
         !--- Biggs (1953) heteorogeneous freezing (e.g., Lin et al., 1983)
         !   - Rescaled mean drop diameter from microns (INDEXR) to mm (DUM) to prevent underflow
         !
                DUM=.001*FLOAT(INDEXR)
                DUM=(EXP(ABFR*TC)-1.)*DUM*DUM*DUM*DUM*DUM*DUM*DUM
                PIACR=MIN(CBFR*N0r*RRHO*DUM, QTRAIN)
                IF (QLICE .GT. EPSQ) THEN
            !
            !--- Freezing of rain by collisions w/ large ice
            !
                  DUM=GAMMAR*VRAIN(INDEXR)
                  DUM1=DUM-VSNOW
            !
            !--- DUM2 - Difference in spectral fall speeds of rain and
            !      large ice, parameterized following eq. (48) on p. 112 of 
            !      Murakami (J. Meteor. Soc. Japan, 1990)
            !
                  DUM2=(DUM1*DUM1+.04*DUM*VSNOW)**.5
                  DUM1=5.E-12*INDEXR*INDEXR+2.E-12*INDEXR*INDEXS        &
     &                 +.5E-12*INDEXS*INDEXS
                  FIR=MIN(1., CIACR*NLICE*DUM1*DUM2)
            !
            !--- Future?  Should COLLECTION BY SMALL ICE SHOULD BE INCLUDED???
            !
                  PIACR=MIN(PIACR+FIR*QTRAIN, QTRAIN)
                ENDIF        ! End IF (QLICE .GT. EPSQ)
                DUM=PREVP-PIACR
                If (DUM .LT. PRLOSS) THEN
                  DUM1=PRLOSS/DUM
                  PREVP=DUM1*PREVP
                  PIACR=DUM1*PIACR
                ENDIF        ! End If (DUM .LT. PRLOSS)
              ENDIF          ! End IF (TC.LT.0. .AND. TCC.LT.0.)
            ENDIF            ! End IF (TC .LT. T_ICE)
          ENDIF              ! End IF (RAIN_logical) 
!
!----------------------------------------------------------------------
!---------------------- Main Budget Equations -------------------------
!----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!--- Update fields, determine characteristics for next lower layer ----
!-----------------------------------------------------------------------
!
!--- Carefully limit sinks of cloud water
!
          DUM1=PIACW+PRAUT+PRACW-MIN(0.,PCOND)
          IF (DUM1 .GT. QW) THEN
            DUM=QW/DUM1
            PIACW=DUM*PIACW
            PIACWI=DUM*PIACWI
            PRAUT=DUM*PRAUT
            PRACW=DUM*PRACW
            IF (PCOND .LT. 0.) PCOND=DUM*PCOND
          ENDIF
          PIACWR=PIACW-PIACWI          ! TC >= 0C
!
!--- QWnew - updated cloud water mixing ratio
!
          DELW=PCOND-PIACW-PRAUT-PRACW
          QWnew=QW+DELW
          IF (QWnew .LE. EPSQ) QWnew=0.
          IF (QW.GT.0. .AND. QWnew.NE.0.) THEN
            DUM=QWnew/QW
            IF (DUM .LT. TOLER) QWnew=0.
          ENDIF
!
!--- Update temperature and water vapor mixing ratios
!
          DELT= XLV1*(PCOND+PIEVP+PICND+PREVP)                          &
     &         +XLS1*PIDEP+XLF1*(PIACWI+PIACR-PIMLT)
          Tnew=TK+DELT
!
          DELV=-PCOND-PIDEP-PIEVP-PICND-PREVP
          WVnew=WV+DELV
!
!--- Update ice mixing ratios
!
!---
!  * TOT_ICEnew - total mass (small & large) ice after microphysics,
!                 which is the sum of the total mass of large ice in the 
!                 current layer and the flux of ice out of the grid box below
!  * RimeF      - Rime Factor, which is the mass ratio of total (unrimed & 
!                 rimed) ice mass to the unrimed ice mass (>=1)
!  * QInew      - updated mixing ratio of total (large & small) ice in layer
!      -> TOT_ICEnew=QInew*THICK+BLDTRH*QLICEnew*VSNOW
!        -> But QLICEnew=QInew*FLIMASS, so
!      -> TOT_ICEnew=QInew*(THICK+BLDTRH*FLIMASS*VSNOW)
!  * ASNOWnew   - updated accumulation of snow at bottom of grid cell
!---
!
          DELI=0.
          RimeF=1.
          IF (ICE_logical) THEN
            DELI=PIDEP+PIEVP+PIACWI+PIACR-PIMLT
            TOT_ICEnew=TOT_ICE+THICK*DELI
            IF (TOT_ICE.GT.0. .AND. TOT_ICEnew.NE.0.) THEN
              DUM=TOT_ICEnew/TOT_ICE
              IF (DUM .LT. TOLER) TOT_ICEnew=0.
            ENDIF
            IF (TOT_ICEnew .LE. CLIMIT) THEN
              TOT_ICEnew=0.
              RimeF=1.
              QInew=0.
              ASNOWnew=0.
            ELSE
      !
      !--- Update rime factor if appropriate
      !
              DUM=PIACWI+PIACR
              IF (DUM.LE.EPSQ .AND. PIDEP.LE.EPSQ) THEN
                RimeF=RimeF1
              ELSE
         !
         !--- Rime Factor, RimeF = (Total ice mass)/(Total unrimed ice mass)
         !      DUM1 - Total ice mass, rimed & unrimed
         !      DUM2 - Estimated mass of *unrimed* ice
         !
                DUM1=TOT_ICE+THICK*(PIDEP+DUM)
                DUM2=TOT_ICE/RimeF1+THICK*PIDEP
                IF (DUM2 .LE. 0.) THEN
                  RimeF=RFmax
                ELSE
                  RimeF=MIN(RFmax, MAX(1., DUM1/DUM2) )
                ENDIF
              ENDIF       ! End IF (DUM.LE.EPSQ .AND. PIDEP.LE.EPSQ)
              QInew=TOT_ICEnew/(THICK+BLDTRH*FLIMASS*VSNOW)
              IF (QInew .LE. EPSQ) QInew=0.
              IF (QI.GT.0. .AND. QInew.NE.0.) THEN
                DUM=QInew/QI
                IF (DUM .LT. TOLER) QInew=0.
              ENDIF
              ASNOWnew=BLDTRH*FLIMASS*VSNOW*QInew
              IF (ASNOW.GT.0. .AND. ASNOWnew.NE.0.) THEN
                DUM=ASNOWnew/ASNOW
                IF (DUM .LT. TOLER) ASNOWnew=0.
              ENDIF
            ENDIF         ! End IF (TOT_ICEnew .LE. CLIMIT)
          ENDIF           ! End IF (ICE_logical)


!
!--- Update rain mixing ratios
!
!---
! * TOT_RAINnew - total mass of rain after microphysics
!                 current layer and the input flux of ice from above
! * VRAIN2      - time-averaged fall speed of rain in grid and below 
!                 (with air resistance correction)
! * QRnew       - updated rain mixing ratio in layer
!      -> TOT_RAINnew=QRnew*(THICK+BLDTRH*VRAIN2)
!  * ARAINnew  - updated accumulation of rain at bottom of grid cell
!---
!
          DELR=PRAUT+PRACW+PIACWR-PIACR+PIMLT+PREVP+PICND
          TOT_RAINnew=TOT_RAIN+THICK*DELR
          IF (TOT_RAIN.GT.0. .AND. TOT_RAINnew.NE.0.) THEN
            DUM=TOT_RAINnew/TOT_RAIN
            IF (DUM .LT. TOLER) TOT_RAINnew=0.
          ENDIF
          IF (TOT_RAINnew .LE. CLIMIT) THEN
            TOT_RAINnew=0.
            VRAIN2=0.
            QRnew=0.
            ARAINnew=0.
          ELSE
   !
   !--- 1st guess time-averaged rain rate at bottom of grid box
   !
            RR=TOT_RAINnew/(DTPH*GAMMAR)
   !
   !--- Use same algorithm as above for calculating mean drop diameter
   !      (IDR, in microns), which is used to estimate the time-averaged
   !      fall speed of rain drops at the bottom of the grid layer.  This
   !      isn't perfect, but the alternative is solving a transcendental 
   !      equation that is numerically inefficient and nasty to program
   !      (coded in earlier versions of GSMCOLUMN prior to 8-22-01).
   !
            IF (RR .LE. RR_DRmin) THEN
              IDR=MDRmin
            ELSE IF (RR .LE. RR_DR1) THEN
              IDR=INT( 1.123E3*RR**.1947 + .5 )
              IDR=MAX( MDRmin, MIN(IDR, MDR1) )
            ELSE IF (RR .LE. RR_DR2) THEN
              IDR=INT( 1.225E3*RR**.2017 + .5 )
              IDR=MAX( MDR1, MIN(IDR, MDR2) )
            ELSE IF (RR .LE. RR_DR3) THEN
              IDR=INT( 1.3006E3*RR**.2083 + .5 )
              IDR=MAX( MDR2, MIN(IDR, MDR3) )
            ELSE IF (RR .LE. RR_DRmax) THEN
              IDR=INT( 1.355E3*RR**.2144 + .5 )
              IDR=MAX( MDR3, MIN(IDR, MDRmax) )
            ELSE 
              IDR=MDRmax
            ENDIF              ! End IF (RR .LE. RR_DRmin)
            VRAIN2=GAMMAR*VRAIN(IDR)
            QRnew=TOT_RAINnew/(THICK+BLDTRH*VRAIN2)
            IF (QRnew .LE. EPSQ) QRnew=0.
            IF (QR.GT.0. .AND. QRnew.NE.0.) THEN
              DUM=QRnew/QR
              IF (DUM .LT. TOLER) QRnew=0.
            ENDIF
            ARAINnew=BLDTRH*VRAIN2*QRnew
            IF (ARAIN.GT.0. .AND. ARAINnew.NE.0.) THEN
              DUM=ARAINnew/ARAIN
              IF (DUM .LT. TOLER) ARAINnew=0.
            ENDIF
          ENDIF
!
          WCnew=QWnew+QRnew+QInew
!
!----------------------------------------------------------------------
!-------------- Begin debugging & verification ------------------------
!----------------------------------------------------------------------
!
!--- QT, QTnew - total water (vapor & condensate) before & after microphysics, resp.
!


          QT=THICK*(WV+WC)+ARAIN+ASNOW
          QTnew=THICK*(WVnew+WCnew)+ARAINnew+ASNOWnew
          BUDGET=QT-QTnew
!
!--- Additional check on budget preservation, accounting for truncation effects
!
          DBG_logical=.FALSE.
!          DUM=ABS(BUDGET)
!          IF (DUM .GT. TOLER) THEN
!            DUM=DUM/MIN(QT, QTnew)
!            IF (DUM .GT. TOLER) DBG_logical=.TRUE.
!          ENDIF
!!
!          DUM=(RHgrd+.001)*QSInew
!          IF ( (QWnew.GT.EPSQ) .OR. QRnew.GT.EPSQ .OR. WVnew.GT.DUM)
!     &        .AND. TC.LT.T_ICE )  DBG_logical=.TRUE.
!
!          IF (TC.GT.5. .AND. QInew.GT.EPSQ) DBG_logical=.TRUE.
!
          IF ((WVnew.LT.EPSQ .OR. DBG_logical) .AND. PRINT_diag) THEN
   !
            WRITE(6,"(/2(a,i4),2(a,i2))") '{} i=',I_index,' j=',J_index,&
     &                                    ' L=',L,' LSFC=',LSFC
   !
            ESW=1000.*FPVS0(Tnew)
            QSWnew=EPS*ESW/(PP-ESW)
            IF (TC.LT.0. .OR. Tnew .LT. 0.) THEN
              ESI=1000.*FPVS(Tnew)
              QSInew=EPS*ESI/(PP-ESI)
            ELSE
              QSI=QSW
              QSInew=QSWnew
            ENDIF
            WSnew=QSInew
            WRITE(6,"(4(a12,g11.4,1x))")                                   &
     & '{} TCold=',TC,'TCnew=',Tnew-T0C,'P=',.01*PP,'RHO=',RHO,            &
     & '{} THICK=',THICK,'RHold=',WV/WS,'RHnew=',WVnew/WSnew,              &
     &   'RHgrd=',RHgrd,                                                   &
     & '{} RHWold=',WV/QSW,'RHWnew=',WVnew/QSWnew,'RHIold=',WV/QSI,        &
     &   'RHInew=',WVnew/QSInew,                                           &
     & '{} QSWold=',QSW,'QSWnew=',QSWnew,'QSIold=',QSI,'QSInew=',QSInew,   &
     & '{} WSold=',WS,'WSnew=',WSnew,'WVold=',WV,'WVnew=',WVnew,           &
     & '{} WCold=',WC,'WCnew=',WCnew,'QWold=',QW,'QWnew=',QWnew,           &
     & '{} QIold=',QI,'QInew=',QInew,'QRold=',QR,'QRnew=',QRnew,           &
     & '{} ARAINold=',ARAIN,'ARAINnew=',ARAINnew,'ASNOWold=',ASNOW,        &
     &   'ASNOWnew=',ASNOWnew,                                             &
     & '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINnew=',TOT_RAINnew,                 &
     &   'TOT_ICE=',TOT_ICE,'TOT_ICEnew=',TOT_ICEnew,                      &
     & '{} BUDGET=',BUDGET,'QTold=',QT,'QTnew=',QTnew                       
   !
            WRITE(6,"(4(a12,g11.4,1x))")                                   &
     & '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,             &
     & '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,       &
     & '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,     &
     & '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',       &
     &    PIMLT,                                                           &
     & '{} PIACR=',PIACR                                                    
   !
            IF (ICE_logical) WRITE(6,"(4(a12,g11.4,1x))")                  &
     & '{} RimeF1=',RimeF1,'GAMMAS=',GAMMAS,'VrimeF=',VrimeF,              &
     &   'VSNOW=',VSNOW,                                                   &
     & '{} INDEXS=',FLOAT(INDEXS),'FLARGE=',FLARGE,'FSMALL=',FSMALL,       &
     &   'FLIMASS=',FLIMASS,                                               &
     & '{} XSIMASS=',XSIMASS,'XLIMASS=',XLIMASS,'QLICE=',QLICE,            &
     &   'QTICE=',QTICE,                                                   &
     & '{} NLICE=',NLICE,'NSmICE=',NSmICE,'PILOSS=',PILOSS,                &
     &   'EMAIRI=',EMAIRI,                                                 &
     & '{} RimeF=',RimeF                                                    
   !
            IF (TOT_RAIN.GT.0. .OR. TOT_RAINnew.GT.0.)                     &
     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
     & '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),               &
     &   'GAMMAR=',GAMMAR,'N0r=',N0r,                                      &
     & '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,   &
     & '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                   &
     &   'VOLR2=',THICK+BLDTRH*VRAIN2
   !
            IF (PRAUT .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} QW0=',QW0
   !
            IF (PRACW .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} FWR=',FWR
   !
            IF (PIACR .GT. 0.) WRITE(6,"(a12,g11.4,1x)") '{} FIR=',FIR
   !
            DUM=PIMLT+PICND-PREVP-PIEVP
            IF (DUM.GT.0. .or. DWVi.NE.0.)                                 &
     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
     & '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                             &
     &   'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
   !
            IF (PREVP .LT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                &
     & '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,     &
     & '{} DWVr=',DWVr,'DENOMW=',DENOMW
   !
            IF (PIDEP.NE.0. .AND. DWVi.NE.0.)                              &
     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
     & '{} DWVi=',DWVi,'DENOMI=',DENOMI,'PIDEP_max=',PIDEP_max,            &
     &   'SFACTOR=',SFACTOR,                                               &
     & '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),           &
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
     & '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
   !
            IF (PIDEP.GT.0. .AND. PCOND.NE.0.)                             &
     &        WRITE(6,"(4(a12,g11.4,1x))")                                 &
     & '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,            &
     &    'DUM2=',PCOND-PIACW
   !
            IF (FWS .GT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                  &
     & '{} FWS=',FWS                     
   !
            DUM=PIMLT+PICND-PIEVP
            IF (DUM.GT. 0.) WRITE(6,"(4(a12,g11.4,1x))")                   &
     & '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),   &
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                                &
     & '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0       
   !
          ENDIF


!
!-----------------------------------------------------------------------
!--------------- Water budget statistics & maximum values --------------
!-----------------------------------------------------------------------
!
          IF (PRINT_diag) THEN
            ITdx=MAX( ITLO, MIN( INT(Tnew-T0C), ITHI ) )
            IF (QInew .GT. EPSQ) NSTATS(ITdx,1)=NSTATS(ITdx,1)+1
            IF (QInew.GT.EPSQ  .AND.  QRnew+QWnew.GT.EPSQ)              &
     &        NSTATS(ITdx,2)=NSTATS(ITdx,2)+1
            IF (QWnew .GT. EPSQ) NSTATS(ITdx,3)=NSTATS(ITdx,3)+1 
            IF (QRnew .GT. EPSQ) NSTATS(ITdx,4)=NSTATS(ITdx,4)+1
  !
            QMAX(ITdx,1)=MAX(QMAX(ITdx,1), QInew)
            QMAX(ITdx,2)=MAX(QMAX(ITdx,2), QWnew)
            QMAX(ITdx,3)=MAX(QMAX(ITdx,3), QRnew)
            QMAX(ITdx,4)=MAX(QMAX(ITdx,4), ASNOWnew)
            QMAX(ITdx,5)=MAX(QMAX(ITdx,5), ARAINnew)
            QTOT(ITdx,1)=QTOT(ITdx,1)+QInew*THICK
            QTOT(ITdx,2)=QTOT(ITdx,2)+QWnew*THICK
            QTOT(ITdx,3)=QTOT(ITdx,3)+QRnew*THICK
  !
            QTOT(ITdx,4)=QTOT(ITdx,4)+PCOND*THICK
            QTOT(ITdx,5)=QTOT(ITdx,5)+PICND*THICK
            QTOT(ITdx,6)=QTOT(ITdx,6)+PIEVP*THICK
            QTOT(ITdx,7)=QTOT(ITdx,7)+PIDEP*THICK
            QTOT(ITdx,8)=QTOT(ITdx,8)+PREVP*THICK
            QTOT(ITdx,9)=QTOT(ITdx,9)+PRAUT*THICK
            QTOT(ITdx,10)=QTOT(ITdx,10)+PRACW*THICK
            QTOT(ITdx,11)=QTOT(ITdx,11)+PIMLT*THICK
            QTOT(ITdx,12)=QTOT(ITdx,12)+PIACW*THICK
            QTOT(ITdx,13)=QTOT(ITdx,13)+PIACWI*THICK
            QTOT(ITdx,14)=QTOT(ITdx,14)+PIACWR*THICK
            QTOT(ITdx,15)=QTOT(ITdx,15)+PIACR*THICK
  !
            QTOT(ITdx,16)=QTOT(ITdx,16)+(WVnew-WV)*THICK
            QTOT(ITdx,17)=QTOT(ITdx,17)+(QWnew-QW)*THICK
            QTOT(ITdx,18)=QTOT(ITdx,18)+(QInew-QI)*THICK
            QTOT(ITdx,19)=QTOT(ITdx,19)+(QRnew-QR)*THICK
            QTOT(ITdx,20)=QTOT(ITdx,20)+(ARAINnew-ARAIN)
            QTOT(ITdx,21)=QTOT(ITdx,21)+(ASNOWnew-ASNOW)
            IF (QInew .GT. 0.)                                          &
     &        QTOT(ITdx,22)=QTOT(ITdx,22)+QInew*THICK/RimeF
  !
          ENDIF
!
!----------------------------------------------------------------------
!------------------------- Update arrays ------------------------------
!----------------------------------------------------------------------
!


          T_col(L)=Tnew                           ! Updated temperature
!
          QV_col(L)=max(EPSQ, WVnew/(1.+WVnew))   ! Updated specific humidity
          WC_col(L)=max(EPSQ, WCnew)              ! Updated total condensate mixing ratio
          QI_col(L)=max(EPSQ, QInew)              ! Updated ice mixing ratio
          QR_col(L)=max(EPSQ, QRnew)              ! Updated rain mixing ratio
          QW_col(L)=max(EPSQ, QWnew)              ! Updated cloud water mixing ratio
          RimeF_col(L)=RimeF                      ! Updated rime factor
          ASNOW=ASNOWnew                          ! Updated accumulated snow
          ARAIN=ARAINnew                          ! Updated accumulated rain
!
!#######################################################################
!
10      CONTINUE         ! ##### End "L" loop through model levels #####


!
!#######################################################################
!
!-----------------------------------------------------------------------
!--------------------------- Return to GSMDRIVE -----------------------
!-----------------------------------------------------------------------
!
        CONTAINS
!#######################################################################
!--------- Produces accurate calculation of cloud condensation ---------
!#######################################################################
!
      REAL FUNCTION CONDENSE (PP, QW, TK, WV)
!
!---------------------------------------------------------------------------------
!------   The Asai (1965) algorithm takes into consideration the release of ------
!------   latent heat in increasing the temperature & in increasing the     ------
!------   saturation mixing ratio (following the Clausius-Clapeyron eqn.).  ------
!---------------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER ::                               &
     & RHLIMIT=.001, RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: COND, SSAT, WCdum
!
      REAL,INTENT(IN) :: QW,PP,WV,TK
      REAL WVdum,Tdum,XLV2,DWV,WS,ESW,XLV1,XLV
integer nsteps
!
!-----------------------------------------------------------------------
!
!--- LV (T) is from Bolton (JAS, 1980)
!
      XLV=3.148E6-2370.*TK
      XLV1=XLV*RCP
      XLV2=XLV*XLV*RCPRV
      Tdum=TK
      WVdum=WV
      WCdum=QW
      ESW=1000.*FPVS0(Tdum)                     ! Saturation vapor press w/r/t water
      WS=RHgrd*EPS*ESW/(PP-ESW)                 ! Saturation mixing ratio w/r/t water
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      CONDENSE=0.
nsteps = 0
      DO WHILE ((SSAT.LT.RHLIMIT1 .AND. WCdum.GT.EPSQ)                  &
     &           .OR. SSAT.GT.RHLIMIT)
        nsteps = nsteps + 1
        COND=DWV/(1.+XLV2*WS/(Tdum*Tdum))       ! Asai (1965, J. Japan)
        COND=MAX(COND, -WCdum)                  ! Limit cloud water evaporation
        Tdum=Tdum+XLV1*COND                     ! Updated temperature
        WVdum=WVdum-COND                        ! Updated water vapor mixing ratio
        WCdum=WCdum+COND                        ! Updated cloud water mixing ratio
        CONDENSE=CONDENSE+COND                  ! Total cloud water condensation
        ESW=1000.*FPVS0(Tdum)                   ! Updated saturation vapor press w/r/t water
        WS=RHgrd*EPS*ESW/(PP-ESW)               ! Updated saturation mixing ratio w/r/t water
        DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
        SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
!
      END FUNCTION CONDENSE
!
!#######################################################################
!---------------- Calculate ice deposition at T<T_ICE ------------------
!#######################################################################
!
      REAL FUNCTION DEPOSIT (PP, Tdum, WVdum)
!
!--- Also uses the Asai (1965) algorithm, but uses a different target
!      vapor pressure for the adjustment
!
      IMPLICIT NONE      
!
      INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT=.001,                 &
     & RHLIMIT1=-RHLIMIT
      REAL (KIND=HIGH_PRES) :: DEP, SSAT
!    
      real,INTENT(IN) ::  PP
      real,INTENT(INOUT) ::  WVdum,Tdum
      real ESI,WS,DWV
!
!-----------------------------------------------------------------------
!
      ESI=1000.*FPVS(Tdum)                      ! Saturation vapor press w/r/t ice
      WS=RHgrd*EPS*ESI/(PP-ESI)                 ! Saturation mixing ratio
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      DEPOSIT=0.
      DO WHILE (SSAT.GT.RHLIMIT .OR. SSAT.LT.RHLIMIT1)
   !
   !--- Note that XLVS2=LS*LV/(CP*RV)=LV*WS/(RV*T*T)*(LS/CP*DEP1), 
   !     where WS is the saturation mixing ratio following Clausius-
   !     Clapeyron (see Asai,1965; Young,1993,p.405) 
   !
        DEP=DWV/(1.+XLS2*WS/(Tdum*Tdum))        ! Asai (1965, J. Japan)
        Tdum=Tdum+XLS1*DEP                      ! Updated temperature
        WVdum=WVdum-DEP                         ! Updated ice mixing ratio
        DEPOSIT=DEPOSIT+DEP                     ! Total ice deposition
        ESI=1000.*FPVS(Tdum)                    ! Updated saturation vapor press w/r/t ice
        WS=RHgrd*EPS*ESI/(PP-ESI)               ! Updated saturation mixing ratio w/r/t ice
        DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
        SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
!
      END FUNCTION DEPOSIT
!
      END SUBROUTINE EGCP01COLUMN 
!#######################################################################
!------- Initialize constants & lookup tables for microphysics ---------
!#######################################################################
!

! SH 0211/2002

!-----------------------------------------------------------------------
!!!   SUBROUTINE FERRIER_INIT (GSMDT,DT,DELX,DELY,LOWLYR,restart,       &
      SUBROUTINE FERRIER_INIT (GSMDT,DT,DELX,DELY,restart,              &
     &   F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,                              &
     &   MP_RESTART_STATE,TBPVS_STATE,TBPVS0_STATE,AVRAIN,              &
     &   ALLOWED_TO_READ,                                               &
     &   IDS,IDE,JDS,JDE,KDS,KDE,                                       &
     &   IMS,IME,JMS,JME,KMS,KME,                                       &
     &   ITS,ITE,JTS,JTE,KTS,KTE                                       )
!-----------------------------------------------------------------------
!-------------------------------------------------------------------------------
!---  SUBPROGRAM DOCUMENTATION BLOCK
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: February 2001
!            Jin             ORG: W/NP22     DATE: January 2002 
!        (Modification for WRF structure)
!                                               
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Reads various microphysical lookup tables used in COLUMN_MICRO
!   * Lookup tables were created "offline" and are read in during execution
!   * Creates lookup tables for saturation vapor pressure w/r/t water & ice
!-------------------------------------------------------------------------------
!     
! USAGE: CALL FERRIER_INIT FROM SUBROUTINE PHYSICS_INITIALIZE
!
!   INPUT ARGUMENT LIST:
!       DTPH - physics time step (s)
!  
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBROUTINES:
!     MY_GROWTH_RATES_NMM - lookup table for growth of nucleated ice
!     GPVS            - lookup tables for saturation vapor pressure (water, ice)
!
!   UNIQUE: NONE
!  
!   LIBRARY: NONE
!  
!   COMMON BLOCKS:
!     CMICRO_CONS - constants used in GSMCOLUMN
!     CMY600       - lookup table for growth of ice crystals in 
!                    water saturated conditions (Miller & Young, 1979)
!     IVENT_TABLES - lookup tables for ventilation effects of ice
!     IACCR_TABLES - lookup tables for accretion rates of ice
!     IMASS_TABLES - lookup tables for mass content of ice
!     IRATE_TABLES - lookup tables for precipitation rates of ice
!     IRIME_TABLES - lookup tables for increase in fall speed of rimed ice
!     MAPOT        - Need lat/lon grid resolution
!     RVENT_TABLES - lookup tables for ventilation effects of rain
!     RACCR_TABLES - lookup tables for accretion rates of rain
!     RMASS_TABLES - lookup tables for mass content of rain
!     RVELR_TABLES - lookup tables for fall speeds of rain
!     RRATE_TABLES - lookup tables for precipitation rates of rain
!   
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!------------------------------------------------------------------------- 
!-------------- Parameters & arrays for lookup tables -------------------- 
!------------------------------------------------------------------------- 
!
!--- Common block of constants used in column microphysics
!
!WRF
!     real DLMD,DPHD
!WRF
!
!-----------------------------------------------------------------------
!--- Parameters & data statement for local calculations
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3
!
!     VARIABLES PASSED IN
      integer,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     & 
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE       
!WRF
!!!!   INTEGER, DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: LOWLYR
!
      real, INTENT(INOUT) ::  AVRAIN
      real, INTENT(IN) ::  DELX,DELY
      real,DIMENSION(*), INTENT(INOUT) :: MP_RESTART_STATE
      real,DIMENSION(NX), INTENT(INOUT) :: TBPVS_STATE,TBPVS0_STATE
      real,DIMENSION(ims:ime, jms:jme, kms:kme),INTENT(OUT) ::          &
     &  F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40
!     integer,DIMENSION(ITLO:ITHI,4),INTENT(INOUT) :: NSTATS
!     real,DIMENSION(ITLO:ITHI,5),INTENT(INOUT) :: QMAX
!     real,DIMENSION(ITLO:ITHI,22),INTENT(INOUT) :: QTOT
!     real,INTENT(INOUT) :: PRECtot(2),PRECmax(2)
      real,INTENT(IN) :: DT,GSMDT
      LOGICAL,INTENT(IN) :: allowed_to_read,restart
!
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      REAL :: BBFR,DTPH,PI,DX,Thour_print
      INTEGER :: I,IM,J,L,K,JTF,KTF,ITF
      INTEGER :: etampnew_unit1
      LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
      LOGICAL :: opened
!!!   LOGICAL , EXTERNAL      :: wrf_dm_on_monitor
      INTEGER :: IRTN,MYPE
      CHARACTER*80 errmess
!
!-----------------------------------------------------------------------
!
      MYPE=MYPE_SHARE
      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(KTE,KDE-1)
      ITF=MIN0(ITE,IDE-1)
!
!     DO J=JTS,JTF
!     DO I=ITS,ITF
!       LOWLYR(I,J)=1
!     ENDDO
!     ENDDO
!    
!     AVRAIN = 0. ! ratko - should we put only if not restart?
!
      IF(.NOT.RESTART)THEN
        DO K = kts,kte
        DO J = jts,jte
        DO I= its,ite
          F_ICE_PHY(i,j,k)=0.
          F_RAIN_PHY(i,j,k)=0.
          F_RIMEF_PHY(i,j,k)=1.
        ENDDO
        ENDDO
        ENDDO
      ENDIF
!    
!-----------------------------------------------------------------------
      IF(ALLOWED_TO_READ)THEN
!-----------------------------------------------------------------------
!
        DX=((DELX)**2+(DELY)**2)**.5/1000.    ! Model resolution at equator (km)
        DX=MIN(100., MAX(5., DX) )
!
!-- Relative humidity threshold for the onset of grid-scale condensation
!!-- 9/1/01:  Assume the following functional dependence for 5 - 100 km resolution:
!!       RHgrd=0.90 for dx=100 km, 0.98 for dx=5 km, where
!        RHgrd=0.90+.08*((100.-DX)/95.)**.5
!
        DTPH=MAX(GSMDT*60.,DT)
        DTPH=NINT(DTPH/DT)*DT
!
!--- Create lookup tables for saturation vapor pressure w/r/t water & ice
!
        CALL GPVS
!
!--- Read in various lookup tables
!
!!!     IF ( wrf_dm_on_monitor() ) THEN
        IF(MYPE==0)THEN
          DO i = 31,99
            INQUIRE ( i , OPENED = opened )
            IF ( .NOT. opened ) THEN
              etampnew_unit1 = i
              GOTO 2061
            ENDIF
          ENDDO
          etampnew_unit1 = -1
 2061     CONTINUE
        ENDIF
!
!!!     CALL wrf_dm_bcast_bytes ( etampnew_unit1 , IWORDSIZE )
        CALL MPI_BCAST(ETAMPNEW_UNIT1,1,MPI_INTEGER,0  &
                      ,MPI_COMM_COMP,IRTN)
!
        IF ( etampnew_unit1 < 0 ) THEN
!!!       CALL wrf_error_fatal ( 'ferrier_init: Can not find unused fortran unit to read in lookup table.' )
          WRITE(0,*)'ferrier_init: Can not find unused fortran unit to read in lookup table.'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
        ENDIF
!
!!!     IF ( wrf_dm_on_monitor() ) THEN
        IF(MYPE==0)THEN
!!was     OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
          OPEN(UNIT=etampnew_unit1,FILE="ETAMPNEW_DATA",                  &
     &        FORM="UNFORMATTED",STATUS="OLD",ERR=9061)
!
          READ(etampnew_unit1) VENTR1
          READ(etampnew_unit1) VENTR2
          READ(etampnew_unit1) ACCRR
          READ(etampnew_unit1) MASSR
          READ(etampnew_unit1) VRAIN
          READ(etampnew_unit1) RRATE
          READ(etampnew_unit1) VENTI1
          READ(etampnew_unit1) VENTI2
          READ(etampnew_unit1) ACCRI
          READ(etampnew_unit1) MASSI
          READ(etampnew_unit1) VSNOWI
          READ(etampnew_unit1) VEL_RF
!        read(etampnew_unit1) my_growth    ! Applicable only for DTPH=180 s
          CLOSE (etampnew_unit1)
        ENDIF
!
!!!     CALL wrf_dm_bcast_bytes ( VENTR1 , size ( VENTR1 ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VENTR2 , size ( VENTR2 ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( ACCRR , size ( ACCRR ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( MASSR , size ( MASSR ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VRAIN , size ( VRAIN ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( RRATE , size ( RRATE ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VENTI1 , size ( VENTI1 ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VENTI2 , size ( VENTI2 ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( ACCRI , size ( ACCRI ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( MASSI , size ( MASSI ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VSNOWI , size ( VSNOWI ) * RWORDSIZE )
!!!     CALL wrf_dm_bcast_bytes ( VEL_RF , size ( VEL_RF ) * RWORDSIZE )
        CALL MPI_BCAST(VENTR1,SIZE(VENTR1),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VENTR2,SIZE(VENTR2),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(ACCRR,SIZE(ACCRR),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(MASSR,SIZE(MASSR),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VRAIN,SIZE(VRAIN),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(RRATE,SIZE(RRATE),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VENTI1,SIZE(VENTI1),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VENTI2,SIZE(VENTI2),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(ACCRI,SIZE(ACCRI),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(MASSI,SIZE(MASSI),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VSNOWI,SIZE(VSNOWI),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
        CALL MPI_BCAST(VEL_RF,SIZE(VEL_RF),MPI_REAL,0  &
                      ,MPI_COMM_COMP,IRTN)
!
!--- Calculates coefficients for growth rates of ice nucleated in water
!    saturated conditions, scaled by physics time step (lookup table)
!
        CALL MY_GROWTH_RATES_NMM (DTPH)
!       CALL MY_GROWTH_RATES_NMM (DTPH,MY_GROWTH_NMM)
!
        PI=ACOS(-1.)
!
!--- Constants associated with Biggs (1953) freezing of rain, as parameterized
!    following Lin et al. (JCAM, 1983) & Reisner et al. (1998, QJRMS).
!
        ABFR=-0.66
        BBFR=100.
        CBFR=20.*PI*PI*BBFR*RHOL*1.E-21
!
!--- CIACW is used in calculating riming rates
!      The assumed effective collection efficiency of cloud water rimed onto
!      ice is =0.5 below:
!
        CIACW=DTPH*0.25*PI*0.5*(1.E5)**C1
!
!--- CIACR is used in calculating freezing of rain colliding with large ice
!      The assumed collection efficiency is 1.0
!
        CIACR=PI*DTPH
!
!--- Based on rain lookup tables for mean diameters from 0.05 to 0.45 mm
!    * Four different functional relationships of mean drop diameter as 
!      a function of rain rate (RR), derived based on simple fits to 
!      mass-weighted fall speeds of rain as functions of mean diameter
!      from the lookup tables.  
!
        RR_DRmin=N0r0*RRATE(MDRmin)     ! RR for mean drop diameter of .05 mm
        RR_DR1=N0r0*RRATE(MDR1)         ! RR for mean drop diameter of .10 mm
        RR_DR2=N0r0*RRATE(MDR2)         ! RR for mean drop diameter of .20 mm
        RR_DR3=N0r0*RRATE(MDR3)         ! RR for mean drop diameter of .32 mm
        RR_DRmax=N0r0*RRATE(MDRmax)     ! RR for mean drop diameter of .45 mm
!
        RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
        RQR_DR1=N0r0*MASSR(MDR1)        ! Rain content for mean drop diameter of .10 mm
        RQR_DR2=N0r0*MASSR(MDR2)        ! Rain content for mean drop diameter of .20 mm
        RQR_DR3=N0r0*MASSR(MDR3)        ! Rain content for mean drop diameter of .32 mm
        RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
        C_N0r0=PI*RHOL*N0r0
        CN0r0=1.E6/C_N0r0**.25
        CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
        CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
!
!--- CRACW is used in calculating collection of cloud water by rain (an
!      assumed collection efficiency of 1.0)
!
        CRACW=DTPH*0.25*PI*1.0
!
        ESW0=1000.*FPVS0(T0C)     ! Saturation vapor pressure at 0C
        RFmax=1.1**Nrime          ! Maximum rime factor allowed
!
!------------------------------------------------------------------------
!--------------- Constants passed through argument list -----------------
!------------------------------------------------------------------------
!
!--- Important parameters for self collection (autoconversion) of 
!    cloud water to rain. 
!
!--- CRAUT is proportional to the rate that cloud water is converted by
!      self collection to rain (autoconversion rate)
!
        CRAUT=1.-EXP(-1.E-3*DTPH)
!
!--- QAUT0 is the threshold cloud content for autoconversion to rain 
!      needed for droplets to reach a diameter of 20 microns (following
!      Manton and Cotton, 1977; Banta and Hanson, 1987, JCAM)
!--- QAUT0=1.2567, 0.8378, or 0.4189 g/m**3 for droplet number concentrations
!          of 300, 200, and 100 cm**-3, respectively
!
        QAUT0=PI*RHOL*NCW*(20.E-6)**3/6.
!
!--- For calculating snow optical depths by considering bulk density of
!      snow based on emails from Q. Fu (6/27-28/01), where optical 
!      depth (T) = 1.5*SWP/(Reff*DENS), SWP is snow water path, Reff 
!      is effective radius, and DENS is the bulk density of snow.
!
!    SWP (kg/m**2)=(1.E-3 kg/g)*SWPrad, SWPrad in g/m**2 used in radiation
!    T = 1.5*1.E3*SWPrad/(Reff*DENS)
!  
!    See derivation for MASSI(INDEXS), note equal to RHO*QSNOW/NSNOW
!
!      SDENS=1.5e3/DENS, DENS=MASSI(INDEXS)/[PI*(1.E-6*INDEXS)**3]
!
        DO I=MDImin,MDImax
          SDENS(I)=PI*1.5E-15*FLOAT(I*I*I)/MASSI(I)
        ENDDO
!
        Thour_print=-DTPH/3600.

! SH 0211/2002
!       IF (PRINT_diag) THEN
       
      !-------- Total and maximum quantities
      !
!         NSTATS=0      ! Microphysical statistics dealing w/ grid-point counts
!         QMAX=0.       ! Microphysical statistics dealing w/ hydrometeor mass
!         QTOT=0.       ! Microphysical statistics dealing w/ hydrometeor mass
!         PRECmax=0.    ! Maximum precip rates (rain, snow) at surface (mm/h)
!         PRECtot=0.    ! Total precipitation (rain, snow) accumulation at surface
!       ENDIF

!wrf
        IF(.NOT.RESTART)THEN
          MP_RESTART_STATE(MY_T1:MY_T2)=MY_GROWTH_NMM(MY_T1:MY_T2)
          MP_RESTART_STATE(MY_T2+1)=C1XPVS0
          MP_RESTART_STATE(MY_T2+2)=C2XPVS0
          MP_RESTART_STATE(MY_T2+3)=C1XPVS
          MP_RESTART_STATE(MY_T2+4)=C2XPVS
          MP_RESTART_STATE(MY_T2+5)=CIACW
          MP_RESTART_STATE(MY_T2+6)=CIACR
          MP_RESTART_STATE(MY_T2+7)=CRACW
          MP_RESTART_STATE(MY_T2+8)=CRAUT
          TBPVS_STATE(1:NX) =TBPVS(1:NX)
          TBPVS0_STATE(1:NX)=TBPVS0(1:NX)
        ENDIF

      ENDIF  ! Allowed_to_read

      RETURN
!
!-----------------------------------------------------------------------
!
9061 CONTINUE
!!!   WRITE( errmess , '(A,I4)' )                                        &
!!!    'module_mp_etanew: error opening ETAMPNEW_DATA on unit '          &
!!!  &, etampnew_unit1
!!!   CALL wrf_error_fatal(errmess)
      WRITE(0,*)' module_mp_etanew: error opening ETAMPNEW_DATA on unit ',etampnew_unit1
!!!   CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
      CALL NMMB_FINALIZE
!
!-----------------------------------------------------------------------
      END SUBROUTINE FERRIER_INIT
!
      SUBROUTINE MY_GROWTH_RATES_NMM (DTPH)
!     SUBROUTINE MY_GROWTH_RATES_NMM (DTPH,MY_GROWTH_NMM)
!
!--- Below are tabulated values for the predicted mass of ice crystals
!    after 600 s of growth in water saturated conditions, based on 
!    calculations from Miller and Young (JAS, 1979).  These values are
!    crudely estimated from tabulated curves at 600 s from Fig. 6.9 of
!    Young (1993).  Values at temperatures colder than -27C were 
!    assumed to be invariant with temperature.  
!
!--- Used to normalize Miller & Young (1979) calculations of ice growth
!    over large time steps using their tabulated values at 600 s.
!    Assumes 3D growth with time**1.5 following eq. (6.3) in Young (1993).
!
      IMPLICIT NONE
!
      REAL,INTENT(IN) :: DTPH
!
      REAL  DT_ICE
      REAL,DIMENSION(35) :: MY_600
!WRF
!
!-----------------------------------------------------------------------
      DATA MY_600 /                                                     &
     & 5.5e-8, 1.4E-7, 2.8E-7, 6.E-7, 3.3E-6,                           & 
     & 2.E-6, 9.E-7, 8.8E-7, 8.2E-7, 9.4e-7,                            & 
     & 1.2E-6, 1.85E-6, 5.5E-6, 1.5E-5, 1.7E-5,                         & 
     & 1.5E-5, 1.E-5, 3.4E-6, 1.85E-6, 1.35E-6,                         & 
     & 1.05E-6, 1.E-6, 9.5E-7, 9.0E-7, 9.5E-7,                          & 
     & 9.5E-7, 9.E-7, 9.E-7, 9.E-7, 9.E-7,                              & 
     & 9.E-7, 9.E-7, 9.E-7, 9.E-7, 9.E-7 /        ! -31 to -35 deg C
!
!-----------------------------------------------------------------------
!
      DT_ICE=(DTPH/600.)**1.5
      MY_GROWTH_NMM=DT_ICE*MY_600
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE MY_GROWTH_RATES_NMM
!
!-----------------------------------------------------------------------
!---------  Old GFS saturation vapor pressure lookup tables  -----------
!-----------------------------------------------------------------------
!
      SUBROUTINE GPVS
!     ******************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    GPVS        COMPUTE SATURATION VAPOR PRESSURE TABLE
!   AUTHOR: N PHILLIPS       W/NP2      DATE: 30 DEC 82
!
! ABSTRACT: COMPUTE SATURATION VAPOR PRESSURE TABLE AS A FUNCTION OF
!   TEMPERATURE FOR THE TABLE LOOKUP FUNCTION FPVS.
!   EXACT SATURATION VAPOR PRESSURES ARE CALCULATED IN SUBPROGRAM FPVSX.
!   THE CURRENT IMPLEMENTATION COMPUTES A TABLE WITH A LENGTH
!   OF 7501 FOR TEMPERATURES RANGING FROM 180.0 TO 330.0 KELVIN.
!
! PROGRAM HISTORY LOG:
!   91-05-07  IREDELL
!   94-12-30  IREDELL             EXPAND TABLE
!   96-02-19  HONG                ICE EFFECT
!   01-11-29  JIN                 MODIFIED FOR WRF
!
! USAGE:  CALL GPVS
!
! SUBPROGRAMS CALLED:
!   (FPVSX)  - INLINABLE FUNCTION TO COMPUTE SATURATION VAPOR PRESSURE
!
! COMMON BLOCKS:
!   COMPVS   - SCALING PARAMETERS AND TABLE FOR FUNCTION FPVS.
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      IMPLICIT NONE
      real :: X,XINC,T
      integer :: JX
!----------------------------------------------------------------------
      XINC=(XMAX-XMIN)/(NX-1)
      C1XPVS=1.-XMIN/XINC
      C2XPVS=1./XINC
      C1XPVS0=1.-XMIN/XINC
      C2XPVS0=1./XINC
!
      DO JX=1,NX
        X=XMIN+(JX-1)*XINC
        T=X
        TBPVS(JX)=FPVSX(T)
        TBPVS0(JX)=FPVSX0(T)
      ENDDO
! 
      END SUBROUTINE GPVS
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
                     REAL   FUNCTION FPVS(T)
!-----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    FPVS        COMPUTE SATURATION VAPOR PRESSURE
!   AUTHOR: N PHILLIPS            W/NP2      DATE: 30 DEC 82
!
! ABSTRACT: COMPUTE SATURATION VAPOR PRESSURE FROM THE TEMPERATURE.
!   A LINEAR INTERPOLATION IS DONE BETWEEN VALUES IN A LOOKUP TABLE
!   COMPUTED IN GPVS. SEE DOCUMENTATION FOR FPVSX FOR DETAILS.
!   INPUT VALUES OUTSIDE TABLE RANGE ARE RESET TO TABLE EXTREMA.
!   THE INTERPOLATION ACCURACY IS ALMOST 6 DECIMAL PLACES.
!   ON THE CRAY, FPVS IS ABOUT 4 TIMES FASTER THAN EXACT CALCULATION.
!   THIS FUNCTION SHOULD BE EXPANDED INLINE IN THE CALLING ROUTINE.
!
! PROGRAM HISTORY LOG:
!   91-05-07  IREDELL             MADE INTO INLINABLE FUNCTION
!   94-12-30  IREDELL             EXPAND TABLE
!   96-02-19  HONG                ICE EFFECT
!   01-11-29  JIN                 MODIFIED FOR WRF
!
! USAGE:   PVS=FPVS(T)
!
!   INPUT ARGUMENT LIST:
!     T        - REAL TEMPERATURE IN KELVIN
!
!   OUTPUT ARGUMENT LIST:
!     FPVS     - REAL SATURATION VAPOR PRESSURE IN KILOPASCALS (CB)
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!$$$
      IMPLICIT NONE
      real,INTENT(IN) :: T
      real XJ
      integer :: JX
!-----------------------------------------------------------------------
      XJ=MIN(MAX(C1XPVS+C2XPVS*T,1.),FLOAT(NX))
      JX=MIN(XJ,NX-1.)
      FPVS=TBPVS(JX)+(XJ-JX)*(TBPVS(JX+1)-TBPVS(JX))
!
      END FUNCTION FPVS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
                     REAL FUNCTION FPVS0(T)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      real,INTENT(IN) :: T
      real :: XJ1
      integer :: JX1
!-----------------------------------------------------------------------
      XJ1=MIN(MAX(C1XPVS0+C2XPVS0*T,1.),FLOAT(NX))
      JX1=MIN(XJ1,NX-1.)
      FPVS0=TBPVS0(JX1)+(XJ1-JX1)*(TBPVS0(JX1+1)-TBPVS0(JX1))
!
      END FUNCTION FPVS0
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
                    REAL FUNCTION FPVSX(T)
!-----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    FPVSX       COMPUTE SATURATION VAPOR PRESSURE
!   AUTHOR: N PHILLIPS            W/NP2      DATE: 30 DEC 82
!
! ABSTRACT: EXACTLY COMPUTE SATURATION VAPOR PRESSURE FROM TEMPERATURE.
!   THE WATER MODEL ASSUMES A PERFECT GAS, CONSTANT SPECIFIC HEATS
!   FOR GAS AND LIQUID, AND NEGLECTS THE VOLUME OF THE LIQUID.
!   THE MODEL DOES ACCOUNT FOR THE VARIATION OF THE LATENT HEAT
!   OF CONDENSATION WITH TEMPERATURE.  THE ICE OPTION IS NOT INCLUDED.
!   THE CLAUSIUS-CLAPEYRON EQUATION IS INTEGRATED FROM THE TRIPLE POINT
!   TO GET THE FORMULA
!       PVS=PSATK*(TR**XA)*EXP(XB*(1.-TR))
!   WHERE TR IS TTP/T AND OTHER VALUES ARE PHYSICAL CONSTANTS
!   THIS FUNCTION SHOULD BE EXPANDED INLINE IN THE CALLING ROUTINE.
!
! PROGRAM HISTORY LOG:
!   91-05-07  IREDELL             MADE INTO INLINABLE FUNCTION
!   94-12-30  IREDELL             EXACT COMPUTATION
!   96-02-19  HONG                ICE EFFECT 
!   01-11-29  JIN                 MODIFIED FOR WRF
!
! USAGE:   PVS=FPVSX(T)
! REFERENCE:   EMANUEL(1994),116-117
!
!   INPUT ARGUMENT LIST:
!     T        - REAL TEMPERATURE IN KELVIN
!
!   OUTPUT ARGUMENT LIST:
!     FPVSX    - REAL SATURATION VAPOR PRESSURE IN KILOPASCALS (CB)
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!$$$
      IMPLICIT NONE
!-----------------------------------------------------------------------
       real, parameter :: TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2   &
      ,         CLIQ=4.1855E+3,CVAP= 1.8460E+3                          &
      ,         CICE=2.1060E+3,HSUB=2.8340E+6
!
      real, parameter :: PSATK=PSAT*1.E-3
      real, parameter :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
      real, parameter :: DLDTI=CVAP-CICE                                &
      ,                  XAI=-DLDTI/RV,XBI=XAI+HSUB/(RV*TTP)
      real T,TR
!-----------------------------------------------------------------------
      TR=TTP/T
!
      IF(T.GE.TTP)THEN
        FPVSX=PSATK*(TR**XA)*EXP(XB*(1.-TR))
      ELSE
        FPVSX=PSATK*(TR**XAI)*EXP(XBI*(1.-TR))
      ENDIF
! 
      END FUNCTION FPVSX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
                 REAL   FUNCTION FPVSX0(T)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      real,parameter :: TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2     &
      ,         CLIQ=4.1855E+3,CVAP=1.8460E+3                           &
      ,         CICE=2.1060E+3,HSUB=2.8340E+6
      real,PARAMETER :: PSATK=PSAT*1.E-3
      real,PARAMETER :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
      real,PARAMETER :: DLDTI=CVAP-CICE                                 &
      ,                 XAI=-DLDT/RV,XBI=XA+HSUB/(RV*TTP)
      real :: T,TR
!-----------------------------------------------------------------------
      TR=TTP/T
      FPVSX0=PSATK*(TR**XA)*EXP(XB*(1.-TR))
!
      END FUNCTION FPVSX0
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
!===================================================================
!
  SUBROUTINE wsm3(th, q, qci, qrs                                  &
                   , w, den, pii, p, delz                          &
                   , delt,g, cpd, cpv, rd, rv, t0c                 &
                   , ep1, ep2, qmin                                &
                   , XLS, XLV0, XLF0, den0, denr                   &
                   , cliq,cice,psat                                &
                   , rain, rainncv                                 &
                   ,snow, snowncv                                  &
                   ,sr                                             &
                   , ids,ide, jds,jde, kds,kde                     &
                   , ims,ime, jms,jme, kms,kme                     &
                   , its,ite, jts,jte, kts,kte                     &
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!
!  This code is a 3-class simple ice microphyiscs scheme (WSM3) of the WRF
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  Production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM3 cloud scheme
!
!  Coded by Song-You Hong (Yonsei Univ.)
!             Jimy Dudhia (NCAR) and Shu-Hua Chen (UC Davis)
!             Summer 2002
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!             Summer 2003
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Dudhia (D89, 1989) J. Atmos. Sci.
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                             qci, &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                       w, &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv

! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte ) ::   t
  INTEGER ::               i,j,k
!-------------------------------------------------------------------
      DO j=jts,jte
         DO k=kts,kte
         DO i=its,ite
            t(i,k)=th(i,k,j)*pii(i,k,j)
!     if(i==101.and.j==163)then
!       write(0,*)' enter wsm3 k=',k,' t=',t(i,k),' th=',th(i,k,j),' pii=',pii(i,k,j)
!     endif
         ENDDO
         ENDDO
         CALL wsm32D(t, q(ims,kms,j), qci(ims,kms,j)               &
                    ,qrs(ims,kms,j),w(ims,kms,j), den(ims,kms,j)   &
                    ,p(ims,kms,j), delz(ims,kms,j)                 &
                    ,delt,g, cpd, cpv, rd, rv, t0c                 &
                    ,ep1, ep2, qmin                                &
                    ,XLS, XLV0, XLF0, den0, denr                   &
                    ,cliq,cice,psat                                &
                    ,j                                             &
                    ,rain(ims,j), rainncv(ims,j)                   &
                    ,sr(ims,j)                                     &
                    ,ids,ide, jds,jde, kds,kde                     &
                    ,ims,ime, jms,jme, kms,kme                     &
                    ,its,ite, jts,jte, kts,kte                     &
                    ,snow(ims,j),snowncv(ims,j)                    &
                                                                   )
         DO K=kts,kte
         DO I=its,ite
            th(i,k,j)=t(i,k)/pii(i,k,j)
!     if(i==101.and.j==163)then
!       write(0,*)' exit wsm3 k=',k,' t=',t(i,k),' pii=',pii(i,k,j)
!     endif
         ENDDO
         ENDDO
      ENDDO
  END SUBROUTINE wsm3
!===================================================================
!
  SUBROUTINE wsm32D(t, q, qci, qrs,w, den, p, delz                &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain, rainncv                                 &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q, &
                                                             qci, &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(IN   ) ::                                       w, &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr

  REAL, DIMENSION( ims:ime ),     OPTIONAL,                       &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
        rh, qs, denfac, rslope, rslope2, rslope3, rslopeb,        &
        pgen, paut, pacr, pisd, pres, pcon, fall, falk,           &
        xl, cpm, work1, work2, xni, qs0, n0sfac
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
              falkc, work1c, work2c, fallc
! variables for optimization
  REAL, DIMENSION( its:ite )           :: tvec1
  INTEGER, DIMENSION( its:ite ) :: mstep, numdt
  LOGICAL, DIMENSION( its:ite ) :: flgcld
  REAL  ::  pi,                                                   &
            cpmcal, xlcal, lamdar, lamdas, diffus,                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            fallsum, fallsum_qsi, vt2i,vt2s,acrfac,               &      
            qdt, pvt, qik, delq, facq, qrsci, frzmlt,             &
            snomlt, hold, holdrs, facqci, supcol, coeres,         &
            supsat, dtcld, xmi, qciik, delqci, eacrs, satdt,      &
            qimax, diameter, xni0, roqi0, supice
  REAL  :: holdc, holdci
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, kk, n
! Temporaries used for inlining fpvs function
  REAL  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
!
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
!
! Optimizatin : A**B => exp(log(A)*(B))
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
!      venfac(a,b,c) = (viscos(b,c)/diffus(b,a))**(.3333333)       &
!                      /viscos(b,c)**(.5)*(den0/c)**0.25
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))    &
                     /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!
      pi = 4. * atan(1.)
!
!----------------------------------------------------------------
!     paddint 0 for negative values generated by dynamics
!
      do k = kts, kte
        do i = its, ite
          qci(i,k) = max(qci(i,k),0.0)
          qrs(i,k) = max(qrs(i,k),0.0)
        enddo
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = kts, kte
        do i = its, ite
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo
!
      do k = kts, kte
        do i = its, ite
          denfac(i,k) = sqrt(den0/den(i,k))
        enddo
      enddo
!!!   do k = kts, kte
!!!     CALL VSREC( tvec1(its), den(its,k), ite-its+1)
!!!     do i = its, ite
!!!       tvec1(i) = tvec1(i)*den0
!!!     enddo
!!!     CALL VSSQRT( denfac(its,k), tvec1(its), ite-its+1)
!!!   enddo
!
! Inline expansion for fpvs
!         qs(i,k) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs0(i,k) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      cvap = cpv
      hvap=xlv0
      hsub=xls
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite
!         tr=ttp/t(i,k)
!         if(t(i,k).lt.ttp) then
!           qs(i,k) =psat*(tr**xai)*exp(xbi*(1.-tr))
!         else
!           qs(i,k) =psat*(tr**xa)*exp(xb*(1.-tr))
!         endif
!         qs0(i,k)  =psat*(tr**xa)*exp(xb*(1.-tr))
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k) =psat*(exp(log(tr)*(xai)))*exp(xbi*(1.-tr))
          else
            qs(i,k) =psat*(exp(log(tr)*(xa)))*exp(xb*(1.-tr))
          endif
          qs0(i,k)  =psat*(exp(log(tr)*(xa)))*exp(xb*(1.-tr))
          qs0(i,k) = (qs0(i,k)-qs(i,k))/qs(i,k)
          qs(i,k) = ep2 * qs(i,k) / (p(i,k) - qs(i,k))
          qs(i,k) = max(qs(i,k),qmin)
          rh(i,k) = max(q(i,k) / qs(i,k),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
        do i = its, ite
          pres(i,k) = 0.
          paut(i,k) = 0.
          pacr(i,k) = 0.
          pgen(i,k) = 0.
          pisd(i,k) = 0.
          pcon(i,k) = 0.
          fall(i,k) = 0.
          falk(i,k) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).ge.t0c) then
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopermax
              rslopeb(i,k) = rsloperbmax
              rslope2(i,k) = rsloper2max
              rslope3(i,k) = rsloper3max
            else
              rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
!             rslopeb(i,k) = rslope(i,k)**bvtr
              rslopeb(i,k) = exp(log(rslope(i,k))*(bvtr))
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          else
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopesmax
              rslopeb(i,k) = rslopesbmax
              rslope2(i,k) = rslopes2max
              rslope3(i,k) = rslopes3max
            else
              rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
!             rslopeb(i,k) = rslope(i,k)**bvts
              rslopeb(i,k) = exp(log(rslope(i,k))*(bvts))
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          endif
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                           &
!                   *max(qci(i,k),qmin))**0.75,1.e3),1.e6)
          xni(i,k) = min(max(5.38e7*exp(log((den(i,k)*max(qci(i,k),qmin)))*(0.75)),1.e3),1.e6)
        enddo
      enddo
!
      mstepmax = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          if(t(i,k).lt.t0c) then
            pvt = pvts
          else
            pvt = pvtr
          endif
          work1(i,k) = pvt*rslopeb(i,k)*denfac(i,k)
          work2(i,k) = work1(i,k)/delz(i,k)
          numdt(i) = max(nint(work2(i,k)*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = its, ite
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = kte
        do i = its, ite
          if(n.le.mstep(i)) then
            falk(i,k) = den(i,k)*qrs(i,k)*work2(i,k)/mstep(i)
            hold = falk(i,k)
            fall(i,k) = fall(i,k)+falk(i,k)
            holdrs = qrs(i,k)
            qrs(i,k) = max(qrs(i,k)-falk(i,k)*dtcld/den(i,k),0.)
          endif
        enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              falk(i,k) = den(i,k)*qrs(i,k)*work2(i,k)/mstep(i)
              hold = falk(i,k)
              fall(i,k) = fall(i,k)+falk(i,k)
              holdrs = qrs(i,k)
              qrs(i,k) = max(qrs(i,k)-(falk(i,k)                        &
                        -falk(i,k+1)*delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
      enddo
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      mstepmax = 1
      mstep = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          if(t(i,k).lt.t0c.and.qci(i,k).gt.0.) then
            xmi = den(i,k)*qci(i,k)/xni(i,k)
!           diameter  = dicon * sqrt(xmi)
!           work1c(i,k) = 1.49e4*diameter**1.31
            diameter  = max(dicon * sqrt(xmi), 1.e-25)
            work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          else
            work1c(i,k) = 0.
          endif
          if(qci(i,k).le.0.) then
            work2c(i,k) = 0.
          else
            work2c(i,k) = work1c(i,k)/delz(i,k)
          endif
          numdt(i) = max(nint(work2c(i,k)*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = its, ite
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = kte
        do i = its, ite
          if (n.le.mstep(i)) then
            falkc(i,k) = den(i,k)*qci(i,k)*work2c(i,k)/mstep(i)
            holdc = falkc(i,k)
            fallc(i,k) = fallc(i,k)+falkc(i,k)
            holdci = qci(i,k)
            qci(i,k) = max(qci(i,k)-falkc(i,k)*dtcld/den(i,k),0.)
          endif
        enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if (n.le.mstep(i)) then
              falkc(i,k) = den(i,k)*qci(i,k)*work2c(i,k)/mstep(i)
              holdc = falkc(i,k)
              fallc(i,k) = fallc(i,k)+falkc(i,k)
              holdci = qci(i,k)
              qci(i,k) = max(qci(i,k)-(falkc(i,k)                       &
                        -falkc(i,k+1)*delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the freezing/melting term. [D89 B16-B17]
!     freezing occurs one layer above the melting level
!
      do i = its, ite
        mstep(i) = 0
      enddo
      do k = kts, kte
!
        do i = its, ite
          if(t(i,k).ge.t0c) then
            mstep(i) = k
          endif
        enddo
      enddo
!
      do i = its, ite
        if(mstep(i).ne.0.and.w(i,mstep(i)).gt.0.) then
          work1(i,1) = float(mstep(i) + 1)
          work1(i,2) = float(mstep(i))
        else
          work1(i,1) = float(mstep(i))
          work1(i,2) = float(mstep(i))
        endif
      enddo
!
      do i = its, ite
        k  = nint(work1(i,1))
        kk = nint(work1(i,2))
        if(k*kk.ge.1) then
          qrsci = qrs(i,k) + qci(i,k)
          if(qrsci.gt.0..or.fall(i,kk).gt.0.) then
            frzmlt = min(max(-w(i,k)*qrsci/delz(i,k),-qrsci/dtcld),    &
                     qrsci/dtcld)
            snomlt = min(max(fall(i,kk)/den(i,kk),-qrs(i,k)/dtcld),    &
                     qrs(i,k)/dtcld)
            if(k.eq.kk) then
              t(i,k) = t(i,k) - xlf0/cpm(i,k)*(frzmlt+snomlt)*dtcld
            else
              t(i,k) = t(i,k) - xlf0/cpm(i,k)*frzmlt*dtcld
              t(i,kk) = t(i,kk) - xlf0/cpm(i,kk)*snomlt*dtcld
            endif
          endif
        endif
      enddo
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
      do i = its, ite
        fallsum = fall(i,1)
        fallsum_qsi = 0.
        if((t0c-t(i,1)).gt.0) then
        fallsum = fallsum+fallc(i,1)
        fallsum_qsi = fall(i,1)+fallc(i,1)
        endif
        rainncv(i) = 0.
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,1)/denr*dtcld*1000.
          rain(i) = fallsum*delz(i,1)/denr*dtcld*1000.                &
                  + rain(i)
        endif
        IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
        snowncv(i) = 0.
        if(fallsum_qsi.gt.0.) then
          snowncv(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.
          snow(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i)
        endif
        ENDIF
        sr(i) = 0.
        if(fallsum.gt.0.)sr(i)=fallsum_qsi*delz(i,kts)/denr*dtcld*1000./(rainncv(i)+1.e-12)
      enddo
!
!----------------------------------------------------------------
!     rsloper: reverse of the slope parameter of the rain(m)
!     xka:    thermal conductivity of air(jm-1s-1k-1)
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
      do k = kts, kte
        do i = its, ite
          if(t(i,k).ge.t0c) then
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopermax
              rslopeb(i,k) = rsloperbmax
              rslope2(i,k) = rsloper2max
              rslope3(i,k) = rsloper3max
            else
              rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
              rslopeb(i,k) = exp(log(rslope(i,k))*(bvtr))
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          else
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopesmax
              rslopeb(i,k) = rslopesbmax
              rslope2(i,k) = rslopes2max
              rslope3(i,k) = rslopes3max
            else
              rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
              rslopeb(i,k) = exp(log(rslope(i,k))*(bvts))
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          if(t(i,k).ge.t0c) then
            work1(i,k) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k))
          else
            work1(i,k) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k))
          endif
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k)
          satdt = supsat/dtcld
          if(t(i,k).ge.t0c) then
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
            if(qci(i,k).gt.qc0) then
!             paut(i,k) = qck1*qci(i,k)**(7./3.)
              paut(i,k) = qck1*exp(log(qci(i,k))*((7./3.)))
              paut(i,k) = min(paut(i,k),qci(i,k)/dtcld)
            endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [D89 B15]
!        (C->R)
!---------------------------------------------------------------
            if(qrs(i,k).gt.qcrmin.and.qci(i,k).gt.qmin) then
                pacr(i,k) = min(pacrr*rslope3(i,k)*rslopeb(i,k)          &
                     *qci(i,k)*denfac(i,k),qci(i,k)/dtcld)
            endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
            if(qrs(i,k).gt.0.) then
                coeres = rslope2(i,k)*sqrt(rslope(i,k)*rslopeb(i,k))
                pres(i,k) = (rh(i,k)-1.)*(precr1*rslope2(i,k)            &
                         +precr2*work2(i,k)*coeres)/work1(i,k)
              if(pres(i,k).lt.0.) then
                pres(i,k) = max(pres(i,k),-qrs(i,k)/dtcld)
                pres(i,k) = max(pres(i,k),satdt/2)
              else
                pres(i,k) = min(pres(i,k),satdt/2)
              endif
            endif
          else
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
            supcol = t0c-t(i,k)
            ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!           xni(i,k) = min(max(5.38e7*(den(i,k)                         &
!                     *max(qci(i,k),qmin))**0.75,1.e3),1.e6)
            xni(i,k) = min(max(5.38e7*exp(log((den(i,k)*max(qci(i,k),qmin)))*(0.75)),1.e3),1.e6)
            eacrs = exp(0.07*(-supcol))
!
            if(qrs(i,k).gt.qcrmin.and.qci(i,k).gt.qmin) then
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter  = min(dicon * sqrt(xmi),dimax)
              vt2i = 1.49e4*diameter**1.31
              vt2s = pvts*rslopeb(i,k)*denfac(i,k)
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k)+2.*diameter*rslope2(i,k)          &
                      +diameter**2*rslope(i,k)
              pacr(i,k) = min(pi*qci(i,k)*eacrs*n0s*n0sfac(i,k)          &
                       *abs(vt2s-vt2i)*acrfac/4.,qci(i,k)/dtcld)
            endif
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter = dicon * sqrt(xmi)
              pisd(i,k) = 4.*diameter*xni(i,k)*(rh(i,k)-1.)/work1(i,k)
              if(pisd(i,k).lt.0.) then
                pisd(i,k) = max(pisd(i,k),satdt/2)
                pisd(i,k) = max(pisd(i,k),-qci(i,k)/dtcld)
              else
                pisd(i,k) = min(pisd(i,k),satdt/2)
              endif
              if(abs(pisd(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k)*sqrt(rslope(i,k)*rslopeb(i,k))
              pres(i,k) = (rh(i,k)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k)   &
                        +precs2*work2(i,k)*coeres)/work1(i,k)
              supice = satdt-pisd(i,k)
              if(pres(i,k).lt.0.) then
                pres(i,k) = max(pres(i,k),-qrs(i,k)/dtcld)
                pres(i,k) = max(max(pres(i,k),satdt/2),supice)
              else
                pres(i,k) = min(min(pres(i,k),satdt/2),supice)
              endif
              if(abs(pisd(i,k)+pres(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-pisd(i,k)-pres(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
!             roqi0 = 4.92e-11*xni0**1.33
              roqi0 = 4.92e-11*exp(log(xni0)*(1.33))
              pgen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k),0.))/dtcld)
              pgen(i,k) = min(min(pgen(i,k),satdt),supice)
            endif
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
              qimax = roqimax/den(i,k)
              paut(i,k) = max(0.,(qci(i,k)-qimax)/dtcld)
            endif
          endif
        enddo
      enddo
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = kts, kte
        do i = its, ite
          qciik = max(qmin,qci(i,k))
          delqci = (paut(i,k)+pacr(i,k)-pgen(i,k)-pisd(i,k))*dtcld
          if(delqci.ge.qciik) then
            facqci = qciik/delqci
            paut(i,k) = paut(i,k)*facqci
            pacr(i,k) = pacr(i,k)*facqci
            pgen(i,k) = pgen(i,k)*facqci
            pisd(i,k) = pisd(i,k)*facqci
          endif
          qik = max(qmin,q(i,k))
          delq = (pres(i,k)+pgen(i,k)+pisd(i,k))*dtcld
          if(delq.ge.qik) then
            facq = qik/delq
            pres(i,k) = pres(i,k)*facq
            pgen(i,k) = pgen(i,k)*facq
            pisd(i,k) = pisd(i,k)*facq
          endif
          work2(i,k) = -pres(i,k)-pgen(i,k)-pisd(i,k)
          q(i,k) = q(i,k)+work2(i,k)*dtcld
          qci(i,k) = max(qci(i,k)-(paut(i,k)+pacr(i,k)-pgen(i,k)     &
                   -pisd(i,k))*dtcld,0.)
          qrs(i,k) = max(qrs(i,k)+(paut(i,k)+pacr(i,k)               &
                   +pres(i,k))*dtcld,0.)
          if(t(i,k).lt.t0c) then
            t(i,k) = t(i,k)-xls*work2(i,k)/cpm(i,k)*dtcld
          else
            t(i,k) = t(i,k)-xl(i,k)*work2(i,k)/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
      cvap = cpv
      hvap = xlv0
      hsub = xls
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite
          tr=ttp/t(i,k)
!         qs(i,k)=psat*(tr**xa)*exp(xb*(1.-tr))
          qs(i,k)=psat*(exp(log(tr)*(xa)))*exp(xb*(1.-tr))
          qs(i,k) = ep2 * qs(i,k) / (p(i,k) - qs(i,k))
          qs(i,k) = max(qs(i,k),qmin)
          denfac(i,k) = sqrt(den0/den(i,k))
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
      do k = kts, kte
        do i = its, ite
          work1(i,k) = conden(t(i,k),q(i,k),qs(i,k),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k)+work1(i,k)
          pcon(i,k) = min(max(work1(i,k),0.),max(q(i,k),0.))/dtcld
          if(qci(i,k).gt.0..and.work1(i,k).lt.0.and.t(i,k).gt.t0c)      &
            pcon(i,k) = max(work1(i,k),-qci(i,k))/dtcld
          q(i,k) = q(i,k)-pcon(i,k)*dtcld
          qci(i,k) = max(qci(i,k)+pcon(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcon(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = kts, kte
        do i = its, ite
          if(qci(i,k).le.qmin) qci(i,k) = 0.0
        enddo
      enddo
!
      enddo                  ! big loops
  END SUBROUTINE wsm32D
! ...................................................................
      REAL FUNCTION rgmma(x)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!     rgmma function:  use infinite product form
      REAL :: euler
      PARAMETER (euler=0.577215664901532)
      REAL :: x, y
      INTEGER :: i
      if(x.eq.1.)then
        rgmma=0.
          else
        rgmma=x*exp(euler*x)
        do i=1,10000
          y=float(i)
          rgmma=rgmma*(1.000+x/y)*exp(-x/y)
        enddo
        rgmma=1./rgmma
      endif
      END FUNCTION rgmma
!
!--------------------------------------------------------------------------
      REAL FUNCTION fpvs_1(t,ice,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c)
!--------------------------------------------------------------------------
      IMPLICIT NONE
!--------------------------------------------------------------------------
      REAL t,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c,dldt,xa,xb,dldti,   &
           xai,xbi,ttp,tr
      INTEGER ice
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      tr=ttp/t
      if(t.lt.ttp.and.ice.eq.1) then
        fpvs_1=psat*(tr**xai)*exp(xbi*(1.-tr))
      else
        fpvs_1=psat*(tr**xa)*exp(xb*(1.-tr))
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END FUNCTION fpvs_1
!-------------------------------------------------------------------
  SUBROUTINE wsm3init(den0,denr,dens,cl,cpv,allowed_to_read)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!.... constants which may not be tunable
   REAL, INTENT(IN) :: den0,denr,dens,cl,cpv
   LOGICAL, INTENT(IN) :: allowed_to_read
   REAL :: pi
!
   pi = 4.*atan(1.)
   xlv1 = cl-cpv
!
   qc0  = 4./3.*pi*denr*r0**3*xncr/den0  ! 0.419e-3 -- .61e-3
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
!
   bvtr1 = 1.+bvtr
   bvtr2 = 2.5+.5*bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            ! 17.837825
   g5pbro2 = rgmma(bvtr2)          ! 1.8273
   pvtr = avtr*g4pbr/6.
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
   xm0  = (di0/dicon)**2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r =  pi*denr*n0r
   pidn0s =  pi*dens*n0s
!
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
!
  END SUBROUTINE wsm3init
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      END MODULE MODULE_MICROPHYSICS_NMM

!-----------------------------------------------------------------------
