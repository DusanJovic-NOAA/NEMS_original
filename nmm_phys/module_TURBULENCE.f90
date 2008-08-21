!-----------------------------------------------------------------------
!
      MODULE MODULE_TURBULENCE
!
!-----------------------------------------------------------------------
!
!***  THE OUTER DRIVER FOR THE SFC LAYER, LSM, AND FULL 3-D TURBULENCE
!***  PLUS THE WRF TURBULENCE DRIVER AND THE VARIOUS TURBULENCE SCHEMES.
!
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!   2008-07-28  Vasic - Turned off counters (now computed in
!                       SET_INTERNAL_STATE_PHY).
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
      USE MODULE_DM_PARALLEL,ONLY : ITS_B1,ITE_B1,ITE_B2                &
                                   ,ITS_B1_H1,ITE_B1_H1,ITE_B1_H2       &
                                   ,ITS_B1_H2,ITE_H1,ITE_H2             &
                                   ,JTS_B1,JTE_B1,JTE_B2                &
                                   ,JTS_B1_H1,JTE_B1_H1,JTE_B1_H2       &
                                   ,JTS_B1_H2,JTE_H2                    &
                                   ,MPI_COMM_COMP                       &
                                   ,MYPE_SHARE,NUM_TILES
!
      USE MODULE_LANDSURFACE
      USE MODULE_SURFACE_LAYER
      USE MODULE_GWD
!
      USE MODULE_H_TO_V
!
      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,CP,ELIV,ELWV,EP_1,EPSQ,EPSQ2,G  &
                                 ,P608,PI,PQ0,R_D,R_V,RHOWATER,STBOLT
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
      PUBLIC :: MYJPBL_INIT,TURBL
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE TURBULENCE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: YSUSCHEME=1                       &
                                     ,MYJPBLSCHEME=2                    &
                                     ,GFSSCHEME=3                       &
                                     ,MRFSCHEME=99
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR MYJ TURBULENCE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER :: ITRMX=5 ! Iteration count for mixing length computation
!
      REAL,PARAMETER :: VKARMAN=0.4
      REAL,PARAMETER :: CAPA=R_D/CP
      REAL,PARAMETER :: XLS=ELIV,XLV=ELWV
      REAL,PARAMETER :: RLIVWV=XLS/XLV,ELOCP=2.72E6/CP
      REAL,PARAMETER :: EPS1=1.E-12,EPS2=0.
      REAL,PARAMETER :: EPSL=0.32,EPSRU=1.E-7,EPSRS=1.E-7               &
                       ,EPSTRB=1.E-24
      REAL,PARAMETER :: EPSA=1.E-8,EPSIT=1.E-4,EPSU2=1.E-4,EPSUST=0.07  &
                       ,FH=1.01
      REAL,PARAMETER :: ALPH=0.30,BETA=1./273.,EL0MAX=1000.,EL0MIN=1.   &
                       ,ELFC=0.23*0.5,GAM1=0.2222222222222222222        &
                       ,PRT=1.
      REAL,PARAMETER :: A1=0.659888514560862645                         &
                       ,A2x=0.6574209922667784586                       &
                       ,B1=11.87799326209552761                         &
                       ,B2=7.226971804046074028                         &
                       ,C1=0.000830955950095854396
      REAL,PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
      REAL,PARAMETER :: ELZ0=0.,ESQ=5.0,EXCM=0.001                      &
                       ,FHNEU=0.8,GLKBR=10.,GLKBS=30.                   &
                       ,QVISC=2.1E-5,RFC=0.191,RIC=0.505,SMALL=0.35     &
                       ,SQPR=0.84,SQSC=0.84,SQVISC=258.2,TVISC=2.1E-5   &
                       ,USTC=0.7,USTR=0.225,VISC=1.5E-5                 &
                       ,WOLD=0.15,WWST=1.2,ZTMAX=1.,ZTFC=1.,ZTMIN=-5.
!
      REAL,PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC
!
      REAL,PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                     &
!                      ,EP_1=R_V/R_D-1.,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS  &
                       ,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS                  &
                       ,RB1=1./B1,RTVISC=1./TVISC,RVISC=1./VISC         &
                       ,ZQRZT=SQSC/SQPR
!
      REAL,PARAMETER :: ADNH= 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG      &
                       ,ADNM=18.*A1*A1*A2x*(B2-3.*A2x)*BTG              &
                       ,ANMH=-9.*A1*A2x*A2x*BTG*BTG                     &
                       ,ANMM=-3.*A1*A2x*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)  &
                                      *BTG                              &
                       ,BDNH= 3.*A2x*(7.*A1+B2)*BTG                     &
                       ,BDNM= 6.*A1*A1                                  &
                       ,BEQH= A2x*B1*BTG+3.*A2x*(7.*A1+B2)*BTG          &
                       ,BEQM=-A1*B1*(1.-3.*C1)+6.*A1*A1                 &
                       ,BNMH=-A2x*BTG                                   &
                       ,BNMM=A1*(1.-3.*C1)                              &
                       ,BSHH=9.*A1*A2x*A2x*BTG                          &
                       ,BSHM=18.*A1*A1*A2x*C1                           &
                       ,BSMH=-3.*A1*A2x*(3.*A2x+3.*B2*C1+12.*A1*C1-B2)  &
                                      *BTG                              &
                       ,CESH=A2x                                        &
                       ,CESM=A1*(1.-3.*C1)                              &
                       ,CNV=EP_1*G/BTG                                  &
                       ,ELFCS=VKARMAN*BTG                               &
                       ,FZQ1=RTVISC*QVISC*ZQRZT                         &
                       ,FZQ2=RTVISC*QVISC*ZQRZT                         &
                       ,FZT1=RVISC *TVISC*SQPR                          &
                       ,FZT2=CZIV*GRRS*TVISC*SQPR                       &
                       ,FZU1=CZIV*VISC                                  &
                       ,PIHF=0.5*PI                                     &
                       ,RFAC=RIC/(FHNEU*RFC*RFC)                        &
                       ,RQVISC=1./QVISC                                 &
                       ,RRIC=1./RIC                                     &
                       ,USTFC=0.018/G                                   &
                       ,WNEW=1.-WOLD                                    &
                       ,WWST2=WWST*WWST
!
!-----------------------------------------------------------------------
!***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: AEQH=9.*A1*A2x*A2x*B1*BTG*BTG                   &
                            +9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG       &
                       ,AEQM=3.*A1*A2x*B1*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)&
                            *BTG+18.*A1*A1*A2x*(B2-3.*A2x)*BTG
!
!-----------------------------------------------------------------------
!***  FORBIDDEN TURBULENCE AREA
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: REQU=-AEQH/AEQM                                 &
                       ,EPSGH=1.E-9,EPSGM=REQU*EPSGH
!
!-----------------------------------------------------------------------
!***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: UBRYL=(18.*REQU*A1*A1*A2x*B2*C1*BTG             &
                               +9.*A1*A2x*A2x*B2*BTG*BTG)               &
                              /(REQU*ADNM+ADNH)                         &
                       ,UBRY=(1.+EPSRS)*UBRYL,UBRY3=3.*UBRY
!
      REAL,PARAMETER :: AUBH=27.*A1*A2x*A2x*B2*BTG*BTG-ADNH*UBRY3       &
                       ,AUBM=54.*A1*A1*A2x*B2*C1*BTG -ADNM*UBRY3        &
                       ,BUBH=(9.*A1*A2x+3.*A2x*B2)*BTG-BDNH*UBRY3       &
                       ,BUBM=18.*A1*A1*C1           -BDNM*UBRY3         &
                       ,CUBR=1.                     -     UBRY3         &
                       ,RCUBR=1./CUBR
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE TURBL(NTSD,DT,NPHS                                     &
                      ,NUM_WATER,NSOIL,SLDPTH,DZSOIL                    &
                      ,DSG2,SGML2,PDSG1,PSGML1,PDTOP,PT                 &
                      ,SM,CZEN,CZMEAN,SIGT4,RLWIN,RSWIN,RADOT           &
!- RLWIN/RSWIN - downward longwave/shortwave at the surface (also TOTLWDN/TOTSWDN in RADIATION)
                      ,PD,T,Q,CWM,F_ICE,F_RAIN,SR                       &
                      ,Q2,U,V,DUDT,DVDT                                 &
                      ,THS,TSFC,SST,PREC,SNO,WATER                      &
                      ,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG                    &
                      ,F_QV,F_QC,F_QR,F_QI,F_QS,F_QG                    &
                      ,FIS,Z0,Z0BASE,USTAR,PBLH,LPBL,XLEN_MIX,RMOL      &
                      ,EXCH_H,AKHS,AKMS,AKHS_OUT,AKMS_OUT               &
                      ,THZ0,QZ0,UZ0,VZ0,UZ0H,VZ0H,QS,MAVAIL             &
                      ,STC,SMC,CMC,SMSTAV,SMSTOT,SSROFF,BGROFF          &
                      ,IVGTYP,ISLTYP,VEGFRC,SHDMIN,SHDMAX,GRNFLX        &
                      ,SFCEXC,ACSNOW,ACSNOM,SNOPCX,SICE,TG,SOILTB       &
                      ,ALBASE,MXSNAL,ALBEDO,SH2O,SI,EPSR                &
                      ,U10,V10,TH10,Q10,TSHLTR,QSHLTR,PSHLTR            &
                      ,T2,QSG,QVG,QCG,SOILT1,TSNAV,SMFR3D,KEEPFR3DFLAG  &
                      ,TWBS,QWBS,SFCSHX,SFCLHX,SFCEVP                   &
                      ,POTEVP,POTFLX,SUBSHX                             &
                      ,APHTIM,ARDSW,ARDLW,ASRFC                         &
                      ,CROT,SROT                                        &
                      ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW      &
                      ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX      &
                      ,RSWOUT,RSWTOA,RLWTOA                             &
                      ,ASWIN,ASWOUT,ASWTOA,ALWIN,ALWOUT,ALWTOA          &
                      ,RTHBLTEN,RQVBLTEN                                & 
                      ,GWDFLG                                           &
                      ,PCPFLG,DDATA                                     & ! PRECIP ASSIM
                      ,UCMCALL                                          &
                      ,TURBULENCE,SFC_LAYER                             &
                      ,LAND_SURFACE,LONGWAVE                            &
                      ,MICROPHYSICS                                     &
                      ,IDS,IDE,JDS,JDE,LM                               &
                      ,IMS,IME,JMS,JME                                  &
                      ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    TURBL       TURBULENCE OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-04-19       
!     
! ABSTRACT:
!     TURBL DRIVES THE TURBULENCE SCHEMES
!     
! PROGRAM HISTORY LOG (with changes to called routines) :
!   95-03-15  JANJIC     - ORIGINATOR OF THE SUBROUTINES CALLED
!   BLACK & JANJIC       - ORIGINATORS OF THE DRIVER
!   95-03-28  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   96-03-29  BLACK      - ADDED EXTERNAL EDGE; REMOVED SCRCH COMMON
!   96-07-19  MESINGER   - ADDED Z0 EFFECTIVE
!   98-~??  TUCCILLO   - MODIFIED FOR CLASS VIII PARALLELISM
!   98-10-27  BLACK      - PARALLEL CHANGES INTO MOST RECENT CODE
!   02-01-10  JANJIC     - MOIST TURBULENCE (DRIVER, MIXLEN, VDIFH)
!   02-01-10  JANJIC     - VERT. DIF OF Q2 INCREASED (Grenier & Bretherton)
!   02-02-02  JANJIC     - NEW SFCDIF
!   02-04-19  BLACK      - ORIGINATOR OF THIS OUTER DRIVER FOR WRF
!   02-05-03  JANJIC     - REMOVAL OF SUPERSATURATION AT 2m AND 10m
!   04-11-18  BLACK      - THREADED
!   06-10-25  BLACK      - BUILT INTO NMMB PHYSICS COMPONENT
!     
! USAGE: CALL TURBL FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : 1
!$$$  
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE                             &
                           ,NPHS,NSOIL,NTSD,NUM_WATER                   &
                           ,UCMCALL
! 
      INTEGER,INTENT(IN) :: P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
! 
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ISLTYP,IVGTYP
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: LPBL
!
      REAL,INTENT(IN) :: DT,PDTOP,PT
!
      REAL,INTENT(INOUT) :: APHTIM,ARDSW,ARDLW,ASRFC
!
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ALBASE,MXSNAL
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CROT,SROT           &
                      ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW      &
                      ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CZEN,CZMEAN         &
                                                   ,FIS,PD              &
                                                   ,RLWIN,RLWTOA        &
                                                   ,RSWIN,RSWOUT,RSWTOA &
                                                   ,SHDMIN,SHDMAX       &
                                                   ,SICE,SIGT4          &
                                                   ,SST,TG,VEGFRC
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: SM,EPSR,SR         !Bandaid
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: GRNFLX,QWBS,RADOT  &
                                                    ,SFCEXC,SMSTAV      &
                                                    ,SOILTB,TWBS
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACSNOM,ACSNOW    &
                                                      ,AKHS,AKMS        &
                                                      ,ALBEDO           &
                                                      ,MAVAIL           &
                                                      ,BGROFF,CMC       &
                                                      ,PBLH,POTEVP      &
                                                      ,POTFLX,PREC      &
                                                      ,QCG,QS,QSG       &
                                                      ,QVG,QZ0,RMOL     &
                                                      ,SFCEVP           &
                                                      ,SFCLHX,SFCSHX    &
                                                      ,SI,SMSTOT        &
                                                      ,SNO,SNOPCX       &
                                                      ,SOILT1           &
                                                      ,SSROFF,SUBSHX    &
                                                      ,T2,THS,THZ0      &
                                                      ,TSFC,TSNAV       &
                                                      ,USTAR,UZ0,VZ0    &
                                                      ,Z0,Z0BASE
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: AKHS_OUT,AKMS_OUT  &
                                                    ,ALWIN,ALWOUT       &
                                                    ,ALWTOA,ASWIN       &
                                                    ,ASWOUT,ASWTOA      &
                                                    ,PSHLTR,Q10,QSHLTR  &
                                                    ,TH10,TSHLTR        &
                                                    ,U10,UZ0H,V10,VZ0H
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM         &
                                                           ,EXCH_H      &
                                                           ,Q,Q2        &
                                                           ,T,U,V
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,NUM_WATER),INTENT(INOUT) :: WATER
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) ::  F_ICE    &
                                                            ,F_RAIN

      REAL,DIMENSION(IMS:IME,1:LM+1,JMS:JME),INTENT(INOUT) ::  RQVBLTEN &
                                                              ,RTHBLTEN
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT) :: DUDT,DVDT     &
                                                         ,XLEN_MIX
!
      REAL,DIMENSION(NSOIL),INTENT(IN) :: DZSOIL,SLDPTH
!
      REAL,DIMENSION(IMS:IME,JMS:JME,NSOIL),INTENT(INOUT) :: SH2O,SMC,STC
!
      REAL,DIMENSION(IMS:IME,NSOIL,JMS:JME),INTENT(INOUT) :: KEEPFR3DFLAG &
                                                            ,SMFR3D
!
      CHARACTER(99),INTENT(IN) :: LAND_SURFACE,LONGWAVE,MICROPHYSICS    &
                                 ,SFC_LAYER,TURBULENCE
!
!  For precip assimilation:
!
      LOGICAL,INTENT(IN) :: GWDFLG,PCPFLG
      LOGICAL,INTENT(IN) :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DDATA
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
      INTEGER :: I,I_M,IDUMMY,IEND,ISFFLX,ISTR,J,K,KFLIP,KOUNT_ALL      &
                ,LENGTH_ROW,LLYR,N,SST_UPDATE
!
      INTEGER :: LW_PHYSICS,MP_PHYSICS,PBL_PHYSICS,SFCLAY_PHYSICS       &
                ,SURFACE_PHYSICS
!
      INTEGER,DIMENSION(NUM_TILES) :: I_START,I_END,J_START,J_END
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: KPBL,LOWLYR
!
      REAL :: TRESH=0.95
!
      REAL :: ALTITUDE,CWML,DQDT,DTDT,DTPHS,DX,DZHALF,EXNER             &
             ,FACTR,FACTRL,G_INV,PDSL,PLM,PLYR,PSFC,QI,QL,QOLD,QR,QW    &
             ,RATIOMX,RDTPHS,ROG,RWMSK,SDEPTH,SNO_FACTR                 &
             ,TL,TLMH,TLMH4,TNEW,TSFC2,U_FRAME,V_FRAME,WMSK,XLVRW
!
      REAL :: APES,CKLQ,FACTOR,FFS,PQ0X,Q2SAT,QFC1,QLOWX,RLIVWV,THBOT
!
      REAL :: RAPHTIM,RARDSW,RARDLW,RASRFC
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: BR,CHKLOWQ,CT,CWMLOW,ELFLX     &
                                        ,EMISS,EXNSFC,FACTRS,FLHC,FLQC  &
                                        ,GZ1OZ0,ONE,PSFC_OUT,PSIH,PSIM  &
                                        ,Q2X,QLOW,RAIN,RAINBL           &
                                        ,RLW_DN_SFC,RSW_NET_SFC         &
                                        ,RSW_DN_SFC                     &
                                        ,SFCEVPX,SFCZ,SNOW,SNOWC,SNOWH  &
                                        ,TH2X,THLOW,TLOW                &
                                        ,VGFRCK,WSPD,XLAND
!
      REAL,DIMENSION(IMS:IME,1:LM+1,JMS:JME) :: DELP                    &
                                               ,DUDT_GWD,DVDT_GWD       &
                                               ,DUDT_PHY,DVDT_PHY,DZ    &
                                               ,EL_MYJ,EXCH_H_PHY       &
                                               ,P8W,P_PHY,PI_PHY        &
                                               ,RQCBLTEN,RQIBLTEN       &
                                               ,RR                      &
                                               ,T_PHY,TH_PHY,TKE        &
                                               ,U_PHY,V_PHY,Z
!
      REAL,DIMENSION(IMS:IME,NSOIL,JMS:JME) :: SH2O_PHY                 &
                                              ,SMC_PHY                  &
                                              ,STC_PHY
!
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: WATER_TRANS
!
      LOGICAL :: E_BDY,N_BDY,S_BDY,W_BDY,WARM_RAIN
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
!***  TRANSLATE THE PACKAGE OPTIONS IN THE CONFIG FILE NEEDED BY
!***  THE TURBULENCE TO THEIR ANALOGS IN THE WRF REGISTRY SO THAT 
!***  THE WRF SURFACE AND PBL DRIVERS REMAIN UNTOUCHED.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(TURBULENCE))
        CASE ('myj')
          PBL_PHYSICS=2
        CASE ('ysu')
          PBL_PHYSICS=1
        CASE DEFAULT
          WRITE(0,*)' User selected TURBULENCE=',TRIM(TURBULENCE)
          WRITE(0,*)' Improper selection of Turbulence scheme in TURBL'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
      END SELECT
!
      SELECT CASE (TRIM(SFC_LAYER))
        CASE ('myj')
          SFCLAY_PHYSICS=2
        CASE ('mm5')
          SFCLAY_PHYSICS=1
        CASE DEFAULT
          WRITE(0,*)' User selected SFC_LAYER=',TRIM(SFC_LAYER)
          WRITE(0,*)' Improper selection of Surface Layer scheme in TURBL'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
      END SELECT
!
      SELECT CASE (TRIM(LAND_SURFACE))
        CASE ('nmm')
          SURFACE_PHYSICS=99
        CASE ('noah')
          SURFACE_PHYSICS=2
        CASE DEFAULT
          WRITE(0,*)' User selected LAND_SURFACE=',TRIM(LAND_SURFACE)
          WRITE(0,*)' Improper selection of Land Surface scheme in TURBL'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
      END SELECT
!
      SELECT CASE (TRIM(LONGWAVE))
        CASE ('gfdl')
          LW_PHYSICS=99
        CASE ('rrtm')
          LW_PHYSICS=1
        CASE DEFAULT
          WRITE(0,*)' User selected LONGWAVE=',TRIM(LONGWAVE)
          WRITE(0,*)' Improper selection of LW scheme in TURBL'
!!!       CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          CALL NMMB_FINALIZE
      END SELECT
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
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,k,i)
!.......................................................................
      DO J=JMS,JME
      DO K=1,LM+1
      DO I=IMS,IME
	U_PHY(I,K,J)=0.
        V_PHY(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  THE WRF PHYSICS DRIVERS EXPECT THE VERTICAL INDEX TO BE RECKONED
!***  SO AS TO INCREASE UPWARD FROM THE GROUND THEREFORE DATA SENT TO
!***  AND RECEIVED FROM THESE DRIVERS MUST BE FLIPPED.
!-----------------------------------------------------------------------
!
      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      G_INV=1./G
      ROG=R_D*G_INV
      FACTOR=-XLV*RHOWATER/DTPHS
!
      U_FRAME=0.
      V_FRAME=0.
!
      IDUMMY=0
      ISFFLX=1
      DX=0.
      SST_UPDATE=0
!
      DO J=JMS,JME
      DO I=IMS,IME
        UZ0H(I,J)=0.
        VZ0H(I,J)=0.
        ONE(I,J)=1.
        RMOL(I,J)=0.     !Reciprocal of Monin-Obukhov length
        SFCEVPX(I,J)=0.  !Dummy for accumulated latent energy, not flux
      ENDDO
      ENDDO
!
      W_BDY=(ITS==IDS)
      E_BDY=(ITE==IDE)
      S_BDY=(JTS==JDS)
      N_BDY=(JTE==JDE)
!
      IF(S_BDY)THEN
        DO K=1,LM
        DO I=ITS,ITE
          U_PHY(I,K,JTS)=0.
          V_PHY(I,K,JTS)=0.
        ENDDO
        ENDDO
      ENDIF 
!
      IF(N_BDY)THEN
        DO K=1,LM
        DO I=ITS,ITE
          U_PHY(I,K,JTE)=0.
          V_PHY(I,K,JTE)=0.
        ENDDO
        ENDDO
      ENDIF
!
      IF(W_BDY)THEN
        DO J=JTS,JTE
        DO K=1,LM
          U_PHY(ITS,K,J)=0.
          V_PHY(ITS,K,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
      IF(E_BDY)THEN
        DO J=JTS,JTE
        DO K=1,LM
          U_PHY(ITE,K,J)=0.
          V_PHY(ITE,K,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
      IF(SURFACE_PHYSICS==99)THEN
        SNO_FACTR=1.
      ELSE
        SNO_FACTR=0.001
      ENDIF
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
      DO J=JTS,JTE
      DO I=ITS,ITE
        LOWLYR(I,J)=1
        VGFRCK(I,J)=100.*VEGFRC(I,J)
        SNOW(I,J)=SNO(I,J)
        SNOWH(I,J)=SI(I,J)*SNO_FACTR
        XLAND(I,J)=SM(I,J)+1.
        T2(I,J)=TSFC(I,J)
        EMISS(I,J)=EPSR(I,J)
      ENDDO
      ENDDO
!
      IF(NTSD==0)THEN
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
        DO J=JTS,JTE
        DO I=ITS,ITE
          Z0BASE(I,J)=Z0(I,J)
          IF(SM(I,J)>0.5.AND.SICE(I,J)>0.5)THEN  !Bandaid
            SM(I,J)=0.        
          ENDIF              
        ENDDO
        ENDDO
      ENDIF
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
      DO K=1,LM
      DO J=JTS,JTE
      DO I=ITS,ITE
        Q2(I,J,K)=MAX(Q2(I,J,K),EPSQ2)
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
      DO J=JTS,JTE
      DO K=1,LM+1
      DO I=ITS,ITE
        Z(I,K,J)=0.
        DZ(I,K,J)=0.
        EXCH_H_PHY(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
!***  PREPARE NEEDED ARRAYS
!
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                     &
!$omp private(j,i,pdsl,psfc,plm,tlmh,factrl,k,kflip,plyr,ql,tl,cwml   &
!$omp        ,exner),SCHEDULE(dynamic)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)         
        PSFC=PD(I,J)+PDTOP+PT
        P8W(I,1,J)=PSFC
        LOWLYR(I,J)=1
        EXNSFC(I,J)=(1.E5/PSFC)**CAPA
        THS(I,J)=(SST(I,J)*EXNSFC(I,J))*SM(I,J)+THS(I,J)*(1.-SM(I,J))
        TSFC(I,J)=THS(I,J)/EXNSFC(I,J)
        SFCZ(I,J)=FIS(I,J)*G_INV
!YL
!       RAIN(I,J)=PREC(I,J)*RHOWATER
        IF(PCPFLG.AND.DDATA(I,J)<100.)THEN
          RAIN(I,J)=DDATA(I,J)*RHOWATER
        ELSE
          RAIN(I,J)=PREC(I,J)*RHOWATER
        ENDIF
!YL
        RAINBL(I,J)=0.
        IF(SNO(I,J)>0.)SNOWC(I,J)=1.
        PLM=PT+PDTOP+SGML2(LM)*PDSL
        TH2X(I,J)=T(I,J,LM)*(1.E5/PLM)**CAPA
        Q2X(I,J)=Q(I,J,LM)
!
!-----------------------------------------------------------------------
!*** LONG AND SHORTWAVE FLUX AT GROUND SURFACE
!-----------------------------------------------------------------------
!
        IF(CZMEAN(I,J)>0.)THEN
          FACTRS(I,J)=CZEN(I,J)/CZMEAN(I,J)
        ELSE
          FACTRS(I,J)=0.
        ENDIF
!
        IF(SIGT4(I,J)>0.)THEN
          TLMH=T(I,J,LM)
          FACTRL=STBOLT*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
        ELSE
          FACTRL=0.
        ENDIF
!     
!- RLWIN/RSWIN - downward longwave/shortwave at the surface
!
        RLW_DN_SFC(I,J)=RLWIN(I,J)*FACTRL
        RSW_NET_SFC(I,J)=(RSWIN(I,J)-RSWOUT(I,J))*FACTRS(I,J)
!
!- Instantaneous downward solar for NMM_LSM
!
        RSW_DN_SFC(I,J)=RSWIN(I,J)*FACTRS(I,J)
!
!-----------------------------------------------------------------------
!***  FILL THE ARRAYS FOR CALLING THE INNER DRIVER.
!-----------------------------------------------------------------------
!
        Z(I,1,J)=SFCZ(I,J)
!
!-----------------------------------------------------------------------
!***  FILL VERTICAL WORKING ARRAYS.
!***  THE FUNDAMENTAL NMMB ARRAYS NEED TO FLIP THEIR DATA
!***  BEFORE SENDING IT TO THE WRF PHYSICS DRIVERS.
!***  USE THE FLIPPED INDEX IN THE NMMB ARRAYS.
!-----------------------------------------------------------------------
!
        DO K=1,LM
          KFLIP=LM+1-K
!
          IF(DSG2(KFLIP)<1.E-10)THEN
            PLYR=PSGML1(KFLIP)
          ELSE
            PLYR=PT+PDTOP+SGML2(KFLIP)*PDSL
          ENDIF
!
          QL=MAX(Q(I,J,KFLIP),EPSQ)
          TL=T(I,J,KFLIP)
          CWML=CWM(I,J,KFLIP)
!
          RR(I,K,J)=PLYR/(R_D*TL)
          T_PHY(I,K,J)=TL
          EXNER=(1.E5/PLYR)**CAPA
          PI_PHY(I,K,J)=1./EXNER
          TH_PHY(I,K,J)=TL*EXNER
          P8W(I,K+1,J)=P8W(I,K,J)-PDSG1(KFLIP)-DSG2(KFLIP)*PDSL
          P_PHY(I,K,J)=PLYR
          TKE(I,K,J)=0.5*Q2(I,J,KFLIP)
!
          RTHBLTEN(I,K,J)=0.
          RQVBLTEN(I,K,J)=0.
          RQCBLTEN(I,K,J)=0.
          RQIBLTEN(I,K,J)=0.
!
          DZ(I,K,J)=T_PHY(I,K,J)*(P608*QL+1.)*R_D                       &
                    *(P8W(I,K,J)-P8W(I,K+1,J))                          &
                    /(PLYR*G)
          Z(I,K+1,J)=Z(I,K,J)+DZ(I,K,J)
!
          DELP(I,K,J)=P8W(I,K,J)-P8W(I,K+1,J)
!
        ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,qlowx)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
        TWBS(I,J)=0.
        QWBS(I,J)=0.
        THLOW(I,J)=TH_PHY(I,1,J)
        TLOW(I,J)=T_PHY(I,1,J)
        QLOW(I,J)=MAX(Q(I,J,LM),EPSQ)
        QLOWX=QLOW(I,J)/(1.-QLOW(I,J))
        QLOW(I,J)=QLOWX/(1.+QLOWX)
        CWMLOW(I,J)=CWM(I,J,LM)
        PBLH(I,J)=MAX(PBLH(I,J),0.)
        PBLH(I,J)=MIN(PBLH(I,J),Z(I,LM,J))
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  COMPUTE VELOCITY COMPONENTS AT MASS POINTS
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,kflip)
!.......................................................................
      DO J=JTS_B1,JTE_B1
        DO I=ITS_B1,ITE_B1
          UZ0H(I,J)=(UZ0(I,J  )+UZ0(I-1,J  )                            &
                    +UZ0(I,J-1)+UZ0(I-1,J-1))*0.25
          VZ0H(I,J)=(VZ0(I,J  )+VZ0(I-1,J  )                            &
                    +VZ0(I,J-1)+VZ0(I-1,J-1))*0.25
        ENDDO
        DO K=1,LM
          KFLIP=LM+1-K
          DO I=ITS_B1,ITE_B1
            U_PHY(I,K,J)=(U(I,J  ,KFLIP)+U(I-1,J  ,KFLIP)               &
                         +U(I,J-1,KFLIP)+U(I-1,J-1,KFLIP))*0.25
            V_PHY(I,K,J)=(V(I,J  ,KFLIP)+V(I-1,J  ,KFLIP)               &
                         +V(I,J-1,KFLIP)+V(I-1,J-1,KFLIP))*0.25
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
!jaa      DO J=JTS_B1,JTE_B1
!jaa        DO I=ITS_B1,ITE_B1
!jaa          UZ0H(I,J)=(UZ0(I,J  )+UZ0(I-1,J  )                            &
!jaa                    +UZ0(I,J-1)+UZ0(I-1,J-1))*0.25
!jaa          VZ0H(I,J)=(VZ0(I,J  )+VZ0(I-1,J  )                            &
!jaa                    +VZ0(I,J-1)+VZ0(I-1,J-1))*0.25
!jaa        ENDDO
!jaa      ENDDO
!
!-----------------------------------------------------------------------
!***  MOISTURE AVAILABILITY
!-----------------------------------------------------------------------
!
      IF(TRIM(LAND_SURFACE)=='nmm')THEN
        DO J=JTS,JTE 
        DO I=ITS,ITE
          ONE(I,J)=1.
        ENDDO
        ENDDO
      ELSE
        DO J=JTS,JTE 
        DO I=ITS,ITE
          ONE(I,J)=MAVAIL(I,J)
        ENDDO
        ENDDO
      ENDIF 
!
!-----------------------------------------------------------------------
!***  TRANSPOSE THE WATER ARRAY (IJK) FOR THE PHYSICS (IKJ).
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,kflip)
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
!***  TRANSPOSE 3-D SOIL ARRAYS (IJK) FOR PHYSICS (IKJ).
!-----------------------------------------------------------------------
!
      DO N=1,NSOIL
        DO J=JMS,JME
        DO I=IMS,IME
          SH2O_PHY(I,N,J)=SH2O(I,J,N)
          SMC_PHY(I,N,J)=SMC(I,J,N)
          STC_PHY(I,N,J)=STC(I,J,N)
        ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
!***  CALL SURFACE LAYER AND LAND SURFACE PHYSICS
!
!-----------------------------------------------------------------------
!
!!!   CALL SET_TILES(GRID,IDS,IDE-1,JDS+1,JDE-1,ITS,ITE,JTS,JTE)
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
      CALL SURFACE_DRIVER(                                              &
                 ACSNOM=ACSNOM,ACSNOW=ACSNOW,AKHS=AKHS,AKMS=AKMS        &
                ,ALBEDO=ALBEDO,BR=BR,CANWAT=CMC,CHKLOWQ=CHKLOWQ         &
                ,DT=DT,DX=DX,DZ8W=DZ,DZS=DZSOIL,GLW=RLW_DN_SFC          &
                ,GRDFLX=GRNFLX,GSW=RSW_NET_SFC,SWDOWN=RSW_DN_SFC        &
                ,GZ1OZ0=GZ1OZ0,HFX=TWBS                                 &
                ,HT=SFCZ,IFSNOW=IDUMMY,ISFFLX=ISFFLX,ISLTYP=ISLTYP      &
                ,ITIMESTEP=NTSD,IVGTYP=IVGTYP,LOWLYR=LOWLYR             &
                ,MAVAIL=ONE,RMOL=RMOL,NUM_SOIL_LAYERS=NSOIL,P8W=P8W &
                ,PBLH=PBLH,PI_PHY=PI_PHY,PSHLTR=PSHLTR,PSIH=PSIH        &
                ,PSIM=PSIM,P_PHY=P_PHY,Q10=Q10,Q2=Q2X,QFX=QWBS,QSFC=QS  &
                ,QSHLTR=QSHLTR,QZ0=QZ0,RAINCV=RAIN                      &
                ,RHO=RR,SFCEVP=SFCEVPX,SFCEXC=SFCEXC,SFCRUNOFF=SSROFF   &
                ,SMOIS=SMC_PHY,SMSTAV=SMSTAV,SMSTOT=SMSTOT,SNOALB=MXSNAL &
                ,SNOW=SNOW,SNOWC=SNOWC,SNOWH=SNOWH,STEPBL=NPHS          &
                ,SST=SST,SST_UPDATE=SST_UPDATE                          &
                ,TH10=TH10,TH2=TH2X,T2=T2,THZ0=THZ0,TH_PHY=TH_PHY       &
                ,TMN=TG,TSHLTR=TSHLTR,TSK=TSFC,TSLB=STC_PHY,T_PHY=T_PHY &
                ,U10=U10,UDRUNOFF=BGROFF,UST=USTAR,UZ0=UZ0H             &
                ,U_FRAME=U_FRAME,U_PHY=U_PHY,V10=V10,VEGFRA=VGFRCK      &
                ,VZ0=VZ0H,V_FRAME=V_FRAME,V_PHY=V_PHY                   &
                ,WARM_RAIN=WARM_RAIN,WSPD=WSPD,XICE=SICE                &
                ,XLAND=XLAND,Z=Z,ZNT=Z0,ZS=SLDPTH,CT=CT,TKE_MYJ=TKE     &
                ,ALBBCK=ALBASE,LH=ELFLX,SH2O=SH2O_PHY,SHDMAX=SHDMAX     &
                ,SHDMIN=SHDMIN,Z0=Z0BASE,FLQC=FLQC,FLHC=FLHC            &
                ,PSFC=PSFC_OUT,EMISS=EPSR                               &
                ,SF_SFCLAY_PHYSICS=SFCLAY_PHYSICS                       &
                ,SF_SURFACE_PHYSICS=SURFACE_PHYSICS                     &
                ,RA_LW_PHYSICS=LW_PHYSICS                               &
                ,UCMCALL=UCMCALL                                        &
                ,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=1,KDE=LM+1         &
                ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=1,KME=LM+1         &
                ,I_START=I_START,I_END=I_END                            &
                ,J_START=J_START,J_END=J_END                            &
                ,KTS=1,KTE=LM,NUM_TILES=NUM_TILES                       &
           ! Optional args
                ,QV_CURR=WATER_TRANS(IMS,1,JMS,P_QV),F_QV=F_QV          &
                ,QC_CURR=WATER_TRANS(IMS,1,JMS,P_QC),F_QC=F_QC          &
                ,QR_CURR=WATER_TRANS(IMS,1,JMS,P_QR),F_QR=F_QR          &
                ,QI_CURR=WATER_TRANS(IMS,1,JMS,P_QI),F_QI=F_QI          &
                ,QS_CURR=WATER_TRANS(IMS,1,JMS,P_QS),F_QS=F_QS          & 
                ,QG_CURR=WATER_TRANS(IMS,1,JMS,P_QG),F_QG=F_QG          &
                ,RAINBL=RAINBL                                          &
! for RUCLSM
                ,QSG=QSG, QVG=QVG, QCG=QCG, SOILT1=SOILT1               &
                ,TSNAV=TSNAV, SMFR3D=SMFR3D, KEEPFR3DFLAG=KEEPFR3DFLAG  &
                ,POTEVP=POTEVP,SNOPCX=SNOPCX,SOILTB=SOILTB,SR=SR)
!
!-----------------------------------------------------------------------
!
!***  CALL FREE ATMOSPHERE TURBULENCE
!
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
      DO J=JMS,JME
      DO K=1,LM+1
      DO I=IMS,IME
        DUDT_PHY(I,K,J)=0.
        DVDT_PHY(I,K,J)=0.
        DUDT_GWD(I,K,J)=0.
        DVDT_GWD(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!***  THE SURFACE EXCHANGE COEFFICIENTS AKHS AND AKMS ARE ACTUALLY
!***  MULTIPLIED BY HALF THE DEPTH OF THE LOWEST LAYER.  WE MUST RETAIN
!***  THOSE VALUES FOR THE NEXT TIMESTEP SO USE AUXILLIARY ARRAYS FOR
!***  THE OUTPUT.
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(dzhalf,i,j)
      DO J=JTS,JTE
      DO I=ITS,ITE
        DZHALF=0.5*DZ(I,1,J)
        AKHS_OUT(I,J)=AKHS(I,J)*DZHALF
        AKMS_OUT(I,J)=AKMS(I,J)*DZHALF
      ENDDO
      ENDDO

!
      CALL PBL_DRIVER(                                                &
                      ITIMESTEP=NTSD,DT=DT                            &
                     ,U_FRAME=U_FRAME,V_FRAME=V_FRAME                 &
                     ,RUBLTEN=DUDT_PHY,RVBLTEN=DVDT_PHY               &
                     ,RTHBLTEN=RTHBLTEN                               &
                     ,RQVBLTEN=RQVBLTEN,RQCBLTEN=RQCBLTEN             &
                     ,RQSBLTEN=RQIBLTEN                               &
                     ,TSK=TSFC,XLAND=XLAND,ZNT=Z0,HT=SFCZ             &
                     ,UST=USTAR, PBLH=PBLH                            &
                     ,HFX=TWBS,QFX=QWBS,  GRDFLX=GRNFLX               &
                     ,U_PHY=U_PHY,V_PHY=V_PHY,TH_PHY=TH_PHY,RHO=RR    &
                     ,P_PHY=P_PHY,PI_PHY=PI_PHY,P8W=P8W,T_PHY=T_PHY   &
                     ,DZ8W=DZ,Z=Z,TKE_MYJ=TKE,EL_MYJ=EL_MYJ           &
                     ,EXCH_H=EXCH_H_PHY,AKHS=AKHS,AKMS=AKMS           &
                     ,THZ0=THZ0,QZ0=QZ0,UZ0=UZ0H,VZ0=VZ0H             &
                     ,QSFC=QS,LOWLYR=LOWLYR                           &
                     ,PSIM=PSIM,PSIH=PSIH,GZ1OZ0=GZ1OZ0               &
                     ,WSPD=WSPD,BR=BR,CHKLOWQ=CHKLOWQ                 &
                     ,DX=DX,STEPBL=NPHS,WARM_RAIN=WARM_RAIN           &
                     ,KPBL=KPBL,CT=CT,LH=ELFLX,SNOW=SNOW,XICE=SICE    &
                     ,BL_PBL_PHYSICS=PBL_PHYSICS                      &
                     ,RA_LW_PHYSICS=LW_PHYSICS                        &
                     ,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=1,KDE=LM+1  &
                     ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=1,KME=LM+1  &
                     ,I_START=I_START,I_END=I_END                     &
                     ,J_START=J_START,J_END=J_END                     &
                     ,KTS=1,KTE=LM,NUM_TILES=NUM_TILES                &
                ! Optional args
                     ,QV_CURR=WATER_TRANS(IMS,1,JMS,P_QV),F_QV=F_QV   &
                     ,QC_CURR=WATER_TRANS(IMS,1,JMS,P_QC),F_QC=F_QC   &
                     ,QR_CURR=WATER_TRANS(IMS,1,JMS,P_QR),F_QR=F_QR   &
                     ,QI_CURR=WATER_TRANS(IMS,1,JMS,P_QI),F_QI=F_QI   &
                     ,QS_CURR=WATER_TRANS(IMS,1,JMS,P_QS),F_QS=F_QS   &
                     ,QG_CURR=WATER_TRANS(IMS,1,JMS,P_QG),F_QG=F_QG  )
!
!***  NOTE THAT THE EXCHANGE COEFFICIENTS FOR HEAT EXCH_H COMING OUT OF
!***  PBL_DRIVER ARE DEFINED AT THE TOPS OF THE LAYERS KTS TO KTE-1
!***  IF MODULE_BL_MYJPBL WAS INVOKED.
!
!-----------------------------------------------------------------------
!***  UNCOMPUTED LOCATIONS MUST BE FILLED IN FOR THE POST-PROCESSOR
!-----------------------------------------------------------------------
!
!***  EASTERN GLOBAL BOUNDARY
!
      IF(E_BDY)THEN
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
        DO J=JDS,JDE
        IF(J>=JTS.AND.J<=JTE)THEN
          TH10(IDE,J)=TH10(IDE-1,J)
          Q10(IDE,J)=Q10(IDE-1,J)
          U10(IDE,J)=U10(IDE-1,J)
          V10(IDE,J)=V10(IDE-1,J)
          TSHLTR(IDE,J)=TSHLTR(IDE-1,J)
          QSHLTR(IDE,J)=QSHLTR(IDE-1,J)
        ENDIF
        ENDDO
      ENDIF
!
!***  SOUTHERN GLOBAL BOUNDARY
!
      IF(S_BDY)THEN
        DO I=IDS,IDE
          IF(I>=ITS.AND.I<=ITE)THEN
            TH10(I,JDS)=TH10(I,JDS+1)
            Q10(I,JDS)=Q10(I,JDS+1)
            U10(I,JDS)=U10(I,JDS+1)
            V10(I,JDS)=V10(I,JDS+1)
            TSHLTR(I,JDS)=TSHLTR(I,JDS+1)
            QSHLTR(I,JDS)=QSHLTR(I,JDS+1)
          ENDIF
        ENDDO
      ENDIF
!
!***  NORTHERN GLOBAL BOUNDARY
!
      IF(N_BDY)THEN
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
        DO I=IDS,JDE
          IF(I>=ITS.AND.I<=ITE)THEN
            TH10(I,JDE)=TH10(I,JDE-1)
            Q10(I,JDE)=Q10(I,JDE-1)
            U10(I,JDE)=U10(I,JDE-1)
            V10(I,JDE)=V10(I,JDE-1)
            TSHLTR(I,JDE)=TSHLTR(I,JDE-1)
            QSHLTR(I,JDE)=QSHLTR(I,JDE-1)
          ENDIF
        ENDDO
      ENDIF
!
      IF(TRIM(SFC_LAYER)/='myj')THEN
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
        DO J=JTS_B1,JTE_B1
        DO I=ITS_B1,ITE_B1
!         TSHLTR(I,J)=TSHLTR(I,J)*(1.E5/PSHLTR(I,J))**RCP
          IF(TSHLTR(I,J)<200..OR.TSHLTR(I,J)>350.)THEN
            WRITE(0,*)'Troublesome TSHLTR...I,J,TSHLTR,PSHLTR: ',       &
               I,J,TSHLTR(I,J),PSHLTR(I,J)
          ENDIF
	ENDDO
	ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  COMPUTE MODEL LAYER CONTAINING THE TOP OF THE BOUNDARY LAYER
!-----------------------------------------------------------------------
!
      IF(TRIM(TURBULENCE)/='myj')THEN
        LENGTH_ROW=ITE_B1-ITS_B1+1
        DO J=JTS_B1,JTE_B1
        DO I=ITS_B1,ITE_B1
          KPBL(I,J)=-1000
        ENDDO
        ENDDO
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,kount_all,k,i,altitude)
!.......................................................................
        DO J=JTS_B1,JTE_B1
          KOUNT_ALL=0
          find_kpbl : DO K=1,LM
          DO I=ITS_B1,ITE_B1
            ALTITUDE=Z(I,K+1,J)-SFCZ(I,J)
            IF(PBLH(I,J)<=ALTITUDE.AND.KPBL(I,J)<0)THEN
              KPBL(I,J)=K
              KOUNT_ALL=KOUNT_ALL+1
            ENDIF
            IF(KOUNT_ALL==LENGTH_ROW)EXIT find_kpbl
          ENDDO
          ENDDO find_kpbl
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      ENDIF
!
      IF(SURFACE_PHYSICS==99)THEN
        SNO_FACTR=1.
      ELSE
        SNO_FACTR=1000.
      ENDIF
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j)
      DO J=JTS_B1,JTE_B1
      DO I=ITS,ITE
        SNO(I,J)=SNOW(I,J)
        SI(I,J)=SNOWH(I,J)*SNO_FACTR
        LPBL(I,J)=LM+1-KPBL(I,J)      !<--- Layer top of PBL counting downward
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  DIAGNOSTIC RADIATION ACCUMULATION
!-----------------------------------------------------------------------
!
      RARDSW=1./ARDSW
      RARDLW=1./ARDLW
!
!.......................................................................
!$omp parallel do private(j,i,tsfc2)
!.......................................................................
      DO J=JTS_B1,JTE_B1
      DO I=ITS,ITE
        ASWIN (I,J)=ASWIN (I,J)+RSWIN(I,J)*FACTRS(I,J)*RARDSW
        ASWOUT(I,J)=ASWOUT(I,J)-RSWOUT(I,J)*FACTRS(I,J)*RARDSW
        ASWTOA(I,J)=ASWTOA(I,J)+RSWTOA(I,J)*FACTRS(I,J)*RARDSW
        ALWIN (I,J)=ALWIN (I,J)+RLW_DN_SFC(I,J)*RARDLW
        ALWOUT(I,J)=ALWOUT(I,J)-RADOT (I,J)*RARDLW
        ALWTOA(I,J)=ALWTOA(I,J)+RLWTOA(I,J)*RARDLW
!
        TSFC2=TSFC(I,J)*TSFC(I,J)
        RADOT(I,J)=EPSR(I,J)*STBOLT*TSFC2*TSFC2
        THS(I,J)=TSFC(I,J)*EXNSFC(I,J)
        PREC(I,J)=0.
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!=======================================================================
!===  Begin gravity wave drag (GWD) and mountain blocking (MB)  ========
!=======================================================================
!
      IF(GWDFLG) THEN
        CALL GWD_DRIVER(U_PHY,V_PHY,T_PHY           &
                       ,WATER_TRANS(IMS,1,JMS,P_QV) &
                       ,Z,DELP,P8W,P_PHY,PI_PHY     &
                       ,KPBL                        &
                       ,HSTDV,HCNVX,HASYW,HASYS     &
                       ,HASYSW,HASYNW,HLENW         &
                       ,HLENS,HLENSW,HLENNW         &
                       ,HANGL,HANIS,HSLOP,HZMAX     &
                       ,CROT,SROT                   &
                       ,DUDT_GWD,DVDT_GWD           &
                       ,IDS,IDE,JDS,JDE,1,LM        &
                       ,IMS,IME,JMS,JME,1,LM        &
                       ,ITS,ITE,JTS,JTE,1,LM )
      ENDIF
!
!=======================================================================
!=====  End gravity wave drag (GWD) and mountain blocking (MB)  ========
!=======================================================================
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD, AND TKE.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,k,kflip,i,dtdt,dqdt,qold,ratiomx,qw,qi,qr,i_m)
!.......................................................................
      DO J=JTS_B1,JTE_B1
        DO K=1,LM
          KFLIP=LM+1-K
          DO I=ITS_B1,ITE_B1
            DTDT=RTHBLTEN(I,K,J)*PI_PHY(I,K,J)
            DQDT=RQVBLTEN(I,K,J)         !Mixing ratio tendency
            T(I,J,KFLIP)=T(I,J,KFLIP)+DTDT*DTPHS
            QOLD=Q(I,J,KFLIP)
            RATIOMX=QOLD/(1.-QOLD)+DQDT*DTPHS
            Q(I,J,KFLIP)=RATIOMX/(1.+RATIOMX)
!           Q(I,J,KFLIP)=MAX(Q(I,J,KFLIP),EPSQ)
            QW=MAX(0.,WATER_TRANS(I,K,J,P_QC)+RQCBLTEN(I,K,J)*DTPHS )
!
            IF(TRIM(MICROPHYSICS)=='fer')THEN
              QI=MAX(0.,WATER_TRANS(I,K,J,P_QS)+RQIBLTEN(I,K,J)*DTPHS )
            ELSE
              QI=MAX(0.,WATER_TRANS(I,K,J,P_QI)+RQIBLTEN(I,K,J)*DTPHS )
            ENDIF
!
            QR=MAX(0.,WATER_TRANS(I,K,J,P_QR) )
!           CWM(I,K,J)=QW+QI+QR
            CWM(I,J,KFLIP)=0.    ! <---- BEWARE of the this line with GFS physics
!
            DO I_M=2,NUM_WATER
!
              IF(I_M/=P_QV)THEN
                CWM(I,J,KFLIP)=CWM(I,J,KFLIP)+WATER_TRANS(I,K,J,I_M)
              ENDIF
!
              IF(I_M==P_QV)THEN
                WATER_TRANS(I,K,J,P_QV)=MAX(EPSQ,(WATER_TRANS(I,K,J,P_QV)+RQVBLTEN(I,K,J)*DTPHS))
              ELSEIF(I_M==P_QC)THEN
                CWM(I,J,KFLIP)=MAX(0.,(CWM(I,J,KFLIP)+RQCBLTEN(I,K,J)*DTPHS))
              ELSEIF(I_M==P_QI)THEN
                CWM(I,J,KFLIP)=MAX(0.,(CWM(I,J,KFLIP)+RQIBLTEN(I,K,J)*DTPHS))
              ENDIF
!
            ENDDO
!
            WATER_TRANS(I,K,J,P_QC)=QW
            WATER_TRANS(I,K,J,P_QR)=QR
!
            IF(TRIM(MICROPHYSICS)=='fer')THEN
              WATER_TRANS(I,K,J,P_QS)=QI
              IF(QI<=EPSQ)THEN  
                F_ICE(I,J,KFLIP)=0.
              ELSE
                F_ICE(I,J,KFLIP)=MAX(0.,MIN(1.,QI/CWM(I,J,KFLIP)))
              ENDIF
!
              IF(QR<=EPSQ)THEN
                F_RAIN(I,J,KFLIP)=0.
              ELSE
                F_RAIN(I,J,KFLIP)=QR/(QW+QR)
              ENDIF
            ELSE
              WATER_TRANS(I,K,J,P_QI)=QI
            ENDIF
!
            Q2(I,J,KFLIP)=2.*TKE(I,K,J)
          ENDDO
        ENDDO
!
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  TRANSPOSE THE WATER_TRANS ARRAY BACK TO THE PROGNOSTIC WATER ARRAY.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(n,j,k,kflip,i)
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
      DEALLOCATE(WATER_TRANS)
!
!-----------------------------------------------------------------------
!***  TRANSPOSE THE TRANSPOSED SOIL ARRAYS BACK TO IJK.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,k,kflip,i)
!.......................................................................
        DO J=JMS,JME
        DO K=1,NSOIL
        KFLIP=NSOIL+1-K
        DO I=IMS,IME
          SH2O(I,J,K)=SH2O_PHY(I,K,J)
          SMC (I,J,K)=SMC_PHY(I,K,J)
          STC (I,J,K)=STC_PHY(I,K,J)
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  TRANSFER THE WIND TENDENCIES.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,k,kflip,i)
!.......................................................................
      DO J=JMS,JME
        DO K=1,LM
          KFLIP=LM+1-K
          DO I=IMS,IME
            DUDT(I,J,KFLIP)=DUDT_PHY(I,K,J)+DUDT_GWD(I,K,J)
            DVDT(I,J,KFLIP)=DVDT_PHY(I,K,J)+DVDT_GWD(I,K,J)
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  TRANSFER THE MIXING LENGTH.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,k,kflip,i)
!.......................................................................
      DO J=JTS,JTE
        DO K=1,LM
          KFLIP=LM+1-K
          DO I=ITS,ITE
            XLEN_MIX(I,J,KFLIP)=EL_MYJ(I,K,J)
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***
!***  SAVE SURFACE-RELATED FIELDS.
!***
!-----------------------------------------------------------------------
!
      RASRFC=1./ASRFC
      RAPHTIM=1./APHTIM
!
!.......................................................................
!$omp parallel do private(j,i,xlvrw)
!.......................................................................
      DO J=JTS_B1,JTE_B1
      DO I=ITS_B1,ITE_B1
!
!-----------------------------------------------------------------------
!***  INSTANTANEOUS SENSIBLE AND LATENT HEAT FLUX
!-----------------------------------------------------------------------
!
        TWBS(I,J)=-TWBS(I,J)
!
        IF ( SM(I,J) + SICE(I,J) <= 0.5 ) THEN
          QWBS(I,J)= QWBS(I,J)*XLV*CHKLOWQ(I,J) ! land
        ELSE
          QWBS(I,J)=-QWBS(I,J)    *CHKLOWQ(I,J) ! ocean
        ENDIF
!
!-----------------------------------------------------------------------
!***  ACCUMULATED QUANTITIES.
!***  IN OPNL LSM, SFCEVP APPEARS TO BE IN UNITS OF
!***  METERS OF LIQUID WATER.  IT IS COMING FROM
!***  WRF MODULE AS KG/M**2.
!-----------------------------------------------------------------------
!
        SFCSHX(I,J)=SFCSHX(I,J)+TWBS(I,J)*RASRFC
        SFCLHX(I,J)=SFCLHX(I,J)+QWBS(I,J)*RASRFC
        XLVRW=DTPHS/RHOWATER
        SFCEVP(I,J)=SFCEVP(I,J)-QWBS(I,J)*XLVRW*RASRFC
        POTEVP(I,J)=POTEVP(I,J)-QWBS(I,J)*SM(I,J)*XLVRW*RASRFC
        POTFLX(I,J)=POTEVP(I,J)*FACTOR
        SUBSHX(I,J)=SUBSHX(I,J)+GRNFLX(I,J)*RASRFC
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  COUNTERS
!-----------------------------------------------------------------------
!
!     APHTIM=APHTIM+1.
!     ARDSW =ARDSW +1.
!     ARDLW =ARDLW +1.
!     ASRFC =ASRFC +1.
!-----------------------------------------------------------------------
!
      END SUBROUTINE TURBL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
   SUBROUTINE pbl_driver(                                          &
                  itimestep,dt,u_frame,v_frame                     &
                 ,rublten,rvblten,rthblten                         &
                 ,tsk,xland,znt,ht                                 &
                 ,ust,pblh,hfx,qfx,grdflx                          &
                 ,u_phy,v_phy,th_phy,rho                           &
                 ,p_phy,pi_phy,p8w,t_phy,dz8w,z                    &
                 ,tke_myj,el_myj,exch_h,akhs,akms                  &
                 ,thz0,qz0,uz0,vz0,qsfc                            &
                 ,lowlyr                                           &
                 ,psim,psih,gz1oz0, wspd,br,chklowq                &
                 ,bl_pbl_physics, ra_lw_physics, dx                &
                 ,stepbl,warm_rain                                 &
                 ,kpbl,ct,lh,snow,xice                             &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,i_start,i_end, j_start,j_end, kts,kte, num_tiles &
             ! Optional
                 ,hol, mol, regime                                 &
             !  Optional moisture tracers
                 ,qv_curr, qc_curr, qr_curr                        &
                 ,qi_curr, qs_curr, qg_curr                        &
                 ,rqvblten,rqcblten,rqiblten                       &
                 ,rqrblten,rqsblten,rqgblten                       &
             !  Optional moisture tracer flags
                 ,f_qv,f_qc,f_qr                                   &
                 ,f_qi,f_qs,f_qg                                   &
                                                                   )
!------------------------------------------------------------------
!------------------------------------------------------------------
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
!-----------------------------------------------------------------
!-- RUBLTEN       U tendency due to 
!                 PBL parameterization (m/s^2)
!-- RVBLTEN       V tendency due to 
!                 PBL parameterization (m/s^2)
!-- RTHBLTEN      Theta tendency due to 
!                 PBL parameterization (K/s)
!-- RQVBLTEN      Qv tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- RQCBLTEN      Qc tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- RQIBLTEN      Qi tendency due to 
!                 PBL parameterization (kg/kg/s)
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           downward short wave flux at ground surface (W/m^2)
!-- EMISS         surface emissivity (between 0 and 1)
!-- TSK           surface temperature (K)
!-- TMN           soil temperature at lower boundary (K)
!-- XLAND         land mask (1 for land, 2 for water)
!-- ZNT           roughness length (m)
!-- MAVAIL        surface moisture availability (between 0 and 1)
!-- UST           u* in similarity theory (m/s)
!-- MOL           T* (similarity theory) (K)
!-- HOL           PBL height over Monin-Obukhov length
!-- PBLH          PBL height (m)
!-- CAPG          heat capacity for soil (J/K/m^3)
!-- THC           thermal inertia (Cal/cm/K/s^0.5)
!-- SNOWC         flag indicating snow coverage (1 for snow cover)
!-- HFX           upward heat flux at the surface (W/m^2)
!-- QFX           upward moisture flux at the surface (kg/m^2/s)
!-- REGIME        flag indicating PBL regime (stable, unstable, etc.)
!-- tke_myj       turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
!-- el_myj        mixing length from Mellor-Yamada-Janjic (MYJ) (m)
!-- akhs          sfc exchange coefficient of heat/moisture from MYJ
!-- akms          sfc exchange coefficient of momentum from MYJ
!-- thz0          potential temperature at roughness length (K)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- p8w           pressure at full levels (Pa)
!-- t_phy         temperature (K)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- DX            horizontal space interval (m)
!-- DT            time step (second)
!-- n_moist       number of moisture species
!-- PSFC          pressure at the surface (Pa)
!-- TSLB          
!-- ZS
!-- DZS
!-- num_soil_layers number of soil layer
!-- IFSNOW      ifsnow=1 for snow-cover effects
!
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QS          species index for snow
!-- P_QG          species index for graupel
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
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!******************************************************************
!------------------------------------------------------------------ 
!


   INTEGER,    INTENT(IN   )    ::     bl_pbl_physics, ra_lw_physics

   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       kts,kte, num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                   &
  &                                    i_start,i_end,j_start,j_end

   INTEGER,    INTENT(IN   )    ::     itimestep,STEPBL
   INTEGER,    DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN   )    ::                        LOWLYR
!
   LOGICAL,      INTENT(IN   )    ::   warm_rain
!
   REAL,       INTENT(IN   )    ::     DT,DX


!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(IN   )    ::                         p_phy, &
                                                          pi_phy, &
                                                             p8w, &
                                                             rho, &
                                                           t_phy, &
                                                           u_phy, &
                                                           v_phy, &
                                                            dz8w, &
                                                               z, &
                                                          th_phy
!
!
   REAL,       DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN   )    ::                         XLAND, &
                                                              HT, &
                                                            PSIM, &
                                                            PSIH, &
                                                          GZ1OZ0, &
                                                              BR, &
                                                         CHKLOWQ
!
   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               INTENT(INOUT)    ::                           TSK, &
                                                             UST, &
                                                            PBLH, &
                                                             HFX, &
                                                             QFX, &
                                                             ZNT, &
                                                            QSFC, &
                                                            AKHS, &
                                                            AKMS, &
                                                             QZ0, &
                                                            THZ0, &
                                                             UZ0, &
                                                             VZ0, &
                                                              CT, &
                                                          GRDFLX  , &
                                                            WSPD

!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       RUBLTEN, &
                                                         RVBLTEN, &
                                                        RTHBLTEN, &
                                                  EXCH_H,TKE_MYJ
!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(OUT)    ::                          EL_MYJ

   REAL ,                             INTENT(IN   )  ::  u_frame, &
                                                         v_frame
!

   INTEGER,    DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(INOUT) ::                             KPBL

   REAL,       DIMENSION( ims:ime , jms:jme ),                    &
               INTENT(IN)    :: XICE, SNOW, LH

!
! Optional
!
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr                  &
                     ,rqvblten,rqcblten,rqrblten                 &
                     ,rqiblten,rqsblten,rqgblten

   REAL,       DIMENSION( ims:ime, jms:jme )                    , &
               OPTIONAL                                         , &
               INTENT(INOUT)    ::                           HOL, &
                                                             MOL, &
                                                          REGIME

!  LOCAL  VAR

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::v_phytmp
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::u_phytmp

   REAL,       DIMENSION( ims:ime, jms:jme )          ::  TSKOLD, &
                                                          USTOLD, &
                                                          ZNTOLD, &
                                                             ZOL, &
                                                            PSFC

!

   REAL    :: DTMIN,DTBL
!
   INTEGER :: i,J,K,NK,jj,ij,its,ite,jts,jte
   LOGICAL :: radiation
   LOGICAL :: flag_qv, flag_qc, flag_qr, flag_qi, flag_qs, flag_qg
   CHARACTER*256 :: message

!------------------------------------------------------------------
!

  flag_qv = .FALSE. ; IF ( PRESENT( F_QV ) ) flag_qv = F_QV
  flag_qc = .FALSE. ; IF ( PRESENT( F_QC ) ) flag_qv = F_QC
  flag_qr = .FALSE. ; IF ( PRESENT( F_QR ) ) flag_qv = F_QR
  flag_qi = .FALSE. ; IF ( PRESENT( F_QI ) ) flag_qv = F_QI
  flag_qs = .FALSE. ; IF ( PRESENT( F_QS ) ) flag_qv = F_QS
  flag_qg = .FALSE. ; IF ( PRESENT( F_QG ) ) flag_qv = F_QG


  if (bl_pbl_physics .eq. 0) return
! RAINBL in mm (Accumulation between PBL calls)


  IF (itimestep .eq. 1 .or. mod(itimestep,STEPBL) .eq. 0) THEN

  radiation = .false.
  IF (ra_lw_physics .gt. 0) radiation = .true.

!---- 
! CALCULATE CONSTANT
 
   DTMIN=DT/60.
! PBL schemes need PBL time step for updates
   DTBL=DT*STEPBL

! SAVE OLD VALUES

!jaa   !$omp parallel do   &
!jaa   !$omp private ( ij,i,j,k )
   DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         TSKOLD(i,j)=TSK(i,j)
         USTOLD(i,j)=UST(i,j)
         ZNTOLD(i,j)=ZNT(i,j)

! REVERSE ORDER IN THE VERTICAL DIRECTION

! testing change later

         DO k=kts,kte
            v_phytmp(i,k,j)=v_phy(i,k,j)+v_frame
            u_phytmp(i,k,j)=u_phy(i,k,j)+u_frame
         ENDDO

! PSFC : in Pa

         PSFC(I,J)=p8w(I,kms,J)

         DO k=kts,min(kte+1,kde)
            RTHBLTEN(I,K,J)=0.
            RUBLTEN(I,K,J)=0.
            RVBLTEN(I,K,J)=0.
            IF ( PRESENT( RQCBLTEN )) RQCBLTEN(I,K,J)=0.
            IF ( PRESENT( RQVBLTEN )) RQVBLTEN(I,K,J)=0.
         ENDDO

         IF (flag_QI .AND. PRESENT(RQIBLTEN) ) THEN
            DO k=kts,min(kte+1,kde)
               RQIBLTEN(I,K,J)=0.
            ENDDO
         ENDIF
      ENDDO
      ENDDO

   ENDDO
!jaa   !$omp end parallel do
!
!jaa  !$omp parallel do   &
!jaa  !$omp private ( ij, i,j,k, its, ite, jts, jte )
  DO ij = 1 , num_tiles

   its = i_start(ij)
   ite = i_end(ij)
   jts = j_start(ij)
   jte = j_end(ij)

   pbl_select: SELECT CASE(bl_pbl_physics)

      CASE (YSUSCHEME)
!       CALL wrf_debug(100,'in YSU PBL')
!          IF ( PRESENT( qv_curr )  .AND. PRESENT( qc_curr )  .AND. &
!               PRESENT( qi_curr )                            .AND. &
!               PRESENT( rqvblten ) .AND. PRESENT( rqcblten ) .AND. &
!               PRESENT( rqiblten )                           .AND. &
!               PRESENT( hol      )                           .AND. &
!                                                       .TRUE.  ) THEN
!            CALL ysu(                                              &
!              U3D=u_phytmp,V3D=v_phytmp,TH3D=th_phy,T3D=t_phy      &
!             ,QV3D=qv_curr,QC3D=qc_curr,QI3D=qi_curr               &
!             ,P3D=p_phy,PI3D=pi_phy                                &
!             ,RUBLTEN=rublten,RVBLTEN=rvblten                      &
!             ,RTHBLTEN=rthblten,RQVBLTEN=rqvblten                  &
!             ,RQCBLTEN=rqcblten,RQIBLTEN=rqiblten                  &
!             ,CP=cp,G=g,ROVCP=rcp,RD=r_D,ROVG=rovg                 &
!             ,DZ8W=dz8w,Z=z,XLV=XLV,RV=r_v,PSFC=PSFC               &
!             ,ZNT=znt,UST=ust,ZOL=zol,HOL=hol,HPBL=pblh            &
!             ,PSIM=psim,PSIH=psih,XLAND=xland                      &
!             ,HFX=hfx,QFX=qfx,TSK=tskold,GZ1OZ0=gz1oz0             &
!             ,WSPD=wspd,BR=br,DT=dtbl,DTMIN=dtmin,KPBL2D=kpbl      &
!             ,SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0            &
!             ,EP1=ep_1,EP2=ep_2,KARMAN=karman,EOMEG=eomeg          &
!             ,STBOLT=stbolt,EXCH_H=exch_h,REGIME=regime            &
!             ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde      &
!             ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme      &
!             ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte      &
!                                                                   )
!          ELSE
!              CALL wrf_error_fatal('Lack arguments to call YSU pbl')
!          ENDIF

      CASE (MRFSCHEME)
!          IF ( PRESENT( qv_curr )  .AND. PRESENT( qc_curr )  .AND. &
!               PRESENT( rqvblten ) .AND. PRESENT( rqcblten ) .AND. &
!               PRESENT( hol      )                           .AND. &
!                                                       .TRUE.  ) THEN

!            CALL wrf_debug(100,'in MRF')
!            CALL mrf(                                              &
!              U3D=u_phytmp,V3D=v_phytmp,TH3D=th_phy,T3D=t_phy      &
!             ,QV3D=qv_curr                                         &
!             ,QC3D=qc_curr                                         &
!             ,QI3D=qi_curr                                         &
!             ,P3D=p_phy,PI3D=pi_phy                                &
!             ,RUBLTEN=rublten,RVBLTEN=rvblten                      &
!             ,RTHBLTEN=rthblten,RQVBLTEN=rqvblten                  &
!             ,RQCBLTEN=rqcblten,RQIBLTEN=rqiblten                  &
!             ,CP=cp,G=g,ROVCP=rcp,R=r_d,ROVG=rovg                  &
!             ,DZ8W=dz8w,Z=z,XLV=xlv,RV=r_v,PSFC=psfc               &
!             ,ZNT=znt,UST=ust,ZOL=zol,HOL=hol                      &
!             ,PBL=pblh,PSIM=psim,PSIH=psih                         &
!             ,XLAND=xland,HFX=hfx,QFX=qfx,TSK=tskold               &
!             ,GZ1OZ0=gz1oz0,WSPD=wspd,BR=br                        &
!             ,DT=dtbl,DTMIN=dtmin,KPBL2D=kpbl                      &
!             ,SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0            &
!             ,EP1=ep_1,EP2=ep_2,KARMAN=karman,EOMEG=eomeg          &
!             ,STBOLT=stbolt,REGIME=regime                          &
!             ,FLAG_QI=flag_qi                                      &
!             ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde      &
!             ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme      &
!             ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte      &
!                                                                   )
!          ELSE
!              CALL wrf_error_fatal('Lack arguments to call MRF pbl')
!          ENDIF

      CASE (GFSSCHEME)
!          IF ( PRESENT( qv_curr )  .AND. PRESENT( qc_curr )  .AND. &
!               PRESENT( rqvblten ) .AND. PRESENT( rqcblten ) .AND. &
!                                                       .TRUE.  ) THEN
!            CALL wrf_debug(100,'in GFS')
!            CALL bl_gfs(                                           &
!              U3D=u_phytmp,V3D=v_phytmp                            &
!             ,TH3D=th_phy,T3D=t_phy                                &
!             ,QV3D=qv_curr,QC3D=qc_curr,QI3D=qi_curr               &
!             ,P3D=p_phy,PI3D=pi_phy                                &
!             ,RUBLTEN=rublten,RVBLTEN=rvblten,RTHBLTEN=rthblten    &
!             ,RQVBLTEN=rqvblten,RQCBLTEN=rqcblten                  &
!             ,RQIBLTEN=rqiblten                                    &
!             ,CP=cp,G=g,ROVCP=rcp,R=r_d,ROVG=rovg,FLAG_QI=flag_qi  &
!             ,DZ8W=dz8w,z=z,PSFC=psfc                              &
!             ,UST=ust,PBL=pblh,PSIM=psim,PSIH=psih                 &
!             ,HFX=hfx,QFX=qfx,TSK=tskold,GZ1OZ0=gz1oz0             &
!             ,WSPD=wspd,BR=br                                      &
!             ,DT=dtbl,KPBL2D=kpbl,EP1=ep_1,KARMAN=karman           &
!             ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde      &
!             ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme      &
!             ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte      &
!                                                                   )
!          ELSE
!              CALL wrf_error_fatal('Lack arguments to call GFS pbl')
!          ENDIF

      CASE (MYJPBLSCHEME)
           IF ( PRESENT( qv_curr )  .AND. PRESENT( qc_curr )  .AND. &
                PRESENT( rqvblten ) .AND. PRESENT( rqcblten ) .AND. &
                                                        .TRUE.  ) THEN
             CALL myjpbl(                                           &
               DT=dt,STEPBL=stepbl,HT=ht,DZ=dz8w                    &
              ,PMID=p_phy,PINT=p8w,TH=th_phy,T=t_phy,EXNER=pi_phy   &
              ,QV=qv_curr, CWM=qc_curr                               &
              ,U=u_phy,V=v_phy,RHO=rho                              &
              ,TSK=tsk,QSFC=qsfc,CHKLOWQ=chklowq,THZ0=thz0          &
              ,QZ0=qz0,UZ0=uz0,VZ0=vz0                              &
              ,LOWLYR=lowlyr                                        &
              ,XLAND=xland,SICE=xice,SNOW=snow                      &
              ,TKE_MYJ=tke_myj,EXCH_H=exch_h,USTAR=ust,ZNT=znt      &
              ,EL_MYJ=el_myj,PBLH=pblh,KPBL=kpbl,CT=ct              &
              ,AKHS=akhs,AKMS=akms,ELFLX=lh                         &
              ,RUBLTEN=rublten,RVBLTEN=rvblten,RTHBLTEN=rthblten    &
              ,RQVBLTEN=rqvblten,RQCBLTEN=rqcblten                  &
              ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde      &
              ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme      &
              ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte      &
                                                                    )
           ELSE
!nmmb          CALL wrf_error_fatal('Lack arguments to call MYJ pbl')
               WRITE(0,*)'Lack arguments to call MYJ pbl'
           ENDIF

     CASE DEFAULT

!nmmb  WRITE( message , * ) 'The pbl option does not exist: bl_pbl_physics = ', bl_pbl_physics
!nmmb  CALL wrf_error_fatal ( message )
       WRITE(0,*)'The pbl option does not exist: bl_pbl_physics = ', bl_pbl_physics

   END SELECT pbl_select

   ENDDO
!jaa   !$omp end parallel do

   ENDIF
!
   END SUBROUTINE pbl_driver
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
! REFERENCES:  Janjic (2002), NCEP Office Note 437
!              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
!
! ABSTRACT:
!     MYJ UPDATES THE TURBULENT KINETIC ENERGY WITH THE PRODUCTION/
!     DISSIPATION TERM AND THE VERTICAL DIFFUSION TERM
!     (USING AN IMPLICIT FORMULATION) FROM MELLOR-YAMADA
!     LEVEL 2.5 AS EXTENDED BY JANJIC.  EXCHANGE COEFFICIENTS FOR
!     THE SURFACE AND FOR ALL LAYER INTERFACES ARE COMPUTED FROM
!     MONIN-OBUKHOV THEORY.
!     THE TURBULENT VERTICAL EXCHANGE IS THEN EXECUTED.
!
!-----------------------------------------------------------------------
      SUBROUTINE MYJPBL(DT,STEPBL,HT,DZ                                &
     &                 ,PMID,PINT,TH,T,EXNER,QV,CWM,U,V,RHO            &
     &                 ,TSK,QSFC,CHKLOWQ,THZ0,QZ0,UZ0,VZ0              &
     &                 ,LOWLYR,XLAND,SICE,SNOW                         &
     &                 ,TKE_MYJ,EXCH_H,USTAR,ZNT,EL_MYJ,PBLH,KPBL,CT   &
     &                 ,AKHS,AKMS,ELFLX                                &
     &                 ,RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN     &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: STEPBL

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: LOWLYR
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: KPBL
!
      REAL,INTENT(IN) :: DT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HT,SICE,SNOW       &
     &                                             ,TSK,XLAND
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: CWM,DZ     &
     &                                                     ,EXNER      &
     &                                                     ,PMID,PINT  &
     &                                                     ,QV,RHO     &
     &                                                     ,T,TH,U,V   
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: PBLH
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: AKHS,AKMS
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME)                          &
     &    ,INTENT(OUT) ::                      EL_MYJ                  &
     &                                        ,RQCBLTEN,RQVBLTEN       &
     &                                        ,RTHBLTEN                &
     &                                        ,RUBLTEN,RVBLTEN        
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CT,QSFC,QZ0     &
     &                                                ,THZ0,USTAR      &
     &                                                ,UZ0,VZ0,ZNT
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME)                          &
     &    ,INTENT(INOUT) ::                    EXCH_H,TKE_MYJ
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CHKLOWQ,ELFLX
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,KFLIP,LLOW,LMH,LMXL
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE) :: LPBL
!
      REAL :: AKHS_DENS,AKMS_DENS,APEX,DCDT,DELTAZ,DQDT,DTDIF,DTDT     &
     &       ,DTTURBL,DUDT,DVDT,EXNSFC,PSFC,PTOP,QFC1,QLOW,QOLD        &
     &       ,RATIOMX,RDTTURBL,RG,RWMSK,SEAMASK,THNEW,THOLD,TX         &
     &       ,ULOW,VLOW,WMSK
!
      REAL,DIMENSION(KTS:KTE) :: CWMK,PK,Q2K,QK,THEK,TK,UK,VK
!
      REAL,DIMENSION(KTS:KTE-1) :: AKHK,AKMK,EL,GH,GM
!
      REAL,DIMENSION(KTS:KTE+1) :: ZHK
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: THSK
!
      REAL,DIMENSION(KTS:KTE,ITS:ITE) :: RHOK
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: APE,THE
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE-1,JTS:JTE) :: AKH,AKM
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE+1,JTS:JTE) :: ZINT
!
!***  Begin debugging
      REAL :: ZSL_DIAG
      INTEGER :: IMD,JMD,PRINT_DIAG
!***  End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
!***  Begin debugging
      IMD=(IMS+IME)/2
      JMD=(JMS+JME)/2
!***  End debugging
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      DTTURBL=DT*STEPBL
      RDTTURBL=1./DTTURBL
      DTDIF=DTTURBL
      RG=1./G
!
      DO J=JTS,JTE
      DO K=KTS,KTE-1
      DO I=ITS,ITE
        AKM(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTS,KTE+1
      DO I=ITS,ITE
        ZINT(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        ZINT(I,KTE+1,J)=HT(I,J)     ! Z at bottom of lowest sigma layer
!
!!!!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,KTE+1,J)=1.E-4       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(KTE+1)=1.E-4            ! Z of bottom of lowest eta layer
!
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTE,KTS,-1
        KFLIP=KTE+1-K
        DO I=ITS,ITE
          ZINT(I,K,J)=ZINT(I,K+1,J)+DZ(I,KFLIP,J)
          APEX=1./EXNER(I,K,J)
          APE(I,K,J)=APEX
          TX=T(I,K,J)
          THE(I,K,J)=(CWM(I,K,J)*(-ELOCP/TX)+1.)*TH(I,K,J)
        ENDDO
      ENDDO
      ENDDO
!
      EL_MYJ = 0.
!
!----------------------------------------------------------------------
!.......................................................................
!$omp parallel do &
!$omp private(j,i,lmh,ptop,psfc,seamask,k,kflip,tk,thek,ratiomx,qk,cwmk, &
!$omp         pk,uk,vk,q2k,zhk,lmxl,gm,gh,el,akmk,akhk,deltaz),          &
!$omp         SCHEDULE(dynamic)
!.......................................................................
!----------------------------------------------------------------------
      setup_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
          PTOP=PINT(I,KTE+1,J)      ! KTE+1=KME
          PSFC=PINT(I,LOWLYR(I,J),J)
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            TK(K)=T(I,KFLIP,J)
            THEK(K)=THE(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            CWMK(K)=CWM(I,KFLIP,J)
            PK(K)=PMID(I,KFLIP,J)
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
!
!***  TKE=0.5*(q**2) ==> q**2=2.*TKE
!
            Q2K(K)=2.*TKE_MYJ(I,KFLIP,J)
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,K,J)
!
          ENDDO
          ZHK(KTE+1)=HT(I,J)          ! Z at bottom of lowest sigma layer
!
!***  Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=1
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=2
!***  End debugging
!
!----------------------------------------------------------------------
!***
!***  FIND THE MIXING LENGTH
!***
          CALL MIXLEN(LMH,UK,VK,TK,THEK,QK,CWMK                        &
     &               ,Q2K,ZHK,GM,GH,EL                                 &
     &               ,PBLH(I,J),LPBL(I,J),LMXL,CT(I,J)                 &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!
!----------------------------------------------------------------------
!***
!***  SOLVE FOR THE PRODUCTION/DISSIPATION OF
!***  THE TURBULENT KINETIC ENERGY
!***
!
          CALL PRODQ2(LMH,DTTURBL,USTAR(I,J),GM,GH,EL,Q2K              &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!
!----------------------------------------------------------------------
!*** THE MODEL LAYER (COUNTING UPWARD) CONTAINING THE TOP OF THE PBL
!----------------------------------------------------------------------
!
          KPBL(I,J)=KTE-LPBL(I,J)+1
!
!----------------------------------------------------------------------
!***
!***  FIND THE EXCHANGE COEFFICIENTS IN THE FREE ATMOSPHERE
!***
          CALL DIFCOF(LMH,LMXL,GM,GH,EL,TK,Q2K,ZHK,AKMK,AKHK      &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,I,J,PRINT_DIAG)   ! debug
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH 
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS KTS TO KTE-1.  COUNTING 
!***  COUNTING UPWARD FROM THE BOTTOM, THOSE SAME COEFFICIENTS EXCH_H
!***  ARE DEFINED ON THE TOPS OF THE LAYERS KTS TO KTE-1.
!
          DO K=KTS,KTE-1
            KFLIP=KTE-K
            AKH(I,K,J)=AKHK(K)
            AKM(I,K,J)=AKMK(K)
            DELTAZ=0.5*(ZHK(KFLIP)-ZHK(KFLIP+2))
            EXCH_H(I,K,J)=AKHK(KFLIP)*DELTAZ
          ENDDO
!
!----------------------------------------------------------------------
!***
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TURBULENT KINETIC ENERGY
!***
!
          CALL VDIFQ(LMH,DTDIF,Q2K,EL,ZHK                              &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!***  SAVE THE NEW TKE AND MIXING LENGTH.
!
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            Q2K(KFLIP)=AMAX1(Q2K(KFLIP),EPSQ2)
            TKE_MYJ(I,K,J)=0.5*Q2K(KFLIP)
            IF(K<KTE)EL_MYJ(I,K,J)=EL(K)   ! EL IS NOT DEFINED AT KTE
          ENDDO
!
        ENDDO
!
!----------------------------------------------------------------------
!
      ENDDO setup_integration
!
!.......................................................................
!$omp end parallel do
!.......................................................................
!----------------------------------------------------------------------
!
!***  CONVERT SURFACE SENSIBLE TEMPERATURE TO POTENTIAL TEMPERATURE.
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        PSFC=PINT(I,LOWLYR(I,J),J)
        THSK(I,J)=TSK(I,J)*(1.E5/PSFC)**CAPA
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!.......................................................................
!$omp parallel do  private(i,j,k    &
!$omp  & ,kflip,thek,ratiomx,qk,cwmk,zhk,rhok   &
!$omp  & ,akhk,seamask,llow,akhs_dens,qfc1,qlow,psfc,exnsfc,lmh &
!$omp  & ,thold,thnew,dtdt,qold,dqdt,dcdt,zsl_diag,akmk,akms_dens,uk,vk &
!$omp  & ,dudt,dvdt),SCHEDULE(dynamic)
!.......................................................................
!----------------------------------------------------------------------
      main_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            THEK(K)=THE(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            CWMK(K)=CWM(I,KFLIP,J)
            ZHK(K)=ZINT(I,K,J)
            RHOK(K,I)=PMID(I,KFLIP,J)/(R_D*T(I,KFLIP,J)*               &
     &                                (1.+P608*QK(K)-CWMK(K)))
          ENDDO
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS KTS TO KTE-1.  THESE COEFFICIENTS
!***  ARE ALSO MULTIPLIED BY THE DENSITY AT THE BOTTOM INTERFACE LEVEL.
!
          DO K=KTS,KTE-1
            AKHK(K)=AKH(I,K,J)*0.5*(RHOK(K,I)+RHOK(K+1,I))
          ENDDO
!
          ZHK(KTE+1)=ZINT(I,KTE+1,J)
!
          SEAMASK=XLAND(I,J)-1.
          THZ0(I,J)=(1.-SEAMASK)*THSK(I,J)+SEAMASK*THZ0(I,J)
!
          LLOW=LOWLYR(I,J)
          AKHS_DENS=AKHS(I,J)*RHOK(KTE+1-LLOW,I)
!
          IF(SEAMASK<0.5)THEN
            QFC1=XLV*CHKLOWQ(I,J)*AKHS_DENS
!
            IF(SNOW(I,J)>0..OR.SICE(I,J)>0.5)THEN
              QFC1=QFC1*RLIVWV
            ENDIF
!
            IF(QFC1>0.)THEN
              QLOW=QK(KTE+1-LLOW)
              QSFC(I,J)=QLOW+ELFLX(I,J)/QFC1
            ENDIF
!
          ELSE
            PSFC=PINT(I,LOWLYR(I,J),J)
            EXNSFC=(1.E5/PSFC)**CAPA

            QSFC(I,J)=PQ0SEA/PSFC                                      &
     &         *EXP(A2*(THSK(I,J)-A3*EXNSFC)/(THSK(I,J)-A4*EXNSFC))
          ENDIF
!
          QZ0 (I,J)=(1.-SEAMASK)*QSFC(I,J)+SEAMASK*QZ0 (I,J)
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TEMPERATURE AND WATER VAPOR
!----------------------------------------------------------------------
!
          CALL VDIFH(DTDIF,LMH,THZ0(I,J),QZ0(I,J)                      &
     &              ,AKHS_DENS,CHKLOWQ(I,J),CT(I,J)                    &
     &              ,THEK,QK,CWMK,AKHK,ZHK,RHOK(KTS,I)                 &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            THOLD=TH(I,K,J)
            THNEW=THEK(KFLIP)+CWMK(KFLIP)*ELOCP*APE(I,K,J)
            DTDT=(THNEW-THOLD)*RDTTURBL
            QOLD=QV(I,K,J)/(1.+QV(I,K,J))
            DQDT=(QK(KFLIP)-QOLD)*RDTTURBL
            DCDT=(CWMK(KFLIP)-CWM(I,K,J))*RDTTURBL
!
            RTHBLTEN(I,K,J)=DTDT
            RQVBLTEN(I,K,J)=DQDT/(1.-QK(KFLIP))**2
            RQCBLTEN(I,K,J)=DCDT
          ENDDO
!
!*** Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=0
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=0
!*** End debugging
!
        PSFC=.01*PINT(I,LOWLYR(I,J),J)
        ZSL_DIAG=0.5*DZ(I,1,J)
!
!*** Begin debugging
!         IF(PRINT_DIAG==1)THEN
!
!           write(6,"(a, 2i5, 2i3, 2f8.2, f6.2, 2f8.2)") &
!           '{turb4 i,j, Kpbl, Kmxl, Psfc, Zsfc, Zsl, Zpbl, Zmxl = ' &
!           , i, j, KPBL(i,j), KTE-LMXL+1, PSFC, ZHK(LMH+1), ZSL_diag  &
!           , PBLH(i,j), ZHK(LMXL)-ZHK(LMH+1)
!           write(6,"(a, 2f7.2, f7.3, 3e11.4)") &
!           '{turb4 tsk, thsk, qz0, q**2_0, akhs, exch_0 = ' &
!           , tsk(i,j)-273.15, thsk(i,j), 1000.*qz0(i,j) &
!           , 2.*tke_myj(i,1,j), akhs(i,j), akhs(i,j)*ZSL_diag
!           write(6,"(a)") &
!           '{turb5 k, Pmid, Pint_1, Tc, TH, DTH, GH, GM, EL, Q**2, Akh, EXCH_h, Dz, Dp'
!           do k=kts,kte/2
!             KFLIP=KTE-K   !-- Includes the KFLIP-1 in earlier versions
!             write(6,"(a,i3, 2f8.2, 2f8.3, 3e12.4, 4e11.4, f7.2, f6.2)") &
!            '{turb5 ', k, .01*pmid(i,k,j),.01*pint(i,k,j), T(i,k,j)-273.15 &
!            , th(i,k,j), DTTURBL*rthblten(i,k,j), GH(KFLIP), GM(KFLIP) &
!            , el_myj(i,KFLIP,j), 2.*tke_myj(i,k+1,j), akh(i,KFLIP,j) &
!            , exch_h(i,k,j), dz(i,k,j), .01*(pint(i,k,j)-pint(i,k+1,j))
!           enddo
!
!         ELSEIF(PRINT_DIAG==2)THEN
!
!           write(6,"(a, 2i5, 2i3, 2f8.2, f6.2, 2f8.2)") &
!           '}turb4 i,j, Kpbl, Kmxl, Psfc, Zsfc, Zsl, Zpbl, Zmxl = ' &
!           , i, j, KPBL(i,j), KTE-LMXL+1, PSFC, ZHK(LMH+1), ZSL_diag  &
!           , PBLH(i,j), ZHK(LMXL)-ZHK(LMH+1)
!           write(6,"(a, 2f7.2, f7.3, 3e11.4)") &
!           '}turb4 tsk, thsk, qz0, q**2_0, akhs, exch_0 = ' &
!           , tsk(i,j)-273.15, thsk(i,j), 1000.*qz0(i,j) &
!           , 2.*tke_myj(i,1,j), akhs(i,j), akhs(i,j)*ZSL_diag
!           write(6,"(a)") &
!           '}turb5 k, Pmid, Pint_1, Tc, TH, DTH, GH, GM, EL, Q**2, Akh, EXCH_h, Dz, Dp'
!           do k=kts,kte/2
!             KFLIP=KTE-K   !-- Includes the KFLIP-1 in earlier versions
!             write(6,"(a,i3, 2f8.2, 2f8.3, 3e12.4, 4e11.4, f7.2, f6.2)") &
!            '}turb5 ', k, .01*pmid(i,k,j),.01*pint(i,k,j), T(i,k,j)-273.15 &
!            , th(i,k,j), DTTURBL*rthblten(i,k,j), GH(KFLIP), GM(KFLIP) &
!            , el_myj(i,KFLIP,j), 2.*tke_myj(i,k+1,j), akh(i,KFLIP,j) &
!            , exch_h(i,k,j), dz(i,k,j), .01*(pint(i,k,j)-pint(i,k+1,j))
!           enddo
!         ENDIF
!*** End debugging
!
!----------------------------------------------------------------------
        ENDDO
!----------------------------------------------------------------------
        DO I=ITS,ITE
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
!
          DO K=KTS,KTE-1
            AKMK(K)=AKM(I,K,J)
            AKMK(K)=AKMK(K)*(RHOK(K,I)+RHOK(K+1,I))*0.5
          ENDDO
!
          LLOW=LOWLYR(I,J)
          AKMS_DENS=AKMS(I,J)*RHOK(KTE+1-LLOW,I)
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
            ZHK(K)=ZINT(I,K,J)
          ENDDO
          ZHK(KTE+1)=ZINT(I,KTE+1,J)
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  VELOCITY COMPONENTS
!----------------------------------------------------------------------
!
          CALL VDIFV(LMH,DTDIF,UZ0(I,J),VZ0(I,J)                       &
     &              ,AKMS_DENS,UK,VK,AKMK,ZHK,RHOK(KTS,I)              &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            DUDT=(UK(KFLIP)-U(I,K,J))*RDTTURBL
            DVDT=(VK(KFLIP)-V(I,K,J))*RDTTURBL
            RUBLTEN(I,K,J)=DUDT
            RVBLTEN(I,K,J)=DVDT
          ENDDO
!
        ENDDO
!----------------------------------------------------------------------
!
      ENDDO main_integration
!jaa!$omp end parallel do
!
!----------------------------------------------------------------------
!
      END SUBROUTINE MYJPBL
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE MIXLEN                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *                   LEVEL 2.5 MIXING LENGTH                      *
!   *                                                                *
!   ******************************************************************
!
     &(LMH,U,V,T,THE,Q,CWM,Q2,Z,GM,GH,EL,PBLH,LPBL,LMXL,CT             &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      INTEGER,INTENT(OUT) :: LMXL,LPBL
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: CWM,Q,Q2,T,THE,U,V
!
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,INTENT(OUT) :: PBLH
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(OUT) :: EL,GH,GM
!
      REAL,INTENT(INOUT) :: CT
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K,LPBLM
!
      REAL :: A,ADEN,B,BDEN,AUBR,BUBR,BLMX,EL0,ELOQ2X,GHL,GML           &
     &       ,QOL2ST,QOL2UN,QDZL,RDZ,SQ,SREL,SZQ,TEM,THM,VKRMZ
!
      REAL,DIMENSION(KTS:KTE) :: Q1
!
      REAL,DIMENSION(KTS:KTE-1) :: DTH,ELM,REL
!
!----------------------------------------------------------------------
!**********************************************************************
!--------------FIND THE HEIGHT OF THE PBL-------------------------------
      LPBL=LMH
!
      DO K=LMH-1,1,-1
        IF(Q2(K)<=EPSQ2*FH)THEN
          LPBL=K
          GO TO 110
        ENDIF
      ENDDO
!
      LPBL=1
!
!--------------THE HEIGHT OF THE PBL------------------------------------
!
 110  PBLH=Z(LPBL)-Z(LMH+1)
!
!-----------------------------------------------------------------------
      DO K=KTS,LMH
        Q1(K)=0.
      ENDDO
!
      DO K=1,LMH-1
        DTH(K)=THE(K)-THE(K+1)
      ENDDO
!
      DO K=LMH-2,1,-1
        IF(DTH(K)>0..AND.DTH(K+1)<=0.)THEN
          DTH(K)=DTH(K)+CT
          EXIT
        ENDIF
      ENDDO
!
      CT=0.
!----------------------------------------------------------------------
      DO K=KTS,LMH-1
        RDZ=2./(Z(K)-Z(K+2))
        GML=((U(K)-U(K+1))**2+(V(K)-V(K+1))**2)*RDZ*RDZ
        GM(K)=MAX(GML,EPSGM)
!
        TEM=(T(K)+T(K+1))*0.5
        THM=(THE(K)+THE(K+1))*0.5
!
        A=THM*P608
        B=(ELOCP/TEM-1.-P608)*THM
!
        GHL=(DTH(K)*((Q(K)+Q(K+1)+CWM(K)+CWM(K+1))*(0.5*P608)+1.)      &
     &     +(Q(K)-Q(K+1)+CWM(K)-CWM(K+1))*A                            &
     &     +(CWM(K)-CWM(K+1))*B)*RDZ
!
        IF(ABS(GHL)<=EPSGH)GHL=EPSGH
        GH(K)=GHL
      ENDDO
!
!----------------------------------------------------------------------
!***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
!----------------------------------------------------------------------
!
      LMXL=LMH
!
      DO K=KTS,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
        IF(GHL>=EPSGH)THEN
          IF(GML/GHL<=REQU)THEN
            ELM(K)=EPSL
            LMXL=K
          ELSE
            AUBR=(AUBM*GML+AUBH*GHL)*GHL
            BUBR= BUBM*GML+BUBH*GHL
            QOL2ST=(-0.5*BUBR+SQRT(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR
            ELOQ2X=1./QOL2ST
            ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
          ENDIF
        ELSE
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          QOL2UN=-0.5*BDEN+SQRT(BDEN*BDEN*0.25-ADEN)
          ELOQ2X=1./(QOL2UN+EPSRU)       ! repsr1/qol2un
          ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
        ENDIF
      ENDDO
!
      IF(ELM(LMH-1)==EPSL)LMXL=LMH
!
!----------------------------------------------------------------------
!***  THE HEIGHT OF THE MIXED LAYER
!----------------------------------------------------------------------
!
      BLMX=Z(LMXL)-Z(LMH+1)
!
!----------------------------------------------------------------------
      DO K=LPBL,LMH
        Q1(K)=SQRT(Q2(K))
      ENDDO
!----------------------------------------------------------------------
      SZQ=0.
      SQ =0.
!
      DO K=KTS,LMH-1
        QDZL=(Q1(K)+Q1(K+1))*(Z(K+1)-Z(K+2))
        SZQ=(Z(K+1)+Z(K+2)-Z(LMH+1)-Z(LMH+1))*QDZL+SZQ
        SQ=QDZL+SQ
      ENDDO
!
!----------------------------------------------------------------------
!***  COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA
!----------------------------------------------------------------------
!
      EL0=MIN(ALPH*SZQ*0.5/SQ,EL0MAX)
      EL0=MAX(EL0            ,EL0MIN)
!
!----------------------------------------------------------------------
!***  ABOVE THE PBL TOP
!----------------------------------------------------------------------
!
      LPBLM=MAX(LPBL-1,1)
!
      DO K=KTS,LPBLM
        EL(K)=MIN((Z(K)-Z(K+2))*ELFC,ELM(K))
        REL(K)=EL(K)/ELM(K)
      ENDDO
!
!----------------------------------------------------------------------
!***  INSIDE THE PBL
!----------------------------------------------------------------------
!
      IF(LPBL<LMH)THEN
        DO K=LPBL,LMH-1
          VKRMZ=(Z(K+1)-Z(LMH+1))*VKARMAN
          EL(K)=MIN(VKRMZ/(VKRMZ/EL0+1.),ELM(K))
          REL(K)=EL(K)/ELM(K)
        ENDDO
      ENDIF
!
      DO K=LPBL+1,LMH-2
        SREL=MIN(((REL(K-1)+REL(K+1))*0.5+REL(K))*0.5,REL(K))
        EL(K)=MAX(SREL*ELM(K),EPSL)
      ENDDO
!
!----------------------------------------------------------------------
      END SUBROUTINE MIXLEN
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE PRODQ2                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *            LEVEL 2.5 Q2 PRODUCTION/DISSIPATION                 *
!   *                                                                *
!   ******************************************************************
!
     &(LMH,DTTURBL,USTAR,GM,GH,EL,Q2                                   &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTTURBL,USTAR
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: GH,GM
      REAL,DIMENSION(KTS:KTE-1),INTENT(INOUT) :: EL
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: ADEN,AEQU,ANUM,ARHS,BDEN,BEQU,BNUM,BRHS,CDEN,CRHS        &
     &       ,DLOQ1,ELOQ11,ELOQ12,ELOQ13,ELOQ21,ELOQ22,ELOQ31,ELOQ32   &
     &       ,ELOQ41,ELOQ42,ELOQ51,ELOQ52,ELOQN,EQOL2,GHL,GML          &
     &       ,RDEN1,RDEN2,RHS2,RHSP1,RHSP2,RHST2
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      main_integration: DO K=1,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE EQUILIBRIUM EQUATION
!----------------------------------------------------------------------
!
        AEQU=(AEQM*GML+AEQH*GHL)*GHL
        BEQU= BEQM*GML+BEQH*GHL
!
!----------------------------------------------------------------------
!***  EQUILIBRIUM SOLUTION FOR L/Q
!----------------------------------------------------------------------
!
        EQOL2=-0.5*BEQU+SQRT(BEQU*BEQU*0.25-AEQU)
!
!----------------------------------------------------------------------
!***  IS THERE PRODUCTION/DISSIPATION ?
!----------------------------------------------------------------------
!
        IF((GML+GHL*GHL<=EPSTRB)                                       &
     &   .OR.(GHL>=EPSGH.AND.GML/GHL<=REQU)                            &
     &   .OR.(EQOL2<=EPS2))THEN
!
!----------------------------------------------------------------------
!***  NO TURBULENCE
!----------------------------------------------------------------------
!
          Q2(K)=EPSQ2
          EL(K)=EPSL
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  TURBULENCE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE NUMERATOR
!----------------------------------------------------------------------
!
          ANUM=(ANMM*GML+ANMH*GHL)*GHL
          BNUM= BNMM*GML+BNMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE NUMERATOR OF THE LINEARIZED EQ.
!----------------------------------------------------------------------
!
          ARHS=-(ANUM*BDEN-BNUM*ADEN)*2.
          BRHS=- ANUM*4.
          CRHS=- BNUM*2.
!
!----------------------------------------------------------------------
!***  INITIAL VALUE OF L/Q
!----------------------------------------------------------------------
!
          DLOQ1=EL(K)/SQRT(Q2(K))
!
!----------------------------------------------------------------------
!***  FIRST ITERATION FOR L/Q, RHS=0
!----------------------------------------------------------------------
!
          ELOQ21=1./EQOL2
          ELOQ11=SQRT(ELOQ21)
          ELOQ31=ELOQ21*ELOQ11
          ELOQ41=ELOQ21*ELOQ21
          ELOQ51=ELOQ21*ELOQ31
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN1=1./(ADEN*ELOQ41+BDEN*ELOQ21+CDEN)
!
!----------------------------------------------------------------------
!***  D(RHS)/D(L/Q)
!----------------------------------------------------------------------
!
          RHSP1=(ARHS*ELOQ51+BRHS*ELOQ31+CRHS*ELOQ11)*RDEN1*RDEN1
!
!----------------------------------------------------------------------
!***  FIRST-GUESS SOLUTION
!----------------------------------------------------------------------
!
          ELOQ12=ELOQ11+(DLOQ1-ELOQ11)*EXP(RHSP1*DTTURBL)
          ELOQ12=MAX(ELOQ12,EPS1)
!
!----------------------------------------------------------------------
!***  SECOND ITERATION FOR L/Q
!----------------------------------------------------------------------
!
          ELOQ22=ELOQ12*ELOQ12
          ELOQ32=ELOQ22*ELOQ12
          ELOQ42=ELOQ22*ELOQ22
          ELOQ52=ELOQ22*ELOQ32
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN2=1./(ADEN*ELOQ42+BDEN*ELOQ22+CDEN)
          RHS2 =-(ANUM*ELOQ42+BNUM*ELOQ22)*RDEN2+RB1
          RHSP2= (ARHS*ELOQ52+BRHS*ELOQ32+CRHS*ELOQ12)*RDEN2*RDEN2
          RHST2=RHS2/RHSP2
!
!----------------------------------------------------------------------
!***  CORRECTED SOLUTION
!----------------------------------------------------------------------
!
          ELOQ13=ELOQ12-RHST2+(RHST2+DLOQ1-ELOQ12)*EXP(RHSP2*DTTURBL)
          ELOQ13=AMAX1(ELOQ13,EPS1)
!
!----------------------------------------------------------------------
!***  TWO ITERATIONS IS ENOUGH IN MOST CASES ...
!----------------------------------------------------------------------
!
          ELOQN=ELOQ13
!
          IF(ELOQN>EPS1)THEN
            Q2(K)=EL(K)*EL(K)/(ELOQN*ELOQN)
            Q2(K)=AMAX1(Q2(K),EPSQ2)
            IF(Q2(K)==EPSQ2)THEN
              EL(K)=EPSL
            ENDIF
          ELSE
            Q2(K)=EPSQ2
            EL(K)=EPSL
          ENDIF
!
!----------------------------------------------------------------------
!***  END OF TURBULENT BRANCH
!----------------------------------------------------------------------
!
        ENDIF
!----------------------------------------------------------------------
!***  END OF PRODUCTION/DISSIPATION LOOP
!----------------------------------------------------------------------
!
      ENDDO main_integration
!
!----------------------------------------------------------------------
!***  LOWER BOUNDARY CONDITION FOR Q2
!----------------------------------------------------------------------
!
      Q2(LMH)=AMAX1(B1**(2./3.)*USTAR*USTAR,EPSQ2)
!----------------------------------------------------------------------
!
      END SUBROUTINE PRODQ2
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE DIFCOF                           &
!   ******************************************************************
!   *                                                                *
!   *                LEVEL 2.5 DIFFUSION COEFFICIENTS                *
!   *                                                                *
!   ******************************************************************
     &(LMH,LMXL,GM,GH,EL,T,Q2,Z,AKM,AKH                                &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE,I,J,PRINT_DIAG)   ! debug
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH,LMXL
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: Q2,T
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: EL,GH,GM
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(OUT) :: AKH,AKM
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K,KINV
!
      REAL :: ADEN,AKMIN,BDEN,BESH,BESM,CDEN,D2T,ELL,ELOQ2,ELOQ4,ELQDZ &
     &       ,ESH,ESM,GHL,GML,Q1L,RDEN,RDZ
!
!*** Begin debugging
      INTEGER,INTENT(IN) :: PRINT_DIAG
!     REAL :: D2Tmin
!*** End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      DO K=1,LMH-1
        ELL=EL(K)
!
        ELOQ2=ELL*ELL/Q2(K)
        ELOQ4=ELOQ2*ELOQ2
!
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
        ADEN=(ADNM*GML+ADNH*GHL)*GHL
        BDEN= BDNM*GML+BDNH*GHL
        CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SM DETERMINANT
!----------------------------------------------------------------------
!
        BESM=BSMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SH DETERMINANT
!----------------------------------------------------------------------
!
        BESH=BSHM*GML+BSHH*GHL
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
        RDEN=1./(ADEN*ELOQ4+BDEN*ELOQ2+CDEN)
!
!----------------------------------------------------------------------
!***  SM AND SH
!----------------------------------------------------------------------
!
        ESM=(BESM*ELOQ2+CESM)*RDEN
        ESH=(BESH*ELOQ2+CESH)*RDEN
!
!----------------------------------------------------------------------
!***  DIFFUSION COEFFICIENTS
!----------------------------------------------------------------------
!
        RDZ=2./(Z(K)-Z(K+2))
        Q1L=SQRT(Q2(K))
        ELQDZ=ELL*Q1L*RDZ
        AKM(K)=ELQDZ*ESM
        AKH(K)=ELQDZ*ESH
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!***  INVERSIONS
!----------------------------------------------------------------------
!
!     IF(LMXL==LMH)THEN
!       KINV=LMH
!       D2Tmin=0.
!
!       DO K=LMH/2,LMH-1
!         D2T=T(K-1)-2.*T(K)+T(K+1)
!         IF(D2T<D2Tmin)THEN
!           D2Tmin=D2T
!           IF(D2T<0)KINV=K
!         ENDIF
!       ENDDO
!
!       IF(KINV<LMH)THEN
!         DO K=KINV-1,LMH-1
!           RDZ=2./(Z(K)-Z(K+2))
!           AKMIN=0.5*RDZ
!           AKM(K)=MAX(AKM(K),AKMIN)
!           AKH(K)=MAX(AKH(K),AKMIN)
!         ENDDO
!
!*** Begin debugging
!         IF(PRINT_DIAG>0)THEN
!           write(6,"(a,3i3)") '{turb1 lmxl,lmh,kinv=',lmxl,lmh,kinv
!           write(6,"(a,3i3)") '}turb1 lmxl,lmh,kinv=',lmxl,lmh,kinv
!           IF(PRINT_DIAG==1)THEN
!             write(6,"(a)") &
!               '{turb3 k, t, d2t, rdz, z(k), z(k+2), akmin, akh '
!           ELSE
!             write(6,"(a)") &
!               '}turb3 k, t, d2t, rdz, z(k), z(k+2), akmin, akh '
!           ENDIF
!           DO K=LMH-1,KINV-1,-1
!             D2T=T(K-1)-2.*T(K)+T(K+1)
!             RDZ=2./(Z(K)-Z(K+2))
!             AKMIN=0.5*RDZ
!             IF(PRINT_DIAG==1)THEN
!               write(6,"(a,i3,f8.3,2e12.5,2f9.2,2e12.5)") '{turb3 ' &
!               ,k,t(k)-273.15,d2t,rdz,z(k),z(k+2),akmin,akh(k)
!             ELSE
!               write(6,"(a,i3,f8.3,2e12.5,2f9.2,2e12.5)") '}turb3 ' &
!               ,k,t(k)-273.15,d2t,rdz,z(k),z(k+2),akmin,akh(k)
!             ENDIF
!           ENDDO
!         ENDIF     !- IF (print_diag > 0) THEN
!       ENDIF       !- IF(KINV<LMH)THEN
!*** End debugging
!
!     ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE DIFCOF
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE VDIFQ                            &
!   ******************************************************************
!   *                                                                *
!   *               VERTICAL DIFFUSION OF Q2 (TKE)                   *
!   *                                                                *
!   ******************************************************************
     &(LMH,DTDIF,Q2,EL,Z                                               &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTDIF
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: EL
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: ADEN,AKQS,BDEN,BESH,BESM,CDEN,CF,DTOZS,ELL,ELOQ2,ELOQ4   &
     &       ,ELQDZ,ESH,ESM,ESQHF,GHL,GML,Q1L,RDEN,RDZ
!
      REAL,DIMENSION(KTS:KTE-2) :: AKQ,CM,CR,DTOZ,RSQ2
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!***
!***  VERTICAL TURBULENT DIFFUSION
!***
!----------------------------------------------------------------------
      ESQHF=0.5*ESQ
!
      DO K=KTS,LMH-2
        DTOZ(K)=(DTDIF+DTDIF)/(Z(K)-Z(K+2))
        AKQ(K)=SQRT((Q2(K)+Q2(K+1))*0.5)*(EL(K)+EL(K+1))*ESQHF         &
     &        /(Z(K+1)-Z(K+2))
        CR(K)=-DTOZ(K)*AKQ(K)
      ENDDO
!
      CM(1)=DTOZ(1)*AKQ(1)+1.
      RSQ2(1)=Q2(1)
!
      DO K=KTS+1,LMH-2
        CF=-DTOZ(K)*AKQ(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(AKQ(K-1)+AKQ(K))*DTOZ(K)+1.
        RSQ2(K)=-RSQ2(K-1)*CF+Q2(K)
      ENDDO
!
      DTOZS=(DTDIF+DTDIF)/(Z(LMH-1)-Z(LMH+1))
      AKQS=SQRT((Q2(LMH-1)+Q2(LMH))*0.5)*(EL(LMH-1)+ELZ0)*ESQHF        &
     &    /(Z(LMH)-Z(LMH+1))
!
      CF=-DTOZS*AKQ(LMH-2)/CM(LMH-2)
!
      Q2(LMH-1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH-2)*CF+Q2(LMH-1))          &
     &        /((AKQ(LMH-2)+AKQS)*DTOZS-CR(LMH-2)*CF+1.)
!
      DO K=LMH-2,KTS,-1
        Q2(K)=(-CR(K)*Q2(K+1)+RSQ2(K))/CM(K)
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFQ
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFH(DTDIF,LMH,THZ0,QZ0,RKHS,CHKLOWQ,CT             &
     &                ,THE,Q,CWM,RKH,Z,RHO                            &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                ,IMS,IME,JMS,JME,KMS,KME                        &
     &                ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!     ***************************************************************
!     *                                                             *
!     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: CHKLOWQ,CT,DTDIF,QZ0,RKHS,THZ0
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: RKH
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: RHO
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: CWM,Q,THE
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: CF,CMB,CMCB,CMQB,CMTB,CTHF,DTOZL,DTOZS                   &
     &       ,RCML,RKHH,RKQS,RSCB,RSQB,RSTB
!
      REAL,DIMENSION(KTS:KTE-1) :: CM,CR,DTOZ,RKCT,RSC,RSQ,RST
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      CTHF=0.5*CT
!
      DO K=KTS,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKH(K)
        RKCT(K)=RKH(K)*(Z(K)-Z(K+2))*CTHF
      ENDDO
!
      CM(KTS)=DTOZ(KTS)*RKH(KTS)+RHO(KTS)
!----------------------------------------------------------------------
      RST(KTS)=-RKCT(KTS)*DTOZ(KTS)                                    &
     &         +THE(KTS)*RHO(KTS)
      RSQ(KTS)=Q(KTS)  *RHO(KTS)
      RSC(KTS)=CWM(KTS)*RHO(KTS)
!----------------------------------------------------------------------
      DO K=KTS+1,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKH(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(RKH(K-1)+RKH(K))*DTOZL+RHO(K)
        RST(K)=-RST(K-1)*CF+(RKCT(K-1)-RKCT(K))*DTOZL+THE(K)*RHO(K)
        RSQ(K)=-RSQ(K-1)*CF+Q(K)  *RHO(K)
        RSC(K)=-RSC(K-1)*CF+CWM(K)*RHO(K)
      ENDDO
!
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKHH=RKH(LMH-1)
!
      CF=-DTOZS*RKHH/CM(LMH-1)
      RKQS=RKHS*CHKLOWQ
!
      CMB=CR(LMH-1)*CF
      CMTB=-CMB+(RKHH+RKHS)*DTOZS+RHO(LMH)
      CMQB=-CMB+(RKHH+RKQS)*DTOZS+RHO(LMH)
      CMCB=-CMB+(RKHH     )*DTOZS+RHO(LMH)
!
      RSTB=-RST(LMH-1)*CF+RKCT(LMH-1)*DTOZS+THE(LMH)*RHO(LMH)
      RSQB=-RSQ(LMH-1)*CF+Q(LMH)  *RHO(LMH)
      RSCB=-RSC(LMH-1)*CF+CWM(LMH)*RHO(LMH)
!----------------------------------------------------------------------
      THE(LMH)=(DTOZS*RKHS*THZ0+RSTB)/CMTB
      Q(LMH)  =(DTOZS*RKQS*QZ0 +RSQB)/CMQB
      CWM(LMH)=(                RSCB)/CMCB
!----------------------------------------------------------------------
      DO K=LMH-1,KTS,-1
        RCML=1./CM(K)
        THE(K)=(-CR(K)*THE(K+1)+RST(K))*RCML
        Q(K)  =(-CR(K)*  Q(K+1)+RSQ(K))*RCML
        CWM(K)=(-CR(K)*CWM(K+1)+RSC(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFH
!
!---------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFV(LMH,DTDIF,UZ0,VZ0,RKMS,U,V,RKM,Z,RHO           &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                ,IMS,IME,JMS,JME,KMS,KME                        &
                      ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!     ***************************************************************
!     *                                                             *
!     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                   &
     &                     ,IMS,IME,JMS,JME,KMS,KME                   &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: RKMS,DTDIF,UZ0,VZ0
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: RKM
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: RHO
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: U,V
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: CF,DTOZAK,DTOZL,DTOZS,RCML,RCMVB,RHOK,RKMH
!
      REAL,DIMENSION(KTS:KTE-1) :: CM,CR,DTOZ,RSU,RSV
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      DO K=1,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKM(K)
      ENDDO
!
      RHOK=RHO(1)
      CM(1)=DTOZ(1)*RKM(1)+RHOK
      RSU(1)=U(1)*RHOK
      RSV(1)=V(1)*RHOK
!----------------------------------------------------------------------
      DO K=2,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKM(K-1)/CM(K-1)
        RHOK=RHO(K)
        CM(K)=-CR(K-1)*CF+(RKM(K-1)+RKM(K))*DTOZL+RHOK
        RSU(K)=-RSU(K-1)*CF+U(K)*RHOK
        RSV(K)=-RSV(K-1)*CF+V(K)*RHOK
      ENDDO
!----------------------------------------------------------------------
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKMH=RKM(LMH-1)
!
      CF=-DTOZS*RKMH/CM(LMH-1)
      RHOK=RHO(LMH)
      RCMVB=1./((RKMH+RKMS)*DTOZS-CR(LMH-1)*CF+RHOK)
      DTOZAK=DTOZS*RKMS
!----------------------------------------------------------------------
      U(LMH)=(DTOZAK*UZ0-RSU(LMH-1)*CF+U(LMH)*RHOK)*RCMVB
      V(LMH)=(DTOZAK*VZ0-RSV(LMH-1)*CF+V(LMH)*RHOK)*RCMVB
!----------------------------------------------------------------------
      DO K=LMH-1,1,-1
        RCML=1./CM(K)
        U(K)=(-CR(K)*U(K+1)+RSU(K))*RCML
        V(K)=(-CR(K)*V(K+1)+RSV(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFV
!
!-----------------------------------------------------------------------
!
!=======================================================================
!!!!  SUBROUTINE MYJPBL_INIT(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,         &
!!!! &                       TKE_MYJ,EXCH_H,RESTART,ALLOWED_TO_READ,    &
      SUBROUTINE MYJPBL_INIT(EXCH_H,RESTART                             &
     &                      ,IDS,IDE,JDS,JDE,LM                         &
     &                      ,IMS,IME,JMS,JME                            &
     &                      ,ITS,ITE,JTS,JTE                         )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!!!!  LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART
      LOGICAL,INTENT(IN) :: RESTART
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
     &                     ,IMS,IME,JMS,JME                             &
     &                     ,ITS,ITE,JTS,JTE          

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT) ::    EXCH_H
!
      INTEGER :: I,J,K,ITF,JTF,KTF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      JTF=MIN0(JTE,JDE-1)
!!!   KTF=MIN0(KTE,KDE-1)
      KTF=LM
      ITF=MIN0(ITE,IDE-1)
!
      IF(.NOT.RESTART)THEN
        DO K=1,KTF
        DO J=JTS,JTF
        DO I=ITS,ITF
!         TKE_MYJ(I,K,J)=EPSQ2
!         RUBLTEN(I,K,J)=0.
!         RVBLTEN(I,K,J)=0.
!         RTHBLTEN(I,K,J)=0.
!         RQVBLTEN(I,K,J)=0.
          EXCH_H(I,J,K)=0.
        ENDDO
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE MYJPBL_INIT
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_TURBULENCE
!
!-----------------------------------------------------------------------
