!-----------------------------------------------------------------------
!
      MODULE MODULE_CONVECTION_DEV
!
!-----------------------------------------------------------------------
!
!***  THE CONVECTION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
      USE MODULE_MY_DOMAIN_SPECS
!
      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,cappa,CP,ELIV,ELWV,EPSQ,G    &
                                 ,P608,PQ0,R_D,TIW
!
      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE

      USE MODULE_CU_BMJ_DEV
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: CUCNVC_DEV
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE CONVECTION OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: KFETASCHEME=1                     &
                                     ,BMJSCHEME=2                       &
                                     ,GDSCHEME=3                        & 
                                     ,SASSCHEME=4 
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE CUCNVC_DEV(NTSD,DT,NCNVC,NRADS,NRADL,MINUTES_HISTORY   &
                       ,FRES,FR,FSL,FSS                                 &
                       ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP           &
                       ,DYH,RESTRT,HYDRO                                &
                       ,CLDEFI,NUM_WATER                                &
                       ,F_ICE,F_RAIN                                    &
                       ,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG                   &
                       ,F_QV,F_QC,F_QR,F_QI,F_QS,F_QG                   &
                       ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1                &
                       ,PT,PD,T,Q,CWM,TCUCN,WATER                       &
                       ,OMGALF,U,V                                      &
                       ,FIS,W0AVG                                       &
                       ,PREC,ACPREC,CUPREC,CUPPT,CPRATE                 &
                       ,CNVBOT,CNVTOP,SM,LPBL                           &
                       ,HTOP,HTOPD,HTOPS                                &
                       ,HBOT,HBOTD,HBOTS                                &
                       ,AVCNVC,ACUTIM                                   &
                       ,RSWIN,RSWOUT                                    &
                       ,CONVECTION                                      &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CUCNVC      CONVECTIVE PRECIPITATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-21       
!     
! ABSTRACT:
!     CUCVNC DRIVES THE WRF CONVECTION SCHEMES
!     
! PROGRAM HISTORY LOG:
!   02-03-21  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-10-11  BLACK      - BUILT INTO UMO PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!     
! USAGE: CALL CUCNVC FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM 
!$$$  
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      CHARACTER(99),INTENT(IN):: &
       CONVECTION
!
      LOGICAL(kind=klog),INTENT(IN):: &
       HYDRO,RESTRT &
      ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP
!
      LOGICAL(kind=klog),INTENT(IN):: &
       F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
!
      INTEGER(kind=kint),INTENT(IN):: &
       IDS,IDE,JDS,JDE,LM &
      ,IMS,IME,JMS,JME &
      ,ITS,ITE,JTS,JTE &
      ,NCNVC,MINUTES_HISTORY &
      ,NRADS,NRADL,NTSD,NUM_WATER &
      ,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
!
      INTEGER(kind=kint),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       LPBL
!
      REAL(kind=kfpt),INTENT(IN):: &
       DT,DYH,FRES,FR,FSL,FSS,PT
!
      REAL(kind=kfpt),DIMENSION(1:LM),INTENT(IN):: &
       DSG2,PDSG1,PSGML1,SGML2
!
      REAL(kind=kfpt),DIMENSION(1:LM+1),INTENT(IN):: &
       PSG1,SG2 
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       FIS,PD &
      ,RSWIN,RSWOUT,SM
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       ACPREC,CLDEFI &
      ,CNVBOT,CNVTOP &
      ,CUPPT,CUPREC &
      ,HBOT,HTOP &
      ,HBOTD,HTOPD &
      ,HBOTS,HTOPS &
      ,PREC,CPRATE &
      ,ACUTIM,AVCNVC  !<-- Were scalars
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN):: &
       F_ICE &
      ,F_RAIN
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,1:LM+1,JMS:JME),INTENT(INOUT):: &
       W0AVG
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT):: &
       Q,T &
      ,CWM         &
      ,TCUCN
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN):: &
       OMGALF,U,V
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM,NUM_WATER) &
                     ,INTENT(INOUT):: &
       WATER
!
!---------------------
!***  Local Variables
!---------------------
!
      LOGICAL(kind=klog):: &
       RESTART,WARM_RAIN
!
      LOGICAL(kind=klog),DIMENSION(IMS:IME,JMS:JME):: &
       CU_ACT_FLAG
!
      INTEGER(kind=kint):: &
       CU_PHYSICS,I,ICLDCK,IJ,J,K,MNTO &
      ,N,NCUBOT,NCUTOP,N_TIMSTPS_OUTPUT
!
!
      INTEGER(kind=kint),DIMENSION(IMS:IME,JMS:JME):: &
       KPBL,LBOT,LTOP
!
      REAL(kind=kfpt):: &
       CF_HI,DQDT,DTCNVC,DTDT,FICE,FRAIN,G_INV &
      ,PCPCOL,PDSL,PLYR,QI,QL_K,QR,QW,RDTCNVC,WC &
      ,QL,TL
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME):: &
       CUBOT,CUTOP,NCA &
      ,RAINC,RAINCV,SFCZ,XLAND
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM+1):: &
       PINT,pmid,exner &
      ,th &
      ,RQCCUTEN,RQRCUTEN &
      ,RQICUTEN,RQSCUTEN &
      ,RQVCUTEN,RTHCUTEN
!
!-----------------------------------------------------------------------
!***  For temperature change check only.
!-----------------------------------------------------------------------
!zj      REAL(kind=kfpt) :: DTEMP_CHECK=1.0
      REAL(kind=kfpt) :: TCHANGE
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Translate the convection options in the config file to their
!***  analogs in the WRF Registry so that the WRF cumulus driver
!***  remains untouched.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(CONVECTION))
        CASE ('bmj_dev')
          CU_PHYSICS=2
        CASE ('kf')
          CU_PHYSICS=1
        CASE ('sas')
          CU_PHYSICS=4
        CASE ('gd')
          CU_PHYSICS=3
        CASE ('none')
!         WRITE(0,*)' User selected to run without parameterized convection.'
        CASE DEFAULT
          WRITE(0,*)' User selected CONVECTION=',TRIM(CONVECTION)
          WRITE(0,*)' Improper selection of Convection scheme in CUCNVC'
          CALL NMMB_FINALIZE
      END SELECT
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Reset the HBOT/HTOP convective cloud bottom (base) and top arrays
!***  used in radiation.  They store the maximum vertical limits of 
!***  convective cloud between radiation calls.  These arrays are out
!***  of the WRF physics and thus their values increase upward.
!***  CUPPT is the accumulated convective precipitation between
!***  radiation calls.
!-----------------------------------------------------------------------
!
      IF(MOD(NTSD,NRADS)==0.OR.MOD(NTSD,NRADL)==0)THEN
         DO J=JMS,JME
         DO I=IMS,IME
           HTOP(I,J)=0.
           HBOT(I,J)=REAL(LM+1)
           CUPPT(I,J)=0.
         ENDDO
         ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='bmj')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='sas')RETURN
!-----------------------------------------------------------------------
!
      RESTART=RESTRT
!
!-----------------------------------------------------------------------
!
      IF(CONVECTION=='kf')THEN
!
        IF(.NOT.RESTART.AND.NTSD==0)THEN
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,k)
          DO J=JTS,JTE
          DO K=1,LM+1
          DO I=ITS,ITE
            W0AVG(I,K,J)=0.
          ENDDO
          ENDDO
          ENDDO
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,k)
!-----------------------------------------------------------------------
!***  General preparation 
!-----------------------------------------------------------------------
!
!-- AVCNVC,ACUTIM were scalars but changed to 2D arrays to allow for updates in ESMF
!
      DO J=JTS,JTE
      DO I=ITS,ITE
         AVCNVC(I,J)=AVCNVC(I,J)+1.
         ACUTIM(I,J)=ACUTIM(I,J)+1.
      ENDDO
      ENDDO
!
      DTCNVC=NCNVC*DT
      RDTCNVC=1./DTCNVC
      G_INV=1./G
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,pdsl,plyr,ql,tl)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)
        RAINCV(I,J)=0.
        RAINC(I,J)=0.
        PINT(I,J,LM+1)=SG2(LM+1)*PDSL+PSG1(LM+1)
        XLAND(I,J)=SM(I,J)+1.
        NCA(I,J)=0.
        SFCZ(I,J)=FIS(I,J)*G_INV
!
        CUTOP(I,J)=999.
        CUBOT(I,J)=999.
!
!***  LPBL is the model layer containing the PBL top
!***  counting downward from the top of the domain
!***  so KPBL is the same layer counting upward from 
!***  the ground.
!
        KPBL(I,J)=LPBL(I,J)
!
!-----------------------------------------------------------------------
!***  Fill vertical working arrays.
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1
!
          PLYR=SGML2(K)*PDSL+PSGML1(K)

          QL=MAX(Q(I,J,K),EPSQ)
          TL=T(I,J,K)
          t(I,J,K)=TL
!
          th(I,J,K)=TL*(1.E5/PLYR)**cappa

!zj          PINT(I,J,K)=PINT(I,J,K+1)-PDSG1(K)-DSG2(K)*PDSL !negative p at top!

          pint(i,j,k)=sg2(k)*pdsl+psg1(k) !zj
          pmid(I,J,K)=PLYR
          exner(I,J,K)=(PLYR*1.E-5)**cappa
!
          RTHCUTEN(I,J,K)=0.
          RQVCUTEN(I,J,K)=0.
          RQCCUTEN(I,J,K)=0.
          RQRCUTEN(I,J,K)=0.
          RQICUTEN(I,J,K)=0.
          RQSCUTEN(I,J,K)=0.
        ENDDO
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Synchronize mixing ratio in water array with specific humidity.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(i,j,k)
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
!
!***  Single-Column Convection
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      IF (CU_PHYSICS /= 0) THEN

          cps_select: SELECT CASE(cu_physics)

            CASE (BMJSCHEME)

              CALL BMJDRV_DEV(TH=th,T=t ,RAINCV=raincv                  &
                         ,FRES=fres,FR=fr,FSL=fsl,FSS=fss               &
                         ,ENTRAIN=ENTRAIN ,NEWALL=NEWALL                &
                         ,NEWSWAP=NEWSWAP ,NEWUPUP=NEWUPUP              &
                         ,NODEEP=NODEEP                                 &
                         ,DT=dt ,ntsd=NTSD ,NCNVC=NCNVC                 &
                         ,CUTOP=CUTOP, CUBOT=CUBOT, KPBL=kpbl           &
                         ,PINT=PINT, PMID=pmid, exner=exner             &
                         ,CP=cp ,R=r_d ,ELWV=ELWV ,ELIV=ELIV ,G=g       &
                         ,TFRZ=TIW ,D608=P608 ,CLDEFI=cldefi            &
                         ,XLAND=xland                                   &
                         ,CU_ACT_FLAG=cu_act_flag                       &
                         ,QV=WATER(IMS,JMS,1,P_QV)                      &
                         ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=1,KDE=lm+1 &
                         ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=1,KME=lm+1 &
                         ,ITS=ITS_B1,ITE=ITE_B1                         &
                         ,JTS=JTS_B1,JTE=JTE_B1                         &
                         ,KTS=1,KTE=lm                                  &
! optionals
                         ,RTHCUTEN=rthcuten ,RQVCUTEN=rqvcuten          &
                                                               )

            CASE DEFAULT

              WRITE( 0 , * ) 'The cumulus option does not exist: cu_physics = ', cu_physics

          END SELECT cps_select


      END IF
!
!-----------------------------------------------------------------------
!***  CNVTOP/CNVBOT hold the maximum vertical limits of convective cloud 
!***  between history output times.  HBOTS/HTOPS store similiar information
!***  for shallow (nonprecipitating) convection, and HBOTD/HTOPD are for
!***  deep (precipitating) convection.  
!-----------------------------------------------------------------------
!
      CF_HI=REAL(MINUTES_HISTORY)/60.
      N_TIMSTPS_OUTPUT=NINT(3600.*CF_HI/DT)
      MNTO=MOD(NTSD,N_TIMSTPS_OUTPUT)
!
      IF(MNTO>0.AND.MNTO<=NCNVC)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          CNVBOT(I,J)=REAL(LM+1.)
          CNVTOP(I,J)=0.
          HBOTD(I,J)=REAL(LM+1.)
          HTOPD(I,J)=0.
          HBOTS(I,J)=REAL(LM+1.)
          HTOPS(I,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,k,i,dqdt,dtdt,tchange,pcpcol,ncubot,ncutop)
!.......................................................................
!-----------------------------------------------------------------------
      DO J=JTS,JTE
      DO I=ITS,ITE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Update temperature, specific humidity, and heating.
!-----------------------------------------------------------------------
!
        DO K=1,LM
!
!-----------------------------------------------------------------------
!***  RQVCUTEN in BMJDRV is the mixing ratio tendency,
!***  so retrieve DQDT by converting to specific humidity.
!-----------------------------------------------------------------------
!
          DQDT=RQVCUTEN(I,J,K)/(1.+WATER(I,J,K,P_QV))**2
!
!-----------------------------------------------------------------------
!***  RTHCUTEN in BMJDRV is DTDT over Exner.
!-----------------------------------------------------------------------
!
          DTDT=RTHCUTEN(I,J,K)*exner(I,J,K)
          T(I,J,K)=T(I,J,K)+DTDT*DTCNVC
          Q(I,J,K)=Q(I,J,K)+DQDT*DTCNVC
          TCUCN(I,J,K)=TCUCN(I,J,K)+DTDT
          WATER(I,J,K,P_QV)=Q(I,J,K)/(1.-Q(I,J,K))       !Convert to mixing ratio
!
!zj          TCHANGE=DTDT*DTCNVC
!zj          IF(ABS(TCHANGE)>DTEMP_CHECK)THEN
!zj            WRITE(0,*)'BIG T CHANGE BY CONVECTION:',TCHANGE,' at (',I,',',J,',',K,')' 
!zj	  ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Update precipitation
!-----------------------------------------------------------------------
!
        PCPCOL=RAINCV(I,J)*1.E-3*NCNVC
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
        CUPREC(I,J)=CUPREC(I,J)+PCPCOL
        CUPPT(I,J)=CUPPT(I,J)+PCPCOL
        CPRATE(I,J)=PCPCOL
!
!-----------------------------------------------------------------------
!***  Save cloud top and bottom for radiation (HTOP/HBOT) and
!***  for output (CNVTOP/CNVBOT, HTOPS/HBOTS, HTOPD/HBOTD) arrays.
!***  Must be treated separately from each other.
!-----------------------------------------------------------------------
!
        NCUTOP=NINT(CUTOP(I,J))
        NCUBOT=NINT(CUBOT(I,J))
!
        IF(NCUTOP>1.AND.NCUTOP<LM+1)THEN
          HTOP(I,J)=MAX(CUTOP(I,J),HTOP(I,J))
          CNVTOP(I,J)=MAX(CUTOP(I,J),CNVTOP(I,J))
          IF(PCPCOL>0.)THEN
            HTOPD(I,J)=MAX(CUTOP(I,J),HTOPD(I,J))
          ELSE
            HTOPS(I,J)=MAX(CUTOP(I,J),HTOPS(I,J))
          ENDIF
        ENDIF
        IF(NCUBOT>0.AND.NCUBOT<LM+1)THEN
          HBOT(I,J)=MIN(CUBOT(I,J),HBOT(I,J))
          CNVBOT(I,J)=MIN(CUBOT(I,J),CNVBOT(I,J))
          IF(PCPCOL>0.)THEN
            HBOTD(I,J)=MIN(CUBOT(I,J),HBOTD(I,J))
          ELSE
            HBOTS(I,J)=MIN(CUBOT(I,J),HBOTS(I,J))
          ENDIF
        ENDIF
!
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CUCNVC_DEV
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CONVECTION_DEV
!
!-----------------------------------------------------------------------
