!-----------------------------------------------------------------------
!
      MODULE MODULE_CU_BMJ_DEV
!
!-----------------------------------------------------------------------
!
!***  THE CONVECTION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,CP,ELIV,ELWV,EPSQ,G          &
                                 ,P608,PQ0,R_D,TIW
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: BMJ_INIT_DEV,BMJDRV_DEV
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR BMJ CONVECTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     REAL,PARAMETER ::                                                  &
     &                  DSPC=-3000.                                     &
     &                 ,DTTOP=0.,EFIFC=5.0,EFIMN=0.20                   &
     &                 ,efmntl=0.70,efmnts=0.70                         &
     &                 ,ELIWV=2.683E6,ENPLO=98500.,ENPUP=95000.         &
     &                 ,EPSDN=1.05,EPSDT=0.                             &
     &                 ,EPSNTP=.0001,EPSNTT=.0001,EPSPR=1.E-7           &
     &                 ,EPSUP=1.00                                      &
     &                 ,FUP=1./200000.                                  &
     &                 ,PBM=13000.,PFRZ=15000.,PNO=1000.                &
     &                 ,PONE=2500.,PQM=20000.                           &
     &                 ,PSH=20000.,PSHU=45000.                          &
     &                 ,RENDP=1./(ENPLO-ENPUP)                          &
     &                 ,RHLSC=0.00,RHHSC=1.10                           &
     &                 ,ROW=1.E3                                        &
     &                 ,STABDF=0.90,STABDS=0.90                         &
     &                 ,STABS=1.0,STRESH=1.10                           &
     &                 ,DTSHAL=-1.0,TREL=2400.
!
      REAL,PARAMETER :: DTtrigr=-0.0                                    &
                       ,DTPtrigr=DTtrigr*PONE      !<-- Average parcel virtual temperature deficit over depth PONE.
                                                   !<-- NOTE: CAPEtrigr is scaled by the cloud base temperature (see below)
!
      real(kind=kfpt),parameter:: &
       dspbfl_base=-3875. &
      ,dsp0fl_base=-5875. &
      ,dsptfl_base=-1875. &
      ,dspbfs_base=-3875. &
      ,dsp0fs_base=-5875. &
      ,dsptfs_base=-1875.
!
      REAL,PARAMETER :: PL=2500.,PLQ=70000.,PH=105000.                  &
     &                 ,THL=210.,THH=365.,THHQ=325.
!
      INTEGER,PARAMETER :: ITB=76,JTB=134,ITBQ=152,JTBQ=440
!
      INTEGER,PARAMETER :: ITREFI_MAX=3
!
!***  ARRAYS FOR LOOKUP TABLES
!
      REAL,DIMENSION(ITB),PRIVATE,SAVE :: STHE,THE0
      REAL,DIMENSION(JTB),PRIVATE,SAVE :: QS0,SQS
      REAL,DIMENSION(ITBQ),PRIVATE,SAVE :: STHEQ,THE0Q
      REAL,DIMENSION(ITB,JTB),PRIVATE,SAVE :: PTBL
      REAL,DIMENSION(JTB,ITB),PRIVATE,SAVE :: TTBL
      REAL,DIMENSION(JTBQ,ITBQ),PRIVATE,SAVE :: TTBLQ
                         
!***  SHARE COPIES FOR MODULE_BL_MYJPBL
!
      REAL,DIMENSION(JTB) :: QS0_EXP,SQS_EXP
      REAL,DIMENSION(ITB,JTB) :: PTBL_EXP
!
      REAL,PARAMETER :: RDP=(ITB-1.)/(PH-PL),RDPQ=(ITBQ-1.)/(PH-PLQ)  &
     &                 ,RDQ=ITB-1,RDTH=(JTB-1.)/(THH-THL)             &
     &                 ,RDTHE=JTB-1.,RDTHEQ=JTBQ-1.                   &
     &                 ,RSFCP=1./101300.
!
      REAL,PARAMETER :: AVGEFI=(EFIMN+1.)*0.5
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE BMJDRV_DEV(                                            &
     &                  IDS,IDE,JDS,JDE,KDS,KDE                         &
     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE                         &
     &                 ,fres,fr,fsl,fss                                 &
     &                 ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP           &
     &                 ,DT,NTSD,NCNVC                                   &
     &                 ,RAINCV,CUTOP,CUBOT,KPBL                         &
     &                 ,TH,T,QV                                         &
     &                 ,PINT,PMID,exner                                 &
     &                 ,CP,R,ELWV,ELIV,G,TFRZ,D608                      &
     &                 ,CLDEFI,XLAND,CU_ACT_FLAG                        &
                      ! optional
     &                 ,RTHCUTEN, RQVCUTEN                              &
     &                                                                  )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
      LOGICAL(kind=klog),INTENT(IN):: &
       ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP
!
      LOGICAL(kind=klog),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       CU_ACT_FLAG
!
      INTEGER(kind=kint),INTENT(IN):: &
       IDS,IDE,JDS,JDE,KDS,KDE &
      ,IMS,IME,JMS,JME,KMS,KME & 
      ,ITS,ITE,JTS,JTE,KTS,KTE &
      ,ntsd,NCNVC
!
      INTEGER(kind=kint),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       KPBL
!
      REAL,INTENT(IN) :: CP,DT,ELIV,ELWV,fres,fr,fsl,fss,G,R,TFRZ,D608
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       XLAND
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN):: &
       PINT &
      ,QV &
      ,exner,PMID,T,TH
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,KMS:KME) &
                     ,OPTIONAL &
                     ,INTENT(INOUT):: &
       RQVCUTEN,RTHCUTEN 
! 
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       CLDEFI,RAINCV
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT):: &
       CUBOT,CUTOP
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
      INTEGER(kind=kint):: &
       i,icldck,ierr,J,K,LBOT,lmh,LPBL,LTOP
!
      REAL(kind=kfpt):: &
       dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
      ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
      ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts &
      ,DTCNVC,seamask,PCPCOL,PSFC,PTOP,qnew,qold
! 
      REAL(kind=kfpt),DIMENSION(KTS:KTE):: &
       DPCOL,DQDT,DTDT,PCOL,QCOL,TCOL
!
!***  Begin debugging convection
      REAL :: DELQ,DELT,PLYR
      INTEGER :: IMD,JMD
      LOGICAL :: PRINT_DIAG
!***  End debugging convection
!
!-----------------------------------------------------------------------
!*********************************************************************** 
!-----------------------------------------------------------------------
!
!***  PREPARE TO CALL BMJ CONVECTION SCHEME
!
!-----------------------------------------------------------------------
!
!***  CHECK TO SEE IF THIS IS A CONVECTION TIMESTEP
!                                                                        
      ICLDCK=MOD(ntsd,NCNVC)                                              
!-----------------------------------------------------------------------
!                                                                      
!***  COMPUTE CONVECTION EVERY NCNVC*DT/60.0 MINUTES
!                                                                     
!***  Begin debugging convection
      IMD=(IMS+IME)/2
      JMD=(JMS+JME)/2
      PRINT_DIAG=.FALSE.
!***  End debugging convection

      IF(ICLDCK==0.OR.ntsd==0)THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          CU_ACT_FLAG(I,J)=.TRUE.
        ENDDO
        ENDDO

!
        DTCNVC=DT*NCNVC
!
!--set up deficit saturation pressures depending on resolution----------
        dspbfl=dspbfl_base*fres*fr
        dsp0fl=dsp0fl_base*fres*fr
        dsptfl=dsptfl_base*fres*fr
        dspbfs=dspbfs_base*fres
        dsp0fs=dsp0fs_base*fres
        dsptfs=dsptfs_base*fres
!
        dspbsl=dspbfl*fsl
        dsp0sl=dsp0fl*fsl
        dsptsl=dsptfl*fsl
        dspbss=dspbfs*fss
        dsp0ss=dsp0fs*fss
        dsptss=dsptfs*fss
!
        slopbl=(dspbfl-dspbsl)/(1.-efimn)
        slop0l=(dsp0fl-dsp0sl)/(1.-efimn)
        sloptl=(dsptfl-dsptsl)/(1.-efimn)
        slopbs=(dspbfs-dspbss)/(1.-efimn)
        slop0s=(dsp0fs-dsp0ss)/(1.-efimn)
        slopts=(dsptfs-dsptss)/(1.-efimn)
!.......................................................................
!$omp parallel do &
!$omp private (j,i,dqdt,dtdt,pcpcol,psfc,ptop,seamask,k,qcol   &
!$omp         ,tcol,pcol,dpcol,lmh,lpbl,delt,delq,plyr,lbot,ltop)
!.......................................................................
        DO J=JTS,JTE  
        DO I=ITS,ITE
!
          DO K=KTS,KTE
            DQDT(K)=0.
            DTDT(K)=0.
          ENDDO
!
          RAINCV(I,J)=0.
          PCPCOL=0.
          PSFC=PINT(I,J,kte+1)
          PTOP=PINT(I,J,1)      ! KTE+1=KME
!
!***  CONVERT TO BMJ LAND MASK (1.0 FOR SEA; 0.0 FOR LAND)
!
          seamask=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS 
!
          DO K=KTS,KTE
!
!***  CONVERT FROM MIXING RATIO TO SPECIFIC HUMIDITY
!
            QCOL(K)=MAX(EPSQ,QV(I,J,K)/(1.+QV(I,J,K)))
            TCOL(K)=T(I,J,K)
            PCOL(K)=PMID(I,J,K)
!zj            DPCOL(K)=RR(I,J,K)*G*DZ(I,J,K)
            dpcol(k)=pint(i,j,k+1)-pint(i,j,k) !zj
          ENDDO
!
!***  LOWEST LAYER ABOVE GROUND MUST ALSO BE FLIPPED
!
          LMH=KTE!!! LOWLYR(I,J)
          LPBL=KPBL(I,J)
!-----------------------------------------------------------------------
!***
!***  CALL CONVECTION
!***
!-----------------------------------------------------------------------
!***  Begin debugging convection
!         PRINT_DIAG=.FALSE.
!         IF(I==IMD.AND.J==JMD)PRINT_DIAG=.TRUE.
!***  End debugging convection
!-----------------------------------------------------------------------
          CALL BMJ_DEV(ntsd,I,J,DTCNVC,LMH,seamask,CLDEFI(I,J)          &
     &            ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs            &
     &            ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss            &
     &            ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts            &
     &            ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP                &
     &            ,DPCOL,PCOL,QCOL,TCOL,PSFC,PTOP                       &
     &            ,DQDT,DTDT,PCPCOL,LBOT,LTOP,LPBL                      &
     &            ,CP,R,ELWV,ELIV,G,TFRZ,D608                           &   
     &            ,PRINT_DIAG                                           &   
     &            ,IDS,IDE,JDS,JDE,KDS,KDE                              &     
     &            ,IMS,IME,JMS,JME,KMS,KME                              &
     &            ,ITS,ITE,JTS,JTE,KTS,KTE)
!-----------------------------------------------------------------------
! 
!***  COMPUTE HEATING AND MOISTENING TENDENCIES
!
          IF(PRESENT(RTHCUTEN).AND.PRESENT(RQVCUTEN))THEN
            DO K=KTS,KTE
              RTHCUTEN(I,J,K)=DTDT(K)/exner(I,J,K)
!
!***  CONVERT FROM SPECIFIC HUMIDTY BACK TO MIXING RATIO
!
              qold=qcol(k)
              qnew=qold+dqdt(k)
!
              rqvcuten(i,j,k)=qnew/(1.-qnew)-qold/(1.-qold)
!
!looks like bug              RQVCUTEN(I,J,K)=DQDT(K)/(1.-QCOL(K))**2
!
            ENDDO
          ENDIF
!
!***  ALL UNITS IN BMJ SCHEME ARE MKS, THUS CONVERT PRECIP FROM METERS
!***  TO MILLIMETERS PER STEP FOR OUTPUT.
!
          RAINCV(I,J)=PCPCOL*1.E3/NCNVC
!
!***  CONVECTIVE CLOUD TOP AND BOTTOM FROM THIS CALL
!
          CUTOP(I,J)=REAL(KTE+1-LTOP)
          CUBOT(I,J)=REAL(KTE+1-LBOT)
!
!-----------------------------------------------------------------------
!***  Begin debugging convection
          IF(PRINT_DIAG)THEN
            DELT=0.
            DELQ=0.
            PLYR=0.
            IF(LBOT>0.AND.LTOP<LBOT)THEN
              DO K=LTOP,LBOT
                PLYR=PLYR+DPCOL(K)
                DELQ=DELQ+DPCOL(K)*DTCNVC*ABS(DQDT(K))
                DELT=DELT+DPCOL(K)*DTCNVC*ABS(DTDT(K))
              ENDDO
              DELQ=DELQ/PLYR
              DELT=DELT/PLYR
            ENDIF
!
            WRITE(6,"(2a,2i4,3e12.4,f7.2,4i3)") &
                 '{cu3 i,j,PCPCOL,DTavg,DQavg,PLYR,'  &
                 ,'LBOT,LTOP,CUBOT,CUTOP = '  &
                 ,i,j, PCPCOL,DELT,1000.*DELQ,.01*PLYR  &
                 ,LBOT,LTOP,NINT(CUBOT(I,J)),NINT(CUTOP(I,J))
!
            IF(PLYR> 0.)THEN
              DO K=LBOT,LTOP,-1
                WRITE(6,"(a,i3,2e12.4,f7.2)") '{cu3a K,DT,DQ,DP = ' &
                     ,K,DTCNVC*DTDT(K),1000.*DTCNVC*DQDT(K),.01*DPCOL(K)
              ENDDO
            ENDIF
          ENDIF
!***  End debugging convection
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      ENDIF
!
      END SUBROUTINE BMJDRV_DEV
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
                          SUBROUTINE BMJ_DEV                            &
!-----------------------------------------------------------------------
     & (ntsd,I,J,DTCNVC,LMH,SM,CLDEFI                                   &
     & ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs                       &
     & ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss                       &
     & ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts                       &
     & ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP                           &
     & ,DPRS,PRSMID,Q,T,PSFC,PT                                         &
     & ,DQDT,DTDT,PCPCOL,LBOT,LTOP,LPBL                                 &
     & ,CP,R,ELWV,ELIV,G,TFRZ,D608                                      &
     & ,PRINT_DIAG                                                      &   
     & ,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     & ,IMS,IME,JMS,JME,KMS,KME                                         &
     & ,ITS,ITE,JTS,JTE,KTS,KTE)
!-----------------------------------------------------------------------
!zj  new shallow cloud added in June 2008 to address swap point
!zj  convection and shallow convection transporting both heat and moisture
!zj  up.  'soft' version of the cloud is implemented.
!zj  reintroduced entrainment on 11/02/2008
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
      LOGICAL(kind=klog),INTENT(IN):: &
       ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
                           ,IMS,IME,JMS,JME,KMS,KME                     &
                           ,ITS,ITE,JTS,JTE,KTS,KTE                     &
                           ,I,J,ntsd
! 
      INTEGER,INTENT(IN) :: LMH,LPBL
!
      INTEGER,INTENT(OUT) :: LBOT,LTOP
!
      REAL(kind=kfpt),INTENT(IN):: &
       CP,D608,DTCNVC,ELIV,ELWV,G,PSFC,PT,R,SM,TFRZ &
      ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
      ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
      ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: DPRS,PRSMID,Q,T
!
      REAL,INTENT(INOUT) :: CLDEFI,PCPCOL
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: DQDT,DTDT
!
!-----------------------------------------------------------------------
!***  DEFINE LOCAL VARIABLES
!-----------------------------------------------------------------------
!                                                            
      REAL,DIMENSION(KTS:KTE) :: APEK,APESK,EL,FPK                      &
                                ,PK,PSK,qbte,qbtk,QK,QREFK,QSATK        &
                                ,THERK,thesp,THEVRF,THSK                &
                                ,THVMOD,THVREF,TK,TREFK
!
      REAL,DIMENSION(KTS:KTE) :: APE,DIFQ,DIFT,THEE,THES,TREF
!
      REAL,DIMENSION(KTS:KTE) :: CPE,CPEcnv,DTV,DTVcnv,THEScnv    !<-- CPE for shallow convection buoyancy check (24 Aug 2006)
!
      real(kind=kfpt),dimension(kts:kte):: &
       dpk,rhk,thmak,thvmk
!
      LOGICAL(kind=klog):: &
       DEEP,SHALLOW
!
!***  Begin debugging convection
      LOGICAL :: PRINT_DIAG
!***  End debugging convection
!
!-----------------------------------------------------------------------
!***
!***  LOCAL SCALARS
!***
!-----------------------------------------------------------------------
      REAL :: APEKL,APEKXX,APEKXY,APES,APESTS                           &
     &            ,AVRGT,AVRGTL,BQ,BQK,BQS00K,BQS10K                    &
     &            ,CAPA,CUP,DEN,DENTPY,DEPMIN,DEPTH                     &
     &            ,DEPWL,DHDT,DIFQL,DIFTL,DP,DPKL,DPLO,DPMIX,DPOT       &
     &            ,DPUP,DQREF,DRHDP,DRHEAT,DSP                          &
     &            ,DSP0,DSP0K,DSPB,DSPBK,DSPT,DSPTK                     &
     &            ,DSQ,DST,DSTQ,DTHEM,DTDP,EFI                          &
     &            ,FEFI,FFUP,FPRS,FPTK,fup,HCORR                        &
     &            ,OTSUM,P,P00K,P01K,P10K,P11K                          &
     &            ,PART1,PART2,PART3,PBOT,PBOTFC,PBTK                   &
     &            ,PK0,PKB,PKL,PKT,PKXXXX,PKXXXY                        &
     &            ,PLMH,PELEVFC,PBTmx,plo,POTSUM,PP1,PPK,PRECK          &
     &            ,PRESK,PSP,PSUM,PTHRS,PTOP,PTPK,PUP                   &
     &            ,QBT,QKL,QNEW,QOTSUM,QQ1,QQK,QRFKL                    &
     &            ,QRFTP,QSP,QSUM,QUP,RDP0T                             &
     &            ,RDPSUM,RDTCNVC,RHH,RHL,RHMAX,ROTSUM,RTBAR,RHAVG      &
     &            ,SM1,SMIX,SQ,SQK,SQS00K,SQS10K,STABDL,SUMDE,SUMDP     &
     &            ,SUMDT,TAUK,TAUKSC,TCORR,THBT,THERKX,THERKY           &
     &            ,THSKL,THTPK,THVMKL,TKL,TNEW                          &
     &            ,TQ,TQK,TREFKX,TRFKL,trmlo,trmup,TSKL,tsp,TTH         &
     &            ,TTHK,TUP                                             &
     &            ,CAPEcnv,PSPcnv,THBTcnv,CAPEtrigr,CAPE                &
     &            ,TLEV2,QSAT1,QSAT2,RHSHmax
!
      INTEGER :: IQ,IQTB,iswap,IT,ITER,ITREFI,ITTB,ITTBK                &
     &          ,KB,KNUMH,KNUML &
     &          ,L,L0,L0M1,LB,LBM1,LCOR,LPT1                            &
     &          ,LQM,LSHU,LTP1,LTP2,LTSH, LBOTcnv,LTOPcnv,LMID
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: ELEVFC=0.6,STEFI=1.
!
      REAL,PARAMETER :: SLOPST=(STABDF-STABDS)/(1.-EFIMN)               &
     &                 ,slopel=(1.-efmntl)/(1.-efimn)                   &
     &                 ,slopes=(1.-efmnts)/(1.-efimn)
!
      real,parameter:: &
!       wdry=0.50 &
       wdry=0.75 &
!       wdry=0.90 &
      ,deftop=.95 !soft2
!zj      ,deftop=1.0 ! hard
!
      REAL :: A23M4L,CPRLG,ELOCP,RCP,QWAT &
      ,a11,a12,a21,a22,ama,aqs,arh,avm &
      ,b1qsat,b1rh,b1thma,b1thvm,b2qsat,b2rh,b2thma,b2thvm &
      ,bma,bqs,brh,bvm &
      ,qcorr,rden,rhenv,rhmean,rhref,sumdq,sumrh,wcld &
      ,adef,fk
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      CAPA=R/CP
      CPRLG=CP/(ROW*G*ELWV)
      ELOCP=ELIWV/CP
      RCP=1./CP
      A23M4L=A2*(A3-A4)*ELWV
      RDTCNVC=1./DTCNVC
      DEPMIN=PSH*PSFC*RSFCP
!--shallow convection switches------------------------------------------
      DEEP=.FALSE.
      SHALLOW=.FALSE.
!-----------------------------------------------------------------------
      DSP0=0.
      DSPB=0.
      DSPT=0.
!-----------------------------------------------------------------------
      TAUK  =DTCNVC/(TREL*01.0)
      TAUKSC=DTCNVC/(TREL*01.0)
!-----------------------------------------------------------------------
!-----------------------------PREPARATIONS------------------------------
!-----------------------------------------------------------------------
      iswap=0
      LBOT=LMH
!
      cup=0.
      dsp0=0.
      dspb=0.
      dspt=0.
      pcpcol=0.
      depmin=psh*psfc*rsfcp
      TREF(KTS)=T(KTS)
!
      DO L=KTS,LMH
        tk     (l)=t(l)
        qk     (l)=q(l)
        ape    (l)=(1.e5/prsmid(l))**capa
        cpecnv (l)=0.
        dtvcnv (l)=0.
        thes   (l)=0.
        thesp  (l)=0.
        thescnv(l)=0.
      ENDDO
!
!-----------------------------------------------------------------------
!----------------SEARCH FOR MAXIMUM BUOYANCY LEVEL----------------------
!-----------------------------------------------------------------------
!
      pelevfc=prsmid(lmh)*elevfc
      pbtmx  =prsmid(lmh)-pone
      capecnv=0.
      pspcnv =0.
      thbtcnv=0.
      lbotcnv=lbot
      ltopcnv=lbot
!
!-----------------------------------------------------------------------
!----------------TRIAL MAXIMUM BUOYANCY LEVEL VARIABLES-----------------
!-----------------------------------------------------------------------
!
      prep_loop: do kb=lmh,1,-1
!-----------------------------------------------------------------------
        if(prsmid(kb).lt.pelevfc.and..not.entrain) exit
!---preparation for search for max cape---------------------------------
        qbt=q(kb)
        thbt=t(kb)*ape(kb)
        tth=(thbt-thl)*rdth
        qq1=tth-aint(tth)
        ittb=int(tth)+1
!---keeping indices within the table------------------------------------
        if(ittb.lt.1)then
          ittb=1
          qq1=0.
        else if(ittb.ge.jtb)then
          ittb=jtb-1
          qq1=0.
        endif
!---base and scaling factor for spec. humidity--------------------------
        ittbk=ittb
        bqs00k=qs0(ittbk)
        sqs00k=sqs(ittbk)
        bqs10k=qs0(ittbk+1)
        sqs10k=sqs(ittbk+1)
!--------------scaling spec. humidity & table index---------------------
        bq=(bqs10k-bqs00k)*qq1+bqs00k
        sq=(sqs10k-sqs00k)*qq1+sqs00k
        tq=(qbt-bq)/sq*rdq
        pp1=tq-aint(tq)
        iqtb=int(tq)+1
!----------------keeping indices within the table-----------------------
        if(iqtb.lt.1)then
          iqtb=1
          pp1=0.
        else if(iqtb.ge.itb)then
          iqtb=itb-1
          pp1=0.
        endif
!--------------saturation pressure at four surrounding table pts.-------
        iq=iqtb
        it=ittb
        p00k=ptbl(iq  ,it  )
        p10k=ptbl(iq+1,it  )
        p01k=ptbl(iq  ,it+1)
        p11k=ptbl(iq+1,it+1)
!--------------saturation point variables at the bottom-----------------
        psp=p00k+(p10k-p00k)*pp1+(p01k-p00k)*qq1 &
     &     +(p00k-p10k-p01k+p11k)*pp1*qq1
        apes=(1.e5/psp)**capa
        thesp(kb)=thbt*exp(elocp*qbt*apes/thbt)
        psk  (kb)=psp
        apesk(kb)=apes
!-----------------------------------------------------------------------
      enddo prep_loop
!-----------------------------------------------------------------------
      max_buoy_loop: do kb=lmh,1,-1
!---choose cloud base as model level just below psp---------------------
!-----------------------------------------------------------------------
        if(prsmid(kb).lt.pelevfc) exit
!
        lbot=lmh
        ltop=lmh
!---search over a scaled depth to find the parcel with the max cape-----
        qbt=q(kb)
        thbt=t(kb)*ape(kb)
        psp=psk(kb)
        apes=apesk(kb)
!
        do l=kts,lmh-1
          p=prsmid(l)
          if(p.lt.psp.and.p.ge.pqm) lbot=l+1
        enddo
!***
!*** warning: lbot must not be .gt. lm-1 in shallow convection
!*** make sure cloud base is at least pone above the surface
!***
        pbot=prsmid(lbot)
        if(pbot.ge.pbtmx.or.lbot.ge.lmh)then
          do l=kts,lmh-1
            p=prsmid(l)
            if(p.lt.pbtmx)lbot=l
          enddo
          pbot=prsmid(lbot)
        endif
!---cloud top computation-----------------------------------------------
        ltop=lbot
        ptop=pbot
!---entrainment during parcel ascent------------------------------------
        if(entrain) then
!-----------------------------------------------------------------------
          do l=lmh,kb,-1
            thes(l)=thesp(kb)
            qbtk(l)=qk   (kb)
          enddo
          do l=kb,1,-1
            thes(l)=thesp(l)
            qbtk(l)=qk   (l)
          enddo
!
          do l=kts,lmh
            thee(l)=thes(l)
            qbte(l)=qbtk(l)
          enddo
!
          if(sm.gt.0.5) then
            fefi=(cldefi-efimn)*slopes+efmnts
          else
            fefi=(cldefi-efimn)*slopel+efmntl
          endif
!
          ffup=fup/(fefi*fefi)
!
          if(pbot.gt.enplo)then
            fprs=1.
          elseif(pbot.gt.enpup)then
            fprs=(pbot-enpup)*rendp
          else
            fprs=0.
          endif
!
          ffup=ffup*fprs*fprs*0.5
          dpup=dprs(kb)
!
          do l=kb-1,1,-1
            dplo=dpup
            dpup=dprs(l)
!
            thes(l)=((-ffup*dplo+1.)*thes(l+1) &
                    +(thee(l)*dpup+thee(l+1)*dplo)*ffup) &
                   /(ffup*dpup+1.)
            qbtk(l)=((-ffup*dplo+1.)*qbtk(l+1) &
                    +(qbte(l)*dpup+qbte(l+1)*dplo)*ffup) &
                   /(ffup*dpup+1.)
          enddo
!---no entrainment------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do l=lmh,1,-1
            thes(l)=thesp(kb)
            qbtk(l)=qk   (kb)
          enddo
!-----------------------------------------------------------------------
        endif ! end of entrainment/no entrainment
!------------------FIRST ENTROPY CHECK----------------------------------
!
        DO L=kts,lmh
          CPE(L)=0.
          DTV(L)=0.
        ENDDO
!-----------------------------------------------------------------------
!-- Begin: Buoyancy check including deep convection (24 Aug 2006)
!-----------------------------------------------------------------------
          DENTPY=0.
          L=KB
          PLO=PRSMID(L)
          TRMLO=0.
          CAPEtrigr=DTPtrigr/T(LBOT)
!
!--- Below cloud base
!
          IF(KB>LBOT) THEN
            DO L=KB-1,LBOT+1,-1
              PUP=PRSMID(L)
              TUP=THBT/APE(L)
              DP=PLO-PUP
              TRMUP=(TUP*(QBT*0.608+1.)                                 &
     &            -T(L)*(Q(L)*0.608+1.))*0.5                            &
     &             /(T(L)*(Q(L)*0.608+1.))
              DTV(L)=TRMLO+TRMUP
              DENTPY=DTV(L)*DP+DENTPY
              CPE(L)=DENTPY
              IF (DENTPY < CAPEtrigr) GO TO 170
              PLO=PUP
              TRMLO=TRMUP
            ENDDO
          ELSE
            L=LBOT+1
            PLO=PRSMID(L)
            TUP=THBT/APE(L)
            TRMLO=(TUP*(QBT*0.608+1.)                                   &
     &            -T(L)*(Q(L)*0.608+1.))*0.5                            &
     &             /(T(L)*(Q(L)*0.608+1.))
          ENDIF  ! IF(KB>LBOT) THEN
!
!--- At cloud base
!
          L=LBOT
          PUP=PSP
          TUP=THBT/APES
          TSP=(T(L+1)-T(L))/(PLO-PBOT)                                  &
     &       *(PUP-PBOT)+T(L)
          QSP=(Q(L+1)-Q(L))/(PLO-PBOT)                                  &
     &       *(PUP-PBOT)+Q(L)
          DP=PLO-PUP
          TRMUP=(TUP*(QBT*0.608+1.)                                     &
     &          -TSP*(QSP*0.608+1.))*0.5                                &
     &         /(TSP*(QSP*0.608+1.))
          DTV(L)=TRMLO+TRMUP
          DENTPY=DTV(L)*DP+DENTPY
          CPE(L)=DENTPY
          DTV(L)=DTV(L)*DP
          PLO=PUP
          TRMLO=TRMUP
          PUP=PRSMID(L)
!
!--- Calculate updraft temperature along moist adiabat (TUP)
!
          IF(PUP<PLQ)THEN
            CALL TTBLEX(ITB,JTB,PL,PUP,RDP,RDTHE                        &
     &                 ,STHE,THE0,THES(L),TTBL,TUP)
          ELSE
            CALL TTBLEX(ITBQ,JTBQ,PLQ,PUP,RDPQ,RDTHEQ                   &
     &                 ,STHEQ,THE0Q,THES(L),TTBLQ,TUP)
          ENDIF
!
          QUP=PQ0/PUP*EXP(A2*(TUP-A3)/(TUP-A4))
          QWAT=QBT-QUP  !-- Water loading effects, reversible adiabat
          DP=PLO-PUP
          TRMUP=(TUP*(QUP*0.608+1.-QWAT)                                &
     &          -T(L)*(Q(L)*0.608+1.))*0.5                              &
     &         /(T(L)*(Q(L)*0.608+1.))
          DENTPY=(TRMLO+TRMUP)*DP+DENTPY
          CPE(L)=DENTPY
          DTV(L)=(DTV(L)+(TRMLO+TRMUP)*DP)/(PRSMID(LBOT+1)-PRSMID(LBOT))
!
          IF (DENTPY < CAPEtrigr) GO TO 170
!
          PLO=PUP
          TRMLO=TRMUP
!
!-----------------------------------------------------------------------
!--- In cloud above cloud base
!-----------------------------------------------------------------------
!
          DO L=LBOT-1,KTS,-1
            PUP=PRSMID(L)
!
!--- Calculate updraft temperature along moist adiabat (TUP)
!
            IF(PUP<PLQ)THEN
              CALL TTBLEX(ITB,JTB,PL,PUP,RDP,RDTHE                      &
       &                 ,STHE,THE0,THES(L),TTBL,TUP)
            ELSE
              CALL TTBLEX(ITBQ,JTBQ,PLQ,PUP,RDPQ,RDTHEQ                 &
       &                 ,STHEQ,THE0Q,THES(L),TTBLQ,TUP)
            ENDIF
!
            QUP=PQ0/PUP*EXP(A2*(TUP-A3)/(TUP-A4))
            QWAT=QBT-QUP  !-- Water loading effects, reversible adiabat
            DP=PLO-PUP
            TRMUP=(TUP*(QUP*0.608+1.-QWAT)                              &
     &            -T(L)*(Q(L)*0.608+1.))*0.5                            &
     &           /(T(L)*(Q(L)*0.608+1.))
            DTV(L)=TRMLO+TRMUP
            DENTPY=DTV(L)*DP+DENTPY
            CPE(L)=DENTPY
!
            IF (DENTPY < CAPEtrigr) GO TO 170
!
            PLO=PUP
            TRMLO=TRMUP
          ENDDO
!
!-----------------------------------------------------------------------
!
170       LTP1=KB
          CAPE=0.
!
!-----------------------------------------------------------------------
!--- Cloud top level (LTOP) is located where CAPE is a maximum
!--- Exit cloud-top search if CAPE < CAPEtrigr
!-----------------------------------------------------------------------
!
          DO L=KB,KTS,-1
            IF (CPE(L) < CAPEtrigr) THEN
              EXIT
            ELSE IF (CPE(L) > CAPE) THEN
              LTP1=L
              CAPE=CPE(L)
            ENDIF
          ENDDO      !-- End DO L=KB,KTS,-1
!
          LTOP=MIN(LTP1,LBOT)
!
!-----------------------------------------------------------------------
!--------------- CHECK FOR MAXIMUM INSTABILITY  ------------------------
!-----------------------------------------------------------------------
          IF(CAPE > CAPEcnv) THEN
            CAPEcnv=CAPE
            PSPcnv=PSP
            THBTcnv=THBT
            LBOTcnv=LBOT
            LTOPcnv=LTOP
            DO L=LMH,KTS,-1
              CPEcnv(L)=CPE(L)
              DTVcnv(L)=DTV(L)
              THEScnv(L)=THES(L)
            ENDDO
          ENDIF    ! End IF(CAPE > CAPEcnv) THEN
!
!-----------------------------------------------------------------------
!
      ENDDO max_buoy_loop
!
!-----------------------------------------------------------------------
!------------------------  MAXIMUM INSTABILITY  ------------------------
!-----------------------------------------------------------------------
!
      IF(CAPEcnv > 0.) THEN
        PSP=PSPcnv
        THBT=THBTcnv
        LBOT=LBOTcnv
        LTOP=LTOPcnv
        PBOT=PRSMID(LBOT)
        PTOP=PRSMID(LTOP)
!
        DO L=LMH,KTS,-1
          CPE(L)=CPEcnv(L)
          DTV(L)=DTVcnv(L)
          THES(L)=THEScnv(L)
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----  Quick exit if cloud is too thin or no CAPE is present  ---------
!-----------------------------------------------------------------------
!

      if(ptop>pbot-pno.or.ltop>lbot-2.or.capeCNV<=0.)then
        cldefi=avgefi*sm+stefi*(1.-sm)
        ltop=kte
      endif

      IF(PTOP>PBOT-PNO.OR.LTOP>LBOT-2.OR.CAPEcnv<=0.)THEN
        LBOT=0
        PBOT=PRSMID(LMH)
        PTOP=PBOT
        return
      ENDIF
!
!***  DEPTH OF CLOUD REQUIRED TO MAKE THE POINT A DEEP CONVECTION POINT
!***  IS A SCALED VALUE OF PSFC.
!
      DEPTH=PBOT-PTOP
!
      IF(DEPTH>=DEPMIN) THEN
        DEEP=.TRUE.
      ELSE
        SHALLOW=.TRUE.
        GO TO 600
      ENDIF
!
!-----------------------------------------------------------------------
!DCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCD
!DCDCDCDCDCDCDCDCDCDCDC    DEEP CONVECTION   DCDCDCDCDCDCDCDCDCDCDCDCDCD
!DCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCD
!-----------------------------------------------------------------------
!
  300 CONTINUE
!
      LB =LBOT
      EFI=CLDEFI
!-----------------------------------------------------------------------
!--------------INITIALIZE VARIABLES IN THE CONVECTIVE COLUMN------------
!-----------------------------------------------------------------------
!***
!***  ONE SHOULD NOTE THAT THE VALUES ASSIGNED TO THE ARRAY TREFK
!***  IN THE FOLLOWING LOOP ARE REALLY ONLY RELEVANT IN ANCHORING THE
!***  REFERENCE TEMPERATURE PROFILE AT LEVEL LB.  WHEN BUILDING THE
!***  REFERENCE PROFILE FROM CLOUD BASE, THEN ASSIGNING THE
!***  AMBIENT TEMPERATURE TO TREFK IS ACCEPTABLE.  HOWEVER, WHEN
!***  BUILDING THE REFERENCE PROFILE FROM SOME OTHER LEVEL (SUCH AS
!***  ONE LEVEL ABOVE THE GROUND), THEN TREFK SHOULD BE FILLED WITH
!***  THE TEMPERATURES IN TREF(L) WHICH ARE THE TEMPERATURES OF
!***  THE MOIST ADIABAT THROUGH CLOUD BASE.  BY THE TIME THE LINE
!***  NUMBERED 450 HAS BEEN REACHED, TREFK ACTUALLY DOES HOLD THE
!***  REFERENCE TEMPERATURE PROFILE.
!***
      DO L=KTS,LMH
        DIFT  (L)=0.
        DIFQ  (L)=0.
        TKL      =T(L)
        TK    (L)=TKL
        TREFK (L)=TKL
        QKL      =Q(L)
        QK    (L)=QKL
        QREFK (L)=QKL
        PKL      =PRSMID(L)
        PK    (L)=PKL
        PSK   (L)=PKL
        APEKL    =APE(L)
        APEK  (L)=APEKL
!
!--- Calculate temperature along moist adiabat (TREF)
!
        IF(PKL<PLQ)THEN
          CALL TTBLEX(ITB,JTB,PL,PKL,RDP,RDTHE                          &
     &               ,STHE,THE0,THES(L),TTBL,TREF(L))
        ELSE
          CALL TTBLEX(ITBQ,JTBQ,PLQ,PKL,RDPQ,RDTHEQ                     &
     &               ,STHEQ,THE0Q,THES(L),TTBLQ,TREF(L))
        ENDIF
        THERK (L)=TREF(L)*APEKL
      ENDDO
!
!------------DEEP CONVECTION REFERENCE TEMPERATURE PROFILE------------
!
      LTP1=LTOP+1
      LBM1=LB-1
      PKB=PK(LB)
      PKT=PK(LTOP)
      STABDL=(EFI-EFIMN)*SLOPST+STABDS
!
!------------TEMPERATURE REFERENCE PROFILE BELOW FREEZING LEVEL-------
!
      EL(LB) = ELWV
      L0=LB
      PK0=PK(LB)
      TREFKX=TREFK(LB)
      THERKX=THERK(LB)
      APEKXX=APEK(LB)
      THERKY=THERK(LBM1)
      APEKXY=APEK(LBM1)
!
      DO L=LBM1,LTOP,-1
        IF(T(L+1)<TFRZ)GO TO 430
        TREFKX=((THERKY-THERKX)*STABDL                                  &
     &          +TREFKX*APEKXX)/APEKXY
        TREFK(L)=TREFKX
        EL(L)=ELWV
        APEKXX=APEKXY
        THERKX=THERKY
        APEKXY=APEK(L-1)
        THERKY=THERK(L-1)
        L0=L
        PK0=PK(L0)
      ENDDO
!
!--------------FREEZING LEVEL AT OR ABOVE THE CLOUD TOP-----------------
!
      GO TO 450
!
!--------------TEMPERATURE REFERENCE PROFILE ABOVE FREEZING LEVEL-------
!
  430 L0M1=L0-1
      RDP0T=1./(PK0-PKT)
      DTHEM=THERK(L0)-TREFK(L0)*APEK(L0)
!
      DO L=LTOP,L0M1
        TREFK(L)=(THERK(L)-(PK(L)-PKT)*DTHEM*RDP0T)/APEK(L)
        EL(L)=ELWV !ELIV
      ENDDO
!
!-----------------------------------------------------------------------
!--------------DEEP CONVECTION REFERENCE HUMIDITY PROFILE---------------
!-----------------------------------------------------------------------
!
!***  DEPWL IS THE PRESSURE DIFFERENCE BETWEEN CLOUD BASE AND
!***  THE FREEZING LEVEL
!
  450 CONTINUE
      DEPWL=PKB-PK0
      DEPTH=PFRZ*PSFC*RSFCP
      SM1=1.-SM
      PBOTFC=1.
!
!-------------FIRST ADJUSTMENT OF TEMPERATURE PROFILE-------------------
!!
!      SUMDT=0.
!      SUMDP=0.
!!
!      DO L=LTOP,LB
!        SUMDT=(TK(L)-TREFK(L))*DPRS(L)+SUMDT
!        SUMDP=SUMDP+DPRS(L)
!      ENDDO
!!
!      TCORR=SUMDT/SUMDP
!!
!      DO L=LTOP,LB
!        TREFK(L)=TREFK(L)+TCORR
!      ENDDO
!!
!-----------------------------------------------------------------------
!--------------- ITERATION LOOP FOR CLOUD EFFICIENCY -------------------
!-----------------------------------------------------------------------
!
      cloud_efficiency : DO ITREFI=1,ITREFI_MAX
!
!-----------------------------------------------------------------------
        DSPBK=((EFI-EFIMN)*SLOPBS+DSPBSS*PBOTFC)*SM                     &
     &       +((EFI-EFIMN)*SLOPBL+DSPBSL*PBOTFC)*SM1
        DSP0K=((EFI-EFIMN)*SLOP0S+DSP0SS*PBOTFC)*SM                     &
     &       +((EFI-EFIMN)*SLOP0L+DSP0SL*PBOTFC)*SM1
        DSPTK=((EFI-EFIMN)*SLOPTS+DSPTSS*PBOTFC)*SM                     &
     &       +((EFI-EFIMN)*SLOPTL+DSPTSL*PBOTFC)*SM1
!
!-----------------------------------------------------------------------
!
        DO L=LTOP,LB
!
!***
!***  SATURATION PRESSURE DIFFERENCE
!***
          IF(DEPWL>=DEPTH)THEN
            IF(L<L0)THEN
              DSP=((PK0-PK(L))*DSPTK+(PK(L)-PKT)*DSP0K)/(PK0-PKT)
            ELSE
              DSP=((PKB-PK(L))*DSP0K+(PK(L)-PK0)*DSPBK)/(PKB-PK0)
            ENDIF
          ELSE
            DSP=DSP0K
            IF(L<L0)THEN
              DSP=((PK0-PK(L))*DSPTK+(PK(L)-PKT)*DSP0K)/(PK0-PKT)
            ENDIF
          ENDIF
!***
!***  HUMIDITY PROFILE
!***
          IF(PK(L)>PQM)THEN
            PSK(L)=PK(L)+DSP
            APESK(L)=(1.E5/PSK(L))**CAPA
            THSK(L)=TREFK(L)*APEK(L)
            QREFK(L)=PQ0/PSK(L)*EXP(A2*(THSK(L)-A3*APESK(L))            &
     &                                /(THSK(L)-A4*APESK(L)))
          ELSE
            QREFK(L)=QK(L)
          ENDIF
!
        ENDDO
!-----------------------------------------------------------------------
!***
!***  ENTHALPY CONSERVATION INTEGRAL
!***
!-----------------------------------------------------------------------
        enthalpy_conservation : DO ITER=1,2
!
          SUMDE=0.
          SUMDP=0.
!
          DO L=LTOP,LB
            SUMDE=((TK(L)-TREFK(L))*CP+(QK(L)-QREFK(L))*EL(L))*DPRS(L)  &
     &            +SUMDE
            SUMDP=SUMDP+DPRS(L)
          ENDDO
!
          HCORR=SUMDE/(SUMDP-DPRS(LTOP))
          LCOR=LTOP+1
!***
!***  FIND LQM
!***
          LQM=1
          DO L=KTS,LB
            IF(PK(L)<=PQM)LQM=L
          ENDDO
!***
!***  ABOVE LQM CORRECT TEMPERATURE ONLY
!***
          IF(LCOR<=LQM)THEN
            DO L=LCOR,LQM
              TREFK(L)=TREFK(L)+HCORR*RCP
            ENDDO
            LCOR=LQM+1
          ENDIF
!***
!***  BELOW LQM CORRECT BOTH TEMPERATURE AND MOISTURE
!***
          DO L=LCOR,LB
            TSKL=TREFK(L)*APEK(L)/APESK(L)
            DHDT=QREFK(L)*A23M4L/(TSKL-A4)**2+CP
            TREFK(L)=HCORR/DHDT+TREFK(L)
            THSKL=TREFK(L)*APEK(L)
            QREFK(L)=PQ0/PSK(L)*EXP(A2*(THSKL-A3*APESK(L))              &
     &                                /(THSKL-A4*APESK(L)))
          ENDDO
!
        ENDDO  enthalpy_conservation
!-----------------------------------------------------------------------
!
!***  HEATING, MOISTENING, PRECIPITATION
!
!-----------------------------------------------------------------------
        AVRGT=0.
        PRECK=0.
        DSQ=0.
        DST=0.
!
        DO L=LTOP,LB
          TKL=TK(L)
          DIFTL=(TREFK(L)-TKL  )*TAUK
          DIFQL=(QREFK(L)-QK(L))*TAUK
          AVRGTL=(TKL+TKL+DIFTL)
          DPOT=DPRS(L)/AVRGTL
          DST=DIFTL*DPOT+DST
          DSQ=DIFQL*EL(L)*DPOT+DSQ
          AVRGT=AVRGTL*DPRS(L)+AVRGT
          PRECK=DIFTL*DPRS(L)+PRECK
          DIFT(L)=DIFTL
          DIFQ(L)=DIFQL
        ENDDO
!
        DST=(DST+DST)*CP
        DSQ=DSQ+DSQ
        DENTPY=DST+DSQ
        AVRGT=AVRGT/(SUMDP+SUMDP)
!
!        DRHEAT=PRECK*CP/AVRGT
        DRHEAT=(PRECK*SM+MAX(1.E-7,PRECK)*(1.-SM))*CP/AVRGT !As in Eta!
        DRHEAT=MAX(DRHEAT,1.E-20)
        EFI=EFIFC*DENTPY/DRHEAT
!-----------------------------------------------------------------------
        EFI=MIN(EFI,1.)
        EFI=MAX(EFI,EFIMN)
!-----------------------------------------------------------------------
!
      ENDDO  cloud_efficiency
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---------------------- DEEP CONVECTION --------------------------------
!-----------------------------------------------------------------------
!
      IF(DENTPY>=EPSNTP.AND.PRECK>EPSPR.and..not.nodeep) THEN
!
        iswap=0 ! deep convection, no swap
        CLDEFI=EFI
!
        if(sm.gt.0.5) then
          fefi=(cldefi-efimn)*slopes+efmnts
        else
          fefi=(cldefi-efimn)*slopel+efmntl
        endif
!
        FEFI=(DENTPY-EPSNTP)*FEFI/DENTPY
        PRECK=PRECK*FEFI
!
!***  UPDATE PRECIPITATION AND TENDENCIES OF TEMPERATURE AND MOISTURE
!
        CUP=PRECK*CPRLG
        PCPCOL=CUP
!
        DO L=LTOP,LB
          DTDT(L)=DIFT(L)*FEFI*RDTCNVC
          DQDT(L)=DIFQ(L)*FEFI*RDTCNVC
        ENDDO
!
      ELSE
!
!-----------------------------------------------------------------------
!***  REDUCE THE CLOUD TOP
!-----------------------------------------------------------------------
!
!        LTOP=LTOP+3           !iterate cloud top
!        PTOP=PRSMID(LTOP)     !iterate cloud top
!        DEPMIN=PSH*PSFC*RSFCP !iterate cloud top
!        DEPTH=PBOT-PTOP       !iterate cloud top
!***
!***  ITERATE DEEP CONVECTION PROCEDURE IF NEEDED
!***
!        IF(DEPTH>=DEPMIN)THEN !iterate cloud top
!          GO TO 300           !iterate cloud top
!        ENDIF                 !iterate cloud top
!
!        CLDEFI=AVGEFI
         CLDEFI=EFIMN*SM+STEFI*(1.-SM)
!***
!***  SEARCH FOR SHALLOW CLOUD TOP
!***
!        LTSH=LBOT
!        LBM1=LBOT-1
!        PBTK=PK(LBOT)
!        DEPMIN=PSH*PSFC*RSFCP
!        PTPK=PBTK-DEPMIN
        PTPK=MAX(PSHU, PK(LBOT)-DEPMIN)
!***
!***  CLOUD TOP IS THE LEVEL JUST BELOW PBTK-PSH or JUST BELOW PSHU
!***
        DO L=KTS,LMH
          IF(PK(L)<=PTPK)LTOP=L+1
        ENDDO
!
!        PTPK=PK(LTOP)
!!***
!!***  HIGHEST LEVEL ALLOWED IS LEVEL JUST BELOW PSHU
!!***
!        IF(PTPK<=PSHU)THEN
!!
!          DO L=KTS,LMH
!            IF(PK(L)<=PSHU)LSHU=L+1
!          ENDDO
!!
!          LTOP=LSHU
!          PTPK=PK(LTOP)
!        ENDIF
!
!        if(ltop>=lbot)then
!!!!!!     lbot=0
!          ltop=lmh
!!!!!!     pbot=pk(lbot)
!          ptop=pk(ltop)
!          pbot=ptop
!          go to 600
!        endif
!
!        LTP1=LTOP+1
!        LTP2=LTOP+2
!!
!        DO L=LTOP,LBOT
!          QSATK(L)=PQ0/PK(L)*EXP(A2*(TK(L)-A3)/(TK(L)-A4))
!        ENDDO
!!
!        RHH=QK(LTOP)/QSATK(LTOP)
!        RHMAX=0.
!        LTSH=LTOP
!!
!        DO L=LTP1,LBM1
!          RHL=QK(L)/QSATK(L)
!          DRHDP=(RHH-RHL)/(PK(L-1)-PK(L))
!!
!          IF(DRHDP>RHMAX)THEN
!            LTSH=L-1
!            RHMAX=DRHDP
!          ENDIF
!!
!          RHH=RHL
!        ENDDO
!
!-----------------------------------------------------------------------
!-- Make shallow cloud top a function of virtual temperature excess (DTV)
!-----------------------------------------------------------------------
!
        LTP1=LBOT
        DO L=LBOT-1,LTOP,-1
          IF (DTV(L) > 0.) THEN
            LTP1=L
          ELSE
            EXIT
          ENDIF
        ENDDO
        LTOP=MIN(LTP1,LBOT)
!***
!***  CLOUD MUST BE AT LEAST TWO LAYERS THICK
!***
!        IF(LBOT-LTOP<2)LTOP=LBOT-2  (eliminate this criterion)
!
!-- End: Buoyancy check (24 Aug 2006)
!
        iswap=1 ! failed deep convection, shallow swap point
        PTOP=PK(LTOP)
        SHALLOW=.TRUE.
        DEEP=.FALSE.
!
      ENDIF
!DCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCD
!DCDCDCDCDCDCDC          END OF DEEP CONVECTION            DCDCDCDCDCDCD
!DCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCD
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  600 CONTINUE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------GATHER SHALLOW CONVECTION POINTS-----------------------
!
!      IF(PTOP<=PBOT-PNO.AND.LTOP<=LBOT-2)THEN
!         DEPMIN=PSH*PSFC*RSFCP
!!
!!        if(lpbl<lbot)lbot=lpbl
!!        if(lbot>lmh-1)lbot=lmh-1
!!        pbot=prsmid(lbot)
!!
!         IF(PTOP+1.>=PBOT-DEPMIN)SHALLOW=.TRUE.
!      ELSE
!         LBOT=0
!         LTOP=KTE
!      ENDIF
!
!***********************************************************************
!-----------------------------------------------------------------------
!***  Begin debugging convection
      IF(PRINT_DIAG)THEN
        WRITE(6,"(a,2i3,L2,3e12.4)")  &
             '{cu2a lbot,ltop,shallow,pbot,ptop,depmin = ' &
             ,lbot,ltop,shallow,pbot,ptop,depmin
      ENDIF
!***  End debugging convection
!-----------------------------------------------------------------------
!
      IF(.NOT.SHALLOW)GO TO 800
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS
!SCSCSCSCSCSCSC         SHALLOW CONVECTION          CSCSCSCSCSCSCSCSCSCS
!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS
!-----------------------------------------------------------------------
      do l=kts,lmh
        pk(l)=prsmid(l)
        tk(l)=t(l)
        qk(l)=q(l)
        trefk(l)=t(l)
        qrefk(l)=q(l)
        qsatk(l)=pq0/pk(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
        apek(l)=ape(l)
        thvref(l)=tk(l)*apek(l)*(qk(l)*d608+1.)
!
!        if(tk(l)>=tfrz)then
          el(l)=elwv
!        else
!          el(l)=eliv
!        endif
      enddo
!
!-----------------------------------------------------------------------
!-- Begin: Raise cloud top if avg RH>RHSHmax and CAPE>0
!   RHSHmax=RH at cloud base associated with a DSP of PONE
!-----------------------------------------------------------------------
!
      TLEV2=T(LBOT)*((PK(LBOT)-PONE)/PK(LBOT))**CAPA
      QSAT1=PQ0/PK(LBOT)*EXP(A2*(T(LBOT)-A3)/(TK(LBOT)-A4))
      QSAT2=PQ0/(PK(LBOT)-PONE)*EXP(A2*(TLEV2-A3)/(TLEV2-A4))
      RHSHmax=QSAT2/QSAT1
      SUMDP=0.
      RHAVG=0.
!
      DO L=LBOT,LTOP,-1
        RHAVG=RHAVG+DPRS(L)*QK(L)/QSATK(L)
        SUMDP=SUMDP+DPRS(L)
      ENDDO
!
      IF (RHAVG/SUMDP > RHSHmax) THEN
        LTSH=LTOP
        DO L=LTOP-1,KTS,-1
          RHAVG=RHAVG+DPRS(L)*QK(L)/QSATK(L)
          SUMDP=SUMDP+DPRS(L)
          IF (CPE(L) > 0.) THEN
            LTSH=L
          ELSE
            EXIT
          ENDIF
          IF (RHAVG/SUMDP <= RHSHmax) EXIT
          IF (PK(L) <= PSHU) EXIT
        ENDDO
        LTOP=LTSH
!swapnomoist        iswap=0 ! old cloud for moist clouds
      ENDIF
!
!-- End: Raise cloud top if avg RH>RHSHmax and CAPE>0
!
!---------------------------SHALLOW CLOUD TOP---------------------------
      LBM1=LBOT-1
      PTPK=PTOP
      LTP1=LTOP-1
      DEPTH=PBOT-PTOP
!-----------------------------------------------------------------------
!***  Begin debugging convection
      IF(PRINT_DIAG)THEN
        WRITE(6,"(a,4e12.4)") '{cu2b PBOT,PTOP,DEPTH,DEPMIN= ' &
             ,PBOT,PTOP,DEPTH,DEPMIN
      ENDIF
!***  End debugging convection
!-----------------------------------------------------------------------
!
!BSF      IF(DEPTH<DEPMIN)THEN
!BSF        GO TO 800
!BSF      ENDIF
!-----------------------------------------------------------------------
      IF(PTOP>PBOT-PNO.OR.LTOP>LBOT-2)THEN
        LBOT=0
        LTOP=KTE
        PTOP=PBOT
        GO TO 800
      ENDIF
!-----------------------------------------------------------------------
!***  New cloud at all shallow points
!-----------------------------------------------------------------------
      if(newall) go to 810 ! new cloud at all shallow points
!-----------------------------------------------------------------------
!***  New cloud only at swap shallow points
!-----------------------------------------------------------------------
      if(newswap.and.iswap.gt.0) go to 810 ! new cloud only at swap pts.
!-----------------------------------------------------------------------
!
!--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX AT TOP-------
!
      THTPK=T(LTP1)*APE(LTP1)
!
      TTHK=(THTPK-THL)*RDTH
      QQK =TTHK-AINT(TTHK)
      IT  =INT(TTHK)+1
!
      IF(IT<1)THEN
        IT=1
        QQK=0.
      ENDIF
!
      IF(IT>=JTB)THEN
        IT=JTB-1
        QQK=0.
      ENDIF
!
!--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY AT TOP--------
!
      BQS00K=QS0(IT)
      SQS00K=SQS(IT)
      BQS10K=QS0(IT+1)
      SQS10K=SQS(IT+1)
!
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX AT TOP--------------
!
      BQK=(BQS10K-BQS00K)*QQK+BQS00K
      SQK=(SQS10K-SQS00K)*QQK+SQS00K
!
!     TQK=(Q(LTOP)-BQK)/SQK*RDQ
      TQK=(Q(LTP1)-BQK)/SQK*RDQ
!
      PPK=TQK-AINT(TQK)
      IQ =INT(TQK)+1
!
      IF(IQ<1)THEN
        IQ=1
        PPK=0.
      ENDIF
!
      IF(IQ>=ITB)THEN
        IQ=ITB-1
        PPK=0.
      ENDIF
!
!----------------CLOUD TOP SATURATION POINT PRESSURE--------------------
!
      PART1=(PTBL(IQ+1,IT)-PTBL(IQ,IT))*PPK
      PART2=(PTBL(IQ,IT+1)-PTBL(IQ,IT))*QQK
      PART3=(PTBL(IQ  ,IT  )-PTBL(IQ+1,IT  )                            &
     &      -PTBL(IQ  ,IT+1)+PTBL(IQ+1,IT+1))*PPK*QQK
      PTPK=PTBL(IQ,IT)+PART1+PART2+PART3
!-----------------------------------------------------------------------
      DPMIX=PTPK-PSP
      IF(ABS(DPMIX).LT.3000.)DPMIX=-3000.
!
!----------------TEMPERATURE PROFILE SLOPE------------------------------
!
      SMIX=(THTPK-THBT)/DPMIX*STABS
!
      TREFKX=TREFK(LBOT+1)
      PKXXXX=PK(LBOT+1)
      PKXXXY=PK(LBOT)
      APEKXX=APEK(LBOT+1)
      APEKXY=APEK(LBOT)
!
      LMID=.5*(LBOT+LTOP)
!
      DO L=LBOT,LTOP,-1
        TREFKX=((PKXXXY-PKXXXX)*SMIX                                    &
     &          +TREFKX*APEKXX)/APEKXY
        TREFK(L)=TREFKX
        IF(L<=LMID) TREFK(L)=MAX(TREFK(L), TK(L)+DTSHAL)
        APEKXX=APEKXY
        PKXXXX=PKXXXY
        APEKXY=APEK(L-1)
        PKXXXY=PK(L-1)
      ENDDO
!
!----------------TEMPERATURE REFERENCE PROFILE CORRECTION---------------
!
      SUMDT=0.
      SUMDP=0.
!
      DO L=LTOP,LBOT
        SUMDT=(TK(L)-TREFK(L))*DPRS(L)+SUMDT
        SUMDP=SUMDP+DPRS(L)
      ENDDO
!
      RDPSUM=1./SUMDP
      FPK(LBOT)=TREFK(LBOT)
!
      TCORR=SUMDT*RDPSUM
!
      DO L=LTOP,LBOT
        TRFKL   =TREFK(L)+TCORR
        TREFK(L)=TRFKL
        FPK  (L)=TRFKL
      ENDDO
!
!----------------HUMIDITY PROFILE EQUATIONS-----------------------------
!
      PSUM  =0.
      QSUM  =0.
      POTSUM=0.
      QOTSUM=0.
      OTSUM =0.
      DST   =0.
      FPTK  =FPK(LTOP)
!
      DO L=LTOP,LBOT
        DPKL  =FPK(L)-FPTK
        PSUM  =DPKL *DPRS(L)+PSUM
        QSUM  =QK(L)*DPRS(L)+QSUM
        RTBAR =2./(TREFK(L)+TK(L))
        OTSUM =DPRS(L)*RTBAR+OTSUM
        POTSUM=DPKL    *RTBAR*DPRS(L)+POTSUM
        QOTSUM=QK(L)   *RTBAR*DPRS(L)+QOTSUM
        DST   =(TREFK(L)-TK(L))*RTBAR*DPRS(L)/EL(L)+DST
      ENDDO
!
      PSUM  =PSUM*RDPSUM
      QSUM  =QSUM*RDPSUM
      ROTSUM=1./OTSUM
      POTSUM=POTSUM*ROTSUM
      QOTSUM=QOTSUM*ROTSUM
      DST   =DST*ROTSUM*CP
!
!-----------------------------------------------------------------------
!***  Begin debugging convection
      IF(PRINT_DIAG)THEN
        WRITE(6,"(a,5e12.4)") '{cu2c DST,PSUM,QSUM,POTSUM,QOTSUM = ' &
             ,DST,PSUM,QSUM,POTSUM,QOTSUM
      ENDIF
!***  End debugging convection
!-----------------------------------------------------------------------
!***  If upward transport of temperature go to new cloud
!-----------------------------------------------------------------------
      if(newupup.and.dst.gt.0.) go to 810 ! new shallow cloud for both heat and moisture upm
!-----------------------------------------------------------------------
!*** otherwise old cloud
!-----------------------------------------------------------------------
      if(dst.gt.0.) then 
        lbot=0          
!!!!!    ltop=lbot    
        ltop=kte     
        ptop=pbot   
        go to 800 
      endif
!-----------------------------------------------------------------------
!***  Otherwise continue with old cloud
!----------------ensure positive entropy change-------------------------
      DSTQ=DST*EPSDN
!----------------CHECK FOR ISOTHERMAL ATMOSPHERE------------------------
!
      DEN=POTSUM-PSUM
!
      IF(-DEN/PSUM<5.E-5)THEN
        LBOT=0
!!!!    LTOP=LBOT
        LTOP=KTE
        PTOP=PBOT
        GO TO 800
!
!----------------SLOPE OF THE REFERENCE HUMIDITY PROFILE----------------
!
      ELSE
        DQREF=(QOTSUM-DSTQ-QSUM)/DEN
      ENDIF
!
!-------------- HUMIDITY DOES NOT INCREASE WITH HEIGHT------------------
!
      IF(DQREF<0.)THEN
        LBOT=0
!!!!    LTOP=LBOT
        LTOP=KTE
        PTOP=PBOT
        GO TO 800
      ENDIF
!
!----------------HUMIDITY AT THE CLOUD TOP------------------------------
!
      QRFTP=QSUM-DQREF*PSUM
!
!----------------HUMIDITY PROFILE---------------------------------------
!
      DO L=LTOP,LBOT
        QRFKL=(FPK(L)-FPTK)*DQREF+QRFTP
!
!***  TOO DRY CLOUDS NOT ALLOWED
!
        TNEW=(TREFK(L)-TK(L))*TAUKSC+TK(L)
        QSATK(L)=PQ0/PK(L)*EXP(A2*(TNEW-A3)/(TNEW-A4))
        QNEW=(QRFKL-QK(L))*TAUKSC+QK(L)
!
        IF(QNEW<QSATK(L)*RHLSC)THEN
          LBOT=0
!!!!      LTOP=LBOT
          LTOP=KTE
          PTOP=PBOT
          GO TO 800
        ENDIF
!
!-------------TOO MOIST CLOUDS NOT ALLOWED------------------------------
!
        IF(QNEW>QSATK(L)*RHHSC)THEN
          LBOT=0
!!!!      LTOP=LBOT
          LTOP=KTE
          PTOP=PBOT
          GO TO 800
        ENDIF

!
        THVREF(L)=TREFK(L)*APEK(L)*(QRFKL*D608+1.)
        QREFK(L)=QRFKL
      ENDDO
!
!------------------ ELIMINATE CLOUDS WITH BOTTOMS TOO DRY --------------
!!
!      qnew=(qrefk(lbot)-qk(lbot))*tauksc+qk(lbot)
!!
!      if(qnew<qk(lbot+1)*stresh)then  !!?? stresh too large!!
!        lbot=0
!!!!!!   ltop=lbot
!        ltop=kte
!        ptop=pbot
!        go to 800
!      endif
!!
!-------------- ELIMINATE IMPOSSIBLE SLOPES (BETTS,DTHETA/DQ)------------
!
      DO L=LTOP,LBOT
        DTDP=(THVREF(L-1)-THVREF(L))/(PRSMID(L)-PRSMID(L-1))
!
        IF(DTDP<EPSDT)THEN
          LBOT=0
!!!!!     LTOP=LBOT
          LTOP=KTE
          PTOP=PBOT
          GO TO 800
        ENDIF
!
      ENDDO
!-----------------------------------------------------------------------
!***  Relaxation to reference profiles
!-----------------------------------------------------------------------
              go to 820 ! relaxation
!-----------------------------------------------------------------------
!***  New cloud starts here
!-----------------------------------------------------------------------
 810  do l=kts,lmh
        dpk(l)=dprs(l)
        rhk(l)=qk(l)/qsatk(l)
        thvmk(l)=tk(l)*apek(l) !zj *(qk(l)*0.608+1.)
!----calculate updraft temperature along moist adiabat tref(l)----------
        if(prsmid(l).lt.plq) then
          call ttblex(itb,jtb,pl,pk(l),rdp,rdthe &
                     ,sthe,the0,thescnv(l),ttbl,tref(l))
        else
          call ttblex(itbq,jtbq,plq,pk(l),rdpq,rdtheq &
                     ,stheq,the0q,thescnv(l),ttblq,tref(l))
        endif
!-----------------------------------------------------------------------
        thmak(l)=tref(l)*apek(l)
      enddo
!-------------mean rh and slopes within cloud---------------------------
      sumdp=0.
      sumrh=0.
      a11=0.
      a12=0.
      b1qsat=0.
      b2qsat=0.
      b1thvm=0.
      b1thma=0.
      b2thvm=0.
      b2thma=0.
      b1rh  =0.
      b2rh  =0.
!
      do l=ltop,lbot
        sumdp=dprs(l)+sumdp
        sumrh=rhk(l)*dprs(l)+sumrh
        a11=prsmid(l)**2*dprs(l)+a11
        a12=prsmid(l)*dprs(l)+a12
        b1qsat=qsatk(l)*prsmid(l)*dprs(l)+b1qsat
        b1thvm=thvmk(l)*prsmid(l)*dprs(l)+b1thvm
        b1thma=thmak(l)*prsmid(l)*dprs(l)+b1thma
        b1rh  =rhk  (l)*prsmid(l)*dprs(l)+b1rh
        b2qsat=qsatk(l)*dprs(l)+b2qsat
        b2thvm=thvmk(l)*dprs(l)+b2thvm
        b2thma=thmak(l)*dprs(l)+b2thma
        b2rh  =rhk  (l)*dprs(l)+b2rh
      enddo
!
      rhmean=sumrh/sumdp
!-------------no shallow convection if the cloud is saturated-----------
      if(rhmean.gt.0.95) then
        lbot=0
        ltop=kte
        ptop=pbot
        go to 800
      endif
!-----------------------------------------------------------------------
      a21=a12
      a22=sumdp
!
      rden=1./(a11*a22-a12*a21)
!
      aqs=(b1qsat*a22-a12*b2qsat)*rden
      avm=(b1thvm*a22-a12*b2thvm)*rden
      ama=(b1thma*a22-a12*b2thma)*rden
      arh=(b1rh  *a22-a12*b2rh  )*rden
!
      bqs=(a11*b2qsat-b1qsat*a21)*rden
      bvm=(a11*b2thvm-b1thvm*a21)*rden
      bma=(a11*b2thma-b1thma*a21)*rden
      brh=(a11*b2rh  -b1rh  *a21)*rden
!-------------no shallow convection if the cloud moister on top---------
      if(arh.lt.0.) then !soft2
        lbot=0           !soft2
        ltop=kte         !soft2
        ptop=pbot        !soft2
        go to 800        !soft2
      endif              !soft2
!-------------first guess t & q profiles--------------------------------
      adef=(1.-deftop)*2./(pk(lbot)-pk(ltop)) !soft2
!
      do l=ltop,lbot
        fk=(pk(l)-pk(ltop))*adef+deftop !soft2
        rhref=rhmean*fk                 !soft2
!
        wcld=(1.-wdry)*rhref/(1.-wdry*rhref)
        trefk(l)=((1.-wcld)*(avm*pk(l)+bvm) &
                 +    wcld *(ama*pk(l)+bma))/apek(l)
        qrefk(l)=rhref*(aqs*prsmid(l)+bqs)
      enddo
!-------------enthalpy conservation-------------------------------------
      sumdp=0.
      sumdt=0.
      sumdq=0.
!
      do l=ltop,lbot
        sumdp=dpk(l)+sumdp
        sumdt=(tk(l)-trefk(l))*dpk(l)+sumdt
        sumdq=(qk(l)-qrefk(l))*dpk(l)+sumdq
      enddo
!
      rdpsum=1./sumdp
!
      tcorr=sumdt*rdpsum
      qcorr=sumdq*rdpsum
!
      do l=ltop,lbot
        trefk(l)=trefk(l)+tcorr
        qrefk(l)=qrefk(l)+qcorr
      enddo
!-----------------------------------------------------------------------
        dsq=0.
        dst=0.
!
        do l=ltop,lbot
          tkl=tk(l)
          diftl=(trefk(l)-tkl  )*tauksc
          difql=(qrefk(l)-qk(l))*tauksc
          dpot=dprs(l)/(tkl+tkl+diftl)
          dst=diftl      *dpot+dst
          dsq=difql*el(l)*dpot+dsq
        enddo
!
        dst=(dst+dst)*cp
        dsq=dsq+dsq
        dentpy=dst+dsq
!
        if(dentpy.lt.0.) then
          lbot=0
!!!!          ltop=lbot
          ltop=kte
          ptop=pbot
          go to 800
        endif
!--------------relaxation towards reference profiles--------------------
 820  do l=ltop,lbot
        DTDT(L)=(TREFK(L)-TK(L))*TAUKSC*RDTCNVC
        DQDT(L)=(QREFK(L)-QK(L))*TAUKSC*RDTCNVC
      ENDDO
!-----------------------------------------------------------------------
!***  Begin debugging convection
      IF(PRINT_DIAG)THEN
        DO L=LBOT,LTOP,-1
          WRITE(6,"(a,i3,4e12.4)") '{cu2 KFLIP,DT,DTDT,DQ,DQDT = ' &
               ,KTE+1-L,TREFK(L)-TK(L),DTDT(L),QREFK(L)-QK(L),DQDT(L)
        ENDDO
      ENDIF
!***  End debugging convection
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS
!SCSCSCSCSCSCSC         END OF SHALLOW CONVECTION        SCSCSCSCSCSCSCS
!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS
!-----------------------------------------------------------------------
  800 CONTINUE
!-----------------------------------------------------------------------
      END SUBROUTINE BMJ_DEV
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
                           SUBROUTINE TTBLEX                            &
     & (ITBX,JTBX,PLX,PRSMID,RDPX,RDTHEX,STHE                           &
     & ,THE0,THESP,TTBL,TREF)
!-----------------------------------------------------------------------
!     ******************************************************************
!     *                                                                *
!     *           EXTRACT TEMPERATURE OF THE MOIST ADIABAT FROM        *
!     *                      THE APPROPRIATE TTBL                      *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: ITBX,JTBX
!
      REAL,INTENT(IN) :: PLX,PRSMID,RDPX,RDTHEX,THESP
!
      REAL,DIMENSION(ITBX),INTENT(IN) :: STHE,THE0
!
      REAL,DIMENSION(JTBX,ITBX),INTENT(IN) :: TTBL
!
      REAL,INTENT(OUT) :: TREF
!-----------------------------------------------------------------------
      REAL :: BTHE00K,BTHE10K,BTHK,PK,PP,QQ,STHE00K,STHE10K,STHK        &
     &       ,T00K,T01K,T10K,T11K,TPK,TTHK
!
      INTEGER :: IPTB,ITHTB
!-----------------------------------------------------------------------
!----------------SCALING PRESSURE & TT TABLE INDEX----------------------
!-----------------------------------------------------------------------
      PK=PRSMID
      TPK=(PK-PLX)*RDPX
      QQ=TPK-AINT(TPK)
      IPTB=INT(TPK)+1
!----------------KEEPING INDICES WITHIN THE TABLE-----------------------
      IF(IPTB<1)THEN
        IPTB=1
        QQ=0.
      ENDIF
!
      IF(IPTB>=ITBX)THEN
        IPTB=ITBX-1
        QQ=0.
      ENDIF
!----------------BASE AND SCALING FACTOR FOR THETAE---------------------
      BTHE00K=THE0(IPTB)
      STHE00K=STHE(IPTB)
      BTHE10K=THE0(IPTB+1)
      STHE10K=STHE(IPTB+1)
!----------------SCALING THE & TT TABLE INDEX---------------------------
      BTHK=(BTHE10K-BTHE00K)*QQ+BTHE00K
      STHK=(STHE10K-STHE00K)*QQ+STHE00K
      TTHK=(THESP-BTHK)/STHK*RDTHEX
      PP=TTHK-AINT(TTHK)
      ITHTB=INT(TTHK)+1
!----------------KEEPING INDICES WITHIN THE TABLE-----------------------
      IF(ITHTB<1)THEN
        ITHTB=1
        PP=0.
      ENDIF
!
      IF(ITHTB>=JTBX)THEN
        ITHTB=JTBX-1
        PP=0.
      ENDIF
!----------------TEMPERATURE AT FOUR SURROUNDING TT TABLE PTS.----------
      T00K=TTBL(ITHTB,IPTB)
      T10K=TTBL(ITHTB+1,IPTB)
      T01K=TTBL(ITHTB,IPTB+1)
      T11K=TTBL(ITHTB+1,IPTB+1)
!-----------------------------------------------------------------------
!----------------PARCEL TEMPERATURE-------------------------------------
!-----------------------------------------------------------------------
      TREF=(T00K+(T10K-T00K)*PP+(T01K-T00K)*QQ                          &
     &    +(T00K-T10K-T01K+T11K)*PP*QQ)
!-----------------------------------------------------------------------
      END SUBROUTINE TTBLEX
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE BMJ_INIT_DEV(CLDEFI,RESTART                            &
     &                   ,IDS,IDE,JDS,JDE,KDS,KDE                       &
     &                   ,IMS,IME,JMS,JME,KMS,KME                       &
     &                   ,ITS,ITE,JTS,JTE,KTS,KTE)
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      LOGICAL,INTENT(IN) :: RESTART
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
!      REAL,INTENT(OUT) :: ACUTIM,AVCNVC
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CLDEFI
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: ELIWV=2.683E6,EPS=1.E-9
!
      REAL, DIMENSION(JTB) :: APP,APT,AQP,AQT,PNEW,POLD,QSNEW,QSOLD     &
     &                       ,THENEW,THEOLD,TNEW,TOLD,Y2P,Y2T
!
      REAL,DIMENSION(JTBQ) :: APTQ,AQTQ,THENEWQ,THEOLDQ                 &
     &                       ,TNEWQ,TOLDQ,Y2TQ
!
      INTEGER :: I,J,K,ITF,JTF,KTF
      INTEGER :: KTH,KTHM,KTHM1,KP,KPM,KPM1
!
      REAL :: APE,DP,DQS,DTH,DTHE,P,QS,QS0K,RD,SQSK,STHEK               &
     &       ,TH,THE0K,DENOM
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(KTE,KDE-1)
      ITF=MIN0(ITE,IDE-1)
! 
      IF(.NOT.RESTART)THEN
        DO J=JTS,JTF
        DO I=ITS,ITF
!!!       CLDEFI(I,J)=1.
          CLDEFI(I,J)=AVGEFI
        ENDDO
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!    FOR ESMF VERSION ONLY   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       ACUTIM=0
!       AVCNVC=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF
!
!-----------------------------------------------------------------------
!----------------COARSE LOOK-UP TABLE FOR SATURATION POINT--------------
!-----------------------------------------------------------------------
!
      KTHM=JTB
      KPM=ITB
      KTHM1=KTHM-1
      KPM1=KPM-1
!
      DTH=(THH-THL)/REAL(KTHM-1)
      DP =(PH -PL )/REAL(KPM -1)
!
      TH=THL-DTH
!-----------------------------------------------------------------------
!
      RD=R_D
!
      DO 100 KTH=1,KTHM
!
      TH=TH+DTH
      P=PL-DP
!
      DO KP=1,KPM
        P=P+DP
        APE=(100000./P)**(RD/CP)
        DENOM=TH-A4*APE
        IF (DENOM>EPS) THEN
           QSOLD(KP)=PQ0/P*EXP(A2*(TH-A3*APE)/DENOM)
        ELSE
           QSOLD(KP)=0.
        ENDIF
        POLD(KP)=P
      ENDDO
!
      QS0K=QSOLD(1)
      SQSK=QSOLD(KPM)-QSOLD(1)
      QSOLD(1  )=0.
      QSOLD(KPM)=1.
!
      DO KP=2,KPM1
        QSOLD(KP)=(QSOLD(KP)-QS0K)/SQSK
        IF((QSOLD(KP)-QSOLD(KP-1))<EPS)QSOLD(KP)=QSOLD(KP-1)+EPS
      ENDDO
!
      QS0(KTH)=QS0K
      QS0_EXP(KTH)=QS0K
      SQS(KTH)=SQSK
      SQS_EXP(KTH)=SQSK
!-----------------------------------------------------------------------
      QSNEW(1  )=0.
      QSNEW(KPM)=1.
      DQS=1./REAL(KPM-1)
!
      DO KP=2,KPM1
        QSNEW(KP)=QSNEW(KP-1)+DQS
      ENDDO
!
      Y2P(1   )=0.
      Y2P(KPM )=0.
!
      CALL SPLINE(JTB,KPM,QSOLD,POLD,Y2P,KPM,QSNEW,PNEW,APP,AQP)
!
      DO KP=1,KPM
        PTBL(KP,KTH)=PNEW(KP)
        PTBL_EXP(KP,KTH)=PNEW(KP)
      ENDDO
!-----------------------------------------------------------------------
  100 CONTINUE
!-----------------------------------------------------------------------
!------------COARSE LOOK-UP TABLE FOR T(P) FROM CONSTANT THE------------
!-----------------------------------------------------------------------
      P=PL-DP
!
      DO 200 KP=1,KPM
!
      P=P+DP
      TH=THL-DTH
!
      DO KTH=1,KTHM
        TH=TH+DTH
        APE=(1.E5/P)**(RD/CP)
        DENOM=TH-A4*APE
        IF (DENOM>EPS) THEN
           QS=PQ0/P*EXP(A2*(TH-A3*APE)/DENOM)
        ELSE
           QS=0.
        ENDIF
!        QS=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
        TOLD(KTH)=TH/APE
        THEOLD(KTH)=TH*EXP(ELIWV*QS/(CP*TOLD(KTH)))
      ENDDO
!
      THE0K=THEOLD(1)
      STHEK=THEOLD(KTHM)-THEOLD(1)
      THEOLD(1   )=0.
      THEOLD(KTHM)=1.
!
      DO KTH=2,KTHM1
        THEOLD(KTH)=(THEOLD(KTH)-THE0K)/STHEK
        IF((THEOLD(KTH)-THEOLD(KTH-1)).LT.EPS)                          &
     &      THEOLD(KTH)=THEOLD(KTH-1)  +  EPS
      ENDDO
!
      THE0(KP)=THE0K
      STHE(KP)=STHEK
!-----------------------------------------------------------------------
      THENEW(1  )=0.
      THENEW(KTHM)=1.
      DTHE=1./REAL(KTHM-1)
!
      DO KTH=2,KTHM1
        THENEW(KTH)=THENEW(KTH-1)+DTHE
      ENDDO
!
      Y2T(1   )=0.
      Y2T(KTHM)=0.
!
      CALL SPLINE(JTB,KTHM,THEOLD,TOLD,Y2T,KTHM,THENEW,TNEW,APT,AQT)
!
      DO KTH=1,KTHM
        TTBL(KTH,KP)=TNEW(KTH)
      ENDDO
!-----------------------------------------------------------------------
  200 CONTINUE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---------------FINE LOOK-UP TABLE FOR SATURATION POINT-----------------
!-----------------------------------------------------------------------
      KTHM=JTBQ
      KPM=ITBQ
      KTHM1=KTHM-1
      KPM1=KPM-1
!
      DTH=(THHQ-THL)/REAL(KTHM-1)
      DP=(PH-PLQ)/REAL(KPM-1)
!
      TH=THL-DTH
      P=PLQ-DP
!-----------------------------------------------------------------------
!---------------FINE LOOK-UP TABLE FOR T(P) FROM CONSTANT THE-----------
!-----------------------------------------------------------------------
      DO 300 KP=1,KPM
!
      P=P+DP
      TH=THL-DTH
!
      DO KTH=1,KTHM
        TH=TH+DTH
        APE=(1.E5/P)**(RD/CP)
        DENOM=TH-A4*APE
        IF (DENOM>EPS) THEN
           QS=PQ0/P*EXP(A2*(TH-A3*APE)/DENOM)
        ELSE
           QS=0.
        ENDIF
!        QS=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
        TOLDQ(KTH)=TH/APE
        THEOLDQ(KTH)=TH*EXP(ELIWV*QS/(CP*TOLDQ(KTH)))
      ENDDO
!
      THE0K=THEOLDQ(1)
      STHEK=THEOLDQ(KTHM)-THEOLDQ(1)
      THEOLDQ(1   )=0.
      THEOLDQ(KTHM)=1.
!
      DO KTH=2,KTHM1
        THEOLDQ(KTH)=(THEOLDQ(KTH)-THE0K)/STHEK
        IF((THEOLDQ(KTH)-THEOLDQ(KTH-1))<EPS)                           &
     &      THEOLDQ(KTH)=THEOLDQ(KTH-1)+EPS
      ENDDO
!
      THE0Q(KP)=THE0K
      STHEQ(KP)=STHEK
!-----------------------------------------------------------------------
      THENEWQ(1  )=0.
      THENEWQ(KTHM)=1.
      DTHE=1./REAL(KTHM-1)
!
      DO KTH=2,KTHM1
        THENEWQ(KTH)=THENEWQ(KTH-1)+DTHE
      ENDDO
!
      Y2TQ(1   )=0.
      Y2TQ(KTHM)=0.
!
      CALL SPLINE(JTBQ,KTHM,THEOLDQ,TOLDQ,Y2TQ,KTHM                     &
     &           ,THENEWQ,TNEWQ,APTQ,AQTQ)
!
      DO KTH=1,KTHM
        TTBLQ(KTH,KP)=TNEWQ(KTH)
      ENDDO
!-----------------------------------------------------------------------
  300 CONTINUE
!-----------------------------------------------------------------------
      END SUBROUTINE BMJ_INIT_DEV
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
      SUBROUTINE SPLINE(JTBX,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
!   ********************************************************************
!   *                                                                  *
!   *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE          *
!   *  PROGRAMED FOR A SMALL SCALAR MACHINE.                           *
!   *                                                                  *
!   *  PROGRAMER Z. JANJIC                                             *
!   *                                                                  *
!   *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3.   *
!   *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE       *
!   *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.         *
!   *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.     *
!   *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL   *
!   *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE        *
!   *         SPECIFIED.                                               *
!   *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.       *
!   *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE       *
!   *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)     *
!   *         AND LE XOLD(NOLD).                                       *
!   *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.             *
!   *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                  *
!   *                                                                  *
!   ********************************************************************
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: JTBX,NNEW,NOLD
      REAL,DIMENSION(JTBX),INTENT(IN) :: XNEW,XOLD,YOLD
      REAL,DIMENSION(JTBX),INTENT(INOUT) :: P,Q,Y2
      REAL,DIMENSION(JTBX),INTENT(OUT) :: YNEW
!
      INTEGER :: K,K1,K2,KOLD,NOLDM1
      REAL :: AK,BK,CK,DEN,DX,DXC,DXL,DXR,DYDXL,DYDXR                   &
     &       ,RDX,RTDXC,X,XK,XSQ,Y2K,Y2KP1
!-----------------------------------------------------------------------
      NOLDM1=NOLD-1
!
      DXL=XOLD(2)-XOLD(1)
      DXR=XOLD(3)-XOLD(2)
      DYDXL=(YOLD(2)-YOLD(1))/DXL
      DYDXR=(YOLD(3)-YOLD(2))/DXR
      RTDXC=0.5/(DXL+DXR)
!
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
      Q(1)=-RTDXC*DXR
!
      IF(NOLD==3)GO TO 150
!-----------------------------------------------------------------------
      K=3
!
  100 DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
!
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
!
      K=K+1
      IF(K<NOLD)GO TO 100
!-----------------------------------------------------------------------
  150 K=NOLDM1
!
  200 Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
!
      K=K-1
      IF(K>1)GO TO 200
!-----------------------------------------------------------------------
      K1=1
!
  300 XK=XNEW(K1)
!
      DO 400 K2=2,NOLD
!
      IF(XOLD(K2)>XK)THEN
        KOLD=K2-1
        GO TO 450
      ENDIF
!
  400 CONTINUE
!
      YNEW(K1)=YOLD(NOLD)
      GO TO 600
!
  450 IF(K1==1)GO TO 500
      IF(K==KOLD)GO TO 550
!
  500 K=KOLD
!
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
!
      AK=.1666667*RDX*(Y2KP1-Y2K)
      BK=0.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)
!
  550 X=XK-XOLD(K)
      XSQ=X*X
!
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
!
  600 K1=K1+1
      IF(K1<=NNEW)GO TO 300
!-----------------------------------------------------------------------
      END SUBROUTINE SPLINE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CU_BMJ_DEV
!
!-----------------------------------------------------------------------
