      SUBROUTINE SHALCVT3(IM,IX,KM,DT,DEL,PRSI,PRSL,PRSLK,KUO,Q,T)
!
      USE MACHINE , ONLY : kind_phys
      USE PHYSCONS, grav => con_g, CP => con_CP, HVAP => con_HVAP
     &,             RD => con_RD
      implicit none
!
!     include 'constant.h'
!
      integer              IM, IX, KM, KUO(IM)
      real(kind=kind_phys) DEL(IX,KM),   PRSI(IX,KM+1), PRSL(IX,KM),
     &                     PRSLK(IX,KM),
     &                     Q(IX,KM),     T(IX,KM),      DT
!
!     Locals
!
      real(kind=kind_phys) ck,    cpdt,   dmse,   dsdz1, dsdz2,
     &                     dsig,  dtodsl, dtodsu, eldq,  g,
     &                     gocp,  rtdls
!
      integer              k,k1,k2,kliftl,kliftu,kt,N2,I,iku,ik1,ik,ii
      integer              INDEX2(IM), KLCL(IM), KBOT(IM), KTOP(IM),kk
cc
C  PHYSICAL PARAMETERS
      PARAMETER(G=GRAV, GOCP=G/CP)
C  BOUNDS OF PARCEL ORIGIN
      PARAMETER(KLIFTL=2,KLIFTU=2)
      LOGICAL   LSHC(IM)
      real(kind=kind_phys) Q2(IM*KM),     T2(IM*KM),
     &                     PRSL2(IM*KM),  PRSLK2(IM*KM),
     &                     AL(IM*(KM-1)), AD(IM*KM), AU(IM*(KM-1))
C-----------------------------------------------------------------------
C  COMPRESS FIELDS TO POINTS WITH NO DEEP CONVECTION
C  AND MOIST STATIC INSTABILITY.
      DO I=1,IM
        LSHC(I)=.FALSE.
      ENDDO
      DO K=1,KM-1
        DO I=1,IM
          IF(KUO(I).EQ.0) THEN
            ELDQ    = HVAP*(Q(I,K)-Q(I,K+1))
            CPDT    = CP*(T(I,K)-T(I,K+1))
            RTDLS   = (PRSL(I,K)-PRSL(I,K+1)) /
     &                 PRSI(I,K+1)*RD*0.5*(T(I,K)+T(I,K+1))
            DMSE    = ELDQ+CPDT-RTDLS
            LSHC(I) = LSHC(I).OR.DMSE.GT.0.
          ENDIF
        ENDDO
      ENDDO
      N2 = 0
      DO I=1,IM
        IF(LSHC(I)) THEN
          N2         = N2 + 1
          INDEX2(N2) = I
        ENDIF
      ENDDO
      IF(N2.EQ.0) RETURN
      DO K=1,KM
        KK = (K-1)*N2
        DO I=1,N2
          IK         = KK + I
          ii         = index2(i)
          Q2(IK)     = Q(II,K)
          T2(IK)     = T(II,K)
          PRSL2(IK)  = PRSL(II,K)
          PRSLK2(IK) = PRSLK(II,K)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
C  COMPUTE MOIST ADIABAT AND DETERMINE LIMITS OF SHALLOW CONVECTION.
C  CHECK FOR MOIST STATIC INSTABILITY AGAIN WITHIN CLOUD.
      CALL MSTADBT3(N2,KM-1,KLIFTL,KLIFTU,PRSL2,PRSLK2,T2,Q2,
     &            KLCL,KBOT,KTOP,AL,AU)
      DO I=1,N2
        KBOT(I) = KLCL(I)-1
        KTOP(I) = KTOP(I)+1
        LSHC(I) = .FALSE.
      ENDDO
      DO K=1,KM-1
        KK = (K-1)*N2
        DO I=1,N2
          IF(K.GE.KBOT(I).AND.K.LT.KTOP(I)) THEN
            IK      = KK + I
            IKU     = IK + N2
            ELDQ    = HVAP * (Q2(IK)-Q2(IKU))
            CPDT    = CP   * (T2(IK)-T2(IKU))
            RTDLS   = (PRSL2(IK)-PRSL2(IKU)) /
     &                 PRSI(index2(i),K+1)*RD*0.5*(T2(IK)+T2(IKU))
            DMSE    = ELDQ + CPDT - RTDLS
            LSHC(I) = LSHC(I).OR.DMSE.GT.0.
            AU(IK)  = G/RTDLS
          ENDIF
        ENDDO
      ENDDO
      K1=KM+1
      K2=0
      DO I=1,N2
        IF(.NOT.LSHC(I)) THEN
          KBOT(I) = KM+1
          KTOP(I) = 0
        ENDIF
        K1 = MIN(K1,KBOT(I))
        K2 = MAX(K2,KTOP(I))
      ENDDO
      KT = K2-K1+1
      IF(KT.LT.2) RETURN
C-----------------------------------------------------------------------
C  SET EDDY VISCOSITY COEFFICIENT CKU AT SIGMA INTERFACES.
C  COMPUTE DIAGONALS AND RHS FOR TRIDIAGONAL MATRIX SOLVER.
C  EXPAND FINAL FIELDS.
      KK = (K1-1) * N2
      DO I=1,N2
        IK     = KK + I
        AD(IK) = 1.
      ENDDO
!
!     DTODSU=DT/DEL(K1)
      DO K=K1,K2-1
!       DTODSL=DTODSU
!       DTODSU=   DT/DEL(K+1)
!       DSIG=SL(K)-SL(K+1)
        KK = (K-1) * N2
        DO I=1,N2
          ii     = index2(i)
          DTODSL = DT/DEL(II,K)
          DTODSU = DT/DEL(II,K+1)
          DSIG   = PRSL(II,K) - PRSL(II,K+1)
          IK     = KK + I
          IKU    = IK + N2
          IF(K.EQ.KBOT(I)) THEN
            CK=1.5
          ELSEIF(K.EQ.KTOP(I)-1) THEN
            CK=1.
          ELSEIF(K.EQ.KTOP(I)-2) THEN
            CK=3.
          ELSEIF(K.GT.KBOT(I).AND.K.LT.KTOP(I)-2) THEN
            CK=5.
          ELSE
            CK=0.
          ENDIF
          DSDZ1   = CK*DSIG*AU(IK)*GOCP
          DSDZ2   = CK*DSIG*AU(IK)*AU(IK)
          AU(IK)  = -DTODSL*DSDZ2
          AL(IK)  = -DTODSU*DSDZ2
          AD(IK)  = AD(IK)-AU(IK)
          AD(IKU) = 1.-AL(IK)
          T2(IK)  = T2(IK)+DTODSL*DSDZ1
          T2(IKU) = T2(IKU)-DTODSU*DSDZ1
        ENDDO
      ENDDO
      IK1=(K1-1)*N2+1
      CALL TRIDI2T3(N2,KT,AL(IK1),AD(IK1),AU(IK1),Q2(IK1),T2(IK1),
     &                                  AU(IK1),Q2(IK1),T2(IK1))
      DO K=K1,K2
        KK = (K-1)*N2
        DO I=1,N2
          IK = KK + I
          Q(INDEX2(I),K) = Q2(IK)
          T(INDEX2(I),K) = T2(IK)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TRIDI2T3(L,N,CL,CM,CU,R1,R2,AU,A1,A2)
cyt      INCLUDE DBTRIDI2;
cc
      USE MACHINE     , ONLY : kind_phys
      implicit none
      integer             k,n,l,i
      real(kind=kind_phys) fk
cc
      real(kind=kind_phys)
     &          CL(L,2:N),CM(L,N),CU(L,N-1),R1(L,N),R2(L,N),
     &          AU(L,N-1),A1(L,N),A2(L,N)
C-----------------------------------------------------------------------
      DO I=1,L
        FK=1./CM(I,1)
        AU(I,1)=FK*CU(I,1)
        A1(I,1)=FK*R1(I,1)
        A2(I,1)=FK*R2(I,1)
      ENDDO
      DO K=2,N-1
        DO I=1,L
          FK=1./(CM(I,K)-CL(I,K)*AU(I,K-1))
          AU(I,K)=FK*CU(I,K)
          A1(I,K)=FK*(R1(I,K)-CL(I,K)*A1(I,K-1))
          A2(I,K)=FK*(R2(I,K)-CL(I,K)*A2(I,K-1))
        ENDDO
      ENDDO
      DO I=1,L
        FK=1./(CM(I,N)-CL(I,N)*AU(I,N-1))
        A1(I,N)=FK*(R1(I,N)-CL(I,N)*A1(I,N-1))
        A2(I,N)=FK*(R2(I,N)-CL(I,N)*A2(I,N-1))
      ENDDO
      DO K=N-1,1,-1
        DO I=1,L
          A1(I,K)=A1(I,K)-AU(I,K)*A1(I,K+1)
          A2(I,K)=A2(I,K)-AU(I,K)*A2(I,K+1)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
