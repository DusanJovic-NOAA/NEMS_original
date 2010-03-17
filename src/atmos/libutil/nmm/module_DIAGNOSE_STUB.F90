!-----------------------------------------------------------------------
!
      MODULE MODULE_DIAGNOSE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE TWR(ARRAY,KK,FIELD,NTSD,MYPE,NPES,MPI_COMM_COMP       &
                    ,IDS,IDE,JDS,JDE                                   &
                    ,IMS,IME,JMS,JME                                   &
                    ,ITS,ITE,JTS,JTE)
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INCLUDE "../../../../inc/kind.inc"
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE                 &
                                      ,IMS,IME,JMS,JME                 &
                                      ,ITS,ITE,JTS,JTE                 &
                                      ,KK,MPI_COMM_COMP,MYPE,NPES,NTSD
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,KK),INTENT(IN) :: ARRAY
!
      CHARACTER(*),INTENT(IN) :: FIELD
!
!*** LOCAL VARIABLES
!
      INTEGER :: IUNIT=23
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
      INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: STATUS_ARRAY
      INTEGER(KIND=KINT),DIMENSION(2) :: IM_REM,JM_REM,IT_REM,JT_REM
!
      INTEGER(KIND=KINT) :: I,IENDX,IER,IPE,IRECV,IRTN,ISEND           &
                           ,J,K,N,NLEN,NSIZE
      INTEGER(KIND=KINT) :: ITS_REM,ITE_REM,JTS_REM,JTE_REM
!
      REAL(KIND=KFPT),DIMENSION(IDS:IDE,JDS:JDE) :: TWRITE
      REAL(KIND=KFPT),ALLOCATABLE,DIMENSION(:) :: VALUES
      CHARACTER(5) :: TIMESTEP
      CHARACTER(6) :: FMT
      CHARACTER(12) :: FILENAME
      integer :: isee=055,jsee=138,lsee=10,pesee=6
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      IF(NTSD<=9)THEN
        FMT='(I1.1)'
        NLEN=1
      ELSEIF(NTSD<=99)THEN
        FMT='(I2.2)'
        NLEN=2
      ELSEIF(NTSD<=999)THEN
        FMT='(I3.3)'
        NLEN=3
      ELSEIF(NTSD<=9999)THEN
        FMT='(I4.4)'
        NLEN=4
      ELSEIF(NTSD<=99999)THEN
        FMT='(I5.5)'
        NLEN=5
      ENDIF
      WRITE(TIMESTEP,FMT)NTSD
      FILENAME=FIELD//'_'//TIMESTEP(1:NLEN)
!
      IF(MYPE==0)THEN
        CLOSE(IUNIT)
        OPEN(UNIT=IUNIT,FILE=FILENAME,FORM='UNFORMATTED'               &
            ,STATUS='REPLACE',IOSTAT=IER)
      ENDIF
!
!----------------------------------------------------------------------
      DO 500 K=1,KK
!----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          TWRITE(I,J)=ARRAY(I,J,K)
        ENDDO
        ENDDO
!
        DO IPE=1,NPES-1
          CALL MPI_RECV(IT_REM,2,MPI_INTEGER,IPE,IPE                    &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
          CALL MPI_RECV(JT_REM,2,MPI_INTEGER,IPE,IPE                    &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
!
          ITS_REM=IT_REM(1)
          ITE_REM=IT_REM(2)
          JTS_REM=JT_REM(1)
          JTE_REM=JT_REM(2)
!
          NSIZE=(ITE_REM-ITS_REM+1)*(JTE_REM-JTS_REM+1)
          ALLOCATE(VALUES(1:NSIZE))
!
          CALL MPI_RECV(VALUES,NSIZE,MPI_REAL,IPE,IPE                   &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
          N=0
          DO J=JTS_REM,JTE_REM
            DO I=ITS_REM,ITE_REM
              N=N+1
              TWRITE(I,J)=VALUES(N)
            ENDDO
          ENDDO
!
          DEALLOCATE(VALUES)
!
        ENDDO
!
!----------------------------------------------------------------------
      ELSE
!
        NSIZE=(ITE-ITS+1)*(JTE-JTS+1)
        ALLOCATE(VALUES(1:NSIZE))
!
        N=0
        DO J=JTS,JTE
        DO I=ITS,ITE
          N=N+1
          VALUES(N)=ARRAY(I,J,K)
        ENDDO
        ENDDO
!
        IT_REM(1)=ITS
        IT_REM(2)=ITE
        JT_REM(1)=JTS
        JT_REM(2)=JTE
!
        CALL MPI_SEND(IT_REM,2,MPI_INTEGER,0,MYPE                       &
                     ,MPI_COMM_COMP,ISEND)
        CALL MPI_SEND(JT_REM,2,MPI_INTEGER,0,MYPE                       &
                     ,MPI_COMM_COMP,ISEND)
!
        CALL MPI_SEND(VALUES,NSIZE,MPI_REAL,0,MYPE                      &
                     ,MPI_COMM_COMP,ISEND)
!
        DEALLOCATE(VALUES)
!
      ENDIF
!----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
      IF(MYPE==0)THEN
!
        DO J=JDS,JDE
          IENDX=IDE
          WRITE(IUNIT)(TWRITE(I,J),I=IDS,IENDX)
        ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------
  500 CONTINUE
!
      IF(MYPE==0)CLOSE(IUNIT)
!----------------------------------------------------------------------
!
      END SUBROUTINE TWR
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
      SUBROUTINE VWR(ARRAY,KK,FIELD,NTSD,MYPE,NPES,MPI_COMM_COMP       &
                    ,IDS,IDE,JDS,JDE                                   &
                    ,IMS,IME,JMS,JME                                   &
                    ,ITS,ITE,JTS,JTE)
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INCLUDE "../../../../inc/kind.inc"
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE                 &
                                      ,IMS,IME,JMS,JME                 &
                                      ,ITS,ITE,JTS,JTE                 &
                                      ,KK,MPI_COMM_COMP,MYPE,NPES,NTSD
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,KK),INTENT(IN) :: ARRAY
!
      CHARACTER(*),INTENT(IN) :: FIELD
!
!*** LOCAL VARIABLES
!
      INTEGER :: IUNIT=23
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
      INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: STATUS_ARRAY
      INTEGER(KIND=KINT),DIMENSION(2) :: IM_REM,JM_REM,IT_REM,JT_REM
!
      INTEGER(KIND=KINT) :: I,IENDX,IER,IPE,IRECV,IRTN,ISEND           &
                           ,J,K,N,NLEN,NSIZE
      INTEGER(KIND=KINT) :: ITS_REM,ITE_REM,JTS_REM,JTE_REM
!
      REAL(KIND=KFPT),DIMENSION(IDS:IDE,JDS:JDE) :: TWRITE
      REAL(KIND=KFPT),ALLOCATABLE,DIMENSION(:) :: VALUES
      CHARACTER(5) :: TIMESTEP
      CHARACTER(6) :: FMT
      CHARACTER(12) :: FILENAME
      integer :: isee=055,jsee=138,lsee=10,pesee=6
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      IF(NTSD<=9)THEN
        FMT='(I1.1)'
        NLEN=1
      ELSEIF(NTSD<=99)THEN
        FMT='(I2.2)'
        NLEN=2
      ELSEIF(NTSD<=999)THEN
        FMT='(I3.3)'
        NLEN=3
      ELSEIF(NTSD<=9999)THEN
        FMT='(I4.4)'
        NLEN=4
      ELSEIF(NTSD<=99999)THEN
        FMT='(I5.5)'
        NLEN=5
      ENDIF
      WRITE(TIMESTEP,FMT)NTSD
      FILENAME=FIELD//'_'//TIMESTEP(1:NLEN)
!
      IF(MYPE==0)THEN
        CLOSE(IUNIT)
        OPEN(UNIT=IUNIT,FILE=FILENAME,FORM='UNFORMATTED',IOSTAT=IER)
      ENDIF
!
!----------------------------------------------------------------------
      DO 500 K=1,KK
!----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          TWRITE(I,J)=ARRAY(I,J,K)
        ENDDO
        ENDDO
!
        DO IPE=1,NPES-1
          CALL MPI_RECV(IT_REM,2,MPI_INTEGER,IPE,IPE                    &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
          CALL MPI_RECV(JT_REM,2,MPI_INTEGER,IPE,IPE                    &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
!
          ITS_REM=IT_REM(1)
          ITE_REM=IT_REM(2)
          JTS_REM=JT_REM(1)
          JTE_REM=JT_REM(2)
!
          NSIZE=(ITE_REM-ITS_REM+1)*(JTE_REM-JTS_REM+1)
          ALLOCATE(VALUES(1:NSIZE))
!
          CALL MPI_RECV(VALUES,NSIZE,MPI_REAL,IPE,IPE                   &
                       ,MPI_COMM_COMP,JSTAT,IRECV)
          N=0
          DO J=JTS_REM,JTE_REM
            DO I=ITS_REM,ITE_REM
              N=N+1
              TWRITE(I,J)=VALUES(N)
            ENDDO
          ENDDO
!
          DEALLOCATE(VALUES)
!
        ENDDO
!
!----------------------------------------------------------------------
      ELSE
!
        NSIZE=(ITE-ITS+1)*(JTE-JTS+1)
        ALLOCATE(VALUES(1:NSIZE))
!
        N=0
        DO J=JTS,JTE
        DO I=ITS,ITE
          N=N+1
          VALUES(N)=ARRAY(I,J,K)
        ENDDO
        ENDDO
!
        IT_REM(1)=ITS
        IT_REM(2)=ITE
        JT_REM(1)=JTS
        JT_REM(2)=JTE
!
        CALL MPI_SEND(IT_REM,2,MPI_INTEGER,0,MYPE                       &
                     ,MPI_COMM_COMP,ISEND)
        CALL MPI_SEND(JT_REM,2,MPI_INTEGER,0,MYPE                       &
                     ,MPI_COMM_COMP,ISEND)
!
        CALL MPI_SEND(VALUES,NSIZE,MPI_REAL,0,MYPE                      &
                     ,MPI_COMM_COMP,ISEND)
!
        DEALLOCATE(VALUES)
!
      ENDIF
!----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
      IF(MYPE==0)THEN
!
        DO J=JDS,JDE-1
          IENDX=IDE-1
          WRITE(IUNIT)(TWRITE(I,J),I=IDS,IENDX)
        ENDDO
!
      ENDIF
!
!----------------------------------------------------------------------
  500 CONTINUE
!
      IF(MYPE==0)CLOSE(IUNIT)
!----------------------------------------------------------------------
!
      END SUBROUTINE VWR
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
      SUBROUTINE EXIT(NAME,PINT,T,Q,U,V,Q2,W                           &
                     ,NTSD,MYPE,MPI_COMM_COMP                          &
                     ,IDS,IDE,JDS,JDE,LM                               &
                     ,IMS,IME,JMS,JME                                  &
                     ,ITS,ITE,JTS,JTE)
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      INCLUDE "mpif.h"
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                         &
                           ,IMS,IME,JMS,JME                            &
                           ,ITS,ITE,JTS,JTE                            &
                           ,MYPE,MPI_COMM_COMP,NTSD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: T,Q,U,V,Q2,W
      REAL,DIMENSION(IMS:IME,JMS:JME,LM+1),INTENT(IN) :: PINT
      CHARACTER(*),INTENT(IN) :: NAME
!
      INTEGER :: I,J,K,IEND,IERR,IRET
      CHARACTER(256) :: ERRMESS
      LOGICAL :: E_BDY,S_BDY
!----------------------------------------------------------------------
      IRET=0
  100 FORMAT(' EXIT ',A,' AT NTSD=',I5)
      IEND=ITE
      S_BDY=(JTS==JDS)
      E_BDY=(ITE==IDE-1)
!
      DO K=1,LM
      DO J=JTS,JTE
      IEND=ITE
      IF(E_BDY.AND.MOD(J,2)==0)IEND=ITE-1
!
      DO I=ITS,IEND
        IF(T(I,J,K)>330..OR.T(I,J,K)<180..OR.T(I,J,K)/=T(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,200)I,J,K,T(I,J,K),MYPE,NTSD
  200     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' T=',E12.5      &
                ,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
!         WRITE(ERRMESS,205)NAME,T(I,J,K),I,J,K,MYPE
  205     FORMAT(' EXIT ',A,' TEMPERATURE=',E12.5                      &
                ,' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ELSEIF(Q(I,J,K)<-1.E-4.OR.Q(I,J,K)>30.E-3                      &
               .OR.Q(I,J,K)/=Q(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,300)I,J,K,Q(I,J,K),MYPE,NTSD
  300     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' Q=',E12.5      &
                ,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
!         WRITE(ERRMESS,305)NAME,Q(I,J,K),I,J,K,MYPE
  305     FORMAT(' EXIT ',A,' SPEC HUMIDITY=',E12.5                    &
                ,' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ELSEIF(PINT(I,J,K)<0..OR.PINT(I,J,K)/=PINT(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,315)I,J,K,PINT(I,J,K),MYPE,NTSD
  315     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' PINT=',E12.5      &
                ,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ELSEIF(W(I,J,K)/=W(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,325)I,J,K,W(I,J,K),MYPE,NTSD
  325     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' W=',E12.5      &
                ,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!
      DO K=1,LM
      DO J=JTS,JTE
      IEND=ITE
      IF(E_BDY.AND.MOD(J,2)==1)IEND=ITE-1
      DO I=ITS,IEND
        IF(ABS(U(I,J,K))>125..OR.ABS(V(I,J,K))>125.                    &
               .OR.U(I,J,K)/=U(I,J,K).OR.V(I,J,K)/=V(I,J,K))THEN
          WRITE(0,100)NAME,NTSD
          WRITE(0,400)I,J,K,U(I,J,K),V(I,J,K),MYPE,NTSD
  400     FORMAT(' BAD VALUE I=',I3,' J=',I3,' K=',I2,' U=',E12.5      &
                ,' V=',E12.5,' MYPE=',I3,' NTSD=',I5)
          IRET=666
          return
!         WRITE(ERRMESS,405)NAME,U(I,J,K),V(I,J,K),I,J,K,MYPE
  405     FORMAT(' EXIT ',A,' U=',E12.5,' V=',E12.5                    &
                ,' AT (',I3,',',I2,',',I3,')',' MYPE=',I3)
!         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!----------------------------------------------------------------------
      END SUBROUTINE EXIT
!----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
      SUBROUTINE FIELD_STATS(FIELD,MYPE,MPI_COMM_COMP,LM &
                            ,ITS,ITE,JTS,JTE &
                            ,IMS,IME,JMS,JME &
                            ,IDS,IDE,JDS,JDE)
!----------------------------------------------------------------------
!***
!***  GENERATE STANDARD STATISTICS FOR THE DESIRED FIELD.
!***
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      INCLUDE "../../../../inc/kind.inc"
      INCLUDE "mpif.h"
!----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),INTENT(IN) :: LM,MPI_COMM_COMP,MYPE
      INTEGER(KIND=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE &
                                      ,IMS,IME,JMS,JME &
                                      ,IDS,IDE,JDS,JDE
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM) &
                                                  ,INTENT(IN) :: FIELD
!
!----------------------------------------------------------------------
!***  LOCAL
!----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: I,IRTN,I_BY_J,J,K,KFLIP
!
      REAL(KIND=KFPT) :: FIJK,FMAXK,FMINK
      REAL(KIND=KDBL) :: F_MEAN,POINTS,RMS,ST_DEV,SUMFK,SUMF2K
      REAL(KIND=KFPT),DIMENSION(1:LM) :: FMAX,FMAX_0,FMIN,FMIN_0
      REAL(KIND=KDBL),DIMENSION(1:LM) :: SUMF,SUMF_0,SUMF2,SUMF2_0
!----------------------------------------------------------------------
!
      I_BY_J=(IDE-IDS+1)*(JDE-JDS+1)
!
      layer_loop:  DO K=1,LM
!
        FMAXK=-1.E10
        FMINK=1.E10
        SUMFK=0.
        SUMF2K=0.
!        write(0,*) 'FIELD(ITS+5,JTS+5,K): ', ITS+5, JTS+5, K, FIELD(ITS+5,JTS+5,K)

!
        DO J=JTS,JTE
          DO I=ITS,ITE
            FIJK=FIELD(I,J,K)
            FMAXK=MAX(FMAXK,FIJK)
            FMINK=MIN(FMINK,FIJK)
            SUMFK=SUMFK+FIJK
            SUMF2K=SUMF2K+FIJK*FIJK
          ENDDO
        ENDDO
!
        FMAX(K)=FMAXK
        FMIN(K)=FMINK
        SUMF(K)=SUMFK
        SUMF2(K)=SUMF2K
!
      ENDDO layer_loop
!
!----------------------------------------------------------------------
!***  GLOBAL STATS
!----------------------------------------------------------------------
!
      CALL MPI_REDUCE(SUMF,SUMF_0,LM,MPI_REAL8,MPI_SUM,0              &
                     ,MPI_COMM_COMP,IRTN)
      CALL MPI_REDUCE(SUMF2,SUMF2_0,LM,MPI_REAL8,MPI_SUM,0            &
                     ,MPI_COMM_COMP,IRTN)
      CALL MPI_REDUCE(FMAX,FMAX_0,LM,MPI_REAL,MPI_MAX,0               &
                     ,MPI_COMM_COMP,IRTN)
      CALL MPI_REDUCE(FMIN,FMIN_0,LM,MPI_REAL,MPI_MIN,0               &
                     ,MPI_COMM_COMP,IRTN)
!
      IF(MYPE==0)THEN
        POINTS=I_BY_J
        DO K=1,LM
          F_MEAN=SUMF_0(K)/POINTS
          ST_DEV=SQRT((POINTS*SUMF2_0(K)-SUMF_0(K)*SUMF_0(K))/         &
                      (POINTS*(POINTS-1)))
          RMS=SQRT(SUMF2_0(K)/POINTS)
          WRITE(0,101)K,FMAX_0(K),FMIN_0(K)
          WRITE(0,102)F_MEAN,ST_DEV,RMS
  101     FORMAT(' LAYER=',I2,' MAX=',E13.6,' MIN=',E13.6)
  102     FORMAT(9X,' MEAN=',E13.6,' STDEV=',E13.6,' RMS=',E13.6)
        ENDDO
      ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE FIELD_STATS
!
!----------------------------------------------------------------------
!
      END MODULE MODULE_DIAGNOSE
!
!----------------------------------------------------------------------
