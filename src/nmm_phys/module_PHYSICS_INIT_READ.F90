!----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_INIT_READ
!
!----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_NEMSIO
      USE MODULE_INCLUDE
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY: PHYSICS_INTERNAL_STATE
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,DSTRB,IDSTRB                        &
                                   ,MYPE_SHARE,NUM_TILES
!
      USE MODULE_MICROPHYSICS_NMM
      USE MODULE_CONSTANTS,ONLY : G
      USE MODULE_EXCHANGE
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      CONTAINS
!
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_GWD(INFILE,NGWD,INT_STATE                &
                                 ,MYPE,MPI_COMM_COMP                   &
                                 ,IDS,IDE,JDS,JDE)
!----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER,INTENT(IN) :: NGWD,MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE                         
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: INFILE
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER,INTENT(IN)  :: INT_STATE        !<-- The physics internal state
!
!---------------------
!***  Local variables
!---------------------
!
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP_GWD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ALLOCATE(TEMP_GWD(IDS:IDE,JDS:JDE))      
!
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        OPEN(unit=NGWD,file=INFILE,status='old',form='unformatted')
      ENDIF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSTDV,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HCNVX,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYSW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYNW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENSW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENNW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANGL,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANIS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSLOP,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HZMAX,1,1,1,1,1)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CLOSE(NGWD)
      ENDIF
!
      DEALLOCATE(TEMP_GWD)
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_GWD 
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_INPUT_BINARY(INFILE,NFCST                 &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IDAT,IHRST,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: INFILE
!
      INTEGER,INTENT(IN) :: NFCST,MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME,NSOIL
!
      INTEGER,INTENT(OUT) :: IHRST,IRTN
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      REAL,INTENT(OUT) :: PT
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                    !<-- The physics internal state
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: IM,JM,LMM,LNSH
      INTEGER :: LPT2
      INTEGER :: N,I,J,L,K,II,JJ
      REAL,DIMENSION(LM+1) :: PSG1
      REAL :: PDTOP
      REAL,DIMENSION(LM) :: DSG1,DSG2,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: SG1,SG2,SGM
      REAL,DIMENSION(LM) :: PDSG1,PSGML1
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ITEMP
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP1
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: TEMPSOIL
      REAL,DIMENSION(NSOIL)                :: SOIL1DIN
      LOGICAL :: RUN
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted')
      ENDIF
!
      IF(MYPE==0)THEN
        READ(NFCST)RUN,IDAT,IHRST
        READ(NFCST)PT,PDTOP,LPT2,SGM,SG1,DSG1,SGML1,SG2,DSG2,SGML2
        READ(NFCST) ! I_PARENT_START,J_PARENT_START
        READ(NFCST)int_state%DLMD,int_state%DPHD                       &
                  ,int_state%WBD,int_state%SBD                         &
                  ,int_state%TLM0D,int_state%TPH0D
        READ(NFCST)IM,JM,LMM,LNSH
      ENDIF
!
      CALL MPI_BCAST(IDAT           ,3    ,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHRST          ,1    ,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!
      CALL MPI_BCAST(PT             ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(PDTOP          ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(LPT2           ,1    ,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGM(1)         ,LM+1 ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG1(1)         ,LM+1 ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG1(1)        ,LM   ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML1(1)       ,LM   ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG2(1)         ,LM+1 ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG2(1)        ,LM   ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML2(1)       ,LM   ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
!
      CALL MPI_BCAST(int_state%SBD  ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%WBD  ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%DPHD ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%DLMD ,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TPH0D,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TLM0D,1    ,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
!
      IF(MYPE==0)THEN
        write(0,11121)int_state%wbd,int_state%sbd
        write(0,11122)int_state%dlmd,int_state%dphd
        write(0,11123)int_state%tlm0d,int_state%tph0d
11121   format(' physics input read in wbd=',e12.5,' sbd=',e12.5)
11122   format(' dlmd=',e12.5,' dphd=',e12.5)
11123   format(' tlm0d=',e12.5,' tph0d=',e12.5)
      ENDIF
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!***  Now that the start data has been read from the input file,
!***  check to see if it agrees with the start time from the
!***  configure file.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        IF(IDAT(2)/=int_state%START_MONTH.OR.                           &
           IDAT(1)/=int_state%START_DAY.OR.                             &
           IDAT(3)/=int_state%START_YEAR.OR.                            &
           IHRST  /=int_state%START_HOUR)THEN
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' DATES IN INPUT FILE AND CONFIGURE FILE DISAGREE!!'
          WRITE(0,*)' INPUT: HOUR=',IHRST,' DAY=',IDAT(1)               &
                    ,' MONTH=',IDAT(2),' YEAR=',IDAT(3)
          WRITE(0,*)' CONFIG: HOUR=',int_state%START_HOUR               &
                            ,' DAY=',int_state%START_DAY                &
                          ,' MONTH=',int_state%START_MONTH              &
                           ,' YEAR=',int_state%START_YEAR
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
!-----------------------------------------------------------------------
!***  Before moving on, transfer values to the internal state.
!-----------------------------------------------------------------------
!
        int_state%PT=PT
        int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%SGML2(L)=SGML2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
      ALLOCATE(TEMP1(IDS:IDE,JDS:JDE),STAT=I)
!
!-----------------------------------------------------------------------
!***  Proceed with getting fields from input file.
!***  NOTE: Two records were already read at the top of this routine.
!-----------------------------------------------------------------------
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!
!-----------------------------------------------------------------------
!***  STDH (Standard Deviation of the Height)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%STDH(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%STDH,1,1,1,1,1)
      CALL HALO_EXCH(int_state%STDH,1,3,3)
!
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, CW
!-----------------------------------------------------------------------
!
      DO K=1,LM

        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! U
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! V
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! T
!        write(0,*) 'min max of T read in: ', K, minval(TEMP1), maxval(TEMP1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
!
      ENDDO
      DO K=1,LM
        JJ=LDIM2-1
        DO J=JMS,JME
          JJ=JJ+1
          II=LDIM1-1
          DO I=IMS,IME
            II=II+1
            int_state%WATER(II,JJ,K,int_state%P_QV)=                      & ! WR F water array uses mixing ratio for vapor
                      int_state%Q(II,JJ,K)/(1.-int_state%Q(II,JJ,K))
          ENDDO
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CWM
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! O3
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
!         int_state%O3(I,J,K)=0.
        ENDDO
        ENDDO
!
!       CALL DSTRB(TEMP1,int_state%O3,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!****  EPSR
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max of EPSR: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!*** SNOW ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max of MXSNAL: ', minval(TEMP1), maxval(TEMP1)
      ENDIF

      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SST/TSK
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1        ! actually NMM_TSK from WRF
!        write(0,*) 'min, max of TSKIN: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1        ! actually SST from WRF
!        write(0,*) 'min, max of SST: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SNO, SICE, STC, SMC, ISLTYP, IVGTYP, VEGFRC
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max of SNO: ', minval(TEMP1), maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max for SR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max for USTAR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
!        write(0,*) 'min, max for Z0: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
      IF(MYPE==0)THEN 
        READ(NFCST)TEMP1 
!        write(0,*) 'min, max for Z0BASE: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF 
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1) 
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3) 
!
      ALLOCATE(TEMPSOIL(1:NSOIL,IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
!        write(0,*) 'min, max for STC: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,4),1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
!        write(0,*) 'min, max for SMC: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SMC(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SMC(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SMC(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SMC(:,:,4),1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST)TEMPSOIL
!        write(0,*) 'min, max for SH2O: ', minval(TEMPSOIL),maxval(TEMPSOIL)
      ENDIF
!
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SH2O(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SH2O(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SH2O(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SH2O(:,:,4),1,1,1,1,1)
!
      DEALLOCATE(TEMPSOIL)
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
!        write(0,*) 'min, max for ISLTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
!        write(0,*) 'min, max for IVGTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
!        write(0,*) 'min, max for VEGFRC: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        READ(NFCST) SOIL1DIN ! DZSOIL
      ENDIF
!
      IF(MYPE==0)THEN
        READ(NFCST) SOIL1DIN ! SLDPTH
      ENDIF
!
!     IF(MYPE==0)THEN        ! here will be 14 orography fields for GWD
!       DO N=1,14
!         READ(NFCST) TEMP1
!       ENDDO
!     ENDIF
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
      IF(MYPE==0)THEN
        CLOSE(NFCST)
      ENDIF
!----------------------------------------------------------------------
!
      DEALLOCATE(ITEMP)
      DEALLOCATE(TEMP1)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_INPUT_BINARY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_RESTT_BINARY(INFILE,NFCST                 &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IYEAR_FCST                   &
                                          ,IMONTH_FCST                  &
                                          ,IDAY_FCST                    &
                                          ,IHOUR_FCST                   &
                                          ,IMINUTE_FCST                 &
                                          ,SECOND_FCST                  &
                                          ,IHRST,IDAT,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument variables
!-----------------------
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: INFILE
      INTEGER,INTENT(IN) :: NFCST,MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME
      INTEGER,INTENT(OUT) :: IYEAR_FCST,IMONTH_FCST,IDAY_FCST           &
                           ,IHOUR_FCST,IMINUTE_FCST,IHRST
      INTEGER,INTENT(OUT) :: NSOIL
      REAL,INTENT(OUT) :: SECOND_FCST
!
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      REAL,INTENT(OUT) :: PT
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                 
!
      INTEGER,INTENT(OUT) :: IRTN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: N,I,J,K,L,LPT2
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ITEMP
!
      REAL :: PDTOP
      REAL,DIMENSION(LM) :: DSG1,DSG2,PDSG1,PSGML1,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM
      REAL,ALLOCATABLE,DIMENSION(:)  :: SLDPTH
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP1
!
      LOGICAL :: RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First we need the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted')
      ENDIF
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST) IYEAR_FCST
        READ(NFCST) IMONTH_FCST
        READ(NFCST) IDAY_FCST
        READ(NFCST) IHOUR_FCST
        READ(NFCST) IMINUTE_FCST
        READ(NFCST) SECOND_FCST
        READ(NFCST) ! NTSD
        READ(NFCST) ! IM
        READ(NFCST) ! JM
        READ(NFCST) ! LM
        READ(NFCST) IHRST
        READ(NFCST) ! I_PAR_STA
        READ(NFCST) ! J_PAR_STA
        READ(NFCST) LPT2
        READ(NFCST) IDAT
      ENDIF
!
      IF(MYPE==0)THEN
        write(0,*)'**** read in physics ****************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',int_state%dt
        write(0,*)' Steps/hour =',3600./int_state%dt
        write(0,*)'*************************************'
      ENDIF
!
      CALL MPI_BCAST(IYEAR_FCST   ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IMONTH_FCST  ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IDAY_FCST    ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHOUR_FCST   ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IMINUTE_FCST ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SECOND_FCST  ,1,MPI_REAL    ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IHRST        ,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(IDAT,3,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST) ! MP_PHYSICS
        READ(NFCST) ! SF_SURFACE_PHYSICS
        READ(NFCST) NSOIL
        READ(NFCST) ! NPHS
        READ(NFCST) ! NCLOD
        READ(NFCST) ! NHEAT
        READ(NFCST) ! NPREC
        READ(NFCST) ! NRDLW
        READ(NFCST) ! NRDSW
        READ(NFCST) ! NSRFC
      ENDIF
!
      CALL MPI_BCAST(NSOIL,1,MPI_INTEGER ,0,MPI_COMM_COMP,IRTN)
      ALLOCATE(SLDPTH(1:NSOIL))
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST) ! DT
        READ(NFCST) ! DYH
        READ(NFCST) PDTOP
      ENDIF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)PT
        int_state%PT=PT
        READ(NFCST)int_state%TLM0D
        READ(NFCST)int_state%TPH0D
        READ(NFCST) ! TSTART
        READ(NFCST)int_state%DPHD
        READ(NFCST)int_state%DLMD
        READ(NFCST)int_state%SBD
        READ(NFCST)int_state%WBD
        write(0,11131)int_state%wbd,int_state%sbd
        write(0,11132)int_state%dlmd,int_state%dphd
        write(0,11133)int_state%tlm0d,int_state%tph0d
11131   format(' physics restart read wbd=',e12.5,' sbd=',e12.5)
11132   format(' dlmd=',e12.5,' dphd=',e12.5)
11133   format(' tlm0d=',e12.5,' tph0d=',e12.5)
      ENDIF
!
      CALL MPI_BCAST(int_state%PT,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%DPHD,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%DLMD,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TPH0D,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TLM0D,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%SBD,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%WBD,1,MPI_REAL,0,MPI_COMM_COMP,IRTN)
!
      PT=int_state%PT
!
!-----------------------------------------------------------------------
!***  Read from the restart file: Real 1D arrays
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        READ(NFCST) ! DXH
        READ(NFCST) SG1
        READ(NFCST) SG2
        READ(NFCST) DSG1
        READ(NFCST) DSG2
        READ(NFCST) SGML1
        READ(NFCST) SGML2
        READ(NFCST) SGM
!
!-- APHTIM,ARDLW,ARDSW,ASRFC,AVCNVC,AVRAIN are now 2D arrays and not scalars (below)
!
        READ(NFCST) SLDPTH
        READ(NFCST) int_state%MP_RESTART_STATE
        READ(NFCST) int_state%TBPVS_STATE
        READ(NFCST) int_state%TBPVS0_STATE
        READ(NFCST) ! ALL_BC_DATA
      ENDIF
!
      CALL MPI_BCAST(SGM(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG1(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG1(1)  ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML1(1) ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SG2(1)   ,LM+1  ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(DSG2(1)  ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SGML2(1) ,LM    ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
!
      CALL MPI_BCAST(PDTOP    ,1     ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(LPT2     ,1     ,MPI_INTEGER,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(SLDPTH   ,NSOIL ,MPI_REAL   ,0,MPI_COMM_COMP,IRTN)
!
      CALL MPI_BCAST(int_state%MP_RESTART_STATE(1) ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TBPVS_STATE(1)      ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
      CALL MPI_BCAST(int_state%TBPVS0_STATE(1)     ,MICRO_RESTART ,MPI_REAL ,0,MPI_COMM_COMP,IRTN)
!
      DO N=1,NSOIL
        int_state%SLDPTH(N)=SLDPTH(N)
      ENDDO
!
!-----------------------------------------------------------------------
!***  Read from restart file: Logical
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! GLOBAL
        READ(NFCST) RUN
        READ(NFCST) ! ADIABATIC
      ENDIF
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
!
!-----------------------------------------------------------------------
!***  Before moving on, transfer values to the internal state.
!-----------------------------------------------------------------------
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
        int_state%SGML2(L)=SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 2D arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRCV)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRST)
!
      DEALLOCATE(ITEMP)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(TEMP1(IDS:IDE,JDS:JDE),STAT=I)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) !GLAT
        READ(NFCST) !GLON
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) ! VLAT
        READ(NFCST) ! VLON
        READ(NFCST) ! PDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  Skip from restart file: Real 3D arrays (only from DYN)
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! W
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! DWDT
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! PINT
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! OMGALF
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RRW
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! DIV
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RTOP
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCU
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCV
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TCT
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TP
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! UP
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! VP
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! E2
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM-1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! PSGDT
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Z
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRCV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKHS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKMS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%BGROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CLDEFI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CPRATE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPPT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZEN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZMEAN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%GRNFLX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PBLH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Q10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RADOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWINC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEXC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCLHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SIGT4,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTAV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOPCX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SOILTB,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SSROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SUBSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TH10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%U10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%UZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%V10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWTOA,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTFLX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RMOL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%T2,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMIN,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMAX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACUTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%APHTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDLW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDSW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASRFC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVRAIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVCNVC,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, Q2, CW, F_ICE, F_RIMEF, F_RAIN
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CLDFRA
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CLDFRA,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! CWM
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO

        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! Q2
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q2,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RLWTT
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RLWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! RSWTT
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RSWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! T
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! TCUCN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TCUCN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1  ! TRAIN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! U
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! V
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! XLEN_MIX
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%XLEN_MIX,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_ICE
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_ICE,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RIMEF
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RIMEF,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RAIN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SH2O: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SH2O,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SMC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SMC,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for STC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%STC,1,1,1,NSOIL,K)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  TRACERS
!-----------------------------------------------------------------------
      DO N=int_state%INDX_RRW+1,int_state%NUM_TRACERS_TOTAL                !<-- The first 'INDX_RRW' arrays are unallocated pointers
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TRACERS
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRACERS(I,J,K,N)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRACERS(:,:,:,N),1,1,1,LM,K)
      ENDDO
      ENDDO
!rv do not use HALO_EXCH
!-----------------------------------------------------------------------
!
      DEALLOCATE(TEMP1)
      DEALLOCATE(SLDPTH)
!
!-----------------------------------------------------------------------
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
         CLOSE(NFCST)
      ENDIF
!----------------------------------------------------------------------
!-----------------------------------------------------------------------

      END SUBROUTINE PHYSICS_READ_RESTT_BINARY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_INPUT_NEMSIO(INFILE                       &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IDAT,IHRST,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: INFILE
      INTEGER,INTENT(IN) :: MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME,NSOIL
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      INTEGER,INTENT(OUT) :: IHRST
      REAL,INTENT(OUT) :: PT
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                 
!
      INTEGER,INTENT(OUT) :: IRTN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: LPT2
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: N,I,J,L,K,II,JJ
      INTEGER,DIMENSION(:),ALLOCATABLE :: ITEMP
!
      REAL :: PDTOP
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM
      REAL,DIMENSION(LM) :: DSG1,DSG2,PDSG1,PSGML1,SGML1,SGML2
!
      REAL,DIMENSION(:),ALLOCATABLE :: TEMP1
!
      LOGICAL :: RUN
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First we need the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT()
!
      CALL NEMSIO_OPEN(gfile,INFILE,'read',iret=irtn)
      if(irtn/=0) write(0,*)'ERROR: open file ',trim(infile),' has failed'
!
!---------------------------------------------------------------------
!***  Vertical layer information is needed in order to send it to
!***  some specific schemes' initialization routines.
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'PT',int_state%PT,iret=irtn)
      PT=int_state%PT
      CALL NEMSIO_GETHEADVAR(gfile,'RUN',run,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'IHRST',ihrst,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'IDAT',idat,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'PDTOP',pdtop,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'LPT2',lpt2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGM',sgm,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG1',sg1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG1',dsg1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML1',sgml1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG2',sg2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG2',dsg2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML2',sgml2,iret=irtn)
!      write(0,*)'in phys,pt=',pt,'run=',run,'ihrst=',ihrst,'idat=',idat,   &
!         'pdtop=',pdtop,'lpt2=',lpt2
!      write(0,*)'in phys,sgm=',sgm(1:10),'sg1=',sg1(1:10),maxval(sg1),minval(sg1),'sg2=',sg2(1:10), &
!         maxval(sg2),minval(sg2),'dsg1=',dsg1(1:10),'dsg2=',dsg2(1:10),'sgml1=',sgml1(1:10),  &
!         'sgml2=',sgml2(1:10)
!
!-----------------------------------------------------------------------
!***  Before moving on, transfer values to the internal state.
!-----------------------------------------------------------------------
!
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%SGML2(L)=SGML2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
      ALLOCATE(TEMP1((IDE-IDS+1)*(JDE-JDS+1)),STAT=I)
!
!-----------------------------------------------------------------------
!***  Proceed with getting fields from input file.
!***  NOTE: Two records were already read at the top of this routine.
!-----------------------------------------------------------------------
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'fis','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,fis=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3) !zj
!
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sm','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,sm=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'dpres','hybrid sig lev',1,temp1,iret=irtn)
!      write(0,*)'in phys,pd=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, CW
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'ugrd','mid layer',k,temp1,iret=irtn)
!      write(0,*)'in phys,ugrd=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'vgrd','mid layer',k,temp1,iret=irtn)
!      write(0,*)'in phys,vgrd=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.

        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'tmp','mid layer',k,temp1,iret=irtn)
!      write(0,*)'in phys,tmp=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'spfh','mid layer',k,temp1,iret=irtn)
!      write(0,*)'in phys,spfh=',maxval(temp1),minval(temp1)
        ENDIF
!

        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)

!
      DO K=1,LM
        JJ=LDIM2-1
        DO J=JMS,JME
          JJ=JJ+1
          II=LDIM1-1
          DO I=IMS,IME
            II=II+1
            int_state%WATER(II,JJ,K,int_state%P_QV)=                      & ! WRF water array uses mixing ratio for vapor
                      int_state%Q(II,JJ,K)/(1.-int_state%Q(II,JJ,K))
          ENDDO
        ENDDO
      ENDDO

!
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'clwmr','mid layer',k,temp1,iret=irtn)
!      write(0,*)'in phys,clwmr=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'albedo','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,albedo=',maxval(temp1),minval(temp1)
      ENDIF
!
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!
! **** EPSR
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'epsr','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,epsr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!*** SNOW ALBEDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'mxsnal','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,mxsnal=',maxval(temp1),minval(temp1)
      ENDIF

      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SST/TSK
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tskin','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,tsk=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tsea','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,sst=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  SNO, SICE, STC, SMC, ISLTYP, IVGTYP, VEGFRC
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sno','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,sno=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'si','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,si=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sice','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,sice=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tg','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,tg=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cmc','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,cmc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sr','sfc',1,temp1,iret=irtn)
!      write(0,*)'in phys,sr=',maxval(temp1),minval(temp1)
!       write(0,*) 'min, max for SR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'ustar','sfc',1,temp1,iret=irtn)
!       write(0,*) 'min, max for USTAR: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'zorl','sfc',1,temp1,iret=irtn)
!        write(0,*) 'min, max for zorl: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
      IF(MYPE==0)THEN !zj
        CALL NEMSIO_READRECV(gfile,'z0base','sfc',1,temp1,iret=irtn)
!        write(0,*) 'min, max for Z0BASE: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF !zj
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1) !zj
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3) !zj
!
      IF(MYPE==0)THEN !zj
        CALL NEMSIO_READRECV(gfile,'stdh','sfc',1,temp1,iret=irtn)
!        write(0,*) 'min, max for STDH: ', minval(TEMP1),maxval(TEMP1) !zj
      ENDIF !zj
      CALL DSTRB(TEMP1,int_state%STDH,1,1,1,1,1) !zj
      CALL HALO_EXCH(int_state%STDH,1,3,3) !zj
!
      DO L=1,NSOIL
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'stc','soil layer',l,temp1,iret=irtn)
!          write(0,*) 'min, max for STC: ', minval(TEMP1),maxval(TEMP1)
        ENDIF
        CALL DSTRB(TEMP1,int_state%STC(IMS:IME,JMS:JME,L),1,1,1,1,1)
      ENDDO
!
!
      DO L=1,NSOIL
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'smc','soil layer',l,temp1,iret=irtn)
!          write(0,*) 'min, max for SMC: ', minval(TEMP1),maxval(TEMP1)
        ENDIF
        CALL DSTRB(TEMP1,int_state%SMC(IMS:IME,JMS:JME,L),1,1,1,1,1)
      ENDDO
!
      DO L=1,NSOIL
        IF(MYPE==0)THEN
          CALL NEMSIO_READRECV(gfile,'sh2o','soil layer',l,temp1,iret=irtn)
!          write(0,*) 'min, max for SH2O: ', minval(TEMP1),maxval(TEMP1)
        ENDIF
        CALL DSTRB(TEMP1,int_state%SH2O(IMS:IME,JMS:JME,L),1,1,1,1,1)
      ENDDO
!
      ALLOCATE(ITEMP((IDE-IDS+1)*(JDE-JDS+1)),STAT=I)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sltyp','sfc',1,temp1,iret=irtn)
        ITEMP=INT(temp1)
!        write(0,*) 'min, max for ISLTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'vgtyp','sfc',1,temp1,iret=irtn)
        ITEMP=INT(temp1)
!        write(0,*) 'min, max for IVGTYP: ', minval(ITEMP),maxval(ITEMP)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'vegfrc','sfc',1,temp1,iret=irtn)
!        write(0,*) 'min, max for vegfrc: ', minval(TEMP1),maxval(TEMP1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!
!----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'SBD',int_state%SBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'WBD',int_state%WBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DPHD',int_state%DPHD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DLMD',int_state%DLMD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TPH0D',int_state%TPH0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TLM0D',int_state%TLM0D,iret=irtn)
!
!----------------------------------------------------------------------
!
      CALL NEMSIO_CLOSE(gfile)
!
      CALL NEMSIO_FINALIZE()
!
!----------------------------------------------------------------------
!
      DEALLOCATE(ITEMP)
      DEALLOCATE(TEMP1)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_INPUT_NEMSIO
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_RESTT_NEMSIO(INFILE                       &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IYEAR_FCST,IMONTH_FCST       &
                                          ,IDAY_FCST,IHOUR_FCST         &
                                          ,IMINUTE_FCST,SECOND_FCST     &
                                          ,IHRST,IDAT,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: INFILE
!
      INTEGER,INTENT(IN) :: MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME
      INTEGER,INTENT(OUT) :: IYEAR_FCST,IMONTH_FCST,IDAY_FCST           &
                           ,IHOUR_FCST,IMINUTE_FCST,IHRST
      INTEGER,INTENT(OUT) :: NSOIL
!
      REAL,INTENT(OUT) :: SECOND_FCST
!
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      REAL,INTENT(OUT) :: PT
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE
!
      INTEGER,INTENT(OUT) :: IRTN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: N,I,J,K,L
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: LPT2
      INTEGER,DIMENSION(7) :: FCSTDATE
      INTEGER,DIMENSION(:),ALLOCATABLE :: ITEMP
!
      REAL :: PDTOP
      REAL,DIMENSION(LM) :: DSG1,DSG2,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: SG1,SG2,SGM
      REAL,ALLOCATABLE,DIMENSION(:) :: SLDPTH
      REAL,DIMENSION(LM) :: PDSG1,PSGML1
      REAL,DIMENSION(LM+1) :: PSG1
      REAL,DIMENSION(:),ALLOCATABLE :: TEMP1
!
      LOGICAL :: RUN
!
      CHARACTER(10) :: VARNAME
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
       CALL NEMSIO_INIT()
!
       CALL NEMSIO_OPEN(gfile,INFILE,'read',iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'FCSTDATE',FCSTDATE,iret=irtn)
      IYEAR_FCST=FCSTDATE(1)
      IMONTH_FCST=FCSTDATE(2)
      IDAY_FCST=FCSTDATE(3)
      IHOUR_FCST=FCSTDATE(4)
      IMINUTE_FCST=FCSTDATE(5)
      SECOND_FCST=0.
      if(FCSTDATE(7)/=0) SECOND_FCST=FCSTDATE(6)/(FCSTDATE(7)*1.)
!
      CALL NEMSIO_GETHEADVAR(gfile,'IHRST',IHRST,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'LPT2',LPT2,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'IDAT',IDAT,iret=irtn)
!
      IF(MYPE==0)THEN
        write(0,*)'**** read in physics ****************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',int_state%dt
        write(0,*)' Steps/hour =',3600./int_state%dt
        write(0,*)'*************************************'
      ENDIF
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'NSOIL',NSOIL,iret=irtn)
!      write(0,*)'in rst, fcstdate=',fcstdate,'ihrst=',ihrst,'lpt2=',lpt2, &
!        'idat=',idat,'nsoil=',nsoil
      ALLOCATE(SLDPTH(1:NSOIL))
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'PDTOP',PDTOP,iret=irtn)
!
!-----------------------------------------------------------------------
      CALL NEMSIO_GETHEADVAR(gfile,'PT',PT,iret=irtn)
      int_state%PT=PT
!
      CALL NEMSIO_GETHEADVAR(gfile,'TLM0D',int_state%TLM0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TPH0D',int_state%TPH0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DPHD',int_state%DPHD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DLMD',int_state%DLMD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SBD',int_state%SBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'WBD',int_state%WBD,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'SG1',SG1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG2',SG2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG1',DSG1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG2',DSG2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML1',SGML1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML2',SGML2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGM',SGM,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SLDPTH',SLDPTH,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'MP_RESTART',int_state%MP_RESTART_STATE,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TBPVS_STAT',int_state%TBPVS_STATE,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TBPVS0_STA',int_state%TBPVS0_STATE,iret=irtn)
!
      DO N=1,NSOIL
        int_state%SLDPTH(N)=SLDPTH(N)
      ENDDO
!
!-----------------------------------------------------------------------
!*** Read from restart file: Logical
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'RUN',RUN,iret=irtn)
!
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!
!jw
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO

!
!-----------------------------------------------------------------------
!***  Transfer values to the internal state.
!-----------------------------------------------------------------------
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
        int_state%SGML2(L)=SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO

!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 2D arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(TEMP1((IDE-IDS+1)*(JDE-JDS+1)),STAT=I)
      ALLOCATE(ITEMP((IDE-IDS+1)*(JDE-JDS+1)),STAT=I)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'SLTYP','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sltyp=',maxval(temp1),minval(temp1)
        itemp=nint(temp1)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'VGTYP','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,vgtyp=',maxval(temp1),minval(temp1)
        itemp=nint(temp1)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cfrcv','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cfrcv=',maxval(temp1),minval(temp1)
        itemp=nint(temp1)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRCV)
!
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cfrst','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cfrst=',maxval(temp1),minval(temp1)
        itemp=nint(temp1)
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRST)
!
      DEALLOCATE(ITEMP)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays
!-----------------------------------------------------------------------
!
!      ALLOCATE(TEMP1(IDS:IDE,JDS:JDE),STAT=I)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'fis','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,fis=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1)
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'dpres','hybrid sig lev',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,pd=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acfrcv','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,acfrcv=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRCV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acfrst','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,acfrst=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acprec','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,acprec=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acsnom','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,acsnom=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acsnow','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,acsnow=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'akhs_out','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,akhs_out=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKHS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------

!***  AKMS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'akms_out','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,akms_out=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKMS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'albase','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,albase=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'albedo','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,albedo=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'alwin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,alwin=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'alwout','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,alwout=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'alwtoa','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,alwtoa=',maxval(temp1),minval(temp1)


      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'aswin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,aswin=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'aswout','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,aswout=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'aswtoa','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,aswtoa=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'bgroff','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,bgroff=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%BGROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cfrach','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cfrach=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cfracl','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cfracl=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cfracm','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cfracm=',maxval(temp1),minval(temp1)


      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cldefi','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cldefi=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CLDEFI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cmc','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cmc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cnvbot','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cnvbot=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cnvtop','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cnvtop=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cprate','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cnprate=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CPRATE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cuppt','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cuppt=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPPT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cuprec','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,cuprec=',maxval(temp1),minval(temp1)

      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'czen','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,czen=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZEN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'czmean','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,czmean=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZMEAN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'epsr','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,epsr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'grnflx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,grnflx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%GRNFLX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'hbotd','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,hbotd=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'hbots','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,hbots=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'htopd','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,htopd=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'htops','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,htops=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'mxsnal','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,mxsnal=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'pblh','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,pblh=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%PBLH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'potevp','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,potevp=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'prec','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,prec=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%PREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'pshltr','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,pshltr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%PSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'q10','10 m above gnd',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,q10=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Q10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'qsh','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,qsh=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'qshltr','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,qshltr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'qwbs','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,qwbs=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%QWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'qz0','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,qz0=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%QZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'radot','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,radot=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RADOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rlwin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rlwin=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rlwtoa','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rlwtoa=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rswin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rswin=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rswinc','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rswinc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWINC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rswout','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rswout=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sfcevp','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sfcevp=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sfcexc','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sfcexc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEXC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sfclhx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sfclhx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCLHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sfcshx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sfcshx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'si','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,si=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sice','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sice=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sigt4','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sigt4=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SIGT4,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sm','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sm=',maxval(temp1),minval(temp1)
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'smstav','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,smstav=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTAV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'smstot','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,smstot=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sno','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sno=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'snopcx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,snopcx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOPCX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'soiltb','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,soiltb=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SOILTB,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sr','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'ssroff','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,ssroff=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SSROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tsea','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,sst=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'subshx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,subshx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%SUBSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tg','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,tg=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'th10','10 m above gnd',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,th10=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TH10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'ths','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,ths=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%THS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'thz0','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,thz0=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%THZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tshltr','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,tshltr=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'twbs','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,twbs=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'u10','10 m above gnd',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,u10=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%U10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'uustar','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,ustar=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'uz0','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,uz0=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%UZ0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%UZ0,1,3,3)
!     must do HALO_EXCH
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'v10','10 m above gnd',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,v10=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%V10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'vegfrc','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,vegfrc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'vz0','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,vz0=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%VZ0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%VZ0,1,3,3)
!     must do HALO_EXCH
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'zorl','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,z0=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tskin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,tskin=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'akhs','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,akhs=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'akms','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,akms=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'hbot','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,hbot=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'htop','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,htop=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rswtoa','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rswtoa=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWTOA,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'potflx','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,potflx=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTFLX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rmol','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,rmol=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%RMOL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'t2','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,t2=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%T2,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'z0base','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,z0base=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tlmin','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,z0base=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tlmax','sfc',1,temp1,iret=irtn)
!        write(0,*)'read rst phys,z0base=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMAX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'acutim','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,acutim=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACUTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'aphtim','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,aphtim=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%APHTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'ardlw','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,ardlw=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDLW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'ardsw','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,ardsw=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDSW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'asrfc','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,asrfc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASRFC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'avrain','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,avrain=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVRAIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'avcnvc','sfc',1,temp1,iret=irtn)
        write(0,*)'read rst phys,avcnvc=',maxval(temp1),minval(temp1)
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVCNVC,1,1,1,1,1)
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: REAL 3D ARRAYS
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, Q2, CW, F_ICE, F_RIMEF, F_RAIN
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'cldfra','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,cldfra=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CLDFRA,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'clwmr','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,clwmr=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'spfh','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,spfh=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=LDIM2,UDIM2
        DO I=LDIM1,UDIM1
          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'q2','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,q2=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%Q2,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rlwtt','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,rlwtt=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RLWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'rswtt','mid layer',k,temp1,iret=irtn)

!        write(0,*)'read rst phys,rswtt=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%RSWTT,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tmp','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,tmp=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'tcucn','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,tcucn=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TCUCN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'train','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,train=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
       CALL NEMSIO_READRECV(gfile,'ugrd','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,ugrd=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'vgrd','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,vgrd=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'xlen_mix','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,xlen_mix=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%XLEN_MIX,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'f_ice','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,f_ice=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_ICE,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'f_rimef','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,f_rimef=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RIMEF,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'f_rain','mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,f_rrain=',maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'sh2o','soil layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,sh2o=',maxval(temp1),minval(temp1)
!          write(0,*) 'lev, min, max for SH2O: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SH2O,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'smc','soil layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,smc=',maxval(temp1),minval(temp1)
!          write(0,*) 'lev, min, max for SMC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SMC,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
        CALL NEMSIO_READRECV(gfile,'stc','soil layer',k,temp1,iret=irtn)
!        write(0,*)'read rst phys,stc=',maxval(temp1),minval(temp1)
!          write(0,*) 'lev, min, max for STC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%STC,1,1,1,NSOIL,K)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  TRACERS
!-----------------------------------------------------------------------
      DO N=int_state%INDX_RRW+1,int_state%NUM_TRACERS_TOTAL                !<-- The first 'INDX_RRW' arrays are unallocated pointers
      DO K=1,LM
        IF(MYPE==0)THEN
        write(varname,'(a8,I2.2)')'tracers_',N
        CALL NEMSIO_READRECV(gfile,varname,'mid layer',k,temp1,iret=irtn)
!        write(0,*)'read rst,name=',varname,maxval(temp1),minval(temp1)
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRACERS(I,J,K,N)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%TRACERS(:,:,:,N),1,1,1,LM,K)
      ENDDO
      ENDDO
!rv do not use HALO_EXCH
!-----------------------------------------------------------------------
!
      DEALLOCATE(TEMP1)
      DEALLOCATE(SLDPTH)
!
!-----------------------------------------------------------------------
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
      CALL NEMSIO_CLOSE(GFILE)
!-----------------------------------------------------------------------
!
      CALL NEMSIO_FINALIZE()
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      END SUBROUTINE PHYSICS_READ_RESTT_NEMSIO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
      END MODULE MODULE_PHYSICS_INIT_READ
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
