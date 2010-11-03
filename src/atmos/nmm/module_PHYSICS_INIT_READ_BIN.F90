!----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_INIT_READ_BIN
!
!----------------------------------------------------------------------
!
      USE ESMF_MOD
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
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER,INTENT(INOUT) :: INT_STATE     !<-- The physics internal state
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: IERR,RC
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
        OPEN(unit=NGWD,file=INFILE,status='old',form='unformatted'      &
            ,iostat=IERR)
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(INFILE)                &
                   ,' in PHYSICS_READ_GWD'
          WRITE(0,*)' ABORTING!'
          CALL ESMF_FINALIZE(terminationflag=ESMF_ABORT                 &
                            ,rc             =RC)
        ENDIF
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
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
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
      INTEGER :: N,I,IERR,J,L,K,II,JJ,RC
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
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted'     &
            ,iostat=IERR)
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(INFILE)                &
                   ,' in PHYSICS_READ_INPUT_BINARY'
          WRITE(0,*)' ABORTING!'
          CALL ESMF_FINALIZE(terminationflag=ESMF_ABORT                 &
                            ,rc             =RC)
        ENDIF
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
      CALL HALO_EXCH(int_state%PD,1,2,2)
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
      CALL HALO_EXCH(int_state%U,LM,2,2)
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
      CALL HALO_EXCH(int_state%V,LM,2,2)
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
      CALL HALO_EXCH(int_state%T,LM,2,2)
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
      CALL HALO_EXCH(int_state%Q,LM,2,2)

      DO K=1,LM
        JJ=LDIM2-1
        DO J=JMS,JME
          JJ=JJ+1
          II=LDIM1-1
          DO I=IMS,IME
            II=II+1
!d            int_state%WATER(II,JJ,K,int_state%P_QV)=                      & ! WR F water array uses mixing ratio for vapor
!d                      int_state%Q(II,JJ,K)/(1.-int_state%Q(II,JJ,K))
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
      CALL HALO_EXCH(int_state%CW,LM,2,2)
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
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
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
      INTEGER :: IERR,RC
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
        OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted'     &
            ,iostat=IERR)
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(INFILE)                &
                   ,' in PHYSICS_READ_RESTT_BINARY'
          WRITE(0,*)' ABORTING!'
          CALL ESMF_FINALIZE(terminationflag=ESMF_ABORT                 &
                            ,rc             =RC)
        ENDIF
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
      CALL HALO_EXCH(int_state%PD,1,2,2)
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
          READ(NFCST)TEMP1   ! OMGALF
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! O3
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
!-----------------------------------------------------------------------
      DO N=1,int_state%INDX_O3
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! TRACERS_PREV
        ENDIF
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
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
!d          int_state%CW(I,J,K)=0.
        ENDDO
        ENDDO

!d        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,K)
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
!d          int_state%Q(I,J,K)=0.
        ENDDO
        ENDDO
!
!d        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,K)
      ENDDO
!d      CALL HALO_EXCH(int_state%Q,LM,2,2)
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
      CALL HALO_EXCH(int_state%Q2,LM,2,2)
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
      DO K=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! PINT
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! DWDT
        ENDIF
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
      CALL HALO_EXCH(int_state%T,LM,2,2)
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
      CALL HALO_EXCH(int_state%U,LM,2,2)
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
      CALL HALO_EXCH(int_state%V,LM,2,2)
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
      DO N=int_state%INDX_O3+1,int_state%NUM_TRACERS_TOTAL                !<-- The first 'INDX_O3' arrays are unallocated pointers
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
      CALL HALO_EXCH(int_state%TRACERS,LM,int_state%NUM_TRACERS_TOTAL,1,2,2)

        do n=1,int_state%num_water
        do l=1,lm
          do j=jms,jme
          do i=ims,ime
!d            int_state%water(i,j,l,n)=int_state%tracers(i,j,l,n+int_state%num_tracers_total-int_state%num_water)
          enddo
          enddo
        enddo
        enddo
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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
      END MODULE MODULE_PHYSICS_INIT_READ_BIN
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
