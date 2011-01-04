!-----------------------------------------------------------------------------
     MODULE MODULE_MP_GFS
!-----------------------------------------------------------------------------
!
!    12-10-2010   Created by Weiguo Wang
       USE MACHINE , ONLY : kind_phys
       USE FUNCPHYS , ONLY : gpvs, fpvs
       IMPLICIT NONE
       REAL, Private, PARAMETER ::                   & 
                                g99=9.80665,         &
                                t0c=273.15,          &
                                t_ice=-40.0+t0c

      PUBLIC :: GFSMP, GFSMP_INIT
      CONTAINS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE GFSMP    (DT,                                         &
                           dz8w,rho_phy,p_phy,pi_phy,th_phy,           &
                           SR,QT,F_ICE_phy,                            &
                           RAINNC,RAINNCV,                             &
                           WATER,P_QV,P_QC,P_QI,NUM_WATER,             &
                           TP1,TP2,QP1,QP2,PSP1,PSP2,                  &
                           ids,ide, jds,jde, kds,kde,                  &
                           ims,ime, jms,jme, kms,kme,                  &
                           its,ite, jts,jte, kts,kte )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
                           ,IMS,IME,JMS,JME,KMS,KME                     &
                           ,ITS,ITE,JTS,JTE,KTS,KTE                     

      REAL, INTENT(IN) 	    :: DT
      REAL, INTENT(IN),     DIMENSION(ims:ime, jms:jme, kms:kme)::      &
                           dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, jms:jme, kms:kme)::      &
                           th_phy,F_ICE_phy, QT
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
                                                         RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR
!
      INTEGER,INTENT(IN) ::  NUM_WATER 
      INTEGER,INTENT(IN) ::  P_QV,P_QC,P_QI
      REAL,DIMENSION(IMS:IME,JMS:JME,1:KTE,NUM_WATER),INTENT(INOUT) :: WATER
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           :: PSP1,PSP2
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme,1:KTE)     :: TP1,TP2,QP1,QP2 

!-----------------------------------------------------------------------
!     LOCAL VARS
!-----------------------------------------------------------------------

       REAL(kind=kind_phys), PARAMETER    :: rh00 = 0.850
       REAL(kind=kind_phys), PARAMETER    ::  psautco = 4.0E-4    &    ! Zhao scheme default opr value
                             ,prautco = 1.0E-4    &    ! Zhao scheme default opr value
                             ,evpco  = 2.0E-5

!     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
!     the microphysics scheme. Instead, they will be used by Eta precip 
!     assimilation.

      REAL,  DIMENSION( ims:ime, jms:jme,kms:kme ) ::                  &
            TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      !REAL,  DIMENSION(its:ite, jts:jte, kts:kte):: t_phy
      REAL,  DIMENSION( ims:ime, jms:jme,kms:kme ) :: t_phy
      INTEGER :: I,J,K, KFLIP, L, KM
!!!LOCAL VARS FOR GFS MICROPHY
       Integer, parameter :: IX=1, IM=1, ipr=1
       REAL(kind=kind_phys), DIMENSION(IX,KTE) :: PRSL, Q_COL, CWM_COL,T_COL, RHC, DELP &
                                 ,TP1_COL,QP1_COL,TP2_COL,QP2_COL
       REAL(kind=kind_phys), DIMENSION(IX) :: PS ,rain1, psp1_1,psp2_1
       REAL(kind=kind_phys) :: dtp,frain,fice
       INTEGER, dimension(ix,kte) :: IW            !! ice flag
       logical lprnt
       REAL(kind=kind_phys), DIMENSION(IX,KTE) :: RAINP            ! not in use
       logical diag
        lprnt = .false.
        diag = .false.  !.true.
!------------------------------------------------------------------------
!**********************************************************************
        KM = KTE
        DTP = 2.0*DT
        frain = DT/DTP
        RHC(:,:) = rh00

               DO j = jts,jte
               DO k = kts,kte
               DO i = its,ite
                t_phy(i,j,k) = th_phy(i,j,k)*pi_phy(i,j,k)
               ENDDO
               ENDDO
               ENDDO

      DO j = jts,jte
       DO i = its,ite
         ACPREC(i,j)=0.
         APREC (i,j)=0.
         PREC  (i,j)=0.
         SR    (i,j)=0.
       ENDDO
       DO k = kts,kte
       DO i = its,ite
	 TLATGS_PHY (i,j,k)=0.
	 TRAIN_PHY  (i,j,k)=0.
       ENDDO
       ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!-- Start of original driver for EGCP01COLUMN
!-----------------------------------------------------------------------
       DO J=JTS,JTE    
        DO I=ITS,ITE  
          DO K=KTS,KTE
            KFLIP = KTE + 1 -K

            DELP(IX,KFLIP)=RHO_PHY(I,J,K)*g99*dz8w(I,J,K)
            PRSL(IX,KFLIP)=P_phy(I,J,K)
            Q_COL(IX,KFLIP) = WATER(I,J,K,P_QV)/(1.0+WATER(I,J,K,P_QV) )   !! to specific humidity        
            CWM_COL(IX,KFLIP)=WATER(I,J,K,P_QC) &
                             +WATER(I,J,K,P_QI) 
            T_COL(IX,KFLIP) = t_phy(i,j,k) 

            TP1_COL(IX,K) = tp1(i,j,k)
            QP1_COL(IX,K) = qp1(i,j,k)
            TP2_COL(IX,K) = tp2(i,j,k)
            QP2_COL(IX,K) = qp2(i,j,k)
            psp1_1(IX)   = psp1(i,j)
            psp2_1(IX)   = psp2(i,j)
            iw(ix,KFLIP) = 0 !F_ICE_phy(I,J,K) 
          ENDDO

            PS(IX) = DELP(IX,1)*0.5+PRSL(IX,1)
            rain1(ix) = 0.0

               IF(DIAG .and. i == 10 .and. j == 10 ) THEN
                  write(0,*)'before calling MICRO'
                  write(0,*)'tp1=',tp1(i,j,:)
                  write(0,*)'qp1=',qp1(i,j,:)
                  write(0,*)'psp1=',psp1(i,j)
                  write(0,*)'qt=',qt(i,j,:)
                  write(0,*)'cwm=',cwm_col(1,:)
                  write(0,*)'water(p_qc)=',water(i,j,:,p_qc)
                  write(0,*)'water(p_qi)=',water(i,j,:,p_qi)
               ENDIF    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CALL MICROPHY                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call gscond(im, ix, KM, dtp, prsl, ps,                       &
                       Q_COL, CWM_COL, T_COL,                           &
                       TP1_COL, QP1_COL,  PSP1_1,                     &
                       TP2_COL, QP2_COL,  PSP2_1,                     &
                       rhc,lprnt, ipr)
           call precpd(im, ix, KM, dtp, delp, prsl, ps,                 &
                       Q_COL, CWM_COL, T_COL, rain1,                    &
                       rainp, rhc, psautco, prautco, evpco,             &
                       lprnt, ipr)
 300      continue 

            DO K=kts,kte
             tp1(i,j,k) = tp1_COL(ix,k)
             qp1(i,j,k) = qp1_COL(ix,k)
             tp2(i,j,k) = tp2_COL(ix,k)
             qp2(i,j,k) = qp2_COL(ix,k)
             psp1(i,j)  = psp1_1(ix)
             psp2(i,j)  = psp2_1(ix)
            ENDDO
!#######################################################################
!
!--- Update storage arrays
!
          DO L=1,KM
            KFLIP = KM + 1 - L
            TRAIN_phy(I,J,L)= (T_col(1,KFLIP)-T_phy(I,J,L))/DTp
            TLATGS_phy(I,J,L)=(T_col(1,KFLIP)-T_phy(I,J,L))*frain
          ENDDO
          DO K=1,KM
            KFLIP=KM + 1 - K 
            T_phy(I,J,K)=T_col(1,KFLIP)
            WATER(I,J,k,P_QV)= Q_col(1,KFLIP)/(1.0-Q_COL(1,KFLIP) )
              fice=1.0
              IF(T_COL(1,KFLIP) .GT. t_ice .and.                   &
                 T_COL(1,KFLIP) .LE. t0c ) THEN
                 fice = 1.0 - (T_COL(1,KFLIP)-t_ice)/(t0c-t_ice)
              ENDIF
              IF(T_COL(1,KFLIP) .GT. t0c ) fice=0.0
      !      fice = float( IW(1,KFLIP) )
            WATER(I,J,K,P_QC) = CWM_COL(1,KFLIP)*(1.0-fice)
            WATER(I,J,K,P_QI) = CWM_COL(1,KFLIP)*fice
            QT(I,J,K) = CWM_COL(1,KFLIP)
            F_ICE_phy(I,J,K) = fice
          ENDDO
!
!--- Update accumulated precipitation statistics
!
!--- Surface precipitation statistics; SR is fraction of surface 
!    precipitation (if >0) associated with snow
!
     !!   APREC(I,J)=(ARAIN+ASNOW)*RRHOL       ! Accumulated surface precip (depth in m)  !<--- Ying
        APREC(I,J)=rain1(1)*frain       ! Accumulated surface precip (depth in m)  !<--- Ying
        PREC(I,J)=PREC(I,J)+APREC(I,J)
        ACPREC(I,J)=ACPREC(I,J)+APREC(I,J)
    !    IF(APREC(I,J) .LT. 1.E-8) THEN
    !      SR(I,J)=0.
    !    ELSE
    !      SR(I,J)=RRHOL*ASNOW/APREC(I,J)
    !    ENDIF
!
!
       IF(DIAG .and. i == 10 .and. j == 10 ) THEN
          write(0,*)'RAIN=',APREC(I,J)
        !i  write(0,*)'FICE=',F_ICE_phy(I,J,:)
        !  write(0,*)'PRSL=',PRSL
          write(0,*)'DELP=',DELP  
        !  write(0,*)'Q=',Q_COL
        !  write(0,*)'CWM=',CWM_COL
        !  write(0,*)'T=',T_COL
          write(0,*)'PS=',PS
          write(0,*)'tp1=',tp1(i,j,:)
          write(0,*)'qp1=',qp1(i,j,:)
          write(0,*)'psp1=',psp1(i,j)
          write(0,*)'tp2=',tp2(i,j,:)
          write(0,*)'qp2=',qp2(i,j,:)
          write(0,*)'psp2=',psp2(i,j)

          write(0,*)'p,cwm,T,Fice,Q'
           do k=1,km
            write(0,100)prsl(1,k),cwm_col(1,k),t_col(1,k),f_ice_phy(i,j,km+1-K),q_col(1,k),fpvs(t_col(1,k))
         !!   write(0,100)prsl(1,k),cwm_col(1,k),t_col(1,k),f_ice_phy(i,j,km+1-K),q_col(1,k)
           enddo
!!           write(0,*)fpvs(t_col(1,k))
          write(0,*)'Max,min T',maxval(tp1),minval(tp1)
       ENDIF
100    format(F10.1,E12.3,F6.1,F4.1,2E12.3)
    enddo                          ! End "I" loop
    enddo                          ! End "J" loop
!.......................................................................
     DO j = jts,jte
        DO k = kts,kte
	DO i = its,ite
           th_phy(i,j,k) = t_phy(i,j,k)/pi_phy(i,j,k)
          ENDDO   !- i
        ENDDO     !- k
     ENDDO        !- j
!.......................................................................
! 
!- Update rain (convert from m to kg/m**2, which is also equivalent to mm depth)
! 
       DO j=jts,jte
       DO i=its,ite
          RAINNC(i,j)=APREC(i,j)*1000.+RAINNC(i,j)
          RAINNCV(i,j)=APREC(i,j)*1000.
       ENDDO
       ENDDO
!
!-----------------------------------------------------------------------
!
  END SUBROUTINE GFSMP
!
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE GFSMP_INIT
          CALL GPVS
        END SUBROUTINE GFSMP_INIT    

!
      END MODULE module_mp_gfs
