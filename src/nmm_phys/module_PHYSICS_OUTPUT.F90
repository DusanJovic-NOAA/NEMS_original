!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
!***  LIST THE QUANTITIES FROM THE PHYSICS INTERNAL STATE THAT
!***  CAN BE SELECTED FOR HISTORY OUTPUT AND ASSOCIATE THEM WITH
!***  UNIQUE INTEGERS.
!***  THE USER WILL PLACE AN 'H' IN THE 3rd ELEMENT OF THE FOLLOWING
!***  LIST OF QUANTITIES THAT ARE IN THE PHYSICS INTERNAL STATE
!***  IF THAT QUANTITY IS TO BE WRITTEN TO THE HISTORY FILES.
!
!***  THE ORDER OF NAMES IN THE CHARACTER SELECTION LIST MUST BE
!***  THE SAME AS THE ORDER OF THE ACTUAL DATA BEING TARGETED
!***  BY POINTERS.
!-----------------------------------------------------------------------
!***  WHEN NEW QUANTITIES ARE ADDED TO THE PHYSICS INTERNAL STATE
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT THEN THEY NEED TO BE
!***  ADDED IN TWO PLACES BELOW: 
!*    (1) THE APPROPRIATE DATA LIST PRECEDING THE 'CONTAINS' STATEMENT
!*    (2) 'THE PHYSICS INTERNAL STATE POINTER BLOCK'
!*        IN SUBROUTINE POINT_PHYSICS_OUPUT
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY: INTERNAL_STATE 
      USE MODULE_ERR_MSG               ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: POINT_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_KOUNT=200
!
!-----------------------------------------------------------------------
!***  LIST THE VARIOUS QUANTITIES IN THE PHYSICS INTERNAL STATE 
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT.
!***  GROUP THEM BY TYPE (VARIOUS SCALARS, VARIOUS ARRAYS) SINCE
!***  WE WILL USE A WORKING POINTER FOR EACH TYPE.
!-----------------------------------------------------------------------
!
!-------------------------
!-------------------------
!***  INTEGER SCALARS  ***
!-------------------------
!-------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_ISCALAR     &
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               'IM        ', '-         ', '-         ' &
                              ,'JM        ', '-         ', '-         ' &
                              ,'LM        ', '-         ', '-         ' &
                              ,'NSOIL     ', '-         ', 'R         ' &
                              ,'NPHS      ', 'H         ', 'R         ' &
                              ,'NCLOD     ', 'H         ', 'R         ' &
                              ,'NHEAT     ', 'H         ', 'R         ' &
                              ,'NPREC     ', 'H         ', 'R         ' &
                              ,'NRDLW     ', 'H         ', 'R         ' &
                              ,'NRDSW     ', 'H         ', 'R         ' &
                              ,'NSRFC     ', 'H         ', 'R         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!----------------------
!----------------------
!***  REAL SCALARS  ***
!----------------------
!----------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_RSCALAR     &  
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               'PDTOP     ', '-         ', '-         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!----------------------------
!----------------------------
!***  INTEGER 1-D ARRAYS  ***
!----------------------------
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_1D_I        & 
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               '-         ', '-         ', '-         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!----------------------------
!----------------------------
!***  INTEGER 2-D ARRAYS  ***
!----------------------------
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_2D_I        & 
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               'ISLTYP    ', 'H         ', 'R         ' &
                              ,'IVGTYP    ', 'H         ', 'R         ' &
                              ,'NCFRCV    ', 'H         ', 'R         ' &
                              ,'NCFRST    ', 'H         ', 'R         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!-------------------------
!-------------------------
!***  REAL 1-D ARRAYS  *** 
!-------------------------
!-------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_1D_R        & 
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               'DSG2      ', '-         ', '-         ' &
                              ,'SGML2     ', '-         ', '-         ' &
                              ,'SLDPTH    ', 'H         ', 'R         ' &
                              ,'MP_RESTART', '-         ', 'R         ' &
                              ,'TBPVS_STAT', '-         ', 'R         ' &
                              ,'TBPVS0_STA', '-         ', 'R         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!-------------------------
!-------------------------
!***  REAL 2-D ARRAYS  ***
!-------------------------
!-------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_2D_R        & 
!
       =RESHAPE((/                                                      &
!                              -----------------------------------------
!
                               'ACFRCV    ', 'H         ', 'R         ' &
                              ,'ACFRST    ', 'H         ', 'R         ' &
                              ,'ACPREC    ', 'H         ', 'R         ' &
                              ,'ACSNOM    ', 'H         ', 'R         ' &
                              ,'ACSNOW    ', 'H         ', 'R         ' &
                              ,'AKHS_OUT  ', 'H         ', 'R         ' &
                              ,'AKMS_OUT  ', 'H         ', 'R         ' &
                              ,'ALBASE    ', 'H         ', 'R         ' &
                              ,'ALBEDO    ', 'H         ', 'R         ' &
                              ,'ALWIN     ', 'H         ', 'R         ' &
                              ,'ALWOUT    ', 'H         ', 'R         ' &
                              ,'ALWTOA    ', 'H         ', 'R         ' &
                              ,'ASWIN     ', 'H         ', 'R         ' &
                              ,'ASWOUT    ', 'H         ', 'R         ' &
                              ,'ASWTOA    ', 'H         ', 'R         ' &
                              ,'BGROFF    ', 'H         ', 'R         ' &
                              ,'CFRACH    ', 'H         ', 'R         ' &
                              ,'CFRACL    ', 'H         ', 'R         ' &
                              ,'CFRACM    ', 'H         ', 'R         ' &
                              ,'CLDEFI    ', 'H         ', 'R         ' &
                              ,'CMC       ', 'H         ', 'R         ' &
                              ,'CNVBOT    ', 'H         ', 'R         ' &
                              ,'CNVTOP    ', 'H         ', 'R         ' &
                              ,'CPRATE    ', 'H         ', 'R         ' &
                              ,'CUPPT     ', 'H         ', 'R         ' &
                              ,'CUPREC    ', 'H         ', 'R         ' &
                              ,'CZEN      ', 'H         ', 'R         ' &
                              ,'CZMEAN    ', 'H         ', 'R         ' &
                              ,'EPSR      ', 'H         ', 'R         ' &
                              ,'FIS       ', '-         ', '-         ' &  !<-- Already turned on in Dynamics output
                              ,'GRNFLX    ', 'H         ', 'R         ' &
                              ,'HBOTD     ', 'H         ', 'R         ' &
                              ,'HBOTS     ', 'H         ', 'R         ' &
                              ,'HTOPD     ', 'H         ', 'R         ' &
                              ,'HTOPS     ', 'H         ', 'R         ' &
                              ,'MXSNAL    ', 'H         ', 'R         ' &
                              ,'PBLH      ', 'H         ', 'R         ' &
                              ,'PD        ', '-         ', '-         ' &  !<-- Already turned on in Dynamics output
                              ,'POTEVP    ', 'H         ', 'R         ' &
                              ,'PREC      ', 'H         ', 'R         ' &
                              ,'PSHLTR    ', 'H         ', 'R         ' &
                              ,'Q10       ', 'H         ', 'R         ' &
                              ,'QSH       ', 'H         ', 'R         ' &
                              ,'QSHLTR    ', 'H         ', 'R         ' &
                              ,'QWBS      ', 'H         ', 'R         ' &
                              ,'QZ0       ', 'H         ', 'R         ' &
                              ,'RADOT     ', 'H         ', 'R         ' &
                              ,'RLWIN     ', 'H         ', 'R         ' &
                              ,'RLWTOA    ', 'H         ', 'R         ' &
                              ,'RSWIN     ', 'H         ', 'R         ' &
                              ,'RSWINC    ', 'H         ', 'R         ' &
                              ,'RSWOUT    ', 'H         ', 'R         ' &
                              ,'SFCEVP    ', 'H         ', 'R         ' &
                              ,'SFCEXC    ', 'H         ', 'R         ' &
                              ,'SFCLHX    ', 'H         ', 'R         ' &
                              ,'SFCSHX    ', 'H         ', 'R         ' &
                              ,'SI        ', 'H         ', 'R         ' &
                              ,'SICE      ', 'H         ', 'R         ' &
                              ,'SIGT4     ', 'H         ', 'R         ' &
                              ,'SM        ', 'H         ', 'R         ' &
                              ,'SMSTAV    ', 'H         ', 'R         ' &
                              ,'SMSTOT    ', 'H         ', 'R         ' &
                              ,'SNO       ', 'H         ', 'R         ' &
                              ,'SNOPCX    ', 'H         ', 'R         ' &
                              ,'SOILTB    ', 'H         ', 'R         ' &
                              ,'SR        ', 'H         ', 'R         ' &
                              ,'SSROFF    ', 'H         ', 'R         ' &
                              ,'SST       ', 'H         ', 'R         ' &
                              ,'SUBSHX    ', 'H         ', 'R         ' &
                              ,'TG        ', 'H         ', 'R         ' &
                              ,'TH10      ', 'H         ', 'R         ' &
                              ,'THS       ', 'H         ', 'R         ' &
                              ,'THZ0      ', 'H         ', 'R         ' &
                              ,'TSHLTR    ', 'H         ', 'R         ' &
                              ,'TWBS      ', 'H         ', 'R         ' &
                              ,'U10       ', 'H         ', 'R         ' &
                              ,'USTAR     ', 'H         ', 'R         ' &
                              ,'UZ0       ', 'H         ', 'R         ' &
                              ,'V10       ', 'H         ', 'R         ' &
                              ,'VEGFRC    ', 'H         ', 'R         ' &
                              ,'VZ0       ', 'H         ', 'R         ' &
                              ,'Z0        ', 'H         ', 'R         ' &
                              ,'TSKIN     ', '-         ', 'R         ' &
                              ,'AKHS      ', '-         ', 'R         ' &
                              ,'AKMS      ', '-         ', 'R         ' &
                              ,'HBOT      ', '-         ', 'R         ' &
                              ,'HTOP      ', '-         ', 'R         ' &
                              ,'RSWTOA    ', '-         ', 'R         ' &
                              ,'POTFLX    ', 'H         ', 'R         ' &
                              ,'RMOL      ', '-         ', 'R         ' &
                              ,'T2        ', '-         ', 'R         ' &
                              ,'Z0BASE    ', '-         ', 'R         ' &
                              ,'PSFC      ', 'H         ', '-         ' &
                              ,'TLMIN     ', 'H         ', 'R         ' &
                              ,'TLMAX     ', 'H         ', 'R         ' &
                              ,'LSPA      ', 'H         ', '-         ' &
                              ,'ACUTIM    ', 'H         ', 'R         ' &
                              ,'APHTIM    ', 'H         ', 'R         ' &
                              ,'ARDLW     ', 'H         ', 'R         ' &
                              ,'ARDSW     ', 'H         ', 'R         ' &
                              ,'ASRFC     ', 'H         ', 'R         ' &
                              ,'AVRAIN    ', 'H         ', 'R         ' &
                              ,'AVCNVC    ', 'H         ', 'R         ' &
!
!                              -----------------------------------------
!
         /)                                                             &
        ,(/3,MAX_KOUNT/)                                                &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!-------------------------
!-------------------------
!***  REAL 3-D ARRAYS  ***
!-------------------------
!-------------------------
!
      CHARACTER(12),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_3D_R          & 
!
       =RESHAPE((/                                                        &
!                          ---------------------------------------------
!                                                                            !<-- Append "_ikj" to any 3D variables with IKJ storage order
                           'CLDFRA      ', 'H           ', 'R           ' &
                          ,'CW          ', 'H           ', 'R           ' &
                          ,'EXCH_H      ', 'H           ', '-           ' &
                          ,'Q           ', 'H           ', 'R           ' &
                          ,'Q2          ', 'H           ', 'R           ' &
                          ,'RLWTT       ', 'H           ', 'R           ' &
                          ,'RSWTT       ', 'H           ', 'R           ' &
                          ,'T           ', 'H           ', 'R           ' &
                          ,'TCUCN       ', 'H           ', 'R           ' &
                          ,'TRAIN       ', 'H           ', 'R           ' &
                          ,'U           ', 'H           ', 'R           ' &
                          ,'V           ', 'H           ', 'R           ' &
                          ,'XLEN_MIX    ', 'H           ', 'R           ' &
!
                          ,'F_ICE       ', 'H           ', 'R           ' &
                          ,'F_RIMEF     ', 'H           ', 'R           ' &
                          ,'F_RAIN      ', 'H           ', 'R           ' &
                          ,'RTHBLTEN    ', '-           ', 'R           ' &
                          ,'RQVBLTEN    ', '-           ', 'R           ' &
                          ,'SH2O        ', 'H           ', 'R           ' &
                          ,'SMC         ', 'H           ', 'R           ' &
                          ,'STC         ', 'H           ', 'R           ' &
!
!                          ---------------------------------------------
!
         /)                                                               &
        ,(/3,MAX_KOUNT/)                                                  &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!
!-------------------------
!-------------------------
!***  REAL 4-D ARRAYS  ***
!-------------------------
!-------------------------
!
      CHARACTER(12),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_4D_R          &
!
       =RESHAPE((/                                                        &
!                          ---------------------------------------------
                           'TRACERS     ', '-           ', 'R           ' &
!
!                          ---------------------------------------------
!
         /)                                                               &
        ,(/3,MAX_KOUNT/)                                                  &
        ,(/'**********', '**********', '**********'/))
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_PHYSICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE TAKES THE USER'S SELECTIONS FOR OUTPUT QUANTITIES,
!***  POINTS AT THEM, AND INSERTS THOSE POINTERS INTO THE IMPORT STATE
!***  OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Grid)     ,INTENT(IN)    :: GRID                         !<-- The ESMF Grid
      TYPE(ESMF_State)    ,INTENT(INOUT) :: IMP_STATE_WRITE              !<-- Import state for the write gridded components
!
      TYPE(INTERNAL_STATE),POINTER       :: INT_STATE                    !<-- The physics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER                      :: ITS,ITE,JTS,JTE                   &
                                     ,IMS,IME,JMS,JME                   &
                                     ,IHALO,JHALO
!
      INTEGER                      :: INDEX_IKJ,K,LENGTH                &
                                     ,MP_PHYSICS,MYPE                   &
                                     ,N,NDIM3,NFIND                     &
                                     ,NUM_2D_FIELDS_I,NUM_2D_FIELDS_R   &
                                     ,RC,RC_PHY_OUT                     &
                                     ,SF_SURFACE_PHYSICS
!
      INTEGER                     :: LDIM1,LDIM2,LDIM3,LDIM4            &
                                    ,UDIM1,UDIM2,UDIM3,UDIM4
!
      INTEGER(ESMF_KIND_I4),DIMENSION(:,:),POINTER :: TEMP_I2D
      REAL(ESMF_KIND_R4)   ,DIMENSION(:,:),POINTER :: TEMP_R2D
!
      CHARACTER(2)           :: MODEL_LEVEL,TRACERS_KIND
      CHARACTER(6)           :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR) :: VBL_NAME,VBL_NAME_X
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE
      TYPE(ESMF_FieldBundle),SAVE :: RESTART_BUNDLE
!
      TYPE(ESMF_Field)       :: FIELD
!
      TYPE(ESMF_CopyFlag)    :: COPYFLAG=ESMF_DATA_REF
!     TYPE(ESMF_CopyFlag)    :: COPYFLAG=ESMF_DATA_COPY
!
!-----------------------------------------------------------------------
!***  FIRST WE MUST PROVIDE POINTERS INTO THE PHYSICS INTERNAL STATE. 
!-----------------------------------------------------------------------
      TYPE PHY_ISC
        INTEGER(ESMF_KIND_I4),POINTER :: NAME                             !<-- Pointer for integer scalars
      END TYPE PHY_ISC
!-----------------------------------------------------------------------
      TYPE PHY_RSC
        REAL(ESMF_KIND_R4),POINTER :: NAME                                !<-- Pointer for real scalars
      END TYPE PHY_RSC
!-----------------------------------------------------------------------
      TYPE PHY_I1D
        INTEGER(ESMF_KIND_I4),DIMENSION(:),POINTER :: NAME                !<-- Pointer for 1D integer arrays
      END TYPE PHY_I1D
!-----------------------------------------------------------------------
      TYPE PHY_I2D
        INTEGER(ESMF_KIND_I4),DIMENSION(:,:),POINTER :: NAME              !<-- Pointer for 2D integer arrays
      END TYPE PHY_I2D
!-----------------------------------------------------------------------
      TYPE PHY_R1D
        REAL(ESMF_KIND_R4),DIMENSION(:),POINTER :: NAME                   !<-- Pointer for 1D real arrays
      END TYPE PHY_R1D
!-----------------------------------------------------------------------
      TYPE PHY_R2D
        REAL(ESMF_KIND_R4),DIMENSION(:,:),POINTER :: NAME                 !<-- Pointer for 2D real arrays
      END TYPE PHY_R2D
!-----------------------------------------------------------------------
      TYPE PHY_R3D
        REAL(ESMF_KIND_R4),DIMENSION(:,:,:),POINTER :: NAME               !<-- Pointer for 3D real arrays
      END TYPE PHY_R3D
!-----------------------------------------------------------------------
      TYPE PHY_R4D
        REAL(ESMF_KIND_R4),DIMENSION(:,:,:,:),POINTER :: NAME             !<-- Pointer for 4D real arrays
      END TYPE PHY_R4D
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARRAYS OF POINTERS OF THE ABOVE TYPES
!-----------------------------------------------------------------------
!
      TYPE(PHY_ISC),DIMENSION(MAX_KOUNT) :: I_SC
      TYPE(PHY_RSC),DIMENSION(MAX_KOUNT) :: R_SC
      TYPE(PHY_I1D),DIMENSION(MAX_KOUNT) :: I_1D
      TYPE(PHY_I2D),DIMENSION(MAX_KOUNT) :: I_2D
      TYPE(PHY_R1D),DIMENSION(MAX_KOUNT) :: R_1D
      TYPE(PHY_R2D),DIMENSION(MAX_KOUNT) :: R_2D
      TYPE(PHY_R3D),DIMENSION(MAX_KOUNT) :: R_3D
      TYPE(PHY_R4D),DIMENSION(MAX_KOUNT) :: R_4D
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE PHYSICS INTERNAL STATE POINTER BLOCK
!-----------------------------------------------------------------------
!***  POINT AT ALL VARIABLES IN THE PHYSICS INTERNAL STATE 
!***  THAT COULD BE WRITTEN TO HISTORY OUTPUT, I.E., THOSE
!***  LISTED AT THE TOP OF THIS MODULE.
!-----------------------------------------------------------------------
!
!---------------------
!***  INTEGER SCALARS
!---------------------
!
      I_SC( 1)%NAME=>int_state%IM
      I_SC( 2)%NAME=>int_state%JM
      I_SC( 3)%NAME=>int_state%LM
      I_SC( 4)%NAME=>int_state%NSOIL
      I_SC( 5)%NAME=>int_state%NPHS
      I_SC( 6)%NAME=>int_state%NCLOD
      I_SC( 7)%NAME=>int_state%NHEAT
      I_SC( 8)%NAME=>int_state%NPREC
      I_SC( 9)%NAME=>int_state%NRDLW
      I_SC(10)%NAME=>int_state%NRDSW
      I_SC(11)%NAME=>int_state%NSRFC
!        
!------------------
!***  REAL SCALARS
!------------------
!
      R_SC(1)%NAME=>int_state%PDTOP
!        
!-----------------------
!***  1D INTEGER ARRAYS
!-----------------------
!
!!!   I_1D(1)%NAME=>
!        
!--------------------
!***  1D REAL ARRAYS
!--------------------
!
      R_1D(1)%NAME=>int_state%DSG2
      R_1D(2)%NAME=>int_state%SGML2
      R_1D(3)%NAME=>int_state%SLDPTH
      R_1D(4)%NAME=>int_state%MP_RESTART_STATE
      R_1D(5)%NAME=>int_state%TBPVS_STATE
      R_1D(6)%NAME=>int_state%TBPVS0_STATE
!        
!-----------------------
!***  2D INTEGER ARRAYS
!-----------------------
!
      I_2D( 1)%NAME=>int_state%ISLTYP
      I_2D( 2)%NAME=>int_state%IVGTYP
      I_2D( 3)%NAME=>int_state%NCFRCV
      I_2D( 4)%NAME=>int_state%NCFRST
!
!--------------------
!***  2D REAL ARRAYS
!--------------------
!
      R_2D( 1)%NAME=>int_state%ACFRCV
      R_2D( 2)%NAME=>int_state%ACFRST
      R_2D( 3)%NAME=>int_state%ACPREC
      R_2D( 4)%NAME=>int_state%ACSNOM
      R_2D( 5)%NAME=>int_state%ACSNOW
      R_2D( 6)%NAME=>int_state%AKHS_OUT
      R_2D( 7)%NAME=>int_state%AKMS_OUT
      R_2D( 8)%NAME=>int_state%ALBASE
      R_2D( 9)%NAME=>int_state%ALBEDO
      R_2D(10)%NAME=>int_state%ALWIN
      R_2D(11)%NAME=>int_state%ALWOUT
      R_2D(12)%NAME=>int_state%ALWTOA
      R_2D(13)%NAME=>int_state%ASWIN
      R_2D(14)%NAME=>int_state%ASWOUT
      R_2D(15)%NAME=>int_state%ASWTOA
      R_2D(16)%NAME=>int_state%BGROFF
      R_2D(17)%NAME=>int_state%CFRACH
      R_2D(18)%NAME=>int_state%CFRACL
      R_2D(19)%NAME=>int_state%CFRACM
      R_2D(20)%NAME=>int_state%CLDEFI
      R_2D(21)%NAME=>int_state%CMC
      R_2D(22)%NAME=>int_state%CNVBOT
      R_2D(23)%NAME=>int_state%CNVTOP
      R_2D(24)%NAME=>int_state%CPRATE
      R_2D(25)%NAME=>int_state%CUPPT
      R_2D(26)%NAME=>int_state%CUPREC
      R_2D(27)%NAME=>int_state%CZEN
      R_2D(28)%NAME=>int_state%CZMEAN
      R_2D(29)%NAME=>int_state%EPSR
      R_2D(30)%NAME=>int_state%FIS
      R_2D(31)%NAME=>int_state%GRNFLX
      R_2D(32)%NAME=>int_state%HBOTD
      R_2D(33)%NAME=>int_state%HBOTS
      R_2D(34)%NAME=>int_state%HTOPD
      R_2D(35)%NAME=>int_state%HTOPS
      R_2D(36)%NAME=>int_state%MXSNAL
      R_2D(37)%NAME=>int_state%PBLH
      R_2D(38)%NAME=>int_state%PD
      R_2D(39)%NAME=>int_state%POTEVP
      R_2D(40)%NAME=>int_state%PREC
      R_2D(41)%NAME=>int_state%PSHLTR
      R_2D(42)%NAME=>int_state%Q10
      R_2D(43)%NAME=>int_state%QSH
      R_2D(44)%NAME=>int_state%QSHLTR
      R_2D(45)%NAME=>int_state%QWBS
      R_2D(46)%NAME=>int_state%QZ0
      R_2D(47)%NAME=>int_state%RADOT
      R_2D(48)%NAME=>int_state%RLWIN
      R_2D(49)%NAME=>int_state%RLWTOA
      R_2D(50)%NAME=>int_state%RSWIN
      R_2D(51)%NAME=>int_state%RSWINC
      R_2D(52)%NAME=>int_state%RSWOUT
      R_2D(53)%NAME=>int_state%SFCEVP
      R_2D(54)%NAME=>int_state%SFCEXC
      R_2D(55)%NAME=>int_state%SFCLHX
      R_2D(56)%NAME=>int_state%SFCSHX
      R_2D(57)%NAME=>int_state%SI
      R_2D(58)%NAME=>int_state%SICE
      R_2D(59)%NAME=>int_state%SIGT4
      R_2D(60)%NAME=>int_state%SM
      R_2D(61)%NAME=>int_state%SMSTAV
      R_2D(62)%NAME=>int_state%SMSTOT
      R_2D(63)%NAME=>int_state%SNO
      R_2D(64)%NAME=>int_state%SNOPCX
      R_2D(65)%NAME=>int_state%SOILTB
      R_2D(66)%NAME=>int_state%SR
      R_2D(67)%NAME=>int_state%SSROFF
      R_2D(68)%NAME=>int_state%SST
      R_2D(69)%NAME=>int_state%SUBSHX
      R_2D(70)%NAME=>int_state%TG
      R_2D(71)%NAME=>int_state%TH10
      R_2D(72)%NAME=>int_state%THS
      R_2D(73)%NAME=>int_state%THZ0
      R_2D(74)%NAME=>int_state%TSHLTR
      R_2D(75)%NAME=>int_state%TWBS
      R_2D(76)%NAME=>int_state%U10
      R_2D(77)%NAME=>int_state%USTAR
      R_2D(78)%NAME=>int_state%UZ0
      R_2D(79)%NAME=>int_state%V10
      R_2D(80)%NAME=>int_state%VEGFRC
      R_2D(81)%NAME=>int_state%VZ0
      R_2D(82)%NAME=>int_state%Z0
      R_2D(83)%NAME=>int_state%TSKIN
      R_2D(84)%NAME=>int_state%AKHS
      R_2D(85)%NAME=>int_state%AKMS
      R_2D(86)%NAME=>int_state%HBOT
      R_2D(87)%NAME=>int_state%HTOP
      R_2D(88)%NAME=>int_state%RSWTOA
      R_2D(89)%NAME=>int_state%POTFLX
      R_2D(90)%NAME=>int_state%RMOL
      R_2D(91)%NAME=>int_state%T2
      R_2D(92)%NAME=>int_state%Z0BASE
      R_2D(93)%NAME=>int_state%PSFC
      R_2D(94)%NAME=>int_state%TLMIN
      R_2D(95)%NAME=>int_state%TLMAX
      R_2D(96)%NAME=>int_state%LSPA
      R_2D(97)%NAME=>int_state%ACUTIM
      R_2D(98)%NAME=>int_state%APHTIM
      R_2D(99)%NAME=>int_state%ARDLW
      R_2D(100)%NAME=>int_state%ARDSW
      R_2D(101)%NAME=>int_state%ASRFC
      R_2D(102)%NAME=>int_state%AVCNVC
      R_2D(103)%NAME=>int_state%AVRAIN
!        
!--------------------
!***  3D REAL ARRAYS
!--------------------
!
      R_3D( 1)%NAME=>int_state%CLDFRA
      R_3D( 2)%NAME=>int_state%CW
      R_3D( 3)%NAME=>int_state%EXCH_H
      R_3D( 4)%NAME=>int_state%Q
      R_3D( 5)%NAME=>int_state%Q2
      R_3D( 6)%NAME=>int_state%RLWTT
      R_3D( 7)%NAME=>int_state%RSWTT
      R_3D( 8)%NAME=>int_state%T
      R_3D( 9)%NAME=>int_state%TCUCN
      R_3D(10)%NAME=>int_state%TRAIN
      R_3D(11)%NAME=>int_state%U
      R_3D(12)%NAME=>int_state%V
      R_3D(13)%NAME=>int_state%XLEN_MIX
!
      R_3D(14)%NAME=>int_state%F_ICE
      R_3D(15)%NAME=>int_state%F_RIMEF
      R_3D(16)%NAME=>int_state%F_RAIN
      R_3D(17)%NAME=>int_state%RTHBLTEN
      R_3D(18)%NAME=>int_state%RQVBLTEN
      R_3D(19)%NAME=>int_state%SH2O
      R_3D(20)%NAME=>int_state%SMC
      R_3D(21)%NAME=>int_state%STC
!
!--------------------
!***  4D REAL ARRAYS
!--------------------
!
      R_4D( 1)%NAME=>int_state%TRACERS
!
!-----------------------------------------------------------------------
!
      ITS=int_state%ITS
      ITE=int_state%ITE
      JTS=int_state%JTS
      JTE=int_state%JTE
!
      IMS=int_state%IMS
      IME=int_state%IME
      JMS=int_state%JMS
      JME=int_state%JME
!
      MYPE=int_state%MYPE
!
!-----------------------------------------------------------------------
!***  EXTRACT THE HISTORY OUTPUT Bundle FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  IT ALREADY CONTAINS OUTPUT VARIABLES FROM
!***  THE DYNAMICS.  WE ARE PREPARING TO ADD HISTORY VARIABLES
!***  FROM THE PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract History Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =IMP_STATE_WRITE                    &  !<-- Take Bundle from the Write component's import state
                        ,itemName   ='Bundle_Output_Data'               &  !<-- The Bundle's name
                        ,fieldbundle=HISTORY_BUNDLE                     &  !<-- The Bundle object
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  EXTRACT THE RESTART OUTPUT Bundle FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  IT ALREADY CONTAINS RESTART VARIABLES FROM
!***  THE DYNAMICS.  WE ARE PREPARING TO ADD RESTART VARIABLES
!***  FROM THE PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =IMP_STATE_WRITE                    &  !<-- Take Bundle from the Write component's import state
                        ,itemName   ='Bundle_Restart_Data'              &  !<-- The Bundle's name
                        ,fieldbundle=RESTART_BUNDLE                     &  !<-- The Bundle object
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE MICROPHYSICS SCHEME SPECIFICATION IS NEEDED IN THE OUTPUT
!***  SO ADD IT DIRECTLY.
!-----------------------------------------------------------------------
!
      IF(int_state%MICROPHYSICS=='fer')THEN
        MP_PHYSICS=5
      ELSEIF(int_state%MICROPHYSICS=='kes')THEN
        MP_PHYSICS=1
      ELSEIF(int_state%MICROPHYSICS=='lin')THEN
        MP_PHYSICS=2
      ELSEIF(int_state%MICROPHYSICS=='tho')THEN
        MP_PHYSICS=8
      ELSEIF(int_state%MICROPHYSICS=='wsm3')THEN
        MP_PHYSICS=3
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Microphysics Scheme Specification into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component history Bundle
                            ,name  ='MP_PHYSICS'                        &  !<-- Name of microphysics scheme variable
                            ,value =MP_PHYSICS                          &  !<-- The microphysics scheme integer specification
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Microphysics Scheme Specification into Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                      &  !<-- The Write component restart Bundle
                            ,name  ='MP_PHYSICS'                        &  !<-- Name of microphysics scheme variable
                            ,value =MP_PHYSICS                          &  !<-- The microphysics scheme integer specification
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE LAND SURFACE SCHEME SPECIFICATION IS NEEDED IN THE OUTPUT
!***  SO ADD IT DIRECTLY.
!-----------------------------------------------------------------------
!
      IF(int_state%LAND_SURFACE=='noah')THEN
        SF_SURFACE_PHYSICS=2
      ELSEIF(int_state%LAND_SURFACE=='slab')THEN
        SF_SURFACE_PHYSICS=1
      ELSEIF(int_state%LAND_SURFACE=='ruc')THEN
        SF_SURFACE_PHYSICS=3
      ELSEIF(int_state%LAND_SURFACE=='nmm')THEN
        SF_SURFACE_PHYSICS=99
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Surface Physics Scheme Specification into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component history Bundle
                            ,name  ='SF_SURFACE_PHYSICS'                &  !<-- Name of land surface scheme variable
                            ,value =SF_SURFACE_PHYSICS                  &  !<-- The land surface scheme integer specification
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Surface Physics Scheme Specification into Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                      &  !<-- The Write component restart Bundle
                            ,name  ='SF_SURFACE_PHYSICS'                &  !<-- Name of land surface scheme variable
                            ,value =SF_SURFACE_PHYSICS                  &  !<-- The land surface scheme integer specification
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER
!***  FOR HISTORY OUTPUT.
!***  PUT THE HISTORY DATA INTO THE ESMF Bundle RESIDING IN THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  BEGIN WITH THE INTEGER SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Physics Integer Scalars into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_ISCALAR(2,NFIND)=='H')THEN                        !<-- Take integer scalar data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_ISCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                  &  !<-- The write component's output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =I_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the output data Bundle
                                ,rc    =RC)
!
        ENDIF
!
        IF(PHY_INT_STATE_ISCALAR(3,NFIND)=='R')THEN                        !<-- Take integer scalar data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_ISCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                  &  !<-- The write component's output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =I_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the output data Bundle
                                ,rc    =RC)
        ENDIF
!
        IF(PHY_INT_STATE_ISCALAR(2,NFIND)=='*' .AND.                    &
           PHY_INT_STATE_ISCALAR(3,NFIND)=='*')THEN                        !<-- End of the integer scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE REAL SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Physics Real Scalars into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_RSCALAR(2,NFIND)=='H')THEN                        !<-- Take real scalar data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_RSCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                  &  !<-- The write component's output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =R_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the output data Bundle
                                ,rc    =RC)
!
        ENDIF
!
        IF(PHY_INT_STATE_RSCALAR(3,NFIND)=='R')THEN                        !<-- Take real scalar data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_RSCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                  &  !<-- The write component's output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =R_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the output data Bundle
                                ,rc    =RC)
        ENDIF
!
        IF(PHY_INT_STATE_RSCALAR(2,NFIND)=='*' .AND.                    &
           PHY_INT_STATE_RSCALAR(3,NFIND)=='*'      )THEN                  !<-- End of the real scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D INTEGER ARRAYS
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Physics 1-D Integer Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_1D_I(2,NFIND)=='H')THEN                           !<-- Take 1D integer array data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_1D_I(1,NFIND))
          LENGTH=SIZE(I_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =HISTORY_BUNDLE               &  !<-- The write component's output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=I_1D(NFIND)%NAME             &  !<-- The 1D integer being inserted into the output data Bundle
                                ,rc       =RC)
!
        ENDIF
!
        IF(PHY_INT_STATE_1D_I(3,NFIND)=='R')THEN                           !<-- Take 1D integer array data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_1D_I(1,NFIND))
          LENGTH=SIZE(I_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =RESTART_BUNDLE               &  !<-- The write component's output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=I_1D(NFIND)%NAME             &  !<-- The 1D integer being inserted into the output data Bundle
                                ,rc       =RC)
        ENDIF
!
        IF(PHY_INT_STATE_1D_I(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_1D_I(3,NFIND)=='*'      )THEN                     !<-- End of the 1D integer array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D REAL ARRAYS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Physics 1-D Real Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_1D_R(2,NFIND)=='H')THEN                           !<-- Take 1D real array data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_1D_R(1,NFIND))
          LENGTH=SIZE(R_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =HISTORY_BUNDLE               &  !<-- The write component's output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=R_1D(NFIND)%NAME             &  !<-- The 1D real being inserted into the output data Bundle
                                ,rc       =RC)
!
        ENDIF
!
        IF(PHY_INT_STATE_1D_R(3,NFIND)=='R')THEN                           !<-- Take 1D real array data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_1D_R(1,NFIND))
          LENGTH=SIZE(R_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =RESTART_BUNDLE               &  !<-- The write component's output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=R_1D(NFIND)%NAME             &  !<-- The 1D real being inserted into the output data Bundle
                                ,rc       =RC)
        ENDIF
!
        IF(PHY_INT_STATE_1D_R(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_1D_R(3,NFIND)=='*')THEN                           !<-- End of the 1D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 2D INTEGER ARRAYS.
!-----------------------------------------------------------------------
!
      IHALO=3
      JHALO=3
!
      NUM_2D_FIELDS_I=0
      NULLIFY(TEMP_I2D)
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_2D_I(2,NFIND)=='H')THEN                           !<-- Take 2D integer array data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_2D_I(1,NFIND))
          TEMP_I2D=>I_2D(NFIND)%NAME(IMS:IME,JMS:JME)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Integer Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_I2D                 &  !<-- The 2D integer array being inserted into output data Bundle
                                ,copyflag     =COPYFLAG                 &
                                ,maxHaloUWidth=(/IHALO,JHALO/)          &
                                ,maxHaloLWidth=(/IHALO,JHALO/)          &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Integer Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE                &  !<-- The write component's output data Bundle
                                  ,field =FIELD                         &  !<-- ESMF Field holding the 2D integer array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NUM_2D_FIELDS_I=NUM_2D_FIELDS_I+1                                !<-- Add upt the number of 2D integer Fields
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_I(3,NFIND)=='R')THEN                           !<-- Take 2D integer array data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_2D_I(1,NFIND))
          TEMP_I2D=>I_2D(NFIND)%NAME(IMS:IME,JMS:JME)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Integer Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_I2D                 &  !<-- The 2D integer array being inserted into output data Bundle
                                ,copyflag     =COPYFLAG                 &
                                ,maxHaloUWidth=(/IHALO,JHALO/)          &
                                ,maxHaloLWidth=(/IHALO,JHALO/)          &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Integer Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(bundle=RESTART_BUNDLE                &  !<-- The write component's restart data Bundle
                                  ,field =FIELD                         &  !<-- ESMF Field holding the 2D integer array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_I(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_2D_I(3,NFIND)=='*'      )THEN                     !<-- End of the 2D integer array list
          EXIT
        ENDIF
!
      ENDDO
      IF(MYPE==0)WRITE(0,*)' PHYSICS_OUTPUT: Number of 2-D Integer Fields=',NUM_2D_FIELDS_I
!
!-----------------------------------------------------------------------
!***  THE 2D REAL ARRAYS.
!-----------------------------------------------------------------------
!
      NUM_2D_FIELDS_R=0
      NULLIFY(TEMP_R2D)
!
      DO NFIND=1,MAX_KOUNT
        IF(PHY_INT_STATE_2D_R(2,NFIND)=='H')THEN                           !<-- Take 2D real array data specified for history output
          VBL_NAME=TRIM(PHY_INT_STATE_2D_R(1,NFIND))
          TEMP_R2D=>R_2D(NFIND)%NAME(IMS:IME,JMS:JME)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Real Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_R2D                 &  !<-- The 2D real array being inserted into the output data Bundle
                                ,copyflag     =COPYFLAG                 &
                                ,maxHaloUWidth=(/IHALO,JHALO/)          &
                                ,maxHaloLWidth=(/IHALO,JHALO/)          &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Real Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE                &  !<-- The write component's output data Bundle
                                  ,field =FIELD                         &  !<-- ESMF Field holding the 2D real array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                                !<-- Add up the number of 2D real Fields
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_R(3,NFIND)=='R')THEN                           !<-- Take 2D real array data specified for restart output
          VBL_NAME=TRIM(PHY_INT_STATE_2D_R(1,NFIND))
          TEMP_R2D=>R_2D(NFIND)%NAME(IMS:IME,JMS:JME)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Real Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_R2D                 &  !<-- The 2D real array being inserted into the output data Bundle
                                ,copyflag     =COPYFLAG                 &
                                ,maxHaloUWidth=(/IHALO,JHALO/)          &
                                ,maxHaloLWidth=(/IHALO,JHALO/)          &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Physics 2-D Real Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(bundle=RESTART_BUNDLE                &  !<-- The write component's restart data Bundle
                                  ,field =FIELD                         &  !<-- ESMF Field holding the 2D real array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_R(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_2D_R(3,NFIND)=='*'      )THEN                     !<-- End of the 2D real array list
          EXIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  THE 3D REAL ARRAYS.
!***  WE ARE WORKING WITH 3D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!
!
      DO NFIND=1,MAX_KOUNT                                                 
!
        IF(PHY_INT_STATE_3D_R(2,NFIND)=='H')THEN                           !<-- Take 3D real array data specified for history output
!
          INDEX_IKJ=INDEX(PHY_INT_STATE_3D_R(1,NFIND),'_ikj')              !<-- If an IKJ array, isolate "_ikj" in the NAME
!
          IF(INDEX_IKJ==0)THEN
            NDIM3=UBOUND(R_3D(NFIND)%NAME,3)                               !<-- Determine # of vertical levels in this IJK variable
            LDIM2=LBOUND(R_3D(NFIND)%NAME,2)
            UDIM2=UBOUND(R_3D(NFIND)%NAME,2)
          ELSE
            NDIM3=UBOUND(R_3D(NFIND)%NAME,2)                               !<-- Determine # of vertical levels in this IKJ variable
            LDIM2=LBOUND(R_3D(NFIND)%NAME,3)
            UDIM2=UBOUND(R_3D(NFIND)%NAME,3)
          ENDIF
!
          LDIM1=LBOUND(R_3D(NFIND)%NAME,1)
          UDIM1=UBOUND(R_3D(NFIND)%NAME,1)
!
          DO K=1,NDIM3                                                     !<-- Loop through the levels of the array
            WRITE(MODEL_LEVEL,FMT)K
!
            IF(INDEX_IKJ==0)THEN                                          
              VBL_NAME=TRIM(PHY_INT_STATE_3D_R(1,NFIND))//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_3D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K)        !<-- Point at appropriate section of this IJK level 
            ELSE
              VBL_NAME_X=PHY_INT_STATE_3D_R(1,NFIND)(1:INDEX_IKJ-1)
              VBL_NAME=TRIM(VBL_NAME_X)//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_3D(NFIND)%NAME(LDIM1:UDIM1,K,LDIM2:UDIM2)        !<-- Point at appropriate section of this IKJ level
            ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Physics 3-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 3D real array being inserted into the data Bundle
                                  ,copyflag     =COPYFLAG               &
                                  ,maxHaloUWidth=(/IHALO,JHALO/)        &
                                  ,maxHaloLWidth=(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 3D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL         &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert Physics 3-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE              &  !<-- The write component's output data Bundle
                                    ,field =FIELD                       &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                              !<-- Continue adding up all levels of 3D Fields
!
          ENDDO
!
        ENDIF
!
        IF(PHY_INT_STATE_3D_R(3,NFIND)=='R')THEN                           !<-- Take 3D real array data specified for restart output
!
          INDEX_IKJ=INDEX(PHY_INT_STATE_3D_R(1,NFIND),'_ikj')              !<-- If an IKJ array, isolate "_ikj" in the NAME
!
          IF(INDEX_IKJ==0)THEN
            NDIM3=UBOUND(R_3D(NFIND)%NAME,3)                               !<-- Determine # of vertical levels in this IJK variable
            LDIM2=LBOUND(R_3D(NFIND)%NAME,2)
            UDIM2=UBOUND(R_3D(NFIND)%NAME,2)
          ELSE
            NDIM3=UBOUND(R_3D(NFIND)%NAME,2)                               !<-- Determine # of vertical levels in this IKJ variable
            LDIM2=LBOUND(R_3D(NFIND)%NAME,3)
            UDIM2=UBOUND(R_3D(NFIND)%NAME,3)
          ENDIF
!
          LDIM1=LBOUND(R_3D(NFIND)%NAME,1)
          UDIM1=UBOUND(R_3D(NFIND)%NAME,1)
!
          DO K=1,NDIM3                                                     !<-- Loop through the levels of the array
            WRITE(MODEL_LEVEL,FMT)K
!
            IF(INDEX_IKJ==0)THEN
              VBL_NAME=TRIM(PHY_INT_STATE_3D_R(1,NFIND))//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_3D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K)        !<-- Point at appropriate section of this IJK level
            ELSE
              VBL_NAME_X=PHY_INT_STATE_3D_R(1,NFIND)(1:INDEX_IKJ-1)
              VBL_NAME=TRIM(VBL_NAME_X)//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_3D(NFIND)%NAME(LDIM1:UDIM1,K,LDIM2:UDIM2)        !<-- Point at appropriate section of this IKJ level
            ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Physics 3-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 3D real array being inserted into the data Bundle
                                  ,copyflag     =COPYFLAG               &
                                  ,maxHaloUWidth=(/IHALO,JHALO/)        &
                                  ,maxHaloLWidth=(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 3D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL         &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert Physics 3-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(bundle=RESTART_BUNDLE              &  !<-- The write component's restart data Bundle
                                    ,field =FIELD                       &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                              !<-- Continue adding up all levels of 3D Fields
!
          ENDDO
!
        ENDIF
!
        IF(PHY_INT_STATE_3D_R(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_3D_R(3,NFIND)=='*'      )THEN                     !<-- End of the 3D real array list
          EXIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  THE 4D REAL ARRAYS.
!***  WE ARE WORKING WITH 4D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!
!
      DO NFIND=1,MAX_KOUNT                                                 
!
        IF(PHY_INT_STATE_4D_R(2,NFIND)=='H')THEN                           !<-- Take 4D real array data specified for history output
!
          LDIM1=LBOUND(R_4D(NFIND)%NAME,1)
          UDIM1=UBOUND(R_4D(NFIND)%NAME,1)
          LDIM2=LBOUND(R_4D(NFIND)%NAME,2)
          UDIM2=UBOUND(R_4D(NFIND)%NAME,2)
          LDIM3=LBOUND(R_4D(NFIND)%NAME,3)
          UDIM3=UBOUND(R_4D(NFIND)%NAME,3)
          LDIM4=LBOUND(R_4D(NFIND)%NAME,4)
          UDIM4=UBOUND(R_4D(NFIND)%NAME,4)
!
          DO N=int_state%INDX_RRW+1,UDIM4                                                  !<-- Loop through the tracers kind
          DO K=LDIM3,UDIM3                                                 !<-- Loop through the levels of the array
            WRITE(TRACERS_KIND,FMT)N
            WRITE(MODEL_LEVEL,FMT)K
!
              VBL_NAME=TRIM(PHY_INT_STATE_4D_R(1,NFIND))//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_4D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K,N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Physics 4-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 4D real array being inserted into the data Bundle
                                  ,copyflag     =COPYFLAG               &
                                  ,maxHaloUWidth=(/IHALO,JHALO/)        &
                                  ,maxHaloLWidth=(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 4D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL         &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert Physics 4-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE              &  !<-- The write component's output data Bundle
                                    ,field =FIELD                       &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                              !<-- Continue adding up all levels of 4D Fields
!
          ENDDO
         ENDDO
!
        ENDIF
!
        IF(PHY_INT_STATE_4D_R(3,NFIND)=='R')THEN                           !<-- Take 4D real array data specified for restart output
!
          LDIM1=LBOUND(R_4D(NFIND)%NAME,1)
          UDIM1=UBOUND(R_4D(NFIND)%NAME,1)
          LDIM2=LBOUND(R_4D(NFIND)%NAME,2)
          UDIM2=UBOUND(R_4D(NFIND)%NAME,2)
          LDIM3=LBOUND(R_4D(NFIND)%NAME,3)
          UDIM3=UBOUND(R_4D(NFIND)%NAME,3)
          LDIM4=LBOUND(R_4D(NFIND)%NAME,4)
          UDIM4=UBOUND(R_4D(NFIND)%NAME,4)
!
!         DO N=LDIM4,UDIM4                                                  !<-- Loop through the tracers kind
          DO N=int_state%INDX_RRW+1,UDIM4 
          DO K=LDIM3,UDIM3                                                 !<-- Loop through the levels of the array
            WRITE(TRACERS_KIND,FMT)N
            WRITE(MODEL_LEVEL,FMT)K
!
              VBL_NAME=TRIM(PHY_INT_STATE_4D_R(1,NFIND))//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D=>R_4D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K,N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Physics 4-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 4D real array being inserted into the data Bundle
                                  ,copyflag     =COPYFLAG               &
                                  ,maxHaloUWidth=(/IHALO,JHALO/)        &
                                  ,maxHaloLWidth=(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 4D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL         &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert Physics 4-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(bundle=RESTART_BUNDLE              &  !<-- The write component's restart data Bundle
                                    ,field =FIELD                       &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                              !<-- Continue adding up all levels of 4D Fields
!
          ENDDO
         ENDDO
!
        ENDIF
!
        IF(PHY_INT_STATE_4D_R(2,NFIND)=='*' .AND.                       &
           PHY_INT_STATE_4D_R(3,NFIND)=='*'      )THEN                     !<-- End of the 4D real array list
          EXIT
        ENDIF
!
      ENDDO
!-----------------------------------------------------------------------
!
      IF(MYPE==0)WRITE(0,*)' PHYSICS_OUTPUT: Number of 2-D Real Fields=',NUM_2D_FIELDS_R
!
!-----------------------------------------------------------------------
!***  INSERT THE OUTPUT DATA Bundle INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phsyics: Insert History Bundle into the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =IMP_STATE_WRITE                    &  !<-- The write component's import state
                        ,fieldbundle=HISTORY_BUNDLE                     &  !<-- The ESMF Bundle holding all Physics output data
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE OUTPUT DATA Bundle INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Phsyics: Insert Restart Bundle into the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =IMP_STATE_WRITE                    &  !<-- The write component's import state
                        ,fieldbundle=RESTART_BUNDLE                     &  !<-- The ESMF Bundle holding all Physics restart data
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
