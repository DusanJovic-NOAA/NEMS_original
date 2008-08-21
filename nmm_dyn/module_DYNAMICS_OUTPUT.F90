!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_OUTPUT
!
!-----------------------------------------------------------------------
!***  LIST THE QUANTITIES FROM THE DYNAMICS INTERNAL STATE THAT
!***  CAN BE SELECTED FOR HISTORY OUTPUT AND ASSOCIATE THEM WITH
!***  UNIQUE INTEGERS.
!***  THE USER WILL PLACE AN 'H' IN THE 3rd ELEMENT OF THE FOLLOWING
!***  LIST OF QUANTITIES THAT ARE IN THE DYNAMICS INTERNAL STATE
!***  IF THAT QUANTITY IS TO BE WRITTEN TO THE HISTORY FILES.
!-----------------------------------------------------------------------
!***  WHEN NEW QUANTITIES ARE ADDED TO THE DYNAMICS INTERNAL STATE
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT THEN THEY NEED TO BE
!***  ADDED IN TWO PLACES BELOW: 
!*    (1) THE APPROPRIATE DATA LIST PRECEDING THE 'CONTAINS' STATEMENT
!*    (2) 'THE DYNAMICS INTERNAL STATE POINTER BLOCK'
!*        IN SUBROUTINE POINT_DYNAMICS_OUPUT
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
      USE MODULE_INCLUDE
      USE MODULE_DYNAMICS_INTERNAL_STATE,ONLY: INTERNAL_STATE 
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: POINT_DYNAMICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_KOUNT=100
!
!-----------------------------------------------------------------------
!***  LIST THE VARIOUS QUANTITIES IN THE DYNAMICS INTERNAL STATE 
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_ISCALAR     &
!
       =RESHAPE((/                                                      &
!                                         ------------------------------
!
                                           'IM        ', 'H         '   &
                                          ,'JM        ', 'H         '   &
                                          ,'LM        ', 'H         '   &
                                          ,'IHRST     ', 'H         '   &
!
!                                         ------------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_RSCALAR     &  
!
       =RESHAPE((/                                                      &
!                                          -----------------------------
!
                                           'DT        ', 'H         '   &
                                          ,'DYH       ', 'H         '   &
                                          ,'PDTOP     ', 'H         '   &
                                          ,'PT        ', 'H         '   &
                                          ,'TLM0D     ', 'H         '   &
                                          ,'TPH0D     ', 'H         '   &
                                          ,'TSTART    ', 'H         '   &
!
!                                          -----------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_1D_I        & 
!
       =RESHAPE((/                                                      &
!                                          -----------------------------
!
                                           'IDAT      ', 'H         '   &
!
!                                          -----------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_2D_I        & 
!
       =RESHAPE((/                                                      &
!                                          -----------------------------
!
                                           '-         ', '-         '   &
!
!                                          -----------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_1D_R        & 
!
       =RESHAPE((/                                                      &
!                                         ------------------------------
!
                                           'DXH       ', 'H         '   &
                                          ,'SG1       ', 'H         '   &
                                          ,'SG2       ', 'H         '   &
                                          ,'DSG1      ', 'H         '   &
                                          ,'DSG2      ', 'H         '   &
                                          ,'SGML1     ', 'H         '   &
                                          ,'SGML2     ', 'H         '   &
!
!                                         ------------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_2D_R        & 
!
       =RESHAPE((/                                                      &
!                                         ------------------------------
!
                                           'FIS       ', 'H         '   &
                                          ,'GLAT      ', 'H         '   &
                                          ,'GLON      ', 'H         '   &
                                          ,'PD        ', 'H         '   &
                                          ,'VLAT      ', 'H         '   &
                                          ,'VLON      ', 'H         '   &
!
!                                         ------------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
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
      CHARACTER(10),DIMENSION(2,MAX_KOUNT) :: DYN_INT_STATE_3D_R        & 
!
       =RESHAPE((/                                                      &
!                                         ------------------------------
!
                                           'T         ', '-         '   &  !<-- The physics counterparts of these variables
                                          ,'Q         ', '-         '   &  !    are being designated for history output.
                                          ,'U         ', '-         '   &  !    This assumes that history output always
                                          ,'V         ', '-         '   &  !    immediately follows a call to the Physics.
                                          ,'Q2        ', '-         '   &  !
                                          ,'CW        ', '-         '   &  !<--
                                          ,'W         ', 'H         '   &
                                          ,'DWDT      ', 'H         '   &
                                          ,'PINT      ', 'H         '   &
                                          ,'OMGALF    ', 'H         '   &
                                          ,'RRW       ', 'H         '   &
!
!                                         ------------------------------
!
         /)                                                             &
        ,(/2,MAX_KOUNT/)                                                &
        ,(/'**********', '**********'/))
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_DYNAMICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
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
      TYPE(ESMF_State)    ,INTENT(INOUT) :: IMP_STATE_WRITE              !<-- Import state for the Write components
!
      TYPE(INTERNAL_STATE),POINTER       :: INT_STATE                    !<-- The dynamics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER                     :: ITS,ITE,JTS,JTE                    &
                                    ,IMS,IME,JMS,JME                    &
                                    ,IHALO,JHALO
!
      INTEGER                     :: K,LENGTH,MYPE                      &
                                    ,N,NDIM3,NFIND,NUM_2D_FIELDS        &
                                    ,RC,RC_DYN_OUT
!
      INTEGER                     :: LDIM1,LDIM2                        &
                                    ,UDIM1,UDIM2
!
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: TEMP_R2D
!
      CHARACTER(2)                :: MODEL_LEVEL
      CHARACTER(6)                :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR)      :: VBL_NAME
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE
!
      TYPE(ESMF_Field)            :: FIELD
!
      TYPE(ESMF_CopyFlag)         :: COPYFLAG=ESMF_DATA_REF
!     TYPE(ESMF_CopyFlag)         :: COPYFLAG=ESMF_DATA_COPY
!
!-----------------------------------------------------------------------
!***  ESMF VERSIONS OF THE LOGICALS IN THE DYNAMICS INTERNAL STATE
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Logical),TARGET :: GLOBAL_ESMF,RUN_ESMF
!
!-----------------------------------------------------------------------
!***  FIRST WE MUST PROVIDE POINTERS INTO THE DYNAMICS INTERNAL STATE. 
!-----------------------------------------------------------------------
      TYPE DYN_ISC
        INTEGER(KIND=KINT),POINTER :: NAME                               !<-- Pointer for integer scalars
      END TYPE DYN_ISC
!-----------------------------------------------------------------------
      TYPE DYN_RSC
        REAL(KIND=KFPT),POINTER :: NAME                                  !<-- Pointer for real scalars
      END TYPE DYN_RSC
!-----------------------------------------------------------------------
      TYPE DYN_I1D
        INTEGER(KIND=KINT),DIMENSION(:),POINTER :: NAME                  !<-- Pointer for 1D integer arrays
      END TYPE DYN_I1D
!-----------------------------------------------------------------------
      TYPE DYN_I2D
        INTEGER(KIND=KINT),DIMENSION(:,:),POINTER :: NAME                !<-- Pointer for 2D integer arrays
      END TYPE DYN_I2D
!-----------------------------------------------------------------------
      TYPE DYN_R1D
        REAL(KIND=KFPT),DIMENSION(:),POINTER :: NAME                     !<-- Pointer for 1D real arrays
      END TYPE DYN_R1D
!-----------------------------------------------------------------------
      TYPE DYN_R2D
        REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: NAME                   !<-- Pointer for 2D real arrays
      END TYPE DYN_R2D
!-----------------------------------------------------------------------
      TYPE DYN_R3D
        REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER :: NAME                 !<-- Pointer for 3D real arrays
      END TYPE DYN_R3D
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARRAYS OF POINTERS OF THE ABOVE TYPES
!-----------------------------------------------------------------------
!
      TYPE(DYN_ISC),DIMENSION(MAX_KOUNT) :: I_SC
      TYPE(DYN_RSC),DIMENSION(MAX_KOUNT) :: R_SC
      TYPE(DYN_I1D),DIMENSION(MAX_KOUNT) :: I_1D
      TYPE(DYN_R1D),DIMENSION(MAX_KOUNT) :: R_1D
      TYPE(DYN_R2D),DIMENSION(MAX_KOUNT) :: R_2D
      TYPE(DYN_R3D),DIMENSION(MAX_KOUNT) :: R_3D
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE DYNAMICS INTERNAL STATE POINTER BLOCK
!-----------------------------------------------------------------------
!***  POINT AT ALL VARIABLES IN THE DYNAMICS INTERNAL STATE 
!***  THAT COULD BE WRITTEN TO HISTORY OUTPUT, I.E., THOSE
!***  LISTED AT THE TOP OF THIS MODULE.
!-----------------------------------------------------------------------
!
!***  INTEGER SCALARS
!
      I_SC(1)%NAME=>int_state%IM
      I_SC(2)%NAME=>int_state%JM
      I_SC(3)%NAME=>int_state%LM
      I_SC(4)%NAME=>int_state%IHRST
!        
!***  REAL SCALARS
!
      R_SC(1)%NAME=>int_state%DT
      R_SC(2)%NAME=>int_state%DYH
      R_SC(3)%NAME=>int_state%PDTOP
      R_SC(4)%NAME=>int_state%PT
      R_SC(5)%NAME=>int_state%TLM0D
      R_SC(6)%NAME=>int_state%TPH0D
      R_SC(7)%NAME=>int_state%TSTART
!        
!***  1D INTEGER ARRAYS
!
      I_1D(1)%NAME=>int_state%IDAT
!        
!***  1D REAL ARRAYS
!
      R_1D(1)%NAME=>int_state%DXH
      R_1D(2)%NAME=>int_state%SG1
      R_1D(3)%NAME=>int_state%SG2
      R_1D(4)%NAME=>int_state%DSG1
      R_1D(5)%NAME=>int_state%DSG2
      R_1D(6)%NAME=>int_state%SGML1
      R_1D(7)%NAME=>int_state%SGML2
!        
!***  2D REAL ARRAYS
!
      R_2D(1)%NAME=>int_state%FIS
      R_2D(2)%NAME=>int_state%GLAT
      R_2D(3)%NAME=>int_state%GLON
      R_2D(4)%NAME=>int_state%PD
      R_2D(5)%NAME=>int_state%VLAT
      R_2D(6)%NAME=>int_state%VLON
!        
!***  3D REAL ARRAYS
!
      R_3D(1)%NAME=>int_state%T
      R_3D(2)%NAME=>int_state%Q
      R_3D(3)%NAME=>int_state%U
      R_3D(4)%NAME=>int_state%V
      R_3D(5)%NAME=>int_state%Q2
      R_3D(6)%NAME=>int_state%CW
      R_3D(7)%NAME=>int_state%W
      R_3D(8)%NAME=>int_state%DWDT
      R_3D(9)%NAME=>int_state%PINT
      R_3D(10)%NAME=>int_state%OMGALF
      R_3D(11)%NAME=>int_state%RRW
!
!-----------------------------------------------------------------------
!***  ESMF VERSION OF LOGICALS NEEDED FOR THEIR INSERTION
!***  INTO THE HISTORY OUTPUT Bundle OF THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
      RUN_ESMF=ESMF_False
      GLOBAL_ESMF=ESMF_False
!
      IF(int_state%RUN)RUN_ESMF=ESMF_True
      IF(int_state%GLOBAL)GLOBAL_ESMF=ESMF_True
!
!-----------------------------------------------------------------------
!
      ITS=int_state%ITS
      ITE=int_state%ITE
      JTS=int_state%JTS
      JTE=int_state%JTE
!
      MYPE=int_state%MYPE
!
!-----------------------------------------------------------------------
!***  CREATE AN ESMF Bundle THAT WILL HOLD HISTORY OUTPUT DATA
!***  AND NOTHING ELSE.  THIS WILL SERVE TO ISOLATE THE OUTPUT
!***  DATA FROM EVERYTHING ELSE INSIDE THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create History Data Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      HISTORY_BUNDLE=ESMF_FieldBundleCreate(grid=GRID                   &  !<-- The ESMF integration Grid
                                           ,name='Bundle_Output_Data'   &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FIRST ADD THE LOCAL SUBDOMAIN LIMITS TO THE WRITE COMPONENT'S
!***  IMPORT STATE AS ATTRIBUTES ALONG WITH THE GLOBAL/REGIONAL MODE.
!***  THIS INFORMATION IS NEEDED FOR QUILTING THE LOCAL DOMAIN DATA
!***  INTO FULL DOMAIN FIELDS.
!***  THE LOCAL DOMAIN LIMITS GO DIRECTLY INTO THE WRITE COMPONENT'S
!***  IMPORT STATE TO KEEP THEM SEPARATE FROM THE HISTORY DATA THAT 
!***  WILL BE INSERTED INTO A Bundle.
!
!***  DO THE SAME WITH THE NUMBER OF FCST TASKS (INPESxJNPES) SINCE
!***  THE WRITE COMPONENT ALSO NEEDS THAT INFORMATION AS WELL AS
!***  THE HALO DEPTHS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Local Subdomain Limits to the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_ISTART'                   &  !<-- Name of the integer array
                            ,count    =int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_ISTART           &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_IEND'                     &  !<-- Name of the integer array
                            ,count    =int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_IEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JSTART'                   &  !<-- Name of the integer array
                            ,count    =int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JSTART           &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JEND'                     &  !<-- Name of the integer array
                            ,count    =int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='INPES'                              &  !<-- Name of the integer scalar
                            ,value=int_state%INPES                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JNPES'                              &  !<-- Name of the integer scalar
                            ,value=int_state%JNPES                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IHALO'                              &  !<-- Name of the integer scalar
                            ,value=int_state%IHALO                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JHALO'                              &  !<-- Name of the integer scalar
                            ,value=int_state%JHALO                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='WRITE_TASKS_PER_GROUP'              &  !<-- Name of the integer scalar
                            ,value=int_state%WRITE_TASKS_PER_GROUP      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='WRITE_GROUPS'                       &  !<-- Name of the integer scalar
                            ,value=int_state%WRITE_GROUPS               &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING LOGICAL VARIABLES ARE TO BE PART OF THE
!***  HISTORY OUTPUT THEREFORE PLACE THEM INTO THE OUTPUT Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component output data Bundle
                            ,name  ='GLOBAL'                            &  !<-- Name of the logical
                            ,value =GLOBAL_ESMF                         &  !<-- The logical being inserted into the import state
                            ,rc    =RC)
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component output data Bundle
                            ,name  ='RUN'                               &  !<-- Name of the logical
                            ,value =RUN_ESMF                             &  !<-- The logical being inserted into the import state
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER FOR
!***  HISTORY OUTPUT.  THE DATA IS PLACED INTO AN ESMF Bundle WHICH
!***  ITSELF WILL BE PLACED INTO THE IMPORT STATE AT THE END OF
!***  THE ROUTINE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  BEGIN WITH THE INTEGER SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Dynamics Integer Scalars into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(DYN_INT_STATE_ISCALAR(2,NFIND)=='H')THEN                        !<-- Take integer scalar data specified for history output
          VBL_NAME=TRIM(DYN_INT_STATE_ISCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                  &  !<-- The Write component output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =I_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the import state
                                ,rc    =RC)
!
        ELSEIF(DYN_INT_STATE_ISCALAR(2,NFIND)=='*')THEN                    !<-- End of the integer scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE REAL SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Dynamics Real Scalars into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(DYN_INT_STATE_RSCALAR(2,NFIND)=='H')THEN                        !<-- Take real scalar data specified for history output
          VBL_NAME=TRIM(DYN_INT_STATE_RSCALAR(1,NFIND))
!
          CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                  &  !<-- The Write component output data Bundle
                                ,name  =VBL_NAME                        &  !<-- Name of the integer scalar
                                ,value =R_SC(NFIND)%NAME                &  !<-- The scalar being inserted into the import state
                                ,rc    =RC)
!
        ELSEIF(DYN_INT_STATE_RSCALAR(2,NFIND)=='*')THEN                    !<-- End of the real scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D INTEGER ARRAYS
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Dynamics 1-D Integer Arrays into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(DYN_INT_STATE_1D_I(2,NFIND)=='H')THEN                           !<-- Take 1D integer array data specified for history output
          VBL_NAME=TRIM(DYN_INT_STATE_1D_I(1,NFIND))
          LENGTH=SIZE(I_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =HISTORY_BUNDLE               &  !<-- The Write component output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=I_1D(NFIND)%NAME             &  !<-- The 1D integer being inserted into the import state
                                ,rc       =RC)
!
        ELSEIF(DYN_INT_STATE_1D_I(2,NFIND)=='*')THEN                       !<-- End of the 1D integer array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D REAL ARRAYS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Dynamics 1-D Real Arrays into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(DYN_INT_STATE_1D_R(2,NFIND)=='H')THEN                           !<-- Take 1D real array data specified for history output
          VBL_NAME=TRIM(DYN_INT_STATE_1D_R(1,NFIND))
          LENGTH=SIZE(R_1D(NFIND)%NAME)
!
          CALL ESMF_AttributeSet(bundle   =HISTORY_BUNDLE               &  !<-- The Write component output data Bundle
                                ,name     =VBL_NAME                     &  !<-- Name of the integer scalar
                                ,count    =LENGTH                       &  !<-- # of elements in this attribute
                                ,valueList=R_1D(NFIND)%NAME             &  !<-- The 1D real being inserted into the import state
                                ,rc       =RC)
!
        ELSEIF(DYN_INT_STATE_1D_R(2,NFIND)=='*')THEN                       !<-- End of the 1D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DEREFERENCE THE LOCAL MEMORY LIMTS.  THESE DETERMINE THE TOTAL
!***  NUMBER OF WORDS LOADED PER SUBDOMAIN.
!-----------------------------------------------------------------------
!
      IMS=int_state%IMS
      IME=int_state%IME
      JMS=int_state%JMS
      JME=int_state%JME
!
!-----------------------------------------------------------------------
!***  THE 2D REAL ARRAYS.
!-----------------------------------------------------------------------
!
      IHALO=int_state%IHALO
      JHALO=int_state%JHALO
!
      NUM_2D_FIELDS=0
      NULLIFY(TEMP_R2D)
!
      DO NFIND=1,MAX_KOUNT
        IF(DYN_INT_STATE_2D_R(2,NFIND)=='H')THEN                           !<-- Take 2D real array data specified for history output
          VBL_NAME=TRIM(DYN_INT_STATE_2D_R(1,NFIND))
!!!       TEMP_R2D=>R_2D(NFIND)%NAME(ITS:ITE,JTS:JTE)
          TEMP_R2D=>R_2D(NFIND)%NAME(IMS:IME,JMS:JME)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Dynamics 2-D Real Data into Field"
          CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_R2D                 &  !<-- The 2D real array being inserted into the import state
                                ,copyflag     =COPYFLAG                 &
                                ,maxHaloUWidth=(/IHALO,JHALO/)          &
                                ,maxHaloLWidth=(/IHALO,JHALO/)          &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Dynamics 2-D Real Field into History Bundle"
          CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE                &  !<-- The Write component output data Bundle
                                  ,field =FIELD                         &  !<-- ESMF Field holding the 2D real array
                                  ,rc    =RC)
!
          NUM_2D_FIELDS=NUM_2D_FIELDS+1
!
        ELSEIF(DYN_INT_STATE_2D_R(2,NFIND)=='*')THEN                       !<-- End of the 2D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 3D REAL ARRAYS.
!***  WE ARE WORKING WITH 3D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!
      NULLIFY(TEMP_R2D)
!
      DO NFIND=1,MAX_KOUNT                                                 
!
        IF(DYN_INT_STATE_3D_R(2,NFIND)=='H')THEN                           !<-- Take 3D real array data specified for history output
          NDIM3=UBOUND(R_3D(NFIND)%NAME,3)                                 !<-- Determine # of vertical levels in this variable
          LDIM1=LBOUND(R_3D(NFIND)%NAME,1)
          UDIM1=UBOUND(R_3D(NFIND)%NAME,1)
          LDIM2=LBOUND(R_3D(NFIND)%NAME,2)
          UDIM2=UBOUND(R_3D(NFIND)%NAME,2)
!
          DO K=1,NDIM3
            WRITE(MODEL_LEVEL,FMT)K
            VBL_NAME=TRIM(DYN_INT_STATE_3D_R(1,NFIND))//'_'//MODEL_LEVEL//'_2D'
            TEMP_R2D=>R_3D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Dynamics 3-D Data"
            CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 3D real array being inserted into the import state
                                  ,copyflag     =COPYFLAG               &
                                  ,maxHaloUWidth=(/IHALO,JHALO/)        &
                                  ,maxHaloLWidth=(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 3D real array
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert Dynamics 3-D Data into History Bundle"
            CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(bundle=HISTORY_BUNDLE              &  !<-- The Write component output data Bundle
                                    ,field =FIELD                       &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS=NUM_2D_FIELDS+1
!
          ENDDO
!
        ELSEIF(DYN_INT_STATE_3D_R(2,NFIND)=='*')THEN                       !<-- End of the 3D real array list
          EXIT
        ENDIF
!
      ENDDO
      IF(MYPE==0)WRITE(0,*)' Exit DYNAMICS_OUTPUT num_2d_fields=',NUM_2D_FIELDS
!
!-----------------------------------------------------------------------
!***  INSERT THE OUTPUT DATA Bundle INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!***  SINCE DYNAMICS IS CALLED BEFORE PHYSICS, WE WILL INSERT THE
!***  BUNDLE NOW AND SIMPLY USE IT IN POINT_PHYSICS_OUTPUT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dynamics: Insert History Bundle into the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =IMP_STATE_WRITE                    &  !<-- The Write component's import state
                        ,fieldbundle=HISTORY_BUNDLE                     &  !<-- The ESMF Bundle holding all Dynamics output data
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_DYNAMICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_OUTPUT
!
!-----------------------------------------------------------------------
