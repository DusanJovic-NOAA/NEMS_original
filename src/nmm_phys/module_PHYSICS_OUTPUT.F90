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
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY: PHYSICS_INTERNAL_STATE 
      USE MODULE_ERR_MSG               ,ONLY: ERR_MSG,MESSAGE_CHECK
      USE MODULE_VARS
      USE MODULE_VARS_STATE
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
      TYPE(ESMF_Grid),INTENT(IN)    :: GRID                                !<-- The ESMF Grid
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE                    !<-- Import state for the Write gridded components
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                    !<-- The Physics internal state
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
                                     ,N,M,NDIM3,NFIND                   &
                                     ,NUM_2D_FIELDS_I,NUM_2D_FIELDS_R   &
                                     ,RC,RC_PHY_OUT                     &
                                     ,SF_SURFACE_PHYSICS
!
      INTEGER                     :: LDIM1,LDIM2,LDIM3,LDIM4            &
                                    ,UDIM1,UDIM2,UDIM3,UDIM4
!
!d      INTEGER(ESMF_KIND_I4),DIMENSION(:,:),POINTER :: TEMP_I2D
!d      REAL(ESMF_KIND_R4)   ,DIMENSION(:,:),POINTER :: TEMP_R2D
!
!d      CHARACTER(2)           :: MODEL_LEVEL,TRACERS_KIND
!d      CHARACTER(6)           :: FMT='(I2.2)'
!d      CHARACTER(ESMF_MAXSTR) :: VBL_NAME,VBL_NAME_X
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE
      TYPE(ESMF_FieldBundle),SAVE :: RESTART_BUNDLE
!
!d      TYPE(ESMF_Field)       :: FIELD
!
!d      TYPE(ESMF_CopyFlag)    :: COPYFLAG=ESMF_DATA_REF
!     TYPE(ESMF_CopyFlag)    :: COPYFLAG=ESMF_DATA_COPY
!
!-----------------------------------------------------------------------
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
      CALL PUT_VARS_IN_BUNDLES(int_state%VARS, int_state%NUM_VARS, GRID, HISTORY_BUNDLE, RESTART_BUNDLE)
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
