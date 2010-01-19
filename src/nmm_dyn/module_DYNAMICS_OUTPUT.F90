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
      USE MODULE_DYNAMICS_INTERNAL_STATE,ONLY: DYNAMICS_INTERNAL_STATE 
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
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
      PUBLIC :: POINT_DYNAMICS_OUTPUT
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
      TYPE(ESMF_Grid) ,INTENT(IN)           :: GRID                        !<-- The ESMF Grid
      TYPE(ESMF_State),INTENT(INOUT)        :: IMP_STATE_WRITE             !<-- Import state for the Write components
      TYPE(DYNAMICS_INTERNAL_STATE),POINTER,INTENT(IN) :: INT_STATE        !<-- The Dynamics internal state
!
!-----------------------------------------------------------------------
!
      INTEGER :: IHALO,JHALO
!
      INTEGER :: K,LENGTH,MYPE                                          &
                ,N,NDIM3,NFIND,NUM_2D_FIELDS                            &
                ,RC,RC_DYN_OUT
!
      INTEGER :: LDIM1,LDIM2                                            &
                ,UDIM1,UDIM2
!
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: TEMP_R2D
!
      CHARACTER(2)                :: MODEL_LEVEL
      CHARACTER(6)                :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR)      :: VBL_NAME
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE
      TYPE(ESMF_FieldBundle),SAVE :: RESTART_BUNDLE
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
      TYPE(ESMF_Logical),TARGET :: ADIABATIC_ESMF                       &
                                  ,GLOBAL_ESMF                          &
                                  ,RUN_ESMF
!
!-----------------------------------------------------------------------
!***  ESMF VERSION OF LOGICALS NEEDED FOR THEIR INSERTION
!***  INTO THE HISTORY OUTPUT Bundle OF THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
      RUN_ESMF      =ESMF_False
      GLOBAL_ESMF   =ESMF_False
      ADIABATIC_ESMF=ESMF_False
!
      IF(int_state%RUN)RUN_ESMF=ESMF_True
      IF(int_state%GLOBAL)GLOBAL_ESMF=ESMF_True
      IF(int_state%ADIABATIC)ADIABATIC_ESMF=ESMF_True
!
!-----------------------------------------------------------------------
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
!***  CREATE AN ESMF Bundle THAT WILL HOLD RESTART DATA
!***  AND NOTHING ELSE.  THIS WILL SERVE TO ISOLATE THE RESTART
!***  DATA FROM EVERYTHING ELSE INSIDE THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Restart Data Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      RESTART_BUNDLE=ESMF_FieldBundleCreate(grid=GRID                   &  !<-- The ESMF integration Grid
                                           ,name='Bundle_Restart_Data'  &  !<-- The Bundle's name
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
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='LNSV'                               &  !<-- Name of the integer scalar
                            ,value=int_state%LNSV                       &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JEND'                     &  !<-- Name of the integer array
                            ,count    =int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IDS'                                &  !<-- Name of the integer scalar
                            ,value=int_state%IDS                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IDE'                                &  !<-- Name of the integer scalar
                            ,value=int_state%IDE                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JDS'                                &  !<-- Name of the integer scalar
                            ,value=int_state%JDS                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JDE'                                &  !<-- Name of the integer scalar
                            ,value=int_state%JDE                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING LOGICAL VARIABLES ARE TO BE PART OF THE
!***  HISTORY OUTPUT THEREFORE PLACE THEM INTO THE HISTORY Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into History Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component output history Bundle
                            ,name  ='GLOBAL'                            &  !<-- Name of the logical
                            ,value =GLOBAL_ESMF                         &  !<-- The logical being inserted into the Bundle
                            ,rc    =RC)
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component output history Bundle
                            ,name  ='RUN'                               &  !<-- Name of the logical
                            ,value =RUN_ESMF                            &  !<-- The logical being inserted into the Bundle
                            ,rc    =RC)
!
      CALL ESMF_AttributeSet(bundle=HISTORY_BUNDLE                      &  !<-- The Write component output history Bundle
                            ,name  ='ADIABATIC'                         &  !<-- Name of the logical
                            ,value =ADIABATIC_ESMF                      &  !<-- The logical being inserted into the Bundle
                            ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING LOGICAL VARIABLES ARE TO BE PART OF THE
!***  RESTART THEREFORE PLACE THEM INTO THE RESTART Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into Restart Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                      &  !<-- The Write component restart Bundle
                            ,name  ='GLOBAL'                            &  !<-- Name of the logical
                            ,value =GLOBAL_ESMF                         &  !<-- The logical being inserted into the Bundle
                            ,rc    =RC)
!
      CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                      &  !<-- The Write component restart Bundle
                            ,name  ='RUN'                               &  !<-- Name of the logical
                            ,value =RUN_ESMF                            &  !<-- The logical being inserted into the Bundle
                            ,rc    =RC)
!
      CALL ESMF_AttributeSet(bundle=RESTART_BUNDLE                      &  !<-- The Write component restart Bundle
                            ,name  ='ADIABATIC'                         &  !<-- Name of the logical
                            ,value =ADIABATIC_ESMF                      &  !<-- The logical being inserted into the Bundle
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
      CALL PUT_VARS_IN_BUNDLES(int_state%VARS, int_state%NUM_VARS, GRID, HISTORY_BUNDLE, RESTART_BUNDLE)
!
!-----------------------------------------------------------------------
!***  INSERT THE HISTORY DATA Bundle INTO THE WRITE COMPONENT'S
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
                        ,fieldbundle=HISTORY_BUNDLE                     &  !<-- The ESMF Bundle holding all Dynamics history data
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE RESTART DATA Bundle INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!***  SINCE DYNAMICS IS CALLED BEFORE PHYSICS, WE WILL INSERT THE
!***  BUNDLE NOW AND SIMPLY USE IT IN POINT_PHYSICS_OUTPUT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dynamics: Insert Restart Bundle into the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =IMP_STATE_WRITE                    &  !<-- The Write component's import state
                        ,fieldbundle=RESTART_BUNDLE                     &  !<-- The ESMF Bundle holding all Dynamics restart data
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
