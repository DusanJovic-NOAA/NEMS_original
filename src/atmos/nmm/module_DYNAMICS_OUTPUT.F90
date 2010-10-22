!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_OUTPUT
!
!-----------------------------------------------------------------------
!***  Insert quantities from the Dynamics internal state into the
!***  Write import state for output.
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
!***  This routine takes the user's selections for output quantities,
!***  points at them, and inserts those pointers into the import state
!***  of the Write components.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Grid) ,INTENT(IN) :: GRID                                  !<-- The ESMF Grid
!
      TYPE(DYNAMICS_INTERNAL_STATE),POINTER,INTENT(IN) :: INT_STATE        !<-- The Dynamics internal state
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE                    !<-- Import state for the Write components
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
      INTEGER :: ITWO=2
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
      TYPE(ESMF_FieldBundle),DIMENSION(1:2) :: BUNDLE_ARRAY
!
      TYPE(ESMF_Field)            :: FIELD
!
      TYPE(ESMF_CopyFlag)         :: COPYFLAG=ESMF_DATA_REF
!     TYPE(ESMF_CopyFlag)         :: COPYFLAG=ESMF_DATA_COPY
!
!-----------------------------------------------------------------------
!***  ESMF versions of the logicals in the Dynamics internal state.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Logical),TARGET :: ADIABATIC_ESMF                       &
                                  ,GLOBAL_ESMF                          &
                                  ,RUN_ESMF
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ESMF version of logicals needed for their insertion
!***  into the history output Bundle of the Write component's
!***  import state.
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
!***  Create an ESMF Bundle that will hold history output data
!***  and nothing else.  This will serve to isolate the output
!***  data from everything else inside the Write component's
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create History Data Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      HISTORY_BUNDLE=ESMF_FieldBundleCreate(grid=GRID                   &  !<-- The ESMF integration Grid
                                           ,name='History Bundle'       &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create an ESMF Bundle that will hold restart data
!***  and nothing else.  This will serve to isolate the restart
!***  data from everything else inside the Write component's
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Restart Data Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      RESTART_BUNDLE=ESMF_FieldBundleCreate(grid=GRID                   &  !<-- The ESMF integration Grid
                                           ,name='Restart Bundle'       &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  First add the local subdomain limits to the Write component's
!***  import state as Attributes along with the global/regional mode.
!***  This information is needed for quilting the local domain data
!***  into full domain fields.
!***  The local domain limits go directly into the Write component's
!***  import state to keep them separate from the history data that 
!***  will be inserted into a Bundle.
!
!***  Do the same with the number of fcst tasks (INPESxJNPES) since
!***  the Write component also needs that information as well as
!***  the halo depths.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Local Subdomain Limits to the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
!***  The following logical variables are to be part of the
!***  history output therefore place them into the history Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
!***  The following logical variables are to be part of the
!***  restart output therefore place them into the restart Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
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
!***  Now insert into the Write components' import state the pointers
!***  of only those quantities that are specified by the user for
!***  history output.  The data is placed into an ESMF Bundle which
!***  itself will be placed into the import state at the end of
!***  the routine.
!-----------------------------------------------------------------------
!
      CALL PUT_VARS_IN_BUNDLES(int_state%VARS                           &
                              ,int_state%NUM_VARS                       &
                              ,GRID                                     &
                              ,HISTORY_BUNDLE                           &
                              ,RESTART_BUNDLE)
!
!-----------------------------------------------------------------------
!***  Load the two output Bundles into the working array which is used
!***  to add them to the Write component's import state.
!-----------------------------------------------------------------------
!
      BUNDLE_ARRAY(1)=HISTORY_BUNDLE
      BUNDLE_ARRAY(2)=RESTART_BUNDLE
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dynamics: Insert Bundle Array into the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state          =IMP_STATE_WRITE                &  !<-- The Write component's import state
                        ,fieldbundlelist=BUNDLE_ARRAY                   &  !<-- Array holding the History/Restart Bundles
                        ,count          =ITWO                           &  !<-- There are two Bundles in the array
                        ,rc             =RC)
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
