!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
      USE MODULE_INCLUDE
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY: PHYSICS_INTERNAL_STATE 
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
      PUBLIC :: POINT_PHYSICS_OUTPUT
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_PHYSICS_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  This routine takes the user's selections for output quantities
!***  and inserts them into ESMF Bundles.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_Grid),INTENT(IN)    :: GRID                                !<-- The ESMF Grid
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE                    !<-- Import state for the Write gridded components
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                    !<-- The Physics internal state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE                             &
                           ,IMS,IME,JMS,JME                             &
                           ,IHALO,JHALO
!
      INTEGER(kind=KINT) :: INDEX_IKJ,K,LENGTH                          &
                           ,MP_PHYSICS,MYPE                             &
                           ,N,M,NDIM3,NFIND                             &
                           ,NUM_2D_FIELDS_I,NUM_2D_FIELDS_R             &
                           ,RC,RC_PHY_OUT                               &
                           ,SF_SURFACE_PHYSICS
!
      INTEGER(kind=KINT) :: LDIM1,LDIM2,LDIM3,LDIM4                     &
                           ,UDIM1,UDIM2,UDIM3,UDIM4
!
      TYPE(ESMF_FieldBundle),SAVE :: HISTORY_BUNDLE
      TYPE(ESMF_FieldBundle),SAVE :: RESTART_BUNDLE
!
!-----------------------------------------------------------------------
!***********************************************************************
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
!***  Extract the history and restart output Bundles from the Write 
!***  component's import state.  It already contains output variables from
!***  from the Dynamics.  We are preparing to add history and restart
!***  variables from the Physics.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract History Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =IMP_STATE_WRITE                    &  !<-- Take Bundle from the Write component's import state
                        ,itemName   ='History Bundle'                   &  !<-- The Bundle's name
                        ,fieldbundle=HISTORY_BUNDLE                     &  !<-- The Bundle object
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =IMP_STATE_WRITE                    &  !<-- Take Bundle from the Write component's import state
                        ,itemName   ='Restart Bundle'                   &  !<-- The Bundle's name
                        ,fieldbundle=RESTART_BUNDLE                     &  !<-- The Bundle object
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The microphysics scheme specification is needed in the output
!***  so add it directly.
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
      ELSEIF(int_state%MICROPHYSICS=='wsm6')THEN
        MP_PHYSICS=6
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
!***  The land surface scheme specification is needed in the output
!***  so add it directly.
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
!***  Put the history and restart data into the ESMF Bundles residing 
!***  in the Write component's import state.
!-----------------------------------------------------------------------
!
      CALL PUT_VARS_IN_BUNDLES(int_state%VARS, int_state%NUM_VARS, GRID, HISTORY_BUNDLE, RESTART_BUNDLE)
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
