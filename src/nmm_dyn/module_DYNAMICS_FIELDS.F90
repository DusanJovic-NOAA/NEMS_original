!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_FIELDS
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE CREATES/DESTROYS ESMF Fields FOR THE DYNAMICS
!***  IMPORT/EXPORT STATES.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_DYNAMICS_INTERNAL_STATE,ONLY : INTERNAL_STATE
      USE MODULE_DM_PARALLEL            ,ONLY : IDS,IDE,JDS,JDE         &
                                               ,IMS,IME,JMS,JME         &
                                               ,ITS,ITE,JTS,JTE         &
                                               ,MYPE_SHARE,LM           &
                                               ,IHALO,JHALO
!
      USE MODULE_ERR_MSG                ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: FIELD_T                                                 &
               ,FIELD_U                                                 &
               ,FIELD_V                                                 &
               ,FIELD_Q2                                                &
               ,FIELD_OMGALF
!
      PUBLIC :: FIELD_PD
!
      PUBLIC :: FIELD_TRACERS
!
      PUBLIC :: ALLOC_FIELDS_DYN
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CopyFlag) :: COPYFLAG=ESMF_DATA_REF
!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!***  LIST ESMF Fields THAT ARE TO BE PART OF
!***  THE DYNAMICS IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Field) ,SAVE :: FIELD_T                                 &
                               ,FIELD_U                                 &
                               ,FIELD_V                                 &
                               ,FIELD_Q2                                &
                               ,FIELD_OMGALF                            &
                               ,FIELD_PD                                &
                               ,FIELD_TRACERS
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ALLOC_FIELDS_DYN(GRID,INT_STATE)
!
!-----------------------------------------------------------------------
!***  CREATE ESMF Fields FOR THE DYNAMICS IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Grid),INTENT(IN)   :: GRID                               !<-- The ESMF grid
!
      TYPE(INTERNAL_STATE),POINTER :: INT_STATE                          !<-- The Dynamics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER       :: I,ISTATUS,J,L,N,RC,RC_FINAL,NUM_TRAC
!
      CHARACTER(20) :: FIELD_NAME
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF FIELDS THAT WILL BE ADDED TO THE IMPORT/EXPORT
!***  STATES AND ASSOCIATE THE APPROPRIATE POINTERS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   T   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_T"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='T'
!
      FIELD_T=ESMF_FieldCreate(grid            =GRID                    &  !<-- The ESMF grid
                              ,farray          =int_state%T             &  !<-- Insert this pointer into the Field
                              ,maxHaloUWidth   =(/IHALO,JHALO/)         &  !<-- Upper bound of halo region
                              ,maxHaloLWidth   =(/IHALO,JHALO/)         &  !<-- Lower bound of halo region
                              ,ungriddedLBound =(/1/)                   &
                              ,ungriddedUBound =(/LM/)                  &
                              ,name            =FIELD_NAME              &  !<-- Name of the 2D real array
                              ,indexFlag       =ESMF_INDEX_DELOCAL      &
                              ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   U   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_U"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='U'
!
      FIELD_U=ESMF_FieldCreate(grid            =GRID                    &  !<-- The ESMF grid
                              ,farray          =int_state%U             &  !<-- Insert this pointer into the Field
                              ,maxHaloUWidth   =(/IHALO,JHALO/)         &  !<-- Upper bound of halo region
                              ,maxHaloLWidth   =(/IHALO,JHALO/)         &  !<-- Lower bound of halo region
                              ,ungriddedLBound =(/1/)                   &
                              ,ungriddedUBound =(/LM/)                  &
                              ,name            =FIELD_NAME              &  !<-- Name of the 2D real array
                              ,indexFlag       =ESMF_INDEX_DELOCAL      &
                              ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   V   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_V"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='V'
!
      FIELD_V=ESMF_FieldCreate(grid            =GRID                    &  !<-- The ESMF grid
                              ,farray          =int_state%V             &  !<-- Insert this pointer into the Field
                              ,maxHaloUWidth   =(/IHALO,JHALO/)         &  !<-- Upper bound of halo region
                              ,maxHaloLWidth   =(/IHALO,JHALO/)         &  !<-- Lower bound of halo region
                              ,ungriddedLBound =(/1/)                   &
                              ,ungriddedUBound =(/LM/)                  &
                              ,name            =FIELD_NAME              &  !<-- Name of the 2D real array
                              ,indexFlag       =ESMF_INDEX_DELOCAL      &
                              ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   Q2  - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_Q2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='Q2'
!
      FIELD_Q2=ESMF_FieldCreate(grid            =GRID                   &  !<-- The ESMF grid
                               ,farray          =int_state%Q2           &  !<-- Insert this pointer into the Field
                               ,maxHaloUWidth   =(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                               ,maxHaloLWidth   =(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                               ,ungriddedLBound =(/1/)                  &
                               ,ungriddedUBound =(/LM/)                 &
                               ,name            =FIELD_NAME             &  !<-- Name of the 2D real array
                               ,indexFlag       =ESMF_INDEX_DELOCAL     &
                               ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - -  OMGALF  - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_OMGALF"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='OMGALF'
!
      FIELD_OMGALF=ESMF_FieldCreate(grid            =GRID               &  !<-- The ESMF grid
                                   ,farray          =int_state%OMGALF   &  !<-- Insert this pointer into the Field
                                   ,maxHaloUWidth   =(/IHALO,JHALO/)    &  !<-- Upper bound of halo region
                                   ,maxHaloLWidth   =(/IHALO,JHALO/)    &  !<-- Lower bound of halo region
                                   ,ungriddedLBound =(/1/)              &
                                   ,ungriddedUBound =(/LM/)             &
                                   ,name            =FIELD_NAME         &  !<-- Name of the 2D real array
                                   ,indexFlag       =ESMF_INDEX_DELOCAL &
                                   ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DO THE 2-D FIELDS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -  PD  - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_PD"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='PD'
!
      FIELD_PD=ESMF_FieldCreate(grid          =GRID                     &  !<-- The ESMF grid
                               ,farray        =int_state%PD             &  !<-- Insert this pointer into the Field
                               ,maxHaloUWidth =(/IHALO,JHALO/)          &  !<-- Upper bound of halo region
                               ,maxHaloLWidth =(/IHALO,JHALO/)          &  !<-- Lower bound of halo region
                               ,name          =FIELD_NAME               &  !<-- Name of the 2D real array
                               ,indexFlag     =ESMF_INDEX_DELOCAL       &
                               ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE AND FILL THE Fields FOR THE 4-D TRACERS FIELD.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - -  TRACERS - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create FIELD_TRACERS"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      FIELD_NAME='TRACERS'
!
      NUM_TRAC=int_state%NUM_TRACERS_TOTAL
!
      FIELD_TRACERS=ESMF_FieldCreate(grid            =GRID              &  !<-- The ESMF grid
                                    ,farray          =int_state%TRACERS &  !<-- Insert this pointer into the Field
                                    ,maxHaloUWidth   =(/IHALO,JHALO/)   &  !<-- Upper bound of halo region
                                    ,maxHaloLWidth   =(/IHALO,JHALO/)   &  !<-- Lower bound of halo region
                                    ,ungriddedLBound =(/1,1/)           &
                                    ,ungriddedUBound =(/LM,NUM_TRAC/)   &
                                    ,name            =FIELD_NAME        &  !<-- Name of the 2D real array
                                    ,indexFlag       =ESMF_INDEX_DELOCAL&
                                    ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END SUBROUTINE ALLOC_FIELDS_DYN
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYNAMICS_FIELDS
!
!-----------------------------------------------------------------------
