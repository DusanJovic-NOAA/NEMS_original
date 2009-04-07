!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_FIELDS
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE CREATES/DESTROYS ESMF Arrays FOR THE DYNAMICS
!***  IMPORT/EXPORT STATES.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_DYNAMICS_INTERNAL_STATE,ONLY : INTERNAL_STATE
      USE MODULE_DM_PARALLEL            ,ONLY : IDS,IDE,JDS,JDE         &
                                               ,IMS,IME,JMS,JME         &
                                               ,ITS,ITE,JTS,JTE         &
                                               ,MYPE_SHARE
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
      PUBLIC :: ARRAY_T                                                 &
               ,ARRAY_U                                                 &
               ,ARRAY_V                                                 &
               ,ARRAY_Q2                                                &
               ,ARRAY_OMGALF
!
      PUBLIC :: ARRAY_PD
!
      PUBLIC :: ARRAY_TRACERS
!
      PUBLIC :: ALLOC_FIELDS_DYN
!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!***  LIST ESMF Arrays THAT ARE TO BE PART OF
!***  THE DYNAMICS IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Array),SAVE :: ARRAY_T                                  &
                              ,ARRAY_U                                  &
                              ,ARRAY_V                                  &
                              ,ARRAY_Q2                                 &
                              ,ARRAY_OMGALF
!
      TYPE(ESMF_Array),SAVE :: ARRAY_PD
!
      TYPE(ESMF_Array),SAVE :: ARRAY_TRACERS
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
!***  CREATE ESMF Arrays FOR THE DYNAMICS IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_DistGrid)          :: DISTGRID                           !<-- Information of the distributed ESMF Grid
!
      TYPE(ESMF_Grid),INTENT(IN)   :: GRID                               !<-- The ESMF grid
!
      TYPE(INTERNAL_STATE),POINTER :: INT_STATE                          !<-- The Dynamics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER       :: I,ISTATUS,J,L,N,RC,RC_FINAL
!
      CHARACTER(20) :: ARRAY_NAME
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DISTRIBUTED GRID INFORMATION.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridGet(grid    =GRID                                   &  !<-- The ESMF Grid
                       ,distgrid=DISTGRID                               &  !<-- ESMF distributed grid information
                       ,rc      =RC)
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF Arrays THAT WILL BE ADDED TO THE IMPORT/EXPORT
!***  STATES AND ASSOCIATE THE APPROPRIATE POINTERS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - -   T   - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create ARRAY_T"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='T'
!
      ARRAY_T=ESMF_ArrayCreate(farray  =int_state%T                     &  !<-- The F90 input array
                              ,distgrid=DISTGRID                        &  !<-- ESMF distributed grid information 
                              ,name    =ARRAY_NAME                      &  !<-- ESMF Array name
                              ,indexFlag=ESMF_INDEX_DELOCAL             &
                              ,rc      =RC)
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
      MESSAGE_CHECK="Create ARRAY_U"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='U'
!
      ARRAY_U=ESMF_ArrayCreate(farray  =int_state%U                     &  !<-- The F90 input array
                              ,distgrid=DISTGRID                        &  !<-- ESMF distributed grid information
                              ,name    =ARRAY_NAME                      &  !<-- ESMF Array name
                              ,indexFlag=ESMF_INDEX_DELOCAL             &
                              ,rc      =RC)
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
      MESSAGE_CHECK="Create ARRAY_V"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='V'
!
      ARRAY_V=ESMF_ArrayCreate(farray    =int_state%V                   &  !<-- The F90 input array
                                ,distgrid=distgrid                      &  !<-- ESMF distributed grid information
                                ,name    =ARRAY_NAME                    &  !<-- ESMF Array name
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc      =RC)
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
      MESSAGE_CHECK="Create ARRAY_Q2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='Q2'
!
      ARRAY_Q2=ESMF_ArrayCreate(farray  =int_state%Q2                   &  !<-- The F90 input array
                               ,distgrid=distgrid                       &  !<-- ESMF distributed grid information
                               ,name    =ARRAY_NAME                     &  !<-- ESMF Array name
                               ,indexFlag=ESMF_INDEX_DELOCAL            &
                               ,rc      =RC)
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
      MESSAGE_CHECK="Create ARRAY_OMGALF"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='OMGALF'
!
      ARRAY_OMGALF=ESMF_ArrayCreate(farray  =int_state%OMGALF           &  !<-- The F90 input array
                                   ,distgrid=distgrid                   &  !<-- ESMF distributed grid information
                                   ,name    =ARRAY_NAME                 &  !<-- ESMF Array name
                                   ,indexFlag=ESMF_INDEX_DELOCAL        &
                                   ,rc      =RC)
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
      MESSAGE_CHECK="Create ARRAY_PD"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='PD'
!
      ARRAY_PD=ESMF_ArrayCreate(farray  =int_state%PD                   &  !<-- The F90 input array
                               ,distgrid=distgrid                       &  !<-- ESMF distributed grid information
                               ,name    =ARRAY_NAME                     &  !<-- ESMF Array name
                               ,indexFlag=ESMF_INDEX_DELOCAL            &
                               ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE AND FILL THE Arrays FOR THE 4-D TRACERS ARRAY.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - -  TRACERS - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create ARRAY_TRACERS"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='TRACERS'
!
      ARRAY_TRACERS=ESMF_ArrayCreate(farray  =int_state%TRACERS           &  !<-- The F90 input array
                                    ,distgrid=distgrid                    &  !<-- ESMF distributed grid information
                                    ,name    =ARRAY_NAME                  &  !<-- ESMF Array name
                                    ,indexFlag=ESMF_INDEX_DELOCAL         &
                                    ,rc      =RC)
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
