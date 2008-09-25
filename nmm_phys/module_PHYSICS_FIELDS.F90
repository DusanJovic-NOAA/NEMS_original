!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_FIELDS
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE CREATES ESMF FIELDS FOR THE PHYSICS
!***  IMPORT/EXPORT STATES.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY : INTERNAL_STATE
      USE MODULE_DM_PARALLEL           ,ONLY : IDS,IDE,JDS,JDE          &
                                              ,IMS,IME,JMS,JME          &
                                              ,ITS,ITE,JTS,JTE          &
                                              ,MYPE_SHARE
!
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
      PUBLIC :: ALLOC_FIELDS_PHY
!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!***  LIST AS FIELDS THE ARRAYS THAT ARE TO BE PART OF
!***  THE PHYSICS IMPORT/EXPORT STATES.
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
      SUBROUTINE ALLOC_FIELDS_PHY(GRID,INT_STATE)
!
!-----------------------------------------------------------------------
!***  CREATE FIELDS FOR THE PHYSICS IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Grid),INTENT(IN)   :: GRID                                !<-- The ESMF grid
!
      TYPE(INTERNAL_STATE),POINTER :: INT_STATE                           !<-- The Physics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER             :: I,ISTATUS,J,L,N,RC,RC_FINAL
!
      CHARACTER(20)       :: ARRAY_NAME
!
      TYPE(ESMF_DistGrid) :: DISTGRID
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DISTRIBUTED GRID INFORMATION.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridGet(grid    =GRID                                   &  !<-- The ESMF Grid
                       ,distgrid=DISTGRID                               &  !<-- Information on the distributed ESMF Grid
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
      MESSAGE_CHECK="Create PHY ARRAY_T"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='T'
!
      ARRAY_T=ESMF_ArrayCreate(farray  =int_state%T                     &  !<-- The F90 input array
                              ,distgrid=DISTGRID                        &  !<-- ESMF distributed grid information
                              ,name    =ARRAY_NAME                      &  !<-- ESMF Array name
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
      MESSAGE_CHECK="Create PHY ARRAY_U"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='U'
!
      ARRAY_U=ESMF_ArrayCreate(farray  =int_state%U                     &  !<-- The F90 input array
                              ,distgrid=DISTGRID                        &  !<-- ESMF distributed grid information
                              ,name    =ARRAY_NAME                      &  !<-- ESMF Array name
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
      MESSAGE_CHECK="Create PHY ARRAY_V"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='V'
!
      ARRAY_V=ESMF_ArrayCreate(farray  =int_state%V                     &  !<-- The F90 input array
                              ,distgrid=DISTGRID                        &  !<-- ESMF distributed grid information
                              ,name    =ARRAY_NAME                      &  !<-- ESMF Array name
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
      MESSAGE_CHECK="Create PHY ARRAY_Q2"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='Q2'
!
      ARRAY_Q2=ESMF_ArrayCreate(farray  =int_state%Q2                   &  !<-- The F90 input array
                               ,distgrid=DISTGRID                       &  !<-- ESMF distributed grid information
                               ,name    =ARRAY_NAME                     &  !<-- ESMF Array name
                               ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - -   OMGALF  - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create PHY ARRAY_OMGALF"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='OMGALF'
!
      ARRAY_OMGALF=ESMF_ArrayCreate(farray  =int_state%OMGALF           &  !<-- The F90 input array
                                   ,distgrid=DISTGRID                   &  !<-- ESMF distributed grid information
                                   ,name    =ARRAY_NAME                 &  !<-- ESMF Array name
                                   ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!***  DO THE 2-D FIELDS.
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - - -  PD  - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create PHY ARRAY_PD"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='PD'
!
      ARRAY_PD=ESMF_ArrayCreate(farray  =int_state%PD                   &  !<-- The F90 input array
                               ,distgrid=DISTGRID                       &  !<-- ESMF distributed grid information
                               ,name    =ARRAY_NAME                     &  !<-- ESMF Array name
                               ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DO THE 4-D TRACERS ARRAY.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - TRACERS - - - - - - - - - - - - - - - - -
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create PHY ARRAY_TRACERS"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ARRAY_NAME='TRACERS'
!
      ARRAY_TRACERS=ESMF_ArrayCreate(farray  =int_state%TRACERS         &  !<-- The F90 input array
                                    ,distgrid=DISTGRID                  &  !<-- ESMF distributed grid information
                                    ,name    =ARRAY_NAME                &  !<-- ESMF Array name
                                    ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ALLOC_FIELDS_PHY
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_FIELDS
!
!-----------------------------------------------------------------------
