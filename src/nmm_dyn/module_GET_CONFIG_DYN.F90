!-----------------------------------------------------------------------
!
      MODULE MODULE_GET_CONFIG_DYN
!
!-----------------------------------------------------------------------
!
!***  EXTRACT DATA FROM THE ESMF CONFIGURATION FILES
!***  AND LOAD IT INTO THE INTERNAL STATES.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
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
      PUBLIC :: GET_CONFIG_DYN
!
!-----------------------------------------------------------------------
!
      INTERFACE GET_CONFIG_DYN
        MODULE PROCEDURE GET_CONFIG_DYN
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_CONFIG_DYN(GRID_COMP,INT_STATE,RC_CONF)
!
!-----------------------------------------------------------------------
      USE MODULE_DYNAMICS_INTERNAL_STATE
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp), INTENT(INOUT) :: GRID_COMP                     !<-- The Dynamics gridded component
      TYPE(INTERNAL_STATE),INTENT(INOUT) :: INT_STATE                     !<-- The Dynamics internal state
      INTEGER             ,INTENT(OUT)   :: RC_CONF                       !<-- Final return code
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config) :: CF
      INTEGER           :: RC
      CHARACTER(50)     :: MODE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_CONF=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ESMF CONFIG OBJECT CF FROM THE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Retrieve Config Object from Dynamics Component"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &   !<--- The Dynamics gridded component
                           ,config  =CF                                 &   !<--- The configure (namelist) object
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  BEGIN EXTRACTION OF THE CONFIGURATION FILE DATA.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract IM from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%IM                  &  !<-- Put extracted quantity here
                                  ,label ='im:'                         &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract JM from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%JM                  &  !<-- Put extracted quantity here
                                  ,label ='jm:'                         &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract LM from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%LM                  &  !<-- Put extracted quantity here
                                  ,label ='lm:'                         &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract TPH0D from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%TPH0D               &  !<-- Put extracted quantity here
                                  ,label ='tph0d:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract TLM0D from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%TLM0D               &  !<-- Put extracted quantity here
                                  ,label ='tlm0d:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_YEAR from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_YEAR          &  !<-- Put extracted quantity here
                                  ,label ='start_year:'                 &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_MONTH from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_MONTH         &  !<-- Put extracted quantity here
                                  ,label ='start_month:'                &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_DAY from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_DAY           &  !<-- Put extracted quantity here
                                  ,label ='start_day:'                  &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_HOUR from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_HOUR          &  !<-- Put extracted quantity here
                                  ,label ='start_hour:'                 &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_MINUTE from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_MINUTE        &  !<-- Put extracted quantity here
                                  ,label ='start_minute:'               &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_SECOND from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%START_SECOND        &  !<-- Put extracted quantity here
                                  ,label ='start_second:'               &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract START_TSTART from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%TSTART              &  !<-- Put extracted quantity here
                                  ,label ='tstart:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract FILTER METHOD from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =int_state%FILTER_METHOD       &
                                  ,label ='filter_method:'              &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract RESTART from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =MODE                          &  !<-- Put extracted quantity here
                                  ,label ='restart:'                    &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      int_state%RESTART=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%RESTART=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract GLOBAL from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =MODE                          &  !<-- Put extracted quantity here
                                  ,label ='global:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      int_state%GLOBAL=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%GLOBAL=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract HYDRO from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =MODE                          &  !<-- Put extracted quantity here
                                  ,label ='hydro:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      int_state%HYDRO=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%HYDRO=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract WBD from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%WBD                 &  !<-- Put extracted quantity here
                                  ,label ='wbd:'                        &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SBD from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SBD                 &  !<-- Put extracted quantity here
                                  ,label ='sbd:'                        &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract NHOURS_FCST from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NHOURS_FCST         &  !<-- Put extracted quantity here
                                  ,label ='nhours_fcst:'                &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract NHOURS_HISTORY from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NHOURS_HISTORY      &  !<-- Put extracted quantity here
                                  ,label ='nhours_history:'             &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract NHOURS_RESTART from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NHOURS_RESTART      &  !<-- Put extracted quantity here
                                  ,label ='nhours_restart:'             &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SECADV from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SECADV              &  !<-- Put extracted quantity here
                                  ,label ='secadv:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SECDIF from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SECDIF              &  !<-- Put extracted quantity here
                                  ,label ='secdif:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SMAG2 from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SMAG2               &  !<-- Put extracted quantity here
                                  ,label ='smag2:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SMAG4 from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SMAG4               &  !<-- Put extracted quantity here
                                  ,label ='smag4:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract CODAMP from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%CODAMP              &  !<-- Put extracted quantity here
                                  ,label ='codamp:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract WCOR from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%WCOR                &  !<-- Put extracted quantity here
                                  ,label ='wcor:'                       &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract PWRC from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%PWRC                &  !<-- Put extracted quantity here
                                  ,label ='pwrc:'                       &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract IDTAD from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%IDTAD               &  !<-- Put extracted quantity here
                                  ,label ='idtad:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract ADVECT_TRACERS from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%ADVECT_TRACERS      &  !<-- Put extracted quantity here
                                  ,label ='advect_tracers:'             &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract IDTADT from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%IDTADT              &  !<-- Put extracted quantity here
                                  ,label ='idtadt:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract NUM_TRACERS_MET from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NUM_TRACERS_MET     &  !<-- Put extracted quantity here
                                  ,label ='num_tracers_met:'            &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract NUM_TRACERS_CHEM from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NUM_TRACERS_CHEM    &  !<-- Put extracted quantity here
                                  ,label ='num_tracers_chem:'           &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract INPES from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%INPES               &  !<-- Put extracted quantity here
                                  ,label ='inpes:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract JNPES from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%JNPES               &  !<-- Put extracted quantity here
                                  ,label ='jnpes:'                      &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract WRITE_TASKS_PER_GROUP from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                              &  !<-- The configure file object
                                  ,value =int_state%WRITE_TASKS_PER_GROUP &  !<-- Put extracted quantity here
                                  ,label ='write_tasks_per_group:'        &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract WRITE_GROUPS from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%WRITE_GROUPS        &  !<-- Put extracted quantity here
                                  ,label ='write_groups:'               &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract READ_GLOBAL_SUMS from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =MODE                          &  !<-- Put extracted quantity here
                                  ,label ='read_global_sums:'           &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      int_state%READ_GLOBAL_SUMS=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%READ_GLOBAL_SUMS=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract WRITE_GLOBAL_SUMS from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =MODE                          &  !<-- Put extracted quantity here
                                  ,label ='write_global_sums:'          &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      int_state%WRITE_GLOBAL_SUMS=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%WRITE_GLOBAL_SUMS=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SHORTWAVE from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SHORTWAVE           &  !<-- Put extracted quantity here
                                  ,label ='shortwave:'                  &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract LONGWAVE from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%LONGWAVE            &  !<-- Put extracted quantity here
                                  ,label ='longwave:'                   &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract SFC_LAYER from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%SFC_LAYER           &  !<-- Put extracted quantity here
                                  ,label ='sfc_layer:'                  &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract TURBULENCE from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%TURBULENCE          &  !<-- Put extracted quantity here
                                  ,label ='turbulence:'                 &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract CONVECTION from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%CONVECTION          &  !<-- Put extracted quantity here
                                  ,label ='convection:'                 &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract MICROPHYSICS from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%MICROPHYSICS        &  !<-- Put extracted quantity here
                                  ,label ='microphysics:'               &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="GET_CONFIG_DYN: Extract ADIABATIC from Config File"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
      int_state%ADIABATIC=.FALSE.
      IF(TRIM(MODE)=='true')THEN
        int_state%ADIABATIC=.TRUE.
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CONF)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_CONF==ESMF_SUCCESS)THEN
!!!     WRITE(0,*)'GET_CONFIG PASSED FOR DYNAMICS'
      ELSE
        WRITE(0,*)'GET_CONFIG FAILED FOR DYNAMICS RC_CONF=',RC_CONF
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GET_CONFIG_DYN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GET_CONFIG_DYN
!
!-----------------------------------------------------------------------
