!-----------------------------------------------------------------------
!
      MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE COUPLER'S REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  IN MODULE_ATM_GRID_COMP.F.
!
!***  THE COUPLER PROVIDES 2-WAY COUPLING BETWEEN THE DYNAMICS AND
!***  PHYSICS GRIDDED COMPONENTS BY TRANSFERING THEIR EXPORT AND
!***  IMPORT STATES BETWEEN THE TWO.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,MYPE_SHARE
!
      USE MODULE_DYNAMICS_GRID_COMP,ONLY : LM
      USE MODULE_DYN_PHY_CPL_DATA
      USE MODULE_CONTROL           ,ONLY : TIMEF
      USE MODULE_ERR_MSG           ,ONLY : ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DYN_PHY_CPL_REGISTER
!
!-----------------------------------------------------------------------
      INCLUDE 'kind.inc'
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT) :: btim,btim0
      REAL(KIND=KFPT),PUBLIC :: cpl_dyn_phy_tim                         &
                               ,get_fld_tim,add_fld_tim
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_PHY_CPL_REGISTER(CPL_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- Coupler component
!
      INTEGER,INTENT(OUT)              :: RC_REG                          !<-- Return code for register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_REG=ESMF_SUCCESS
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
      CALL ESMF_LogWrite("Set Entry Point for Coupler Initialize"       &
                        ,ESMF_LOG_INFO,RC=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Initialize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,CPL_INITIALIZE                     &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,CPL_RUN                            &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,CPL_FINALIZE                       &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CHECK THE ERROR SIGNAL VARIABLE.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)" COUPLER_REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" COUPLER_REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_PHY_CPL_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_INITIALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK      &
                               ,RC_CPL)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SET UP THE COUPLER.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,  INTENT(OUT)   :: RC_CPL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER          :: RC,RC_FINAL
      TYPE(ESMF_VM)    :: VM
      TYPE(ESMF_Field) :: SRC_FIELD,DST_FIELD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      cpl_dyn_phy_tim=0.
      get_fld_tim=0.
      add_fld_tim=0.
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
      RC_CPL  =ESMF_SUCCESS
!
!-----------------------------------------------------------------------
      IF(RC_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL INITIALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL INITIALIZE STEP FAILED"
      ENDIF
!
      IF(PRESENT(RC_CPL))THEN
        RC_CPL=RC_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_RUN(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_CPL)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  RUN THE COUPLER TO TRANSFER DATA BETWEEN THE GRIDDED COMPONENTS.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
      
!
      INTEGER,OPTIONAL,INTENT(OUT)     :: RC_CPL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: L,N,NDATA2,NDATA3,RC,RC_FINAL
      INTEGER :: NDATA1I,NDATA2I,NDATA3I,NDATA1O,NDATA2O,NDATA3O
      INTEGER :: IMP_ITEM,EXP_ITEM
!
      CHARACTER(ESMF_MAXSTR) :: IMPORT_STATENAME,EXPORT_STATENAME
!
      TYPE(ESMF_Array)       :: HOLD_ARRAY
!
      CHARACTER(20)     :: ARRAY_NAME,IMP_ITEM_NAME(20), EXP_ITEM_NAME(20)
!
      LOGICAL,SAVE :: POINT_PHY_AT_DYN=.FALSE.                           &  !<-- Has Physics Import pointed to Dynamics Export yet?
                     ,POINT_DYN_AT_PHY=.FALSE.                              !<-- Has Dynamics Import pointed to Physics Export yet?
!
      LOGICAL,SAVE :: FROM_EXP_DYN_TO_IMP_PHY=.FALSE.                    &
                     ,FROM_EXP_PHY_TO_IMP_DYN=.FALSE.                    &
                     ,FROM_IMP_DYN_TO_EXP_DYN=.FALSE.                    &
                     ,FROM_IMP_PHY_TO_EXP_PHY=.FALSE.         
!
      LOGICAL,SAVE :: CHECK
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  IF THE STATE Field POINTERS HAVE NOT ALREADY BEEN SET, 
!***  THEN PROCEED WITH THE DIRECTING OF THE Field POINTERS. 
!***  AFTER THEY HAVE BEEN SET ONCE, THEY DO NOT NEED TO BE
!***  SET AGAIN.
!-----------------------------------------------------------------------
!
      btim0=timef()
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
      RC_CPL  =ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  DETERMINE THE DIRECTION OF THE TRANSFER BY EXTRACTING
!***  THE STATENAME FROM THE IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve State Name in Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state=IMP_STATE                                &  !<-- The Dyn-Phy Coupler's Import State
                        ,name =IMPORT_STATENAME                         &  !<-- The Import State's Name
                        ,rc   =RC)
!
      CALL ESMF_StateGet(state=EXP_STATE                                &  !<-- The Dyn-Phy Coupler's Imp
                        ,name =EXPORT_STATENAME                         &  !<-- The Import State's Name
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Physics Import' .AND.                 &
         FROM_EXP_DYN_TO_IMP_PHY ) RETURN
!
      IF(TRIM(IMPORT_STATENAME)=='Physics Export' .AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import'.AND.                 &
         FROM_EXP_PHY_TO_IMP_DYN ) RETURN
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Import'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Export'.AND.                 &
         FROM_IMP_DYN_TO_EXP_DYN ) RETURN
!
      IF(TRIM(IMPORT_STATENAME)=='Physics Import'.AND.                  &
         TRIM(EXPORT_STATENAME)=='Physics Export'.AND.                  &
         FROM_IMP_PHY_TO_EXP_PHY ) RETURN
!
      CALL ESMF_StateGet(state       =IMP_STATE                         &
                        ,itemcount   =IMP_ITEM                          &
                        ,itemnamelist=IMP_ITEM_NAME                     &
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve State Name in Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state       =EXP_STATE                         &
                        ,itemcount   =EXP_ITEM                          &
                        ,itemnamelist=EXP_ITEM_NAME                     &
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      RC_FINAL=RC
!
!-----------------------------------------------------------------------
!***  THE NUMBER OF FIELDS TRANSFERRED FROM THE DYNAMICS TO
!***  THE PHYSICS MAY NOT EQUAL THE NUMBER OF FIELDS TRANSFERRED
!***  IN THE OTHER DIRECTION.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!------------------------------------
!***  Get From Dynamics Import State
!------------------------------------
!
      IF(TRIM(IMPORT_STATENAME)=="Dynamics Import")THEN
        NDATA3I=NDATA_3D_PHY_TO_DYN
        NDATA2I=NDATA_2D_PHY_TO_DYN
!
!------------------------------------
!***  Get From Dynamics Export State
!------------------------------------
!
      ELSEIF(TRIM(IMPORT_STATENAME)=="Dynamics Export")THEN
        NDATA3I=NDATA_3D_DYN_TO_PHY
        NDATA2I=NDATA_2D_DYN_TO_PHY
!
!-----------------------------------
!***  Get From Physics Import State
!-----------------------------------
!
      ELSEIF(TRIM(IMPORT_STATENAME)=="Physics Import")THEN
        NDATA3I=NDATA_3D_DYN_TO_PHY
        NDATA2I=NDATA_2D_DYN_TO_PHY
!
!-----------------------------------
!***  Get From Physics Export State
!-----------------------------------
!
      ELSEIF(TRIM(IMPORT_STATENAME)=="Physics Export")THEN
        NDATA3I=NDATA_3D_PHY_TO_DYN
        NDATA2I=NDATA_2D_PHY_TO_DYN
      ENDIF

!-----------------------------------------------------------------------
!
!------------------------------------
!***  Put Into Dynamics Import State
!------------------------------------
!
      IF(TRIM(EXPORT_STATENAME)=="Dynamics Import")THEN
        NDATA3O=NDATA_3D_PHY_TO_DYN
        NDATA2O=NDATA_2D_PHY_TO_DYN
!
!------------------------------------
!***  Put Into Dynamics Export State
!------------------------------------
!
      ELSEIF(TRIM(EXPORT_STATENAME)=="Dynamics Export")THEN
        NDATA3O=NDATA_3D_DYN_TO_PHY
        NDATA2O=NDATA_2D_DYN_TO_PHY
!
!-----------------------------------
!***  Put Into Physics Import State
!-----------------------------------
!
      ELSEIF(TRIM(EXPORT_STATENAME)=="Physics Import")THEN
        NDATA3O=NDATA_3D_DYN_TO_PHY
        NDATA2O=NDATA_2D_DYN_TO_PHY
!
!-----------------------------------
!***  Put Into Physics Export State
!-----------------------------------
!
      ELSEIF(TRIM(EXPORT_STATENAME)=="Physics Export")THEN
        NDATA3O=NDATA_3D_PHY_TO_DYN
        NDATA2O=NDATA_2D_PHY_TO_DYN
!
      ELSE
        WRITE(6,*)' Error: No state name match, state_name='            &
                 ,TRIM(EXPORT_STATENAME)
        WRITE(0,*)' Error: No state name match, state_name='            &
                 ,TRIM(EXPORT_STATENAME)
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(NDATA2O>NDATA2I .OR. NDATA3O>NDATA3I)THEN
        WRITE(6,*)' ERROR: Import data is too few for export data.'
        WRITE(0,*)' ERROR: Import data is too few for export data.'
        CALL ABORT
      ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  LOOP THROUGH THE DATA NAMES, EXTRACT THOSE Arrays FROM THE
!***  IMPORT STATE, AND ADD THEM TO THE EXPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  BEGIN WITH THE 3-D FIELDS.
!-----------------------------------------------------------------------
!
      data_3D: DO N=1,NDATA3I
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        ARRAY_NAME=TRIM(DATANAMES_3D(N))
!       write(0,*)' ARRAY_NAME=', ARRAY_NAME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract 3-D Array from Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<--- Import State that holds the Fields
                          ,itemName=ARRAY_NAME                          &  !<--- Extract Array with this name
                          ,array   =HOLD_ARRAY                          &  !<--- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        get_fld_tim=get_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert 3-D Array into Dyn-Phy Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &  !<--- Insert Array into this Export State
                          ,array=HOLD_ARRAY                             &  !<--- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        add_fld_tim=add_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDDO data_3D
!
!-----------------------------------------------------------------------
!***  NOW TRANSFER THE 2-D FIELDS.
!-----------------------------------------------------------------------
!
      data_2D: DO N=1,NDATA2I
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        ARRAY_NAME=TRIM(DATANAMES_2D(N))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract 2-D Array from Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<--- Import State that holds the Arrays
                          ,itemName=ARRAY_NAME                          &  !<--- Extract Array with this name
                          ,array   =HOLD_ARRAY                          &  !<--- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        get_fld_tim=get_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert 2-D Array into Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &  !<--- Insert Array into this Export State
                          ,array=HOLD_ARRAY                             &  !<--- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        add_fld_tim=add_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDDO data_2D
!
!-----------------------------------------------------------------------
!***  TRANSFER THE 4-D TRACERS ARRAY.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      ARRAY_NAME='TRACERS'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Tracers Array from Dyn-Phy Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<--- The Coupler's Import State
                        ,itemName=ARRAY_NAME                            &  !<--- Extract Array with this name
                        ,array   =HOLD_ARRAY                            &  !<--- Put the extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      get_fld_tim=get_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert TRACERS Array into Dyn-Phy Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE                                &  !<--- Insert Array into the Coupler's Export State
                        ,array=HOLD_ARRAY                               &  !<--- The Array to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      add_fld_tim=add_fld_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  AFTER THE ESMF Arrays IN THE EXPORT STATES ARE POINTED AT THE
!***  CORRECT DATA IN THE IMPORT STATES, THEY REMAIN PROPERLY POINTED
!***  THUS THE StateGet/StateAdd PROCEDURE DOES NOT NEED TO BE REPEATED.
!***  SIGNAL THAT FACT WITH THE FOLLOWING LOGICAL FLAGS.
!-----------------------------------------------------------------------
!
     IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                  &
        TRIM(EXPORT_STATENAME)=='Physics Import')                       &
       FROM_EXP_DYN_TO_IMP_PHY = .TRUE.

     IF(TRIM(IMPORT_STATENAME)=='Physics Export'.AND.                   &
        TRIM(EXPORT_STATENAME)=='Dynamics Import')                      &
       FROM_EXP_PHY_TO_IMP_DYN = .TRUE.
!
     IF(TRIM(IMPORT_STATENAME)=='Dynamics Import'.AND.                  &
        TRIM(EXPORT_STATENAME)=='Dynamics Export')                      &
       FROM_IMP_DYN_TO_EXP_DYN = .TRUE.
!
     IF(TRIM(IMPORT_STATENAME)=='Physics Import'.AND.                   &
        TRIM(EXPORT_STATENAME)=='Physics Export')                       &
       FROM_IMP_PHY_TO_EXP_PHY = .TRUE.

!-----------------------------------------------------------------------
!
      RC_FINAL=RC
!
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL RUN STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL RUN STEP FAILED"
      ENDIF
!
      IF(PRESENT(RC_CPL))THEN
        RC_CPL=RC_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      cpl_dyn_phy_tim=cpl_dyn_phy_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_FINALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_CPL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE COUPLER.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                         !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                        !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                        !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                            !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,   INTENT(OUT)  :: RC_CPL
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC,RC_FINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
      RC_CPL  =ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL FINALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL FINALIZE STEP FAILED"
      ENDIF
!
      IF(PRESENT(RC_CPL))THEN
        RC_CPL=RC_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_FINALIZE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
