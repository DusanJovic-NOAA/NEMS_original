#include "../../ESMFVersionDefine.h"

!  2011-05-12  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!--------------------------------------------------------------------------

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
!***  IMPORT STATES BETWEEN THE TWO.  IT IS ALSO USED TO TRANSFER
!***  THE DYNAMICS EXPORT STATE TO THE DYNAMICS IMPORT STATE WHEN
!***  FORECASTS WITH NO PHYSICS ARE MADE.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,MYPE_SHARE
!
      USE MODULE_DYN_PHY_CPL_DATA,ONLY : DATANAMES_2D                   &
                                        ,DATANAMES_3D                   &
                                        ,NDATA_2D_FROM_DYN              &
                                        ,NDATA_2D_FROM_PHY              &
                                        ,NDATA_3D_FROM_DYN              &
                                        ,NDATA_3D_FROM_PHY  
!
      USE MODULE_CONTROL,ONLY : TIMEF
      USE MODULE_CLOCKTIMES,ONLY : add_fld_tim,cpl_dyn_phy_tim,get_fld_tim
      USE MODULE_ERR_MSG,ONLY : ERR_MSG,MESSAGE_CHECK
      USE MODULE_INCLUDE
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
      REAL(KIND=KDBL) :: btim,btim0
!
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
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                      !<-- Coupler component
!
      INTEGER(kind=KINT),INTENT(OUT) :: RC_REG                            !<-- Return code for register
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_REG=ESMF_SUCCESS
                                                                                                                                              
!-----------------------------------------------------------------------
!***  Register the Coupler Initialize subroutine.  Since it is just one
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN,
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,CPL_INITIALIZE                     &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,CPL_INITIALIZE                     &  !<-- User's subroutineName
                                    ,phase=ESMF_SINGLEPHASE             &  !<-- Phase
                                    ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the coupler Run subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,CPL_RUN                            &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,CPL_RUN                            &  !<-- User's subroutineName
                                    ,phase=ESMF_SINGLEPHASE             &  !<-- Phase
                                    ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Coupler Finalize subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Coupler Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,CPL_FINALIZE                       &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Dyn-Phy Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,CPL_FINALIZE                       &  !<-- User's subroutineName
                                    ,phase=ESMF_SINGLEPHASE             &  !<-- Phase
                                    ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Check the error signal variable.
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
      SUBROUTINE CPL_INITIALIZE(CPL_COMP                                &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK                                   &
                               ,RC_CPL)
!-----------------------------------------------------------------------
!***  Set up the Dynamics-Physics coupler.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Dyn-Phy Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_CPL
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC, RC_FINAL, I, itemCount
!
      CHARACTER(ESMF_MAXSTR) :: IMPORT_STATENAME,EXPORT_STATENAME
      CHARACTER(LEN=10), DIMENSION(100) :: itemNameList
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_StateItemType) :: stateItemType

      LOGICAL,SAVE :: FROM_EXP_DYN_TO_IMP_PHY=.FALSE.                   &
                     ,FROM_EXP_PHY_TO_IMP_DYN=.FALSE.                   &
                     ,FROM_IMP_DYN_TO_EXP_DYN=.FALSE.                   &
                     ,FROM_IMP_PHY_TO_EXP_PHY=.FALSE.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
      RC_CPL  =ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet(state=IMP_STATE ,name =IMPORT_STATENAME ,rc   =RC)
      CALL ESMF_StateGet(state=EXP_STATE ,name =EXPORT_STATENAME ,rc   =RC)

      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Physics Import' .AND.                 &
         FROM_EXP_DYN_TO_IMP_PHY ) RETURN

      IF(TRIM(IMPORT_STATENAME)=='Physics Export' .AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import'.AND.                 &
         FROM_EXP_PHY_TO_IMP_DYN ) RETURN
!
!-----------------------------------------------------------------------
!***  Move all Field pointers from the Import State to the Export State.
!
!***  NOTE:  This is a fundamental step in handling the ownership of
!***         variables between Dynamics and Physics.  In the following
!***         DO loop those variables that are exported and thus owned 
!***         by a given component have their allocated pointers taken
!***         from the export state and moved to the other component's
!***         import state.  Then in the 2nd phases os the Dynamics/
!***         Physics Init steps those allocated pointers are unloaded
!***         from the import states and the unowned variables are
!***         pointed at them.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet(state       =IMP_STATE                         &
                        ,itemCount   =itemCount                         &
                        ,itemNameList=itemNameList                      &
                        ,rc          =RC)
!
      DO i=1,itemCount  
        CALL ESMF_StateGet(              IMP_STATE                      &
                          ,              itemNameList(i)                &
                          ,              stateItemType                  &
                          ,rc           =RC)
!
        IF (stateItemType == ESMF_STATEITEM_FIELD) THEN
          if(mype_share==0)then
            write(0,*) TRIM(IMPORT_STATENAME),' -> ',TRIM(EXPORT_STATENAME) &
                                             ,' Field     : ',i,itemNameList(i)
          endif
!
          CALL ESMF_StateGet(state   =IMP_STATE                         &
                            ,itemName=itemNameList(i)                   &
                            ,field   =HOLD_FIELD                        &
                            ,rc      =RC)
!
          CALL ESMF_StateAdd(state=EXP_STATE                            &
                            ,field=HOLD_FIELD                           &
                            ,rc   =RC)
        END IF
      END DO

! DYN -> PHY
     IF(TRIM(IMPORT_STATENAME)=='Dynamics Export' .AND.                 &
        TRIM(EXPORT_STATENAME)=='Physics Import') THEN
       FROM_EXP_DYN_TO_IMP_PHY = .TRUE.
     ENDIF

! PHY -> DYN
      IF(TRIM(IMPORT_STATENAME)=='Physics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import') THEN
       FROM_EXP_PHY_TO_IMP_DYN = .TRUE.
     ENDIF

      ! Call CPL_RUN to exchange State Attributes
!d      CALL CPL_RUN(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,RC_CPL)


      IF(RC_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL INITIALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL INITIALIZE STEP FAILED"
      ENDIF
!
      RC_CPL=RC_FINAL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_RUN(CPL_COMP                                       &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK                                          &
                        ,RC_CPL)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Run the coupler to transfer data between export and import states
!***  of the Dynamics and Physics.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Dyn-Phy Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_CPL
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: L,N,NDATA_2D,NDATA_3D,NUM_TRACERS_TOTAL     &
                           ,RC,RC_FINAL
!
      CHARACTER(ESMF_MAXSTR) :: IMPORT_STATENAME,EXPORT_STATENAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      CHARACTER(20) :: FIELD_NAME
!
      LOGICAL(kind=KLOG),SAVE :: FROM_DYN_EXP_TO_PHY_IMP=.FALSE.        &  !<-- Has Physics Import pointed to Dynamics Export?
                                ,FROM_PHY_EXP_TO_DYN_IMP=.FALSE.        &  !<-- Has Dynamics Import pointed to Physics Export?
                                ,FROM_DYN_EXP_TO_DYN_IMP=.FALSE.           !<-- Has Dynamics Import pointed to Dynamics Export?
!
      LOGICAL(kind=KLOG),SAVE :: CHECK
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  If the State Field pointers have not already been set, 
!***  then proceed with directing them.  After they have been
!***  set once for a given combination they do not need to be
!***  set again.
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
      RC_CPL  =ESMF_SUCCESS
#if 0
!
!-----------------------------------------------------------------------
!***  Determine the direction of the transfer and the states involved
!***  by extracting the State name from the import and export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve State Names in Dyn-Phy Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state=IMP_STATE                                &  !<-- The Dyn-Phy Coupler's import State
                        ,name =IMPORT_STATENAME                         &  !<-- The import state's Name
                        ,rc   =RC)
!
      CALL ESMF_StateGet(state=EXP_STATE                                &  !<-- The Dyn-Phy Coupler's export state
                        ,name =EXPORT_STATENAME                         &  !<-- The export State's Name
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Physics Import' .AND.                 &
         FROM_DYN_EXP_TO_PHY_IMP ) RETURN
!
      IF(TRIM(IMPORT_STATENAME)=='Physics Export' .AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import'.AND.                 &
         FROM_PHY_EXP_TO_DYN_IMP ) RETURN
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import'.AND.                 &
         FROM_DYN_EXP_TO_DYN_IMP ) RETURN
!
!-----------------------------------------------------------------------
!***  What is the number of Fields to be transferred from 
!***  the import state?
!-----------------------------------------------------------------------
!
!-------------------------------------------------
!***  Number of Fields from Dynamics export state
!-------------------------------------------------
!
!!!   IF(TRIM(EXPORT_STATENAME)=="Dynamics Export")THEN
      IF(TRIM(IMPORT_STATENAME)=="Dynamics Export")THEN
        NDATA_3D=NDATA_3D_FROM_DYN
        NDATA_2D=NDATA_2D_FROM_DYN
      ENDIF
!
!------------------------------------------------
!***  Number of Fields from Physics export state
!------------------------------------------------
!
!!!   IF(TRIM(EXPORT_STATENAME)=="Physics Export")THEN
      IF(TRIM(IMPORT_STATENAME)=="Physics Export")THEN
        NDATA_3D=NDATA_3D_FROM_PHY
        NDATA_2D=NDATA_2D_FROM_PHY
      ENDIF
!
!-----------------------------------------------------------------------
!***  Loop through the data names, extract the Fields from the
!***  import state, and add them to the export state.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Begin with the 3-D Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      data_3d: DO N=1,NDATA_3D
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        FIELD_NAME=TRIM(DATANAMES_3D(N))
!       write(0,*)' Dyn-Phy Cpl Run FIELD_NAME=', FIELD_NAME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract 3-D Field from Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- Import State that holds the Fields
                          ,itemName=FIELD_NAME                          &  !<-- Extract Field with this name
                          ,field   =HOLD_FIELD                          &  !<-- Put the extracted Field here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        get_fld_tim=get_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert 3-D Field into Dyn-Phy Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &  !<-- Insert Field into this Export State
                          ,field=HOLD_FIELD                             &  !<-- The Field to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        add_fld_tim=add_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDDO data_3d
!
!-----------------------------------------------------------------------
!***  Now transfer the 2-D Fields.
!-----------------------------------------------------------------------
!
      data_2d: DO N=1,NDATA_2D
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        FIELD_NAME=TRIM(DATANAMES_2D(N))
!       write(0,*)' Dyn-Phy Cpl Run FIELD_NAME=', FIELD_NAME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract 2-D Field from Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<--- Import State that holds the Fields
                          ,itemName=FIELD_NAME                          &  !<--- Extract Field with this name
                          ,field   =HOLD_FIELD                          &  !<--- Put the extracted Field here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        get_fld_tim=get_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert 2-D Field into Dyn-Phy Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE                              &  !<--- Insert Field into this Export State
                          ,field=HOLD_FIELD                             &  !<--- The Field to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        add_fld_tim=add_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDDO data_2d
!
!-----------------------------------------------------------------------
!***  Transfer the 4-D Tracers Field.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      FIELD_NAME='TRACERS'
!     write(0,*)' Dyn-Phy Cpl Run FIELD_NAME=', FIELD_NAME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Tracers Field from Dyn-Phy Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<--- The Coupler's Import State
                        ,itemName=FIELD_NAME                            &  !<--- Extract Field with this name
                        ,field   =HOLD_FIELD                            &  !<--- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      get_fld_tim=get_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert TRACERS Field into Dyn-Phy Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE                                &  !<--- Insert Field into the Coupler's Export State
                        ,field=HOLD_FIELD                               &  !<--- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      add_fld_tim=add_fld_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  After the ESMF Fields in the export states are pointed at the
!***  correct data in the import states, they remain properly pointed
!***  thus the StateGet/StateAdd procedure does not need to be repeated.
!***  Signal that fact with the following logical flags.
!-----------------------------------------------------------------------
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Physics Import')                      &
        FROM_DYN_EXP_TO_PHY_IMP = .TRUE.

      IF(TRIM(IMPORT_STATENAME)=='Physics Export'.AND.                  &
         TRIM(EXPORT_STATENAME)=='Dynamics Import')                     &
        FROM_PHY_EXP_TO_DYN_IMP = .TRUE.
!
      IF(TRIM(IMPORT_STATENAME)=='Dynamics Export'.AND.                 &
         TRIM(EXPORT_STATENAME)=='Dynamics Import')                     &
        FROM_DYN_EXP_TO_DYN_IMP = .TRUE.
!
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
#endif
      RC_CPL=RC_FINAL
!
!-----------------------------------------------------------------------
!
      cpl_dyn_phy_tim=cpl_dyn_phy_tim+(timef()-btim0)
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
!***  Finalize the Dynamics-Physics coupler.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Dyn-Phy Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock),INTENT(IN) :: CLOCK                                 !<-- The ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_CPL
!      
!---------------------
!***  Local Variables
!---------------------
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
      RC_CPL=RC_FINAL
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
