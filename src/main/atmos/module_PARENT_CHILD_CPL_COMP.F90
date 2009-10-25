!-----------------------------------------------------------------------
!
      MODULE MODULE_PARENT_CHILD_CPL_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE CONTAINS THE COUPLER THAT EXCHANGES DATA BETWEEN
!***  PARENT DOMAINS AND THEIR CHILDREN.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!
!   2008-06-12  Black - Module created.
!   2009-02-19  Black - Hydrostatic update of nest boundaries.
!
!-----------------------------------------------------------------------
!
! USAGE: 
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_INCLUDE
!
      USE MODULE_CONTROL,ONLY: TIMEF
!
!!!   USE MODULE_CLOCKTIMES
!
      USE MODULE_CONSTANTS,ONLY: P608,R_D
!
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
      PUBLIC :: PARENT_CHILD_CPL_REGISTER                               &
               ,PARENT_CHILD_COUPLER_SETUP
!
!-----------------------------------------------------------------------
!
      TYPE BOUNDARY_SIDES                                                 !<-- Hold the boundary blending region along each side of a domain
        REAL,DIMENSION(:,:,:),POINTER :: SOUTH
        REAL,DIMENSION(:,:,:),POINTER :: NORTH
        REAL,DIMENSION(:,:,:),POINTER :: WEST 
        REAL,DIMENSION(:,:,:),POINTER :: EAST
      END TYPE BOUNDARY_SIDES
!
      TYPE REAL_DATA     
        REAL,DIMENSION(:),POINTER :: DATA
      END TYPE REAL_DATA    
!
      TYPE REAL_DATA_TASKS  
        TYPE(REAL_DATA),DIMENSION(:),POINTER :: TASKS  
      END TYPE REAL_DATA_TASKS   
!
      TYPE INTEGER_DATA_TASKS 
        INTEGER,DIMENSION(:),POINTER :: TASKS        
      END TYPE INTEGER_DATA_TASKS 
!
      TYPE CHILD_TASKS
        INTEGER,DIMENSION(:,:),POINTER :: LIMITS
      END TYPE CHILD_TASKS
!
      TYPE SIDES_0D  
        INTEGER :: SOUTH           
        INTEGER :: NORTH
        INTEGER :: WEST
        INTEGER :: EAST
      END TYPE SIDES_0D
!
      TYPE SIDES_1D_INT  
        INTEGER,DIMENSION(:),POINTER :: SOUTH           
        INTEGER,DIMENSION(:),POINTER :: NORTH
        INTEGER,DIMENSION(:),POINTER :: WEST
        INTEGER,DIMENSION(:),POINTER :: EAST
      END TYPE SIDES_1D_INT
!
      TYPE SAVE_TASK_IJ                                
        INTEGER,DIMENSION(:),POINTER :: I_LO_SOUTH    
        INTEGER,DIMENSION(:),POINTER :: I_HI_SOUTH
        INTEGER,DIMENSION(:),POINTER :: I_HI_SOUTH_TRANSFER
        INTEGER,DIMENSION(:),POINTER :: I_LO_NORTH
        INTEGER,DIMENSION(:),POINTER :: I_HI_NORTH
        INTEGER,DIMENSION(:),POINTER :: I_HI_NORTH_TRANSFER
        INTEGER,DIMENSION(:),POINTER :: J_LO_WEST
        INTEGER,DIMENSION(:),POINTER :: J_HI_WEST
        INTEGER,DIMENSION(:),POINTER :: J_HI_WEST_TRANSFER
        INTEGER,DIMENSION(:),POINTER :: J_LO_EAST
        INTEGER,DIMENSION(:),POINTER :: J_HI_EAST
        INTEGER,DIMENSION(:),POINTER :: J_HI_EAST_TRANSFER
      END TYPE SAVE_TASK_IJ
!
      TYPE HANDLE_SEND         
        INTEGER,DIMENSION(:),POINTER :: NTASKS_TO_RECV                    !<-- Parent MPI handles used when ISend'ing to each child task 
      END TYPE HANDLE_SEND       
!
      TYPE DATA_INFO
        REAL,DIMENSION(:),POINTER :: STRING
        INTEGER                   :: LENGTH
        INTEGER                   :: ID_SOURCE
        INTEGER                   :: INDX_START
        INTEGER                   :: INDX_END
        INTEGER                   :: INDX_END_EXP
      END TYPE DATA_INFO
!
      TYPE PARENT_DATA
        TYPE(DATA_INFO) :: SOUTH_H
        TYPE(DATA_INFO) :: SOUTH_V
        TYPE(DATA_INFO) :: NORTH_H
        TYPE(DATA_INFO) :: NORTH_V
        TYPE(DATA_INFO) :: WEST_H
        TYPE(DATA_INFO) :: WEST_V
        TYPE(DATA_INFO) :: EAST_H
        TYPE(DATA_INFO) :: EAST_V
      END TYPE PARENT_DATA
!
      TYPE PARENT_POINTS_SURROUND_H                                       !<-- Indices of parent points around each child bndry H point
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_NBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_SBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_EBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_WBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_NBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_SBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_EBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_WBND
      END TYPE PARENT_POINTS_SURROUND_H
!
      TYPE PARENT_POINTS_SURROUND_V                                        !<-- Indices of parent points around each child bndry V point
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_NBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_SBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_EBND
        INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_WBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_NBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_SBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_EBND
        INTEGER,DIMENSION(:,:,:),POINTER :: J_INDX_WBND
      END TYPE PARENT_POINTS_SURROUND_V
!
      TYPE PARENT_WEIGHTS_SURROUND_H                                       !<-- Bilinear interpolation weights of the 4 parent points
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_NBND                      !     around each child boundary H point
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_SBND
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_EBND
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_WBND
      END TYPE PARENT_WEIGHTS_SURROUND_H
!
      TYPE PARENT_WEIGHTS_SURROUND_V                                       !<-- Bilinear interpolation weights of the 4 parent points
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_NBND                      !     around each child boundary H point
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_SBND
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_EBND
        REAL,DIMENSION(:,:,:),POINTER :: WEIGHTS_WBND
      END TYPE PARENT_WEIGHTS_SURROUND_V
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: INDX_SW=1                                    &
                          ,INDX_SE=2                                    &
                          ,INDX_NW=3                                    &
                          ,INDX_NE=4
!
      INTEGER,SAVE :: COMM_MY_DOMAIN                                    &
                     ,COMM_TO_MY_PARENT                                 &
                     ,ITS,ITE,JTS,JTE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,IDS,IDE,JDS,JDE                                   &
                     ,INDX_CW,INDX_Q                                    &
                     ,MY_DOMAIN_ID                                      &
                     ,MY_LOCAL_RANK_CHILD                               &
                     ,MY_LOCAL_RANK_PARENT                              &
                     ,N_BLEND_H,N_BLEND_V                               &
                     ,NHALO                                             &
                     ,NUM_CHILDREN
!
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: IM_CHILD                 &
                                              ,JM_CHILD                 &
                                              ,NSTEP_CHILD_RECV         &
!
                                              ,NUM_TASKS_SEND_H_S       &
                                              ,NUM_TASKS_SEND_H_N       &
                                              ,NUM_TASKS_SEND_H_W       &
                                              ,NUM_TASKS_SEND_H_E       &
                                              ,NUM_TASKS_SEND_V_S       &
                                              ,NUM_TASKS_SEND_V_N       &
                                              ,NUM_TASKS_SEND_V_W       &
                                              ,NUM_TASKS_SEND_V_E 
!
      INTEGER,DIMENSION(:),POINTER,SAVE :: COMM_TO_MY_CHILDREN          &
                                          ,FTASKS_DOMAIN                &
                                          ,ID_PARENTS                   &
                                          ,INC_FIX                      &
                                          ,MY_CHILDREN_ID               &
                                          ,PARENT_CHILD_SPACE_RATIO     &
                                          ,PARENT_CHILD_TIME_RATIO 
!
      REAL,SAVE :: PDTOP,PT
!
      REAL,DIMENSION(:),POINTER,SAVE :: PSGML1                          &
                                       ,SG1                             &
                                       ,SG2                             &
                                       ,SGML2
!
      REAL,DIMENSION(:,:),POINTER,SAVE :: FIS,PD
!
      REAL,DIMENSION(:,:,:),POINTER,SAVE :: T,U,V
!
      REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: PDB_S,PDB_N,PDB_W,PDB_E
!
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: TB_S,TB_N,TB_W,TB_E     &
                                               ,QB_S,QB_N,QB_W,QB_E     &
                                               ,CWB_S,CWB_N,CWB_W,CWB_E &
                                               ,UB_S,UB_N,UB_W,UB_E     &
                                               ,VB_S,VB_N,VB_W,VB_E
!
      REAL,DIMENSION(:),POINTER,SAVE :: BOUND_1D_SOUTH_H                &
                                       ,BOUND_1D_SOUTH_V                &
                                       ,BOUND_1D_NORTH_H                &
                                       ,BOUND_1D_NORTH_V                &
                                       ,BOUND_1D_WEST_H                 &
                                       ,BOUND_1D_WEST_V                 &
                                       ,BOUND_1D_EAST_H                 &
                                       ,BOUND_1D_EAST_V   
!
      REAL,DIMENSION(:,:,:,:),POINTER,SAVE :: TRACERS
!
      LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: SEND_CHILD_DATA
!
      TYPE(ESMF_Logical) :: I_AM_A_FCST_TASK,I_AM_A_PARENT
!
      TYPE(REAL_DATA_TASKS),DIMENSION(:),POINTER,SAVE ::                &
                                                 CHILD_BOUND_H_SOUTH    & 
                                                ,CHILD_BOUND_H_NORTH    & 
                                                ,CHILD_BOUND_H_WEST     & 
                                                ,CHILD_BOUND_H_EAST     &
                                                ,CHILD_BOUND_V_SOUTH    &
                                                ,CHILD_BOUND_V_NORTH    &
                                                ,CHILD_BOUND_V_WEST     &
                                                ,CHILD_BOUND_V_EAST     
!
      TYPE(REAL_DATA_TASKS),DIMENSION(:),POINTER,SAVE ::                &
                                                 PARENT_BOUND_H_SOUTH   & 
                                                ,PARENT_BOUND_H_NORTH   & 
                                                ,PARENT_BOUND_H_WEST    & 
                                                ,PARENT_BOUND_H_EAST    &
                                                ,PARENT_BOUND_V_SOUTH   &
                                                ,PARENT_BOUND_V_NORTH   &
                                                ,PARENT_BOUND_V_WEST    &
                                                ,PARENT_BOUND_V_EAST
!
      TYPE(REAL_DATA_TASKS),DIMENSION(:),POINTER,SAVE :: PD_B_SOUTH     &
                                                        ,PD_B_NORTH     &
                                                        ,PD_B_WEST      &
                                                        ,PD_B_EAST      &
!
                                                        ,PD_B_SOUTH_V   &
                                                        ,PD_B_NORTH_V   &
                                                        ,PD_B_WEST_V    &
                                                        ,PD_B_EAST_V    &
!
                                                        ,T_B_SOUTH      &
                                                        ,T_B_NORTH      &
                                                        ,T_B_WEST       &
                                                        ,T_B_EAST       &
!
                                                        ,Q_B_SOUTH      &
                                                        ,Q_B_NORTH      &
                                                        ,Q_B_WEST       &
                                                        ,Q_B_EAST       &
!
                                                        ,CW_B_SOUTH     &
                                                        ,CW_B_NORTH     &
                                                        ,CW_B_WEST      &
                                                        ,CW_B_EAST      &
!
                                                        ,U_B_SOUTH      &
                                                        ,U_B_NORTH      &
                                                        ,U_B_WEST       &
                                                        ,U_B_EAST       &
!
                                                        ,V_B_SOUTH      &
                                                        ,V_B_NORTH      &
                                                        ,V_B_WEST       &
                                                        ,V_B_EAST       &
!
                                                        ,FIS_CHILD_SOUTH &
                                                        ,FIS_CHILD_NORTH &
                                                        ,FIS_CHILD_WEST  &
                                                        ,FIS_CHILD_EAST
!
      TYPE(INTEGER_DATA_TASKS),DIMENSION(:),POINTER,SAVE ::             &
                                                 WORDS_BOUND_H_SOUTH    & 
                                                ,WORDS_BOUND_H_NORTH    & 
                                                ,WORDS_BOUND_H_WEST     & 
                                                ,WORDS_BOUND_H_EAST     &
                                                ,WORDS_BOUND_V_SOUTH    &
                                                ,WORDS_BOUND_V_NORTH    &
                                                ,WORDS_BOUND_V_WEST     &
                                                ,WORDS_BOUND_V_EAST
!
      TYPE(CHILD_TASKS),DIMENSION(:),POINTER,SAVE :: CTASK_LIMITS
!
      TYPE(SIDES_0D),SAVE :: INDX_MAX_H                                 &
                            ,INDX_MAX_V                                 &
                            ,INDX_MIN_H                                 &
                            ,INDX_MIN_V                                 &
!!!                         ,LENGTH_BND_SEG_H                           &
!!!                         ,LENGTH_BND_SEG_V                           &
                            ,NUM_PARENT_TASKS_SENDING_H                 &
                            ,NUM_PARENT_TASKS_SENDING_V
!
      TYPE(SIDES_1D_INT),DIMENSION(:),POINTER,SAVE ::                   &
                                             CHILDTASK_BNDRY_H_RANKS    &
                                            ,CHILDTASK_BNDRY_V_RANKS
!
      TYPE(SAVE_TASK_IJ),DIMENSION(:),POINTER,SAVE :: CHILDTASK_H_SAVE  &
                                                     ,CHILDTASK_V_SAVE     !<-- 'TRANSFER' subcomponents of SAVE_TASK_IJ are irrelevant
!
      TYPE(PARENT_POINTS_SURROUND_H),DIMENSION(:),POINTER,SAVE ::       &
                                                   PARENT_4_INDICES_H
!
      TYPE(PARENT_POINTS_SURROUND_V),DIMENSION(:),POINTER,SAVE ::       &
                                                   PARENT_4_INDICES_V
!
      TYPE(PARENT_WEIGHTS_SURROUND_H),DIMENSION(:),POINTER,SAVE ::      &
                                                    PARENT_4_WEIGHTS_H
!
      TYPE(PARENT_WEIGHTS_SURROUND_V),DIMENSION(:),POINTER,SAVE ::      &
                                                    PARENT_4_WEIGHTS_V
!
      TYPE(PARENT_DATA),DIMENSION(:),POINTER,SAVE :: PARENT_TASK
!
      TYPE(HANDLE_SEND),DIMENSION(:),POINTER,SAVE :: HANDLE_H_SOUTH     &
                                                    ,HANDLE_H_NORTH     &
                                                    ,HANDLE_H_WEST      &
                                                    ,HANDLE_H_EAST      &
                                                    ,HANDLE_V_SOUTH     &
                                                    ,HANDLE_V_NORTH     &
                                                    ,HANDLE_V_WEST      &
                                                    ,HANDLE_V_EAST
!
!-----------------------------------------------------------------------
!
      REAL(KIND=KDBL) :: btim,btim0
!
      REAL(KIND=KDBL),SAVE :: cpl1_prelim_tim                           &
                             ,cpl1_south_h_tim,cpl1_south_v_tim         &
                             ,cpl1_north_h_tim,cpl1_north_v_tim         &
                             ,cpl1_west_h_tim, cpl1_west_v_tim          &
                             ,cpl1_east_h_tim, cpl1_east_v_tim          &
                             ,cpl1_recv_tim
!
      REAL(KIND=KDBL),SAVE :: cpl1_south_h_recv_tim                     &
                             ,cpl1_south_h_undo_tim                     &
                             ,cpl1_south_h_exp_tim                      &
                             ,cpl1_south_v_recv_tim                     &
                             ,cpl1_south_v_undo_tim                     &
                             ,cpl1_south_v_exp_tim
!
      REAL(KIND=KDBL),SAVE :: cpl2_comp_tim                             &
                             ,cpl2_send_tim                             &
                             ,cpl2_wait_tim
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_CPL_REGISTER(CPL_COMP,RC_NEST_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE NESTING COUPLER COMPONENT'S INITIALIZE, RUN, AND 
!***  FINALIZE ROUTINES.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- Coupler component
!
      INTEGER(kind=KINT),INTENT(OUT)   :: RC_NEST_REG                     !<-- Return code for register
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
      RC         =ESMF_SUCCESS
      RC_NEST_REG=ESMF_SUCCESS

!-----------------------------------------------------------------------
!***  Register the coupler Initialize subroutine.  Since it is just one
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN,
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Nesting Coupler Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_INITIALIZE        &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the coupler Run subroutine.
!***  The Parent-Child Run step of the coupler has two distinct parts.
!***  The first occurs at the top of the timesteps where a child might
!***  need to receive data from its parent.  The second occurs at the
!***  end of every timestep where a parent sends to its children.
!***  So register the coupler's Run step with those two phases.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Phase 1 of Nesting Coupler Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_RECV          &  !<-- User's subroutineName
                                    ,1                                  &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Phase 2 of Nesting Coupler Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_SEND          &  !<-- User's subroutineName
                                    ,2                                  &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the coupler Finalize subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Nesting Coupler Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_FINALIZE          &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Check the error signal variable.
!-----------------------------------------------------------------------
!
      IF(RC_NEST_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)" NESTING COUPLER REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" NESTING COUPLER REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_CPL_INITIALIZE(CPL_COMP                   &
                                            ,IMP_STATE                  &
                                            ,EXP_STATE                  &
                                            ,CLOCK                      &
                                            ,RC_FINAL)   
!-----------------------------------------------------------------------
!***  PERFORM INITIAL WORK NEEDED BY THE PARENT-CHILD COUPLER.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Dyn-Phy Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,  INTENT(OUT)   :: RC_FINAL
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L                                       &
                           ,I_END_TRANSFER,J_END_TRANSFER               &
                           ,NBASE,NBASE_3D,NBASE_EXP,NT
!
      INTEGER(kind=KINT) :: CHILD_ID,IERR,ITE_CHILD_X,JTE_CHILD_X       &
                           ,LIM1_H,LIM1_V,LIM2_H,LIM2_V                 &
                           ,N,NN,N_START,N_END                          &
                           ,N_H_EAST_WEST,N_H_NORTH_SOUTH               &
                           ,N_V_EAST_WEST,N_V_NORTH_SOUTH               &
                           ,N_TASK,NCHILD_TASKS                         &
                           ,NLOC_1,NLOC_2,NLOC_2_EXP                    &
                           ,NTIMESTEP                                   &
                           ,NUM_BOUNDARY_WORDS,NUM_DOMAINS              &
                           ,NUM_TASKS_MINE,NUM_TASKS_PARENT             &
                           ,NWORDS
!
      INTEGER(kind=KINT) :: RC,RC_CPL_INIT
!
      REAL(kind=KFPT) :: DIST_NESTV_SOUTH_TO_PARENTV_SOUTH
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT='(I2.2)'
      CHARACTER(99) :: CONFIG_FILE_NAME
!
      TYPE(ESMF_Array) :: HOLD_ARRAY
!
      TYPE(ESMF_Config),DIMENSION(:),ALLOCATABLE :: CF
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      integer :: mype
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINAL   =ESMF_SUCCESS
      RC_CPL_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Extract key variables from the import state.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------
!***  Current Domain ID
!-----------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Current Domain ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='MY_DOMAIN_ID'                       &  !<-- Name of the attribute to extract
                            ,value=MY_DOMAIN_ID                         &  !<-- Current domain's ID
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------
!***  Number of Children
!------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Number of Children on This Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='NUM_CHILDREN'                       &  !<-- Name of the attribute to extract
                            ,value=NUM_CHILDREN                         &  !<-- How many children does this ATM Component have?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------
!***  Communicator for each domain
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain Intercommunicator in Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='Single Domain Comm'                 &  !<-- Name of the attribute to extract
                            ,value=COMM_MY_DOMAIN                       &  !<-- MPI communicator for each individual domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------
!***  Child-to-Parent Communicator
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Child-to-Parent Comm in Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='Child-to-Parent Comm'               &  !<-- Name of the attribute to extract
                            ,value=COMM_TO_MY_PARENT                    &  !<-- MPI communicator to my parent
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  Total Number of Domains
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Total Number of Domains in Init Step Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='NUM_DOMAINS'                        &  !<-- Name of the attribute to extract
                            ,value=NUM_DOMAINS                          &  !<-- Total number of domains
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------
!***  Forecast Tasks On Each Domain
!-----------------------------------
!
      ALLOCATE(FTASKS_DOMAIN(1:NUM_DOMAINS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Number of Fcst Tasks on Each Domain in Init Step Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Name of the attribute to extract
                            ,count    =NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valuelist=FTASKS_DOMAIN                    &  !<-- # of forecast tasks on each domain
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------
!***  Domain IDs of Parents
!---------------------------
!
      ALLOCATE(ID_PARENTS(1:NUM_DOMAINS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain IDs of Parents in Init Step Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='ID_PARENTS'                     &  !<-- Name of the attribute to extract
                            ,count    =NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valuelist=ID_PARENTS                       &  !<-- Domain IDs of parents 
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  The Tasks' Subdomain Integration Limits
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Integration Subdomain Limits in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='ITS'                                &  !<-- Name of the attribute to extract
                            ,value=ITS                                  &  !<-- This task's integration limit: Starting I
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='ITE'                                &  !<-- Name of the attribute to extract
                            ,value=ITE                                  &  !<-- This task's integration limit: Ending I
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='JTS'                                &  !<-- Name of the attribute to extract
                            ,value=JTS                                  &  !<-- This task's integration limit: Starting J
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='JTE'                                &  !<-- Name of the attribute to extract
                            ,value=JTE                                  &  !<-- This task's integration limit: Ending J
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='LM'                                 &  !<-- Name of the attribute to extract
                            ,value=LM                                   &  !<-- This task's integration limit: # of layers in vertical
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='NHALO'                              &  !<-- Name of the attribute to extract
                            ,value=NHALO                                &  !<-- Width of the task subdomain haloes
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMS=ITS-NHALO
      IME=ITE+NHALO
      JMS=JTS-NHALO
      JME=JTE+NHALO
!
!--------------------------------
!***  The Full Domain Dimensions
!--------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Full Domain Dimensions in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='IDS'                                &  !<-- Name of the attribute to extract
                            ,value=IDS                                  &  !<-- This task's integration limit: Starting I
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='IDE'                                &  !<-- Name of the attribute to extract
                            ,value=IDE                                  &  !<-- This task's integration limit: Ending I
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='JDS'                                &  !<-- Name of the attribute to extract
                            ,value=JDS                                  &  !<-- This task's integration limit: Starting J
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='JDE'                                &  !<-- Name of the attribute to extract
                            ,value=JDE                                  &  !<-- This task's integration limit: Ending J
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  The Widths of the Boundary Blending Region
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Widths of Bndry Blending Region in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='N_BLEND_H'                          &  !<-- Name of the attribute to extract
                            ,value=N_BLEND_H                            &  !<-- # of boundary blending rows for H points
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='N_BLEND_V'                          &  !<-- Name of the attribute to extract
                            ,value=N_BLEND_V                            &  !<-- # of boundary blending rows for V points
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(N_BLEND_V>N_BLEND_H)THEN
        WRITE(0,*)' N_BLEND_V CANNOT EXCEED N_BLEND_H DUE TO PD AVERAGING!!!'
        WRITE(0,*)' ABORTING'
        CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
      ENDIF
!
!-----------------------------------------------------------------------
!
      child_block: IF(NUM_CHILDREN>0)THEN                                  !<-- Tasks proceed if their domain has children
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST THE PARENT TASKS EXTRACT RELEVANT DATA ABOUT THEMSELVES
!***  FROM THE COUPLER IMPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------
!***  Parent-to-Child Communicators
!-----------------------------------
!
        ALLOCATE(COMM_TO_MY_CHILDREN(1:NUM_CHILDREN))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-to-Child Comm in Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-to-Child Comms'        &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=COMM_TO_MY_CHILDREN            &  !<-- MPI communicators to my children
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  The IDs of the Children
!-----------------------------
!
        ALLOCATE(MY_CHILDREN_ID(1:NUM_CHILDREN))

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract IDs of Children in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valuelist=MY_CHILDREN_ID                 &  !<-- The domain IDs of the current domain's children
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  UNLOAD THE PARENTS'S OWN PROGNOSTIC DATA SO IT CAN BE USED
!***  TO INTERPOLATE TO ITS CHILDREN.
!-----------------------------------------------------------------------
!
!--------
!***  PD
!--------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PD from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='PD'                                &  !<-- Extract PD
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PD from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=PD                                 &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------
!***  Temperature
!-----------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract T from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='T'                                 &  !<-- Extract temperature
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract T from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=T                                  &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------
!***  U Wind
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract U from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='U'                                 &  !<-- Extract U wind
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract U from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=U                                  &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------
!***  V Wind
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract V from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='V'                                 &  !<-- Extract V wind
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract V from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=V                                  &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  Tracers
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Tracers from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='TRACERS'                           &  !<-- Extract tracers
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Tracers from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=TRACERS                            &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------
!***  Sfc Geopotential
!----------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract FIS from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='FIS'                               &  !<-- Extract PD
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract FIS from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=FIS                                &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------
!***  PT,PDTOP,PSGML1,SG1,SG2,SGML2
!-----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PT from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The parent-child coupler import state
                              ,name ='PT'                               &  !<-- Name of Attribute to extract
                              ,value=PT                                 &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PDTOP from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The parent-child coupler import state
                              ,name ='PDTOP'                            &  !<-- Name of Attribute to extract
                              ,value=PDTOP                              &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALLOCATE(PSGML1(1:LM))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PSGML1 from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='PSGML1'                       &  !<-- Name of Attribute to extract
                              ,count    =LM                             &  !<-- # of words in data list
                              ,valueList=PSGML1                         &  !<-- Put the extracted Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALLOCATE(SG1(1:LM+1))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract SG1 from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='SG1'                          &  !<-- Name of Attribute to extract
                              ,count    =LM+1                           &  !<-- # of words in data list
                              ,valueList=SG1                            &  !<-- Put the extracted Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALLOCATE(SG2(1:LM+1))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract SG2 from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='SG2'                          &  !<-- Name of Attribute to extract
                              ,count    =LM+1                           &  !<-- # of words in data list
                              ,valueList=SG2                            &  !<-- Put the extracted Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALLOCATE(SGML2(1:LM))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract SGML2 from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='SGML2'                        &  !<-- Name of Attribute to extract
                              ,count    =LM                             &  !<-- # of words in data list
                              ,valueList=SGML2                          &  !<-- Put the extracted Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------
!***  INDX_Q
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract INDX_Q from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The parent-child coupler import state
                              ,name ='INDX_Q'                           &  !<-- Name of Attribute to extract
                              ,value=INDX_Q                             &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  INDX_CW
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract INDX_CW from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The parent-child coupler import state
                              ,name ='INDX_CW'                          &  !<-- Name of Attribute to extract
                              ,value=INDX_CW                            &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The parents create and load the configure objects
!***  of their children.
!-----------------------------------------------------------------------
!
        ALLOCATE(CF(1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
          CF(N)=ESMF_ConfigCreate(rc=RC)
!
          CHILD_ID=MY_CHILDREN_ID(N)
          WRITE(INT_TO_CHAR,FMT)CHILD_ID 
          CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                  !<-- Prepare the config file names
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Load Configure Files"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(N)                       &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate the arrays holding the children's domain sizes 
!***  and extract those values from the configure files.
!-----------------------------------------------------------------------
!
        ALLOCATE(IM_CHILD(1:NUM_CHILDREN))
        ALLOCATE(JM_CHILD(1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Extract Global IM,JM of Child"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- The child's config object
                                      ,value =IM_CHILD(N)               &  !<-- The variable filled (IM of child domain)
                                      ,label ='im:'                     &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- The child's config object
                                      ,value =JM_CHILD(N)               &  !<-- The variable filled (IM of child domain)
                                      ,label ='jm:'                     &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Initialize the timestep of the children at which they receive
!***  data from their parent.  This must be known by the parents
!***  since it will provide the proper tag to the MPI data sent.
!-----------------------------------------------------------------------
!
        ALLOCATE(PARENT_CHILD_TIME_RATIO(1:NUM_CHILDREN))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-to-Child DT Ratio from Imp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=PARENT_CHILD_TIME_RATIO        &  !<-- Ratio of parent to child DTs 
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the timestep from the Clock to set the initial value for
!***  the child's timestep at which it will next receive parent data.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK                           &
                          ,advanceCount=NTIMESTEP_ESMF                  &
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        ALLOCATE(NSTEP_CHILD_RECV(1:NUM_CHILDREN))                         !<-- Children's timesteps at which they recv data
!
        DO N=1,NUM_CHILDREN
          NSTEP_CHILD_RECV(N)=(NTIMESTEP-1)*PARENT_CHILD_TIME_RATIO(N)
        ENDDO
!
!-----------------------------------------------------------------------
!***  EXTRACT THE RATIO OF THE PARENT GRID INCREMENT TO THE CHILDREN'S
!***  FROM THE CONFIGURE FILES.
!-----------------------------------------------------------------------
!
        ALLOCATE(PARENT_CHILD_SPACE_RATIO(1:NUM_CHILDREN))
!
        ALLOCATE(INC_FIX(1:NUM_CHILDREN))                                  !<-- See immediately below
!
        DO N=1,NUM_CHILDREN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Parent-to-Child Space Ratio"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                      ,value =PARENT_CHILD_SPACE_RATIO(N) &  !<-- The variable filled (# of child grid inc's to parent's)
                                      ,label ='parent_child_space_ratio:' &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  USE THIS RATIO TO COMPUTE AN INCREMENT THAT IS NEEDED FOR
!***  SELECTING THE APPROPRIATE NEST TASKS AS MASS VALUES ON THE
!***  NEST BOUNDARIES ARE AVERAGED TO THE V POINTS.  ITS VALUES
!***  ARE BASED ON THE NEST GRID INCREMENT DISTANCE FROM THE
!***  SOUTHERNMOST V POINT ON A NEST's SOUTHERNMOST TASKS TO
!***  THE NEAREST PARENT V POINT TO THE NORTH.  FRACTIONAL VALUES
!***  ARE INCREASED TO THE NEXT INTEGER.
!-----------------------------------------------------------------------
!
          DIST_NESTV_SOUTH_TO_PARENTV_SOUTH=                            &
                                     (PARENT_CHILD_SPACE_RATIO(N)-1)*0.5
          INC_FIX(N)=INT(DIST_NESTV_SOUTH_TO_PARENTV_SOUTH+0.9)
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Pre-compute various indices and weights needed by the parents
!***  to compute boundary data for their children.
!***  Only the parents are calling this routine.
!-----------------------------------------------------------------------
!
        CALL PARENT_CHILD_INTERP_SETUP(NUM_CHILDREN                     &
                                      ,MY_CHILDREN_ID                   &
                                      ,IM_CHILD                         &
                                      ,JM_CHILD                         &
                                      ,FTASKS_DOMAIN                    &
                                      ,CTASK_LIMITS                     &
                                      ,N_BLEND_H                        &
                                      ,N_BLEND_V                        &
                                      ,CF                               &
                                      ,ITS,ITE,JTS,JTE                  &
                                      ,IDS,IDE,JDS,JDE )
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that will hold all of the interpolated
!***  boundary data for the child tasks if the parent task contains
!***  child boundary points on the four sides.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The data pointer within the CHILD_BOUND_* arrays will
!***  hold each of the children's tasks' boundary data in a 
!***  1-D string that will be used to send the data from the parent
!***  to the child tasks.  The *_B_* arrays are what will be sent
!***  into subroutine PARENT_TO_CHILD_BNDRY_COMPUTE where the
!***  actual computations of the boundary data are carried out.
!***  The unallocated subcomponents of the *_B_* arrays are filled
!***  by the routine and thus the 1-D string will automatically
!***  be filled and ready for sending to each child task that contains
!***  boundary points.
!
!***  NOTE the heirarchy of the derived data variables' pointers
!***  that hold the boundary data: 
!
!***  (1) Primary variable dimensioned 1-D over the children.
!***  (2) For each child, TASKS is dimensioned 1-D over the given
!***      child's tasks that contain segments of boundary on the
!***      parent task.
!***  (3) For each child task, DATA is dimensioned 1-D since 
!***      1-D strings are required for MPI Send/Recv.
!***      (a) The 1-D CHILD_BOUND_* DATA pointers are allocated
!***          and are filled with the child boundary data destined
!***          for each individual child task that contains any 
!***          segment of a child's boundary on a given parent
!***          task.
!***      (b) The specific boundary variables (T_B_*, Q_B_*, etc.)
!***          DATA subcomponent pointers are declared but never
!***          allocated but instead are simply pointed into the
!***          allocated 1-D CHILD_BOUND_* DATA subcomponent pointer.
!***          These specific boundary variables are sent into the
!***          subroutine where the child boundary data is computed
!***          which then leads to the allocated CHILD_BOUND_* 1-D
!***          DATA pointer being filled automatically.  Thus the
!***          allocated 1-D pointer is immediately ready for subsequent
!***          sending to child tasks.
!***          The subcomponent pointer for PD_B_* though must be
!***          allocated.  That is because it contains values one row
!***          beyond those actually needed to be sent to the child
!***          tasks to update their boundary values of PD.  That 
!***          extra row is used to do 4-pt averaging to obtain PD
!***          on V points within the nests' boundary regions which
!***          is needed to do the hydrostatic updating of the winds
!***          there.  So the PD_B sections of the full 1-D DATA
!***          pointer CHILD_BOUND_* that is actually sent from
!***          parents to children must be filled explicitly inside
!***          subroutine PARENT_UPDATE_CHILD_PSFC.
!-----------------------------------------------------------------------
!
        ALLOCATE(CHILD_BOUND_H_SOUTH(1:NUM_CHILDREN))                       !<-- 1-D bndry data string for child tasks with Sbndry H points
        ALLOCATE(CHILD_BOUND_V_SOUTH(1:NUM_CHILDREN))                       !<-- 1-D bndry data string for child tasks with Sbndry V points
        ALLOCATE(WORDS_BOUND_H_SOUTH(1:NUM_CHILDREN))                       !<-- # of words in Sbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_SOUTH(1:NUM_CHILDREN))                       !<-- # of words in Sbndry V point 1-D data string
!
        ALLOCATE(PD_B_SOUTH(1:NUM_CHILDREN))                                !<-- South boundary PD
        ALLOCATE(T_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary temperature
        ALLOCATE(Q_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary specific humidity
        ALLOCATE(CW_B_SOUTH(1:NUM_CHILDREN))                                !<-- South boundary cloud condensate
        ALLOCATE(U_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary U wind component
        ALLOCATE(V_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary V wind component
!
        ALLOCATE(PD_B_SOUTH_V(1:NUM_CHILDREN))                              !<-- South boundary PD on V points
!
        ALLOCATE(CHILD_BOUND_H_NORTH(1:NUM_CHILDREN))                       !<-- 1-D bndry data string for child tasks with Nbndry H points
        ALLOCATE(CHILD_BOUND_V_NORTH(1:NUM_CHILDREN))                       !<-- 1-D bndry data string for child tasks with Nbndry V points
        ALLOCATE(WORDS_BOUND_H_NORTH(1:NUM_CHILDREN))                       !<-- # of words in Nbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_NORTH(1:NUM_CHILDREN))                       !<-- # of words in Nbndry V point 1-D data string
!
        ALLOCATE(PD_B_NORTH(1:NUM_CHILDREN))                                !<-- North boundary PD
        ALLOCATE(T_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary temperature
        ALLOCATE(Q_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary specific humidity
        ALLOCATE(CW_B_NORTH(1:NUM_CHILDREN))                                !<-- North boundary cloud condensate
        ALLOCATE(U_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary U wind component
        ALLOCATE(V_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary V wind component
!
        ALLOCATE(PD_B_NORTH_V(1:NUM_CHILDREN))                              !<-- North boundary PD on V points
!
        ALLOCATE(CHILD_BOUND_H_WEST(1:NUM_CHILDREN))                        !<-- 1-D bndry data string for child tasks with Wbndry H points
        ALLOCATE(CHILD_BOUND_V_WEST(1:NUM_CHILDREN))                        !<-- 1-D bndry data string for child tasks with Wbndry V points
        ALLOCATE(WORDS_BOUND_H_WEST(1:NUM_CHILDREN))                        !<-- # of words in Wbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_WEST(1:NUM_CHILDREN))                        !<-- # of words in Wbndry V point 1-D data string
!
        ALLOCATE(PD_B_WEST(1:NUM_CHILDREN))                                 !<-- West boundary PD
        ALLOCATE(T_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary temperature
        ALLOCATE(Q_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary specific humidity
        ALLOCATE(CW_B_WEST(1:NUM_CHILDREN))                                 !<-- West boundary cloud condensate
        ALLOCATE(U_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary U wind component
        ALLOCATE(V_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary V wind component
!
        ALLOCATE(PD_B_WEST_V(1:NUM_CHILDREN))                               !<-- West boundary PD on V points
!
        ALLOCATE(CHILD_BOUND_H_EAST(1:NUM_CHILDREN))                        !<-- 1-D bndry data string for child tasks with Ebndry H points
        ALLOCATE(CHILD_BOUND_V_EAST(1:NUM_CHILDREN))                        !<-- 1-D bndry data string for child tasks with Ebndry V points
        ALLOCATE(WORDS_BOUND_H_EAST(1:NUM_CHILDREN))                        !<-- # of words in Ebndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_EAST(1:NUM_CHILDREN))                        !<-- # of words in Ebndry V point 1-D data string
!
        ALLOCATE(PD_B_EAST(1:NUM_CHILDREN))                                 !<-- East boundary PD
        ALLOCATE(T_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary temperature
        ALLOCATE(Q_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary specific humidity
        ALLOCATE(CW_B_EAST(1:NUM_CHILDREN))                                 !<-- East boundary cloud condensate
        ALLOCATE(U_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary U wind component
        ALLOCATE(V_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary V wind component
!
        ALLOCATE(PD_B_EAST_V(1:NUM_CHILDREN))                               !<-- East boundary PD on V points
!
!-----------------------------------------------------------------------
        point_loop: DO N=1,NUM_CHILDREN
!-----------------------------------------------------------------------
!
!-----------
!***  South
!-----------
!
          south_h: IF(NUM_TASKS_SEND_H_S(N)>0)THEN                          !<-- Parent task has child south boundary H points?
            NCHILD_TASKS=NUM_TASKS_SEND_H_S(N)
            ALLOCATE(CHILD_BOUND_H_SOUTH(N)%TASKS(1:NCHILD_TASKS))          !<-- 1-D bndry data string for child tasks with Sbndry H points
            ALLOCATE(WORDS_BOUND_H_SOUTH(N)%TASKS(1:NCHILD_TASKS))          !<-- # of words in Sbndry H point 1-D data string
!
            ALLOCATE(PD_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                   !<-- PD_B_SOUTH for each child task
            ALLOCATE(T_B_SOUTH (N)%TASKS(1:NCHILD_TASKS))                   !<-- T_B_SOUTH for each child task
            ALLOCATE(Q_B_SOUTH (N)%TASKS(1:NCHILD_TASKS))                   !<-- Q_B_SOUTH for each child task
            ALLOCATE(CW_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                   !<-- CW_B_SOUTH for each child task
!
            DO NT=1,NCHILD_TASKS
              N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NT)+1
              ITE_CHILD_X=CTASK_LIMITS(N)%LIMITS(2,N_TASK)
              I_END_TRANSFER=MIN(ITE_CHILD_X+2,IM_CHILD(N))
!
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further east than the easternmost V in the segment in order to
!***  be able to do 4-pt averaging of PDB onto the V points in order to
!***  do hydrostatic updating of V by the parent.  Now indicate that
!***  reduction of the points to be transferred.
!
!!!           CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)=I_END_TRANSFER    !<-- Sbndry I limit for transfer to child
              CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)=              &   !<-- Sbndry I limit for transfer to child
                CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)-1
              IF(CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)==ITE_CHILD_X)       &   !<-- We do not reduce the area for H data   
                CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)=            &   !    transfer if the bndry segment reaches
                CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)+1               !    the physical limit of that bndry
!
              NBASE=CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)         &
                   -CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NT)+1          
!
              NBASE_3D=LM*NBASE*N_BLEND_H
              NWORDS  =(3*LM+1)*NBASE*N_BLEND_H                             !<-- # of Sbndry words to transfer from parent to child
              WORDS_BOUND_H_SOUTH(N)%TASKS(NT)=NWORDS                       !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA(1:NWORDS))     !<-- 1-D bndry data string for child tasks
!                                                                                with south boundary H points
              NLOC_1    =1
              NLOC_2    =NLOC_1+NBASE*N_BLEND_H-1
!
              NBASE_EXP=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)             &
                       -CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NT)+1
!
              NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H+1)-1                   !<-- Extend PD_B_* to allow 4-pt averaging to V pts
!!!!!         PD_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- PD_B_SOUTH storage location
              ALLOCATE(PD_B_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP))
!
              NLOC_1=NLOC_2+1                                               !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
              NLOC_2=NLOC_1+NBASE_3D-1
              T_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_SOUTH storage location
      
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              Q_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_SOUTH storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              CW_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_SOUTH storage location
            ENDDO
!
          ELSE  south_h                                                     !<-- Dummy nonzero length 
!
            ALLOCATE(PD_B_SOUTH(N)%TASKS(1:1))              
            ALLOCATE(T_B_SOUTH (N)%TASKS(1:1))             
            ALLOCATE(Q_B_SOUTH (N)%TASKS(1:1))            
            ALLOCATE(CW_B_SOUTH(N)%TASKS(1:1))           
!
          ENDIF south_h
!
          south_v: IF(NUM_TASKS_SEND_V_S(N)>0)THEN                          !<-- Parent task has child south boundary V points?
            NCHILD_TASKS=NUM_TASKS_SEND_V_S(N)
            ALLOCATE(CHILD_BOUND_V_SOUTH(N)%TASKS(1:NCHILD_TASKS))          !<-- 1-D bndry data string for child tasks with Sbndry V points
            ALLOCATE(WORDS_BOUND_V_SOUTH(N)%TASKS(1:NCHILD_TASKS))          !<-- # of words in Sbndry V point 1-D data string
!
            ALLOCATE(PD_B_SOUTH_V(N)%TASKS(1:NCHILD_TASKS))                 !<-- PD_B_SOUTH_V for each child task
            ALLOCATE(U_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                    !<-- U_B_SOUTH for each child task
            ALLOCATE(V_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                    !<-- V_B_SOUTH for each child task
!
            DO NT=1,NCHILD_TASKS
              NBASE   =CHILDTASK_V_SAVE(N)%I_HI_SOUTH(NT)                &
                      -CHILDTASK_V_SAVE(N)%I_LO_SOUTH(NT)+1
              NBASE_3D=LM*NBASE*N_BLEND_V
              NWORDS  =2*NBASE_3D                                           !<-- Total number of V boundary words for child task's segment
              WORDS_BOUND_V_SOUTH(N)%TASKS(NT)=NWORDS                       !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_V_SOUTH(N)%TASKS(NT)%DATA(1:NWORDS))     !<-- 1-D bndry data string for child tasks with Sbndry V points
!
              ALLOCATE(PD_B_SOUTH_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V))
!
              NLOC_1=1
              NLOC_2=NLOC_1+NBASE_3D-1
              U_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_SOUTH storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              V_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_SOUTH storage location
!
            ENDDO
!
          ELSE  south_v                                                     !<-- Dummy nonzero length
!
            ALLOCATE(U_B_SOUTH(N)%TASKS(1:1))
            ALLOCATE(V_B_SOUTH(N)%TASKS(1:1))
!
          ENDIF south_v
!
!-----------
!***  North
!-----------
!
          north_h: IF(NUM_TASKS_SEND_H_N(N)>0)THEN                          !<-- Parent task has child north boundary H points?
            NCHILD_TASKS=NUM_TASKS_SEND_H_N(N)
            ALLOCATE(CHILD_BOUND_H_NORTH(N)%TASKS(1:NCHILD_TASKS))          !<-- 1-D bndry data string for child tasks with Nbndry H points
            ALLOCATE(WORDS_BOUND_H_NORTH(N)%TASKS(1:NCHILD_TASKS))          !<-- # of words in Nbndry H point 1-D data string
!
            ALLOCATE(PD_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                   !<-- PD_B_NORTH for each child task
            ALLOCATE(T_B_NORTH (N)%TASKS(1:NCHILD_TASKS))                   !<-- T_B_NORTH for each child task
            ALLOCATE(Q_B_NORTH (N)%TASKS(1:NCHILD_TASKS))                   !<-- Q_B_NORTH for each child task
            ALLOCATE(CW_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                   !<-- CW_B_NORTH for each child task
!
            DO NT=1,NCHILD_TASKS
              N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NT)+1
              ITE_CHILD_X=CTASK_LIMITS(N)%LIMITS(2,N_TASK)
              I_END_TRANSFER=MIN(ITE_CHILD_X+2,IM_CHILD(N))
!
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further east than the easternmost V in the segment in order to
!***  be able to do 4-pt averaging of PDB onto the V points in order to
!***  do hydrostatic updating of V by the parent.  Now indicate that
!***  reduction of the points to be transferred.
!
!!!           CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)=I_END_TRANSFER    !<-- Nbndry I limit for transfer to child
              CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)=              &   !<-- Nbndry I limit for transfer to child
                CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)-1
!
              IF(CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)==ITE_CHILD_X)       &   !<-- We do not reduce the area for H data
                CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)=            &   !    transfer if the bndry segment reaches
                CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)+1               !    the physical limit of that bndry
!
              NBASE=CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)         &
                   -CHILDTASK_H_SAVE(N)%I_LO_NORTH(NT)+1
!
              NBASE_3D=LM*NBASE*N_BLEND_H
              NWORDS  =(3*LM+1)*NBASE*N_BLEND_H                             !<-- # of Nbndry words to transfer from parent to child
              WORDS_BOUND_H_NORTH(N)%TASKS(NT)=NWORDS                       !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA(1:NWORDS))     !<-- 1-D bndry data string for child tasks
!                                                                                with north boundary H points
              NLOC_1=1
              NLOC_2=NLOC_1+NBASE*N_BLEND_H-1
!
              NBASE_EXP=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)              &
                       -CHILDTASK_H_SAVE(N)%I_LO_NORTH(NT)+1
!
              NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H+1)-1                   !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!!!!!         PD_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- PD_B_NORTH storage location
              ALLOCATE(PD_B_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP))
!
              NLOC_1=NLOC_2+1                                               !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
              NLOC_2=NLOC_1+NBASE_3D-1
              T_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_NORTH storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              Q_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_NORTH storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              CW_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_NORTH storage location
            ENDDO
!
          ELSE  north_h                                                     !<-- Dummy nonzero length
!
            ALLOCATE(PD_B_NORTH(N)%TASKS(1:1))
            ALLOCATE(T_B_NORTH (N)%TASKS(1:1))
            ALLOCATE(Q_B_NORTH (N)%TASKS(1:1))
            ALLOCATE(CW_B_NORTH(N)%TASKS(1:1))
!
          ENDIF north_h
!
          north_v: IF(NUM_TASKS_SEND_V_N(N)>0)THEN                          !<-- Parent task has child north boundary V points?
            NCHILD_TASKS=NUM_TASKS_SEND_V_N(N)
            ALLOCATE(CHILD_BOUND_V_NORTH(N)%TASKS(1:NCHILD_TASKS))          !<-- 1-D bndry data string for child tasks with Nbndry V points
            ALLOCATE(WORDS_BOUND_V_NORTH(N)%TASKS(1:NCHILD_TASKS))          !<-- # of words in Nbndry V point 1-D data string
!
            ALLOCATE(PD_B_NORTH_V(N)%TASKS(1:NCHILD_TASKS))                 !<-- PD_B_NORTH_V for each child task
            ALLOCATE(U_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                    !<-- U_B_NORTH for each child task
            ALLOCATE(V_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                    !<-- V_B_NORTH for each child task
!
            DO NT=1,NCHILD_TASKS
              NBASE=CHILDTASK_V_SAVE(N)%I_HI_NORTH(NT)                   &
                   -CHILDTASK_V_SAVE(N)%I_LO_NORTH(NT)+1
              NWORDS=2*LM*NBASE*N_BLEND_V                                   !<-- Total number of V boundary words for child task's segment
              WORDS_BOUND_V_NORTH(N)%TASKS(NT)=NWORDS                       !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_V_NORTH(N)%TASKS(NT)%DATA(1:NWORDS))     !<-- 1-D bndry data string for child tasks with Nbndry V points
!
              ALLOCATE(PD_B_NORTH_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V))
!
              NLOC_1=1
              NLOC_2=NLOC_1+LM*NBASE-1
              U_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_NORTH storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+LM*NBASE-1
              V_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_NORTH storage location
!
            ENDDO
!
          ELSE  north_v                                                     !<-- Dummy nonzero length
!
            ALLOCATE(U_B_NORTH(N)%TASKS(1:1))
            ALLOCATE(V_B_NORTH(N)%TASKS(1:1))
!
          ENDIF north_v
!
!----------
!***  West
!----------
!
          west_h: IF(NUM_TASKS_SEND_H_W(N)>0)THEN                           !<-- Parent task has child west boundary H points?
            NCHILD_TASKS=NUM_TASKS_SEND_H_W(N)
            ALLOCATE(CHILD_BOUND_H_WEST(N)%TASKS(1:NCHILD_TASKS))           !<-- 1-D bndry data string for child tasks with Wbndry H points
            ALLOCATE(WORDS_BOUND_H_WEST(N)%TASKS(1:NCHILD_TASKS))           !<-- # of words in Wbndry H point 1-D data string
!
            ALLOCATE(PD_B_WEST(N)%TASKS(1:NCHILD_TASKS))                    !<-- PD_B_WEST for each child task
            ALLOCATE(T_B_WEST (N)%TASKS(1:NCHILD_TASKS))                    !<-- T_B_WEST for each child task
            ALLOCATE(Q_B_WEST (N)%TASKS(1:NCHILD_TASKS))                    !<-- Q_B_WEST for each child task
            ALLOCATE(CW_B_WEST(N)%TASKS(1:NCHILD_TASKS))                    !<-- CW_B_WEST for each child task
!
            DO NT=1,NCHILD_TASKS
              N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%WEST(NT)+1
              JTE_CHILD_X=CTASK_LIMITS(N)%LIMITS(4,N_TASK)
              J_END_TRANSFER=MIN(JTE_CHILD_X+2,JM_CHILD(N))
!
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further north than the northernmost V in the segment in order to
!***  be able to do 4-pt averaging of PDB onto the V points in order to
!***  do hydrostatic updating of V by the parent.  Now indicate that
!***  reduction of the points to be transferred.
!
!!!           CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)=J_END_TRANSFER     !<-- Wbndry J limit for transfer to child
              CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)=                &  !<-- Wbndry J limit for transfer to child
                CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)-1
!
              IF(CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)==JTE_CHILD_X)        &   !<-- We do not reduce the area for H data
                CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)=             &   !    transfer if the bndry segment reaches
                CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)+1                !    the physical limit of that bndry
!
              NBASE=CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)           &
                   -CHILDTASK_H_SAVE(N)%J_LO_WEST(NT)+1
!
              NBASE_3D=LM*NBASE*N_BLEND_H
              NWORDS  =(3*LM+1)*NBASE*N_BLEND_H                             !<-- # of Wbndry words to transfer from parent to child
              WORDS_BOUND_H_WEST(N)%TASKS(NT)=NWORDS                        !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA(1:NWORDS))      !<-- 1-D bndry data string for child tasks
!                                                                                with west boundary H points
              NLOC_1=1
              NLOC_2=NLOC_1+NBASE*N_BLEND_H-1
!
              NBASE_EXP=CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)                &
                       -CHILDTASK_H_SAVE(N)%J_LO_WEST(NT)+1
!
              NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H+1)-1                   !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!!!!!         PD_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- PD_B_WEST storage location
              ALLOCATE(PD_B_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP))
!
              NLOC_1=NLOC_2+1                                               !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
              NLOC_2=NLOC_1+NBASE_3D-1
              T_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_WEST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              Q_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_WEST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              CW_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_WEST storage location
            ENDDO
!
          ELSE  west_h                                                      !<-- Dummy nonzero length
!
            ALLOCATE(PD_B_WEST(N)%TASKS(1:1))
            ALLOCATE(T_B_WEST (N)%TASKS(1:1))
            ALLOCATE(Q_B_WEST (N)%TASKS(1:1))
            ALLOCATE(CW_B_WEST(N)%TASKS(1:1))
!
          ENDIF west_h
!
          west_v: IF(NUM_TASKS_SEND_V_W(N)>0)THEN                           !<-- Parent task has child west boundary V points?
            NCHILD_TASKS=NUM_TASKS_SEND_V_W(N)
            ALLOCATE(CHILD_BOUND_V_WEST(N)%TASKS(1:NCHILD_TASKS))           !<-- 1-D bndry data string for child tasks with Wbndry V points
            ALLOCATE(WORDS_BOUND_V_WEST(N)%TASKS(1:NCHILD_TASKS))           !<-- # of words in Wbndry V point 1-D data string
!
            ALLOCATE(PD_B_WEST_V(N)%TASKS(1:NCHILD_TASKS))                  !<-- PD_B_WEST_V for each child task
            ALLOCATE(U_B_WEST(N)%TASKS(1:NCHILD_TASKS))                     !<-- U_B_WEST for each child task
            ALLOCATE(V_B_WEST(N)%TASKS(1:NCHILD_TASKS))                     !<-- V_B_WEST for each child task
!
            DO NT=1,NCHILD_TASKS
              NBASE=CHILDTASK_V_SAVE(N)%J_HI_WEST(NT)                   &
                   -CHILDTASK_V_SAVE(N)%J_LO_WEST(NT)+1
              NWORDS=2*LM*NBASE*N_BLEND_V                                   !<-- Total number of V boundary words for child task's segment
              WORDS_BOUND_V_WEST(N)%TASKS(NT)=NWORDS                        !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_V_WEST(N)%TASKS(NT)%DATA(1:NWORDS))      !<-- 1-D bndry data string for child tasks with Wbndry V points
!
              ALLOCATE(PD_B_WEST_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V))
!
              NLOC_1=1
              NLOC_2=NLOC_1+LM*NBASE-1
              U_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_WEST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+LM*NBASE-1
              V_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_WEST storage location
!
            ENDDO
!
          ELSE  west_v                                                      !<-- Dummy nonzero length
!
            ALLOCATE(U_B_WEST(N)%TASKS(1:1))
            ALLOCATE(V_B_WEST(N)%TASKS(1:1))
!
          ENDIF west_v
!
!----------
!***  East
!----------
!
          east_h: IF(NUM_TASKS_SEND_H_E(N)>0)THEN                           !<-- Parent task has child east boundary H points?
            NCHILD_TASKS=NUM_TASKS_SEND_H_E(N)
            ALLOCATE(CHILD_BOUND_H_EAST(N)%TASKS(1:NCHILD_TASKS))           !<-- 1-D bndry data string for child tasks with Ebndry H points
            ALLOCATE(WORDS_BOUND_H_EAST(N)%TASKS(1:NCHILD_TASKS))           !<-- # of words in Ebndry H point 1-D data string
!
            ALLOCATE(PD_B_EAST(N)%TASKS(1:NCHILD_TASKS))                    !<-- PD_B_EAST for each child task
            ALLOCATE(T_B_EAST (N)%TASKS(1:NCHILD_TASKS))                    !<-- T_B_EAST for each child task
            ALLOCATE(Q_B_EAST (N)%TASKS(1:NCHILD_TASKS))                    !<-- Q_B_EAST for each child task
            ALLOCATE(CW_B_EAST(N)%TASKS(1:NCHILD_TASKS))                    !<-- CW_B_EAST for each child task
!
            DO NT=1,NCHILD_TASKS
              N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%EAST(NT)+1
              JTE_CHILD_X=CTASK_LIMITS(N)%LIMITS(4,N_TASK)
              J_END_TRANSFER=MIN(JTE_CHILD_X+2,JM_CHILD(N))
!
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further north than the northernmost V in the segment in order to
!***  be able to do 4-pt averaging of PDB onto the V points in order to
!***  do hydrostatic updating of V by the parent.  Now indicate that
!***  reduction of the points to be transferred.
!
!!!           CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)=J_END_TRANSFER     !<-- Ebndry J limit for transfer to child
              CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)=                &  !<-- Ebndry J limit for transfer to child
                CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)-1
!
              IF(CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)==JTE_CHILD_X)        &   !<-- We do not reduce the area for H data
                CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)=             &   !    transfer if the bndry segment reaches
                CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)+1                !    the physical limit of that bndry
!
              NBASE=CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)           &
                   -CHILDTASK_H_SAVE(N)%J_LO_EAST(NT)+1
!
              NBASE_3D=LM*NBASE*N_BLEND_H
              NWORDS  =(3*LM+1)*NBASE*N_BLEND_H                             !<-- # of Ebndry words to transfer from parent to child
              WORDS_BOUND_H_EAST(N)%TASKS(NT)=NWORDS                        !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA(1:NWORDS))      !<-- 1-D bndry data string for child tasks
!                                                                                with east boundary H points
              NLOC_1=1
              NLOC_2=NLOC_1+NBASE*N_BLEND_H-1
!
              NBASE_EXP=CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)                &
                       -CHILDTASK_H_SAVE(N)%J_LO_EAST(NT)+1
!
              NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H+1)-1                    !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!!!!!         PD_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- PD_B_EAST storage location
              ALLOCATE(PD_B_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP))
!
              NLOC_1=NLOC_2+1                                               !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
              NLOC_2=NLOC_1+NBASE_3D-1
              T_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_EAST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              Q_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_EAST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+NBASE_3D-1
              CW_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_EAST storage location
            ENDDO
!
          ELSE east_h                                                       !<-- Dummy nonzero length
!
            ALLOCATE(PD_B_EAST(N)%TASKS(1:1))
            ALLOCATE(T_B_EAST (N)%TASKS(1:1))
            ALLOCATE(Q_B_EAST (N)%TASKS(1:1))
            ALLOCATE(CW_B_EAST(N)%TASKS(1:1))
!
          ENDIF east_h
!
          east_v: IF(NUM_TASKS_SEND_V_E(N)>0)THEN                           !<-- Parent task has child east boundary V points?
            NCHILD_TASKS=NUM_TASKS_SEND_V_E(N)
            ALLOCATE(CHILD_BOUND_V_EAST(N)%TASKS(1:NCHILD_TASKS))           !<-- 1-D bndry data string for child tasks with Ebndry V points
            ALLOCATE(WORDS_BOUND_V_EAST(N)%TASKS(1:NCHILD_TASKS))           !<-- # of words in Ebndry V point 1-D data string
!
            ALLOCATE(PD_B_EAST_V(N)%TASKS(1:NCHILD_TASKS))                  !<-- PD_B_EAST_V for each child task
            ALLOCATE(U_B_EAST(N)%TASKS(1:NCHILD_TASKS))                     !<-- U_B_EAST for each child task
            ALLOCATE(V_B_EAST(N)%TASKS(1:NCHILD_TASKS))                     !<-- V_B_EAST for each child task
!
            DO NT=1,NCHILD_TASKS
              NBASE=CHILDTASK_V_SAVE(N)%J_HI_EAST(NT)                   &
                   -CHILDTASK_V_SAVE(N)%J_LO_EAST(NT)+1
              NWORDS=2*LM*NBASE*N_BLEND_V                                   !<-- Total number of V boundary words for child task's segment
              WORDS_BOUND_V_EAST(N)%TASKS(NT)=NWORDS                        !<-- Save total number of words
!
              ALLOCATE(CHILD_BOUND_V_EAST(N)%TASKS(NT)%DATA(1:NWORDS))      !<-- 1-D bndry data string for child tasks with Ebndry V points
!
              ALLOCATE(PD_B_EAST_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V))
!
              NLOC_1=1
              NLOC_2=NLOC_1+LM*NBASE-1
              U_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_EAST storage location
!
              NLOC_1=NLOC_2+1
              NLOC_2=NLOC_1+LM*NBASE-1
              V_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_EAST storage location
!
            ENDDO
!
          ELSE  east_v                                                      !<-- Dummy nonzero length
!
            ALLOCATE(U_B_EAST(N)%TASKS(1:1))
            ALLOCATE(V_B_EAST(N)%TASKS(1:1))
!
          ENDIF east_v
!
!-----------------------------------------------------------------------
!
        ENDDO point_loop
!
!-----------------------------------------------------------------------
!***  Set logical flag indicating if parent task holds any child
!***  boundary points for the purpose of sending that data to
!***  pertinent child tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(SEND_CHILD_DATA(1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
!
          IF(NUM_TASKS_SEND_H_S(N)>0.OR.                                &
             NUM_TASKS_SEND_H_N(N)>0.OR.                                &
             NUM_TASKS_SEND_H_W(N)>0.OR.                                &
             NUM_TASKS_SEND_H_E(N)>0.OR.                                &
             NUM_TASKS_SEND_V_S(N)>0.OR.                                &
             NUM_TASKS_SEND_V_N(N)>0.OR.                                &
             NUM_TASKS_SEND_V_W(N)>0.OR.                                &
             NUM_TASKS_SEND_V_E(N)>0)THEN
!
            SEND_CHILD_DATA(N)=.TRUE.
!
          ELSE
            SEND_CHILD_DATA(N)=.FALSE.
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate the handles to be used by the parent task when it ISends
!***  data to the individual child tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(HANDLE_H_SOUTH(1:NUM_CHILDREN))
        ALLOCATE(HANDLE_H_NORTH(1:NUM_CHILDREN))
        ALLOCATE(HANDLE_H_WEST (1:NUM_CHILDREN))
        ALLOCATE(HANDLE_H_EAST (1:NUM_CHILDREN))
!
        ALLOCATE(HANDLE_V_SOUTH(1:NUM_CHILDREN))
        ALLOCATE(HANDLE_V_NORTH(1:NUM_CHILDREN))
        ALLOCATE(HANDLE_V_WEST (1:NUM_CHILDREN))
        ALLOCATE(HANDLE_V_EAST (1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
!
!-------------------------------
!***  For child boundary, south
!-------------------------------
!
          IF(NUM_TASKS_SEND_H_S(N)>0)THEN
            ALLOCATE(HANDLE_H_SOUTH(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_S(N)))
            DO NN=1,NUM_TASKS_SEND_H_S(N)
              HANDLE_H_SOUTH(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
          IF(NUM_TASKS_SEND_V_S(N)>0)THEN
            ALLOCATE(HANDLE_V_SOUTH(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_S(N)))
            DO NN=1,NUM_TASKS_SEND_V_S(N)
              HANDLE_V_SOUTH(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
!-------------------------------
!***  For child boundary, north
!-------------------------------
!
          IF(NUM_TASKS_SEND_H_N(N)>0)THEN
            ALLOCATE(HANDLE_H_NORTH(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_N(N)))
            DO NN=1,NUM_TASKS_SEND_H_N(N)
              HANDLE_H_NORTH(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
          IF(NUM_TASKS_SEND_V_N(N)>0)THEN
            ALLOCATE(HANDLE_V_NORTH(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_N(N)))
            DO NN=1,NUM_TASKS_SEND_V_N(N)
              HANDLE_V_NORTH(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
!------------------------------
!***  For child boundary, west
!------------------------------
!
          IF(NUM_TASKS_SEND_H_W(N)>0)THEN
            ALLOCATE(HANDLE_H_WEST(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_W(N)))
            DO NN=1,NUM_TASKS_SEND_H_W(N)
              HANDLE_H_WEST(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
          IF(NUM_TASKS_SEND_V_W(N)>0)THEN
            ALLOCATE(HANDLE_V_WEST(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_W(N)))
            DO NN=1,NUM_TASKS_SEND_V_W(N)
              HANDLE_V_WEST(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
!------------------------------
!***  For child boundary, east
!------------------------------
!
          IF(NUM_TASKS_SEND_H_E(N)>0)THEN
            ALLOCATE(HANDLE_H_EAST(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_E(N)))
            DO NN=1,NUM_TASKS_SEND_H_E(N)
              HANDLE_H_EAST(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
          IF(NUM_TASKS_SEND_V_E(N)>0)THEN
            ALLOCATE(HANDLE_V_EAST(N)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_E(N)))
            DO NN=1,NUM_TASKS_SEND_V_E(N)
              HANDLE_V_EAST(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        DEALLOCATE(CF)
!
!-----------------------------------------------------------------------
!
      ENDIF child_block
!
!-----------------------------------------------------------------------
!***  Now the children need to obtain a few pieces of information
!***  from their parents in order to be prepared to receive boundary
!***  data from parent tasks during the integration.  With that
!***  information the separate boundary data arrays for each child's
!***  segment of the boundary are allocated.  Those arrays are
!***  loaded into the Parent-Child coupler export state.  
!
!***  In addition the parents need to obtain the sfc geopotential
!***  of the child boundary points they cover in order to maintain
!***  balance when parent data is interpolated to child boundaries
!***  where the terrain is different.
!
!***  All domains need to call this prep routine due to the exchange
!***  of the information required between parents and children.
!-----------------------------------------------------------------------
!
      CALL PRELIM_CHILD_INFO(IMP_STATE,EXP_STATE)
!
!-----------------------------------------------------------------------
!
      cpl1_prelim_tim=0.
      cpl1_south_h_tim=0.
      cpl1_south_v_tim=0.
      cpl1_north_h_tim=0.
      cpl1_north_v_tim=0.
      cpl1_west_h_tim=0.
      cpl1_west_v_tim=0.
      cpl1_east_h_tim=0.
      cpl1_east_v_tim=0.
      cpl1_recv_tim=0.
!
      cpl1_south_h_recv_tim=0.
      cpl1_south_h_undo_tim=0.
      cpl1_south_h_exp_tim=0.
      cpl1_south_v_recv_tim=0.
      cpl1_south_v_undo_tim=0.
      cpl1_south_v_exp_tim=0.
!
      cpl2_comp_tim=0.
      cpl2_wait_tim=0.
      cpl2_send_tim=0.
!
!-----------------------------------------------------------------------
!
      IF(RC_CPL_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)"PARENT_CHILD_CPL INITIALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"PARENT_CHILD_CPL INITIALIZE STEP FAILED"
      ENDIF
!
      IF(PRESENT(RC_FINAL))THEN
        RC_FINAL=RC_CPL_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_CPL_RUN_RECV(CPL_COMP                     &
                                          ,IMP_STATE                    &
                                          ,EXP_STATE                    &
                                          ,CLOCK                        &
                                          ,RC_FINAL)
!
!-----------------------------------------------------------------------
!***  RUN THE COUPLER STEP WHERE CHILDREN RECEIVE DATA FROM PARENTS.
!***  THIS IS PHASE 1 OF THE COUPLER SINCE IT OCCURS AT THE BEGINNING
!***  OF THE TIMESTEPS.  ONLY CHILD TASKS ENTER THE ROUTINE.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Parent-Child Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,INTENT(OUT)     :: RC_FINAL
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: IERR,N,NP_H,NP_V,NTAG,NTIMESTEP
      INTEGER :: RC,RC_CPL_RUN
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      integer,dimension(8) :: values
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!     call date_and_time(values=values)
!     write(0,12345)values(5),values(6),values(7),values(8)
12345 format('Enter coupler receive at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
      btim =timef()
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_FINAL  =ESMF_SUCCESS
      RC_CPL_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  EACH CHILD TASK THAT HOLDS PART OF THE DOMAIN BOUNDARY RECEIVES
!***  DATA FROM THE PARENT TASKS WITH WHICH IT SHARES THAT BOUNDARY. 
!***  THE PARENT IS SENDING THE DATA FROM THE END OF ITS TIMESTEP. 
!***  THE CHILD WILL RECEIVE THE DATA FROM THE FUTURE, COMPUTE THE 
!***  TIME TENDENCIES FOR THE BOUNDARY VARIABLES AND THEN PROCEED AS
!***  ANY STAND ALONE DOMAIN WOULD UNTIL IT REACHES THE SAME POINT
!***  IN TIME FROM WHICH THE PARENT SENT THE DATA.  THE BOUNDARY
!***  DATA EXCHANGE PROCESS THEN REPEATS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  NOW FOR EACH SIDE OF THE NESTS' BOUNDARIES:
!
!***   (a) Each child task receives all boundary data from the
!***       relevant parent task(s).  Note that more than one 
!***       parent task might send to each child task and there may
!***       be overlap due to haloes.
!***   (b) The child tasks separate the data received from each of
!***       the parent tasks and combines them into unified segments
!***       on the boundary for each variable.
!***   (c) All boundary data is loaded into the Parent-Child Coupler's
!***       export state.
!
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet(clock       =CLOCK                             &
                        ,advanceCount=NTIMESTEP_ESMF                    &
                        ,rc          =RC)
!
      NTIMESTEP=NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!
      cpl1_prelim_tim=cpl1_prelim_tim+(timef()-btim0)
!
!--------------------
!***  South H Points
!--------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%SOUTH                                ! # of parent tasks sending south boundary H data
!
      IF(NP_H>0)THEN
!
        NTAG=NTIMESTEP+101                                                 !<-- Add 101 to obtain a unique south H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Sboundary H data
!         call date_and_time(values=values)     
!         write(0,123)n,parent_task(n)%south_h%id_source,values(5),values(6),values(7),values(8)
! 123     format(' Ready to recv South_H from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)

          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%SOUTH_H%STRING                   &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%SOUTH_H%LENGTH                   &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%SOUTH_H%ID_SOURCE                &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_south_h_recv_tim=cpl1_south_h_recv_tim+(timef()-btim)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!         write(0,*)' south_h recv time=',((timef()-btim))*1.e-3

!         call date_and_time(values=values)     
!         write(0,124)n,parent_task(n)%south_h%id_source,values(5),values(6),values(7),values(8)
! 124     format(' Recvd South_H from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%SOUTH_H%LENGTH     &  !<-- Length of parent datastring
                                     ,DATASTRING =PARENT_TASK(N)%SOUTH_H%STRING     &  !<-- Parent datastring of child task bndry segment 
                                     ,ILIM_LO    =INDX_MIN_H%SOUTH                  &  !<-- Lower I limit of child's segment of boundary
                                     ,ILIM_HI    =INDX_MAX_H%SOUTH                  &  !<-- Upper I limit of child's segment of boundary
                                     ,JLIM_LO    =1                                 &  !<-- Lower J limit of child's segment of boundary
                                     ,JLIM_HI    =N_BLEND_H                         &  !<-- Upper J limit of child's segment of boundary
                                     ,I_START    =PARENT_TASK(N)%SOUTH_H%INDX_START &  !<-- Child's segment Istart on each parent task
                                     ,I_END      =PARENT_TASK(N)%SOUTH_H%INDX_END   &  !<-- Child's segment Iend on each parent task
                                     ,J_START    =1                                 &  !<-- Child's segment Jstart on each parent task
                                     ,J_END      =N_BLEND_H                         &  !<-- Child's segment Jend on each parent task
                                     ,PDB        =PDB_S                             &  !<-- Child's 1-D segment of PD on this side of bndry
                                     ,TB         =TB_S                              &  !<-- Child's 1-D segment of T on this side of bndry
                                     ,QB         =QB_S                              &  !<-- Child's 1-D segment of Q on this side of bndry
                                     ,CWB        =CWB_S )                              !<-- Child's 1-D segment of CW on this side of bndry
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_H%SOUTH )             !<-- Length of the 1-D segment in I and J
!
          cpl1_south_h_undo_tim=cpl1_south_h_undo_tim+(timef()-btim)
!         write(0,*)' south_h recv time=',((timef()-btim))*1.e-3
!
!         call date_and_time(values=values)     
!         write(0,224)n,parent_task(n)%south_h%id_source,values(5),values(6),values(7),values(8)
! 224     format(' After South_H Recv data_from_string task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
        ENDDO
!
        btim=timef()
        CALL EXPORT_CHILD_BOUNDARY(PDB         =PDB_S                               &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,TB          =TB_S                                &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,QB          =QB_S                                &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,CWB         =CWB_S                               &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ILIM_LO     =INDX_MIN_H%SOUTH                    &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =INDX_MAX_H%SOUTH                    &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =1                                   &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =N_BLEND_H                           &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='SOUTH_H'                           &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_SOUTH_H                    &  !<-- Combined boundary segment H data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
          cpl1_south_h_exp_tim=cpl1_south_h_exp_tim+(timef()-btim)

!         call date_and_time(values=values)     
!         write(0,225)values(5),values(6),values(7),values(8)
! 225     format(' After South_H Recv export at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
      ENDIF
!
      cpl1_south_h_tim=cpl1_south_h_tim+(timef()-btim0)
!
!--------------------
!***  South V Points
!--------------------
!
      btim0=timef()
!
      NP_V=NUM_PARENT_TASKS_SENDING_V%SOUTH                                ! # of parent tasks sending south boundary V data
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+102                                                 !<-- Add 102 to obtain a unique south H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Sboundary V data
!         call date_and_time(values=values)     
!         write(0,125)n,parent_task(n)%south_v%id_source,values(5),values(6),values(7),values(8)
! 125     format(' Ready to recv South_V from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%SOUTH_V%STRING                   &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%SOUTH_V%LENGTH                   &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%SOUTH_V%ID_SOURCE                &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_south_v_recv_tim=cpl1_south_v_recv_tim+(timef()-btim)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!         write(0,*)' south_v recv time=',((timef()-btim))*1.e-3
!
!         call date_and_time(values=values)     
!         write(0,126)n,parent_task(n)%south_v%id_source,values(5),values(6),values(7),values(8)
! 126     format(' Recvd South_V from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%SOUTH_V%LENGTH     &  !<-- Length of parent datastring
                                     ,DATASTRING =PARENT_TASK(N)%SOUTH_V%STRING     &  !<-- Parent datastring of child task bndry segment
                                     ,ILIM_LO    =INDX_MIN_V%SOUTH                  &  !<-- Lower I limit of child's segment of boundary
                                     ,ILIM_HI    =INDX_MAX_V%SOUTH                  &  !<-- Upper I limit of child's segment of boundary
                                     ,JLIM_LO    =1                                 &  !<-- Lower J limit of child's segment of boundary
                                     ,JLIM_HI    =N_BLEND_V                         &  !<-- Upper J limit of child's segment of boundary
                                     ,I_START    =PARENT_TASK(N)%SOUTH_V%INDX_START &  !<-- Child's segment Istart on each parent task
                                     ,I_END      =PARENT_TASK(N)%SOUTH_V%INDX_END   &  !<-- Child's segment Iend on each parent task
                                     ,J_START    =1                                 &  !<-- Child's segment Jstart on each parent task
                                     ,J_END      =N_BLEND_V                         &  !<-- Child's segment Jend on each parent task
                                     ,UB         =UB_S                              &  !<-- Child's 1-D segment of U on this side of bndry
                                     ,VB         =VB_S )                               !<-- Child's 1-D segment of V on this side of bndry
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_V%SOUTH )             !<-- Length of the 1-D segment in I and J
          cpl1_south_v_undo_tim=cpl1_south_v_undo_tim+(timef()-btim)

!         call date_and_time(values=values)     
!         write(0,324)n,parent_task(n)%south_v%id_source,values(5),values(6),values(7),values(8)
! 324     format(' After South_V Recv data_from_string task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO
!
        btim=timef()
        CALL EXPORT_CHILD_BOUNDARY(UB          =UB_S                                &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,VB          =VB_S                                &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ILIM_LO     =INDX_MIN_V%SOUTH                    &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =INDX_MAX_V%SOUTH                    &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =1                                   &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =N_BLEND_V                           &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='SOUTH_V'                           &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_SOUTH_V                    &  !<-- Combined boundary segment V data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
        cpl1_south_v_exp_tim=cpl1_south_v_exp_tim+(timef()-btim)
!
!       call date_and_time(values=values)     
!       write(0,325)values(5),values(6),values(7),values(8)
! 325   format(' After South_V Recv export at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
      ENDIF
!
      cpl1_south_v_tim=cpl1_south_v_tim+(timef()-btim0)
!
!--------------------
!***  North H Points
!--------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%NORTH                                ! # of parent tasks sending north boundary H data
!
      IF(NP_H>0)THEN
        NTAG=NTIMESTEP+103                                                 !<-- Add 103 to obtain a unique north H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Nboundary H data
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%NORTH_H%STRING                   &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%NORTH_H%LENGTH                   &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%NORTH_H%ID_SOURCE                &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%NORTH_H%LENGTH     &
                                     ,DATASTRING =PARENT_TASK(N)%NORTH_H%STRING     & 
                                     ,ILIM_LO    =INDX_MIN_H%NORTH                  &
                                     ,ILIM_HI    =INDX_MAX_H%NORTH                  &
                                     ,JLIM_LO    =1                                 &
                                     ,JLIM_HI    =N_BLEND_H                         &
                                     ,I_START    =PARENT_TASK(N)%NORTH_H%INDX_START &
                                     ,I_END      =PARENT_TASK(N)%NORTH_H%INDX_END   &
                                     ,J_START    =1                                 &
                                     ,J_END      =N_BLEND_H                         &
                                     ,PDB        =PDB_N                             &
                                     ,TB         =TB_N                              &
                                     ,QB         =QB_N                              &
                                     ,CWB        =CWB_N )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_H%NORTH )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(PDB         =PDB_N                               &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,TB          =TB_N                                &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,QB          =QB_N                                &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,CWB         =CWB_N                               &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ILIM_LO     =INDX_MIN_H%NORTH                    &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =INDX_MAX_H%NORTH                    &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =1                                   &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =N_BLEND_H                           &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='NORTH_H'                           &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_NORTH_H                    &  !<-- Combined boundary segment H data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_north_h_tim=cpl1_north_h_tim+(timef()-btim0)
!
!--------------------
!***  North V Points
!--------------------
!
      btim0=timef()
!
      NP_V=NUM_PARENT_TASKS_SENDING_V%NORTH                                ! # of parent tasks sending north boundary V data
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+104                                                 !<-- Add 104 to obtain a unique north H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Nboundary V data
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%NORTH_V%STRING                   &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%NORTH_V%LENGTH                   &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%NORTH_V%ID_SOURCE                &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%NORTH_V%LENGTH     &
                                     ,DATASTRING =PARENT_TASK(N)%NORTH_V%STRING     & 
                                     ,ILIM_LO    =INDX_MIN_V%NORTH                  &
                                     ,ILIM_HI    =INDX_MAX_V%NORTH                  &
                                     ,JLIM_LO    =1                                 &
                                     ,JLIM_HI    =N_BLEND_V                         &
                                     ,I_START    =PARENT_TASK(N)%NORTH_V%INDX_START &
                                     ,I_END      =PARENT_TASK(N)%NORTH_V%INDX_END   &
                                     ,J_START    =1                                 &
                                     ,J_END      =N_BLEND_V                         &
                                     ,UB         =UB_N                              &
                                     ,VB         =VB_N )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_V%NORTH )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(UB          =UB_N                                &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,VB          =VB_N                                &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ILIM_LO     =INDX_MIN_V%NORTH                    &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =INDX_MAX_V%NORTH                    &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =1                                   &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =N_BLEND_V                           &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='NORTH_V'                           &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_NORTH_V                    &  !<-- Combined boundary segment V data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_north_v_tim=cpl1_north_v_tim+(timef()-btim0)
!
!-------------------
!***  West H Points
!-------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%WEST                                 ! # of parent tasks sending west boundary H data
!
      IF(NP_H>0)THEN
        NTAG=NTIMESTEP+105                                                 !<-- Add 105 to obtain a unique west H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Wboundary H data
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%WEST_H%STRING                    &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%WEST_H%LENGTH                    &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%WEST_H%ID_SOURCE                 &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%WEST_H%LENGTH      &
                                     ,DATASTRING =PARENT_TASK(N)%WEST_H%STRING      & 
                                     ,ILIM_LO    =1                                 &
                                     ,ILIM_HI    =N_BLEND_H                         &
                                     ,JLIM_LO    =INDX_MIN_H%WEST                   &
                                     ,JLIM_HI    =INDX_MAX_H%WEST                   &
                                     ,I_START    =1                                 &
                                     ,I_END      =N_BLEND_H                         &
                                     ,J_START    =PARENT_TASK(N)%WEST_H%INDX_START  &
                                     ,J_END      =PARENT_TASK(N)%WEST_H%INDX_END    &
                                     ,PDB        =PDB_W                             &
                                     ,TB         =TB_W                              &
                                     ,QB         =QB_W                              &
                                     ,CWB        =CWB_W )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_H%WEST )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(PDB         =PDB_W                               &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,TB          =TB_W                                &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,QB          =QB_W                                &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,CWB         =CWB_W                               &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ILIM_LO     =1                                   &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =N_BLEND_H                           &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =INDX_MIN_H%WEST                     &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =INDX_MAX_H%WEST                     &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='WEST_H'                            &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_WEST_H                     &  !<-- Combined boundary segment H data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_west_h_tim=cpl1_west_h_tim+(timef()-btim0)
!
!-------------------
!***  West V Points
!-------------------
!
      btim0=timef()
!
      NP_V=NUM_PARENT_TASKS_SENDING_V%WEST                                 ! # of parent tasks sending west boundary V data
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+106                                                 !<-- Add 106 to obtain a unique west H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Sboundary V data
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%WEST_V%STRING                    &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%WEST_V%LENGTH                    &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%WEST_V%ID_SOURCE                 &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%WEST_V%LENGTH      &
                                     ,DATASTRING =PARENT_TASK(N)%WEST_V%STRING      & 
                                     ,ILIM_LO    =1                                 &
                                     ,ILIM_HI    =N_BLEND_V                         &
                                     ,JLIM_LO    =INDX_MIN_V%WEST                   &
                                     ,JLIM_HI    =INDX_MAX_V%WEST                   &
                                     ,I_START    =1                                 &
                                     ,I_END      =N_BLEND_V                         &
                                     ,J_START    =PARENT_TASK(N)%WEST_V%INDX_START  &
                                     ,J_END      =PARENT_TASK(N)%WEST_V%INDX_END    &
                                     ,UB         =UB_W                              &
                                     ,VB         =VB_W )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_V%WEST )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(UB          =UB_W                                &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,VB          =VB_W                                &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ILIM_LO     =1                                   &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =N_BLEND_V                           &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =INDX_MIN_V%WEST                     &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =INDX_MAX_V%WEST                     &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='WEST_V'                            &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_WEST_V                     &  !<-- Combined boundary segment V data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_west_v_tim=cpl1_west_v_tim+(timef()-btim0)
!
!-------------------
!***  East H Points
!-------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%EAST                                 ! # of parent tasks sending east boundary H data
!
      IF(NP_H>0)THEN
        NTAG=NTIMESTEP+107                                                 !<-- Add 107 to obtain a unique east H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Eboundary H data
          btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%EAST_H%STRING                    &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%EAST_H%LENGTH                    &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%EAST_H%ID_SOURCE                 &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
          cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%EAST_H%LENGTH      &
                                     ,DATASTRING =PARENT_TASK(N)%EAST_H%STRING      & 
                                     ,ILIM_LO    =1                                 &
                                     ,ILIM_HI    =N_BLEND_H                         &
                                     ,JLIM_LO    =INDX_MIN_H%EAST                   &
                                     ,JLIM_HI    =INDX_MAX_H%EAST                   &
                                     ,I_START    =1                                 &
                                     ,I_END      =N_BLEND_H                         &
                                     ,J_START    =PARENT_TASK(N)%EAST_H%INDX_START  &
                                     ,J_END      =PARENT_TASK(N)%EAST_H%INDX_END    &
                                     ,PDB        =PDB_E                             &
                                     ,TB         =TB_E                              &
                                     ,QB         =QB_E                              &
                                     ,CWB        =CWB_E )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_H%EAST )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(PDB         =PDB_E                               &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,TB          =TB_E                                &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,QB          =QB_E                                &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,CWB         =CWB_E                               &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ILIM_LO     =1                                   &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =N_BLEND_H                           &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =INDX_MIN_H%EAST                     &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =INDX_MAX_H%EAST                     &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='EAST_H'                            &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_EAST_H                     &  !<-- Combined boundary segment H data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_east_h_tim=cpl1_east_h_tim+(timef()-btim0)
!
!-------------------
!***  East V Points
!-------------------
!
      btim0=timef()
!
      NP_V=NUM_PARENT_TASKS_SENDING_V%EAST                                 ! # of parent tasks sending east boundary V data
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+108                                                 !<-- Add 108 to obtain a unique east H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Eboundary V data
      btim=timef()
          CALL MPI_RECV(PARENT_TASK(N)%EAST_V%STRING                    &  !<-- 1-D boundary datastring from parent task
                       ,PARENT_TASK(N)%EAST_V%LENGTH                    &  !<-- # of words in the datastring
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,PARENT_TASK(N)%EAST_V%ID_SOURCE                 &  !<-- Local rank of the parent task sending the datastring
                       ,NTAG                                            &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
      cpl1_recv_tim=cpl1_recv_tim+(timef()-btim)
!
          CALL CHILD_DATA_FROM_STRING(LENGTH_DATA=PARENT_TASK(N)%EAST_V%LENGTH      &
                                     ,DATASTRING =PARENT_TASK(N)%EAST_V%STRING      & 
                                     ,ILIM_LO    =1                                 &
                                     ,ILIM_HI    =N_BLEND_V                         &
                                     ,JLIM_LO    =INDX_MIN_V%EAST                   &
                                     ,JLIM_HI    =INDX_MAX_V%EAST                   &
                                     ,I_START    =1                                 &
                                     ,I_END      =N_BLEND_V                         &
                                     ,J_START    =PARENT_TASK(N)%EAST_V%INDX_START  &
                                     ,J_END      =PARENT_TASK(N)%EAST_V%INDX_END    &
                                     ,UB         =UB_E                              &
                                     ,VB         =VB_E )
!!!                                  ,LEN_2D     =LENGTH_BND_SEG_V%EAST )
!
        ENDDO
!
        CALL EXPORT_CHILD_BOUNDARY(UB          =UB_E                                &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,VB          =VB_E                                &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ILIM_LO     =1                                   &  !<-- Lower I limit of child's segment of boundary
                                  ,ILIM_HI     =N_BLEND_V                           &  !<-- Upper I limit of child's segment of boundary
                                  ,JLIM_LO     =INDX_MIN_V%EAST                     &  !<-- Lower J limit of child's segment of boundary
                                  ,JLIM_HI     =INDX_MAX_V%EAST                     &  !<-- Upper J limit of child's segment of boundary
                                  ,DATA_NAME   ='EAST_V'                            &  !<-- Name attached to the combined exported data
                                  ,DATA_EXP    =BOUND_1D_EAST_V                     &  !<-- Combined boundary segment V data for child task
                                  ,EXPORT_STATE=EXP_STATE )                            !<-- The Parent-Child Coupler export state
!
      ENDIF
!
      cpl1_east_v_tim=cpl1_east_v_tim+(timef()-btim0)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Clocktime for Recv in Phase1 into Cpl Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='Cpl1_Recv_Time'                     &  !<-- Name of the attribute to insert
                            ,value=cpl1_recv_tim                        &  !<-- Phase 1 Recv time
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_RUN_RECV
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_CPL_RUN_SEND(CPL_COMP                     &
                                          ,IMP_STATE                    &
                                          ,EXP_STATE                    &
                                          ,CLOCK                        &
                                          ,RC_FINAL)
!
!-----------------------------------------------------------------------
!***  RUN THE COUPLER STEP WHERE PARENTS SEND DATA TO THEIR CHILDREN.
!***  THIS IS PHASE 2 OF THE COUPLER SINCE IT OCCURS AT THE END OF
!***  THE TIMESTEPS.  THE CHILDREN RECEIVE THE DATA IN PHASE 1 AT
!***  THE BEGINNING OF THEIR APPROPRIATE TIMESTEPS.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                        !<-- The Parent-Child Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                       !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                       !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                           !<-- The ATM Driver Clock for this parent domain
!
      INTEGER,OPTIONAL,INTENT(OUT)     :: RC_FINAL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
!     type(boundary_data),dimension(:),pointer :: ctasks_work_s         &
!                                                ,ctasks_work_n         &
!                                                ,ctasks_work_w         &
!                                                ,ctasks_work_e
!
      INTEGER :: ID_DOM,IERR,N,NT,NTAG,NUM_CHILD_TASKS
      INTEGER :: RC,RC_CPL_RUN
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: ISTAT
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: PD_V
!
      integer :: nnnn
      integer,dimension(8) :: values
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!     call date_and_time(values=values)
!     write(0,12345)values(5),values(6),values(7),values(8)
12345 format('Enter parent send at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_FINAL  =ESMF_SUCCESS
      RC_CPL_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  PARENT GENERATES BOUNDARY DATA FOR ITS CHILDREN. 
!-----------------------------------------------------------------------
!
      child_loop: DO N=1,NUM_CHILDREN
!
!-----------------------------------------------------------------------
!***  FIRST THE PARENTS COMPUTE THE NEW SURFACE PRESSURE ON THE
!***  NESTS' BOUNDARY POINTS.
!-----------------------------------------------------------------------
!
        NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N))
!
!--------
!***  PD
!--------
!
        btim=timef()
        CALL PARENT_UPDATE_CHILD_PSFC(FIS,PD,T,TRACERS(:,:,:,INDX_Q)    &  !<-- Native parent values
                                     ,PT,PDTOP                          &  !<-- Domain PT and PDTOP
                                     ,SG1,SG2                           &  !<-- General vertical structure (shared by all domains)
                                     ,IMS,IME,JMS,JME                   &  !<-- Parent task subdomain lateral memory dimensions
                                     ,LM                                &  !<-- # of model layers
!
                                     ,NUM_CHILD_TASKS                   &  !<-- # of fcst tasks on child N
                                     ,CTASK_LIMITS(N)%LIMITS            &  !<-- Integration limits on each task of child N
!
                                     ,FIS_CHILD_SOUTH(N)%TASKS          &  !<-- Sfc geopotential on Sbndry points on child tasks
                                     ,FIS_CHILD_NORTH(N)%TASKS          &  !<-- Sfc geopotential on Nbndry points on child tasks
                                     ,FIS_CHILD_WEST(N)%TASKS           &  !<-- Sfc geopotential on Wbndry points on child tasks
                                     ,FIS_CHILD_EAST(N)%TASKS           &  !<-- Sfc geopotential on Ebndry points on child tasks
! 
                                     ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH  &  !<-- Ranks of child N's fcst tasks on its Sbndry
                                     ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH  &  !<-- Ranks of child N's fcst tasks on its Nbndry
                                     ,CHILDTASK_BNDRY_H_RANKS(N)%WEST   &  !<-- Ranks of child N's fcst tasks on its Wbndry
                                     ,CHILDTASK_BNDRY_H_RANKS(N)%EAST   &  !<-- Ranks of child N's fcst tasks on its Ebndry
!
                                     ,NUM_TASKS_SEND_H_S(N)             &  !<-- # of child tasks with south boundary segments
                                     ,NUM_TASKS_SEND_H_N(N)             &  !<-- # of child tasks with north boundary segments
                                     ,NUM_TASKS_SEND_H_W(N)             &  !<-- # of child tasks with west boundary segments
                                     ,NUM_TASKS_SEND_H_E(N)             &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_H(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child Sbndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child Nbndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child Wbndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child Ebndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child Sbndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child Nbndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child Wbndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child Ebndry point
!
                               ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Ending J on each east boundary child task
!
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND &  !    on each side of the child boundary.
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND &  !
!
                                     ,N_BLEND_H                         &  !<-- Width of boundary blending region for mass points
                                     ,IM_CHILD(N)                       &  !<-- East-west points on child domain
                                     ,JM_CHILD(N)                       &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                     ,CHILD_BOUND_H_SOUTH(N)%TASKS      &  !<-- 1-D H-pt Sbndry datastring to be sent by parent to child
                                     ,CHILD_BOUND_H_NORTH(N)%TASKS      &  !<-- 1-D H-pt Nbndry datastring to be sent by parent to child
                                     ,CHILD_BOUND_H_WEST(N)%TASKS       &  !<-- 1-D H-pt Wbndry datastring to be sent by parent to child
                                     ,CHILD_BOUND_H_EAST(N)%TASKS       &  !<-- 1-D H-pt Ebndry datastring to be sent by parent to child
!
                                     ,PD_B_SOUTH(N)%TASKS               &  !<-- Updated sigma domain pressure (Pa) on nest bndry points
                                     ,PD_B_NORTH(N)%TASKS               &  !    for all four sides of nest N's boundary.
                                     ,PD_B_WEST(N)%TASKS                &  !
                                     ,PD_B_EAST(N)%TASKS )                 !
!
!-----------------------------------------------------------------------
!***  NOW COMPUTE THE NEW MASS VARIABLE VALUES IN THE COLUMNS ABOVE
!***  THE NEST BOUNDARY POINTS.
!-----------------------------------------------------------------------
!
!-----------------
!***  Temperature
!-----------------
!
        CALL PARENT_UPDATE_CHILD_BNDRY(T                                &  !<-- Parent sensible temperature (K)
                                      ,PD,PT,PDTOP                      &  !<-- Parent PD; domain PT and PDTOP
                                      ,PSGML1,SGML2,SG1,SG2             &  !<-- General vertical structure (shared by all domains)
!
                                      ,PD_B_SOUTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                      ,PD_B_NORTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                      ,PD_B_WEST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                      ,PD_B_EAST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task subdomain lateral memory dimensions
                                      ,LM                               &  !<-- # of model layers
                                      ,0                                &  !<-- # of rows to ignore on north/east nest boundaries
!
                                      ,NUM_TASKS_SEND_H_S(N)            &  !<-- # of child tasks with south boundary segments
                                      ,NUM_TASKS_SEND_H_N(N)            &  !<-- # of child tasks with north boundary segments
                                      ,NUM_TASKS_SEND_H_W(N)            &  !<-- # of child tasks with west boundary segments
                                      ,NUM_TASKS_SEND_H_E(N)            &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_H(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child E bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child E bndry point
!
                               ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Ending I for transfer to child on each Sbndry child task
                               ,CHILDTASK_H_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Ending I for transfer to child on each Nbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Ending J for transfer to child on each Wbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point.
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND &  !
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND &  !<--
!
                                      ,N_BLEND_H                        &  !<-- Width of boundary blending region
                                      ,IM_CHILD(N)                      &  !<-- East-west points on child domain
                                      ,JM_CHILD(N)                      &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                      ,T_B_SOUTH(N)%TASKS               &  !<-- Updated sensible temperature (K) on nest bndry points.
                                      ,T_B_NORTH(N)%TASKS               &  !
                                      ,T_B_WEST(N)%TASKS                &  !
                                      ,T_B_EAST(N)%TASKS )                 !<-- 
!
!-----------------
!***  Specific humidity
!-----------------
!
        CALL PARENT_UPDATE_CHILD_BNDRY(TRACERS(:,:,:,INDX_Q)            &  !<-- Parent specific humidity (kg/kg) 
                                      ,PD,PT,PDTOP                      &  !<-- Parent PD; domain PT and PDTOP
                                      ,PSGML1,SGML2,SG1,SG2             &  !<-- General vertical structure (shared by all domains)
!
                                      ,PD_B_SOUTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                      ,PD_B_NORTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                      ,PD_B_WEST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                      ,PD_B_EAST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task subdomain lateral memory dimensions
                                      ,LM                               &  !<-- # of model layers
                                      ,0                                &  !<-- # of rows to ignore on north/east nest boundaries
!
                                      ,NUM_TASKS_SEND_H_S(N)            &  !<-- # of child tasks with south boundary segments
                                      ,NUM_TASKS_SEND_H_N(N)            &  !<-- # of child tasks with north boundary segments
                                      ,NUM_TASKS_SEND_H_W(N)            &  !<-- # of child tasks with west boundary segments
                                      ,NUM_TASKS_SEND_H_E(N)            &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_H(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child E bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child E bndry point
!
                               ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Ending I for transfer to child on each Sbndry child task
                               ,CHILDTASK_H_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Ending I for transfer to child on each Nbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Ending J for transfer to child on each Wbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point.
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND &  !
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND &  !<--
!
                                      ,N_BLEND_H                        &  !<-- Width of boundary blending region
                                      ,IM_CHILD(N)                      &  !<-- East-west points on child domain
                                      ,JM_CHILD(N)                      &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                      ,Q_B_SOUTH(N)%TASKS               &  !<-- Updated specific humidity (kg/kg) on nest bndry points.
                                      ,Q_B_NORTH(N)%TASKS               &  !
                                      ,Q_B_WEST(N)%TASKS                &  !
                                      ,Q_B_EAST(N)%TASKS )                 !<-- 
!
!-----------------------------------------------------------------------
!
!----------------------
!***  Cloud Condensate 
!----------------------
!
        CALL PARENT_UPDATE_CHILD_BNDRY(TRACERS(:,:,:,INDX_CW)           &  !<-- Parent cloud condensate (kg/kg)
                                      ,PD,PT,PDTOP                      &  !<-- Parent PD; domain PT and PDTOP
                                      ,PSGML1,SGML2,SG1,SG2             &  !<-- General vertical structure (shared by all domains)
!
                                      ,PD_B_SOUTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                      ,PD_B_NORTH(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                      ,PD_B_WEST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                      ,PD_B_EAST(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task subdomain lateral memory dimensions
                                      ,LM                               &  !<-- # of model layers
                                      ,0                                &  !<-- # of rows to ignore on north/east nest boundaries
!
                                      ,NUM_TASKS_SEND_H_S(N)            &  !<-- # of child tasks with south boundary segments
                                      ,NUM_TASKS_SEND_H_N(N)            &  !<-- # of child tasks with north boundary segments
                                      ,NUM_TASKS_SEND_H_W(N)            &  !<-- # of child tasks with west boundary segments
                                      ,NUM_TASKS_SEND_H_E(N)            &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_H(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child E bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child S bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child N bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child W bndry point
                                     ,PARENT_4_INDICES_H(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child E bndry point
!
                               ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Ending I for transfer to child on each Sbndry child task
                               ,CHILDTASK_H_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Ending I for transfer to child on each Nbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Ending J for transfer to child on each Wbndry child task
                               ,CHILDTASK_H_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point.
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND &  !
                                    ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND &  !<--
!
                                      ,N_BLEND_H                        &  !<-- Width of boundary blending region
                                      ,IM_CHILD(N)                      &  !<-- East-west points on child domain
                                      ,JM_CHILD(N)                      &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                      ,CW_B_SOUTH(N)%TASKS              &  !<-- Updated cloud condensate (kg/kg) on nest bndry points.
                                      ,CW_B_NORTH(N)%TASKS              &  !
                                      ,CW_B_WEST(N)%TASKS               &  !
                                      ,CW_B_EAST(N)%TASKS )                !<-- 
!
!-----------------------------------------------------------------------
!***  BEFORE WE CAN LET THE PARENT PROCEED TO COMPUTE WIND COMPONENT
!***  UPDATES ON THE NESTS' BOUNDARIES, WE MUST GENERATE THE PRESSURE
!***  ON THE PARENT'S V POINTS AND AT THE NEST BOUNDARY V POINTS.
!***  THIS IS BEGUN BY SIMPLE 4-PT AVERAGING.  
!-----------------------------------------------------------------------
!
        CALL PRESSURE_ON_NEST_BNDRY_V(PD                                &  !<-- Sigma domain pressure (Pa) on parent mass points
                                     ,IMS,IME,JMS,JME                   &  !<-- Memory dimensions of PD
!
                                     ,PD_B_SOUTH(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Sbndry mass points
                                     ,PD_B_NORTH(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Nbndry mass points
                                     ,PD_B_WEST (N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry mass points
                                     ,PD_B_EAST (N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry mass points
! 
                                     ,NUM_TASKS_SEND_V_S(N)             &  !<-- # of child tasks with south boundary segments on V
                                     ,NUM_TASKS_SEND_V_N(N)             &  !<-- # of child tasks with north boundary segments on V
                                     ,NUM_TASKS_SEND_V_W(N)             &  !<-- # of child tasks with west boundary segments on V
                                     ,NUM_TASKS_SEND_V_E(N)             &  !<-- # of child tasks with east boundary segments on V
!
                                     ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH    &  !<-- Starting I on each south V boundary child task 
                                     ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH    &  !<-- Ending I on each south V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%I_LO_NORTH    &  !<-- Starting I on each north V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%I_HI_NORTH    &  !<-- Ending I on each north V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%J_LO_WEST     &  !<-- Starting J on each west V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%J_HI_WEST     &  !<-- Ending J on each west V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%J_LO_EAST     &  !<-- Starting J on each east V boundary child task
                                     ,CHILDTASK_V_SAVE(N)%J_HI_EAST     &  !<-- Ending J on each east V boundary child task
!
                                     ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH    &  !<-- Starting I on each south H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH    &  !<-- Ending I on each south H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%I_LO_NORTH    &  !<-- Starting I on each north H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%I_HI_NORTH    &  !<-- Ending I on each north H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%J_LO_WEST     &  !<-- Starting J on each west H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%J_HI_WEST     &  !<-- Ending J on each west H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%J_LO_EAST     &  !<-- Starting J on each east H boundary child task
                                     ,CHILDTASK_H_SAVE(N)%J_HI_EAST     &  !<-- Ending J on each east H boundary child task
!
                                     ,N_BLEND_V                         &  !<-- V rows in nests' boundary regions
                                     ,IM_CHILD(N)                       &  !<-- East-west points on child domain
                                     ,JM_CHILD(N)                       &  !<-- North-south points on child domain
!
                                     ,INC_FIX(N)                        &  !<-- Increment used to select nest tasks for averaging H to V
!
                                     ,PD_V                              &  !<-- Sigma domain pressure (Pa) on parent V points
                                     ,PD_B_SOUTH_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Sbndry V points
                                     ,PD_B_NORTH_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Nbndry V points
                                     ,PD_B_WEST_V(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Wbndry V points
                                     ,PD_B_EAST_V(N)%TASKS )               !<-- Sigma domain pressure (Pa) on nest Ebndry V points
!
!-----------------------------------------------------------------------
!
!----------------------
!***  U Wind Component 
!----------------------
!
        CALL PARENT_UPDATE_CHILD_BNDRY(U                                &  !<-- Parent U wind component (m/s)
                                      ,PD_V,PT,PDTOP                    &  !<-- Parent PD; domain PT and PDTOP
                                      ,PSGML1,SGML2,SG1,SG2             &  !<-- General vertical structure (shared by all domains)
!
                                      ,PD_B_SOUTH_V(N)%TASKS            &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                      ,PD_B_NORTH_V(N)%TASKS            &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                      ,PD_B_WEST_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                      ,PD_B_EAST_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task subdomain lateral memory dimensions
                                      ,LM                               &  !<-- # of model layers
                                      ,1                                &  !<-- # of rows to ignore on north/east nest boundaries
!
                                      ,NUM_TASKS_SEND_V_S(N)            &  !<-- # of child tasks with south boundary segments
                                      ,NUM_TASKS_SEND_V_N(N)            &  !<-- # of child tasks with north boundary segments
                                      ,NUM_TASKS_SEND_V_W(N)            &  !<-- # of child tasks with west boundary segments
                                      ,NUM_TASKS_SEND_V_E(N)            &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_V(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child S bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child N bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child W bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child E bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child S bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child N bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child W bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child E bndry point
!
                               ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Not relevant for V points
!
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point.
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND &  !
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND &  !<--
!
                                      ,N_BLEND_V                        &  !<-- Width of boundary blending region
                                      ,IM_CHILD(N)                      &  !<-- East-west points on child domain
                                      ,JM_CHILD(N)                      &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                      ,U_B_SOUTH(N)%TASKS               &  !<-- Updated U wind component (m/s) on nest bndry points.
                                      ,U_B_NORTH(N)%TASKS               &  !
                                      ,U_B_WEST(N)%TASKS                &  !
                                      ,U_B_EAST(N)%TASKS )                 !<-- 
!
!------------
!***  V Wind
!------------
!
        CALL PARENT_UPDATE_CHILD_BNDRY(V                                &  !<-- Parent V wind component (m/s)
                                      ,PD_V,PT,PDTOP                    &  !<-- Parent PD; domain PT and PDTOP
                                      ,PSGML1,SGML2,SG1,SG2             &  !<-- General vertical structure (shared by all domains)
!
                                      ,PD_B_SOUTH_V(N)%TASKS            &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                      ,PD_B_NORTH_V(N)%TASKS            &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                      ,PD_B_WEST_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                      ,PD_B_EAST_V(N)%TASKS             &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task subdomain lateral memory dimensions
                                      ,LM                               &  !<-- # of model layers
                                      ,1                                &  !<-- # of rows to ignore on north/east nest boundaries
!
                                      ,NUM_TASKS_SEND_V_S(N)            &  !<-- # of child tasks with south boundary segments
                                      ,NUM_TASKS_SEND_V_N(N)            &  !<-- # of child tasks with north boundary segments
                                      ,NUM_TASKS_SEND_V_W(N)            &  !<-- # of child tasks with west boundary segments
                                      ,NUM_TASKS_SEND_V_E(N)            &  !<-- # of child tasks with east boundary segments
!
                                     ,PARENT_4_INDICES_V(N)%I_INDX_SBND &  !<-- Parent I's west and east of each child S bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_NBND &  !<-- Parent I's west and east of each child N bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_WBND &  !<-- Parent I's west and east of each child W bndry point
                                     ,PARENT_4_INDICES_V(N)%I_INDX_EBND &  !<-- Parent I's west and east of each child E bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_SBND &  !<-- Parent J's south and north of each child S bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_NBND &  !<-- Parent J's south and north of each child N bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_WBND &  !<-- Parent J's south and north of each child W bndry point
                                     ,PARENT_4_INDICES_V(N)%J_INDX_EBND &  !<-- Parent J's south and north of each child E bndry point
!
                               ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH          &  !<-- Starting I on each south boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH          &  !<-- Ending I on each south boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%I_LO_NORTH          &  !<-- Starting I on each north boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_NORTH          &  !<-- Ending I on each north boundary child task
                               ,CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%J_LO_WEST           &  !<-- Starting J on each west boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_WEST           &  !<-- Ending J on each west boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER  &  !<-- Not relevant for V points
                               ,CHILDTASK_V_SAVE(N)%J_LO_EAST           &  !<-- Starting J on each east boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_EAST           &  !<-- Ending J on each east boundary child task
                               ,CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER  &  !<-- Not relevant for V points
!
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND &  !<-- Bilinear interpolation wgts of the four parent
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND &  !    points surrounding each child bndry point.
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND &  !
                                    ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND &  !<--
!
                                      ,N_BLEND_V                        &  !<-- Width of boundary blending region
                                      ,IM_CHILD(N)                      &  !<-- East-west points on child domain
                                      ,JM_CHILD(N)                      &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             INPUT
!                                                                          --------------
!                                                                             OUTPUT
!                                                                               |   
!                                                                               v   
                                      ,V_B_SOUTH(N)%TASKS               &  !<-- Updated V wind component (m/s) on nest bndry points.
                                      ,V_B_NORTH(N)%TASKS               &  !
                                      ,V_B_WEST(N)%TASKS                &  !
                                      ,V_B_EAST(N)%TASKS )                 !<-- 
!
        cpl2_comp_tim=cpl2_comp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  PARENT TASKS SEND DATA DIRECTLY TO CHILD TASKS WHOSE BOUNDARY
!***  POINTS THE PARENT TASKS CONTAIN. 
!***
!***  MAKE SURE THE BOUNDARY DATA FROM THE PREVIOUS STEP HAS BEEN
!***  RECEIVED BY THE CHILDREN.  ONLY THEN CAN WE OVERWRITE IT WITH
!***  NEW DATA.  
!-----------------------------------------------------------------------
!
!
        NSTEP_CHILD_RECV(N)=NSTEP_CHILD_RECV(N)+PARENT_CHILD_TIME_RATIO(N) !<-- Child "N" is waiting at this timestep to recv its data
!
!-----------
!***  South
!-----------
!
        NTAG=NSTEP_CHILD_RECV(N)+101                                       !<-- Add 101 to obtain a unique South H tag
!
        IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                    !<-- Parent task has Sbndry H data to send to child tasks?     
!
          DO NT=1,NUM_TASKS_SEND_H_S(N)                                    !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,221)n,nt,childtask_bndry_h_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 221       format(' Ready to Send South_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_H_SOUTH(N)%TASKS(NT)%DATA         &  !<-- Child south boundary data on child task NT
                         ,WORDS_BOUND_H_SOUTH(N)%TASKS(NT)              &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NT)          &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,124)n,nt,childtask_bndry_h_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 124       format(' Sent South_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
          ENDDO 
        ENDIF
!
!---------------
!
        NTAG=NSTEP_CHILD_RECV(N)+102                                       !<-- Add 102 to obtain a unique South V tag
!
        IF(NUM_TASKS_SEND_V_S(N)>0)THEN                                    !<-- Parent task has Sbndry V data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_V_S(N)                                    !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,125)n,nt,childtask_bndry_v_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 125       format(' Ready to send South_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_V_SOUTH(N)%TASKS(NT)%DATA         &  !<-- Child south boundary data on child task NT
                         ,WORDS_BOUND_V_SOUTH(N)%TASKS(NT)              &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_V_RANKS(N)%SOUTH(NT)          &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,126)n,nt,childtask_bndry_v_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 126       format(' Sent South_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          ENDDO 
        ENDIF
!
!-----------
!***  North
!-----------
!
        NTAG=NSTEP_CHILD_RECV(N)+103                                       !<-- Add 103 to obtain a unique North H tag
!
        IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                    !<-- Parent task has Nbndry H data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_H_N(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_H_NORTH(N)%TASKS(NT)%DATA         &  !<-- Child north boundary data on child task NT
                         ,WORDS_BOUND_H_NORTH(N)%TASKS(NT)              &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NT)          &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!---------------
!
        NTAG=NSTEP_CHILD_RECV(N)+104                                       !<-- Add 104 to obtain a unique North V tag
!
        IF(NUM_TASKS_SEND_V_N(N)>0)THEN                                    !<-- Parent task has Nbndry V data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_V_N(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_V_NORTH(N)%TASKS(NT)%DATA         &  !<-- Child north boundary data on child task NT
                         ,WORDS_BOUND_V_NORTH(N)%TASKS(NT)              &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_V_RANKS(N)%NORTH(NT)          &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!----------
!***  West
!----------
!
        NTAG=NSTEP_CHILD_RECV(N)+105                                       !<-- Add 105 to obtain a unique West H tag
!
        IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                    !<-- Parent task has Wbndry H data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_H_W(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_H_WEST(N)%TASKS(NT)%DATA          &  !<-- Child west boundary data on child task NT
                         ,WORDS_BOUND_H_WEST(N)%TASKS(NT)               &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NT)           &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!---------------
!
        NTAG=NSTEP_CHILD_RECV(N)+106                                       !<-- Add 106 to obtain a unique West V tag
!
        IF(NUM_TASKS_SEND_V_W(N)>0)THEN                                    !<-- Parent task has Wbndry V data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_V_W(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_V_WEST(N)%TASKS(NT)%DATA          &  !<-- Child west boundary data on child task NT
                         ,WORDS_BOUND_V_WEST(N)%TASKS(NT)               &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_V_RANKS(N)%WEST(NT)           &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!----------
!***  East
!----------
!
        NTAG=NSTEP_CHILD_RECV(N)+107                                       !<-- Add 107 to obtain a unique East H tag
!
        IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                    !<-- Parent task has Ebndry H data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_H_E(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_H_EAST(N)%TASKS(NT)%DATA          &  !<-- Child east boundary data on child task NT
                         ,WORDS_BOUND_H_EAST(N)%TASKS(NT)               &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NT)           &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!---------------
!
        NTAG=NSTEP_CHILD_RECV(N)+108                                       !<-- Add 108 to obtain a unique East V tag
!
        IF(NUM_TASKS_SEND_V_E(N)>0)THEN                                    !<-- Parent task has Ebndry V data to send to child tasks?     
          DO NT=1,NUM_TASKS_SEND_V_E(N)                                    !<-- Send to all appropriate child tasks
!
            btim=timef()
            CALL MPI_SEND(CHILD_BOUND_V_EAST(N)%TASKS(NT)%DATA          &  !<-- Child east boundary data on child task NT
                         ,WORDS_BOUND_V_EAST(N)%TASKS(NT)               &  !<-- # of words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK_BNDRY_V_RANKS(N)%EAST(NT)           &  !<-- Local rank of child to recv data 
                         ,NTAG                                          &  !<-- MPI tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator
                         ,IERR )
            cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
          ENDDO 
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO child_loop
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Clocktime for Comp in Phase2 into Cpl Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='Cpl2_Comp_Time'                     &  !<-- Name of the attribute to insert
                            ,value=cpl2_comp_tim                        &  !<-- Phase 2 Compute time
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Clocktime for Wait in Phase2 into Cpl Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='Cpl2_Wait_Time'                     &  !<-- Name of the attribute to insert
                            ,value=cpl2_wait_tim                        &  !<-- Phase 2 Wait time
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Clocktime for Send in Phase2 into Cpl Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='Cpl2_Send_Time'                     &  !<-- Name of the attribute to insert
                            ,value=cpl2_send_tim                        &  !<-- Phase 2 Send time
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_RUN_SEND
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_CPL_FINALIZE(CPL_COMP                     &
                                          ,IMP_STATE                    &
                                          ,EXP_STATE                    &
                                          ,CLOCK                        &
                                          ,RC_FINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE COUPLER.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP                         !<-- The Parent-Child Coupler Component
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE                        !<-- The Coupler's Import State
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE                        !<-- The Coupler's Export State
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK                            !<-- The ESMF Clock
!
      INTEGER,OPTIONAL,   INTENT(OUT)  :: RC_FINAL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC,RC_CPL_FINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_FINAL    =ESMF_SUCCESS
      RC_CPL_FINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      WRITE(0,*)' Clocktime Parent-Child Coupler'
      WRITE(0,*)'   Cpl1 Prelim=',cpl1_prelim_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_H=',cpl1_south_h_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_V=',cpl1_south_v_tim*1.e-3
      WRITE(0,*)'   Cpl1 North_H=',cpl1_north_h_tim*1.e-3
      WRITE(0,*)'   Cpl1 North_V=',cpl1_north_v_tim*1.e-3
      WRITE(0,*)'   Cpl1 West_H=',cpl1_west_h_tim*1.e-3
      WRITE(0,*)'   Cpl1 West_V=',cpl1_west_v_tim*1.e-3
      WRITE(0,*)'   Cpl1 East_H=',cpl1_east_h_tim*1.e-3
      WRITE(0,*)'   Cpl1 East_V=',cpl1_east_v_tim*1.e-3
      WRITE(0,*)' '
      WRITE(0,*)'   Cpl1 South_H_Recv=',cpl1_south_h_recv_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_H_Undo=',cpl1_south_h_undo_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_H_Exp =',cpl1_south_h_exp_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_V_Recv=',cpl1_south_v_recv_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_V_Undo=',cpl1_south_v_undo_tim*1.e-3
      WRITE(0,*)'   Cpl1 South_V_Exp =',cpl1_south_v_exp_tim*1.e-3
! 
!-----------------------------------------------------------------------
!
      IF(RC_CPL_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL FINALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL FINALIZE STEP FAILED"
      ENDIF
!
      IF(PRESENT(RC_FINAL))THEN
        RC_FINAL=RC_CPL_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_COUPLER_SETUP(NUM_DOMAINS                 &  !      
                                           ,MY_DOMAIN_ID                &  !      
                                           ,NUM_CHILDREN                &  !      
                                           ,COMM_TO_MY_CHILDREN         &  !      
                                           ,COMM_TO_MY_PARENT           &  !      
                                           ,COMM_MY_DOMAIN              &  !      
                                           ,DT                          &  !
                                           ,CHILD_ID                    &  !     ^
                                           ,EXP_STATE_ATM               &  !     |
                                           ,FTASKS_DOMAIN               &  !     |  
                                           ,ID_PARENTS                  &  !   INPUT 
!                                                                           -----------
                                           ,IMP_STATE_CPL_NEST          &  !   OUTPUT
                                           ,EXP_STATE_CPL_NEST          &  !     |
                                           ,PARENT_CHILD_COUPLER_COMP )    !     v
!
!-----------------------------------------------------------------------
!***  CREATE THE PARENT-CHILD COUPLER THROUGH WHICH THEY WILL
!***  COMMUNICATE.  THIS COUPLER IS CALLED BY THE ATM_DRIVER
!***  COMPONENT.  MOVE DATA FROM THE ATM EXPORT STATE INTO THE
!***  PARENT-CHILD COUPLER IMPORT STATE THAT THE COUPLER WILL
!***  NEED IN ORDER TO GENERATE BOUDARY DATA FOR ITS CHILDREN.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: COMM_MY_DOMAIN                              &  !<-- MPI communicator for each individual domain
                           ,COMM_TO_MY_PARENT                           &  !<-- Current domain's MPI communicator to its parent
                           ,MY_DOMAIN_ID                                &  !<-- ID of current domain
                           ,NUM_CHILDREN                                &  !<-- Current domain's number of children
                           ,NUM_DOMAINS                                    !<-- Total number of domains
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: CHILD_ID               &  !<-- Domain IDs of current domain's children
                                                ,COMM_TO_MY_CHILDREN    &  !<-- Current domain's MPI communicators to its children
                                                ,FTASKS_DOMAIN          &  !<-- # of forecast tasks on each domain
                                                ,ID_PARENTS                !<-- IDs of parents of nested domains
!
      REAL,DIMENSION(1:NUM_DOMAINS),INTENT(IN) :: DT                       !<-- Timesteps for all domains (ATM Components)
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_ATM                      !<-- Export state of the current ATM Component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_CPL_NEST              &  !<-- Parent-Child Coupler import state
                                       ,EXP_STATE_CPL_NEST                 !<-- Parent-Child Coupler export state
!
      TYPE(ESMF_CplComp),INTENT(OUT) :: PARENT_CHILD_COUPLER_COMP          !<-- Parent-Child Coupler Component
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
!!    TYPE CHILD_TASKS
!!      INTEGER,DIMENSION(:,:),POINTER :: LIMITS
!!    END TYPE CHILD_TASKS
!
!!    TYPE HANDLES_MULTI     
!!      INTEGER,DIMENSION(:),POINTER :: RECV
!!    END TYPE HANDLES_MULTI
!
      INTEGER :: ITS,ITE,JTS,JTE,LM                                     &
                ,IDS,IDE,JDS,JDE                                        &
                ,ID_DOM,IERR                                            &
                ,IHANDLE_RECV,IHANDLE_SEND                              &
                ,INDX_CW,INDX_Q                                         &
                ,MYPE                                                   &
                ,N,N_BLEND_H,N_BLEND_V                                  &
                ,NN,NX
!
      INTEGER :: RC,RC_NESTSET
!
      INTEGER,DIMENSION(4) :: LIMITS
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: ISTAT,JSTAT
!
      INTEGER,DIMENSION(:),POINTER,SAVE :: PARENT_CHILD_RATIO 
!
!!!   INTEGER,DIMENSION(:,:),ALLOCATABLE :: ISTAT_ALL
!
      REAL,DIMENSION(:),ALLOCATABLE :: ARRAY_1D
!
!!!!  TYPE(HANDLES_RECV),DIMENSION(:),POINTER :: HANDLES
!
      TYPE(ESMF_Array) :: HOLD_ARRAY
!
      integer :: nsize
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_NESTSET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  CREATE IMPORT/EXPORT STATES FOR THE PARENT-CHILD COUPLER
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Nesting"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_STATE_CPL_NEST=ESMF_StateCreate(stateName="Nesting Coupler Import" &  !<-- The Nesting Coupler import state name
                                         ,statetype=ESMF_STATE_IMPORT        &
                                         ,rc       =RC)
!
      EXP_STATE_CPL_NEST=ESMF_StateCreate(stateName="Nesting Coupler Export" &  !<-- The Nesting Coupler export state name
                                         ,statetype=ESMF_STATE_EXPORT        &
                                         ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE PARENT-CHILD COUPLER
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Nesting Coupler Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      PARENT_CHILD_COUPLER_COMP=ESMF_CplCompCreate(name='Parent_Child Coupler' &
                                                  ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER'S INIT, RUN, AND FINALIZE STEPS
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the Nesting Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetServices(comp          =PARENT_CHILD_COUPLER_COMP &  ! <-- The Nesting coupler component
                                  ,subroutineName=PARENT_CHILD_CPL_REGISTER &  ! <-- The user's subroutineName
                                  ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  WE WANT THE WRITE TASKS TO SEE THE PARENT-CHILD COUPLER BUT THEY
!***  SHOULD NEVER ACTUALLY USE IT THEREFORE THEY MAY RETURN NOW.
!***  THE SAME IS TRUE FOR DOMAINS THAT HAVE NO CHILDREN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Parent-Child CPL Setup: Extract Fcst-or-Write Flag from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM component export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_FCST_TASK                     &  !<-- Am I a forecast task?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Parent/Not-a-Parent Flag from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM component export state
                            ,name ='I-Am-A-Parent Flag'                 &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_PARENT                        &  !<-- Am I on a nested domain?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!   IF(I_AM_A_FCST_TASK==ESMF_FALSE.OR.I_AM_A_PARENT==ESMF_FALSE)RETURN
      IF(I_AM_A_FCST_TASK==ESMF_FALSE)RETURN
!
!-----------------------------------------------------------------------
!***  LOAD KEY VARIABLES INTO THE COUPLER'S IMPORT STATE.
!-----------------------------------------------------------------------
!
!-------------------------
!***  Current Domain's ID 
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Current Domain ID to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='MY_DOMAIN_ID'                       &  !<-- Current domain's ID
                            ,value=MY_DOMAIN_ID                         &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
!-----------------------------
!***  Total Number of Domains
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Number of Domains to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='NUM_DOMAINS'                        &  !<-- Total number of domains
                            ,value=NUM_DOMAINS                          &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------
!***  Number of Fcst Tasks on Domains
!-------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Number of Fcst Tasks Per Domain to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Number of forecast tasks on each domain
                            ,count    =NUM_DOMAINS                      &  !<-- Number of domains
                            ,valuelist=FTASKS_DOMAIN                    &  !<-- Insert this into the import state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------
!***  Domain IDs of Parents
!---------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain IDs of Parents to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='ID_PARENTS'                     &  !<-- IDs of parent domain
                            ,count    =NUM_DOMAINS                      &  !<-- Number of domains
                            ,valuelist=ID_PARENTS                       &  !<-- Insert this into the import state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------
!***  Number of Children
!------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Number of Children to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='NUM_CHILDREN'                       &  !<-- This ATM Component's # of children
                            ,value=NUM_CHILDREN                         &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------
!***  Communicators to Children
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Parent-to-Child Communicators to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='Parent-to-Child Comms'          &  !<-- Name of Attribute
                            ,count    =NUM_CHILDREN                     &  !<-- Length of inserted array
                            ,valueList=COMM_TO_MY_CHILDREN              &  !<-- Communicators to my children
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Communicator to Parent
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Child-to-Parent Communicator to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='Child-to-Parent Comm'               &  !<-- Name of Attribute
                            ,value=COMM_TO_MY_PARENT                    &  !<-- The communicator to my parent
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------
!***  Communicator for each domain
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain Intercommunicator to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='Single Domain Comm'                 &  !<-- Name of Attribute
                            ,value=COMM_MY_DOMAIN                       &  !<-- The communicator to my parent
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------
!***  Subdomain Integration Limits
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Integration Limits from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='ITS'                                &  !<-- The name of the Attribute 
                            ,value=ITS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='ITE'                                &  !<-- The name of the Attribute
                            ,value=ITE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='JTS'                                &  !<-- The name of the Attribute
                            ,value=JTS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='JTE'                                &  !<-- The name of the Attribute
                            ,value=JTE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='LM'                                 &  !<-- The name of the Attribute
                            ,value=LM                                   &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute
                            ,value=NHALO                                &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Subdomain Integration Limits to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='ITS'                                &  !<-- The name of the Attribute
                            ,value=ITS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='ITE'                                &  !<-- The name of the Attribute
                            ,value=ITE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='JTS'                                &  !<-- The name of the Attribute
                            ,value=JTS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='JTE'                                &  !<-- The name of the Attribute
                            ,value=JTE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='LM'                                 &  !<-- The name of the Attribute
                            ,value=LM                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute
                            ,value=NHALO                                &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Full Domain Dimensions
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Full Domain Dimensions from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='IDS'                                &  !<-- The name of the Attribute 
                            ,value=IDS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='IDE'                                &  !<-- The name of the Attribute 
                            ,value=IDE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='JDS'                                &  !<-- The name of the Attribute 
                            ,value=JDS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='JDE'                                &  !<-- The name of the Attribute 
                            ,value=JDE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Full Domain Dimensions to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='IDS'                                &  !<-- The name of the Attribute
                            ,value=IDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='IDE'                                &  !<-- The name of the Attribute
                            ,value=IDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='JDS'                                &  !<-- The name of the Attribute
                            ,value=JDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='JDE'                                &  !<-- The name of the Attribute
                            ,value=JDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------
!***  Width of Boundary Blending Region
!---------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Boundary Blending Region Widths from ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='LNSH'                               &  !<-- The name of the Attribute 
                            ,value=N_BLEND_H                            &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The ATM export state
                            ,name ='LNSV'                               &  !<-- The name of the Attribute
                            ,value=N_BLEND_V                            &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Boundary Blending Region Widths to the Nesting CPL Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='N_BLEND_H'                          &  !<-- The name of the Attribute
                            ,value=N_BLEND_H                            &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='N_BLEND_V'                          &  !<-- The name of the Attribute
                            ,value=N_BLEND_V                            &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------
!***  Transfer Sfc Geopotential
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS from Parent ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_ATM                         &  !<-- The Parent's ATM export state
                        ,itemName='FIS'                                 &  !<-- Extract FIS array
                        ,array   =HOLD_ARRAY                            &  !<-- Put the extracted Array here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert FIS into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Nesting Coupler's import state
                        ,array=HOLD_ARRAY                               &  !<-- The Array to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------
!***  Transfer PT,PDTOP,PSGML1,SG1,SG2,SGML2
!--------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PT from Parent ATM Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The Parent's ATM export state
                            ,name ='PT'                                 &  !<-- Extract PT
                            ,value=PT                                   &  !<-- Put the extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='PT'                                 &  !<-- Insert PT
                            ,value=PT                                   &  !<-- Insert this value
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The Parent's ATM export state
                            ,name ='PDTOP'                              &  !<-- Extract PDTOP
                            ,value=PDTOP                                &  !<-- Put the extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Nesting Coupler's import state
                            ,name ='PDTOP'                              &  !<-- Insert PDTOP
                            ,value=PDTOP                                &  !<-- Insert this value
                            ,rc   =RC)
!
      ALLOCATE(ARRAY_1D(1:LM))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_ATM                    &  !<-- The Parent's ATM export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='PSGML1'                         &  !<-- Insert PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_ATM                    &  !<-- The Parent's ATM export state
                            ,name     ='SGML2'                          &  !<-- Extract SGML2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SGML2'                          &  !<-- Insert SGML2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
      ALLOCATE(ARRAY_1D(1:LM+1))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_ATM                    &  !<-- The Parent's ATM export state
                            ,name     ='SG1'                            &  !<-- Extract SG1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SG1'                            &  !<-- Insert SG1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_ATM                    &  !<-- The Parent's ATM export state
                            ,name     ='SG2'                            &  !<-- Extract SG2
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SG2'                            &  !<-- Insert SG2
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      child_block: IF(NUM_CHILDREN>0)THEN             
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------------------------------------
!***  Ratio of Domain's Timestep to Children's
!----------------------------------------------
!
        ALLOCATE(PARENT_CHILD_RATIO(1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
          PARENT_CHILD_RATIO(N)=NINT(DT(MY_DOMAIN_ID)/DT(CHILD_ID(N)))
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Parent-Child DT Ratios to the Nesting CPL Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of Attribute
                              ,count    =NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=PARENT_CHILD_RATIO             &  !<-- The communicator to my parent
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DEALLOCATE(PARENT_CHILD_RATIO)
!
!-------------------------------
!***  The Children's Domain IDs
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Domain IDs of Children to Nesting CPL Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of Attribute
                              ,count    =NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=CHILD_ID                       &  !<-- The children's IDs of this ATM Component
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW TRANSFER THE PARENT'S PROGNOSTIC ARRAYS FROM THE ATM EXPORT
!***  STATE TO THE PARENT-CHILD COUPLER IMPORT STATE THAT WILL BE
!***  REQUIRED FOR THE CHILDREN'S BOUNDARY DATA.
!-----------------------------------------------------------------------
!
!-----------------
!***  Transfer PD
!-----------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PD from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =EXP_STATE_ATM                       &  !<-- The Parent's ATM export state
                          ,itemName='PD'                                &  !<-- Extract PD array
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert PD into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                     &  !<-- The Nesting Coupler's import state
                          ,array=HOLD_ARRAY                             &  !<-- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------
!***  Transfer Temperature
!--------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract T from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =EXP_STATE_ATM                       &  !<-- The Parent's ATM export state
                          ,itemName='T'                                 &  !<-- Extract T array
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert T into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                     &  !<-- The Nesting Coupler's import state
                          ,array=HOLD_ARRAY                             &  !<-- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------
!***  Transfer U Wind
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract U Wind from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =EXP_STATE_ATM                       &  !<-- The Parent's ATM export state
                          ,itemName='U'                                 &  !<-- Extract U array
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert U into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                     &  !<-- The Nesting Coupler's import state
                          ,array=HOLD_ARRAY                             &  !<-- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------
!***  Transfer V Wind
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract V Wind from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =EXP_STATE_ATM                       &  !<-- The Parent's ATM export state
                          ,itemName='V'                                 &  !<-- Extract V array
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert V into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                     &  !<-- The Nesting Coupler's import state
                          ,array=HOLD_ARRAY                             &  !<-- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Transfer TRACERS Array
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Tracers from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =EXP_STATE_ATM                       &  !<-- The Parent's ATM export state
                          ,itemName='TRACERS'                           &  !<-- Extract TRACERS array
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Tracers into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                     &  !<-- The Nesting Coupler's import state
                          ,array=HOLD_ARRAY                             &  !<-- The Array to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------
!***  Transfer INDX_Q
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract INDX_Q from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=EXP_STATE_ATM                      &  !<-- The Parent's ATM export state
                              ,name ='INDX_Q'                           &  !<-- Name of Attribute to extract
                              ,value=INDX_Q                             &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_Q into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                 &  !<-- The Nesting Coupler's import state
                              ,name ='INDX_Q'                           &  !<-- The name of the Attribute to insert
                              ,value=INDX_Q                             &  !<-- The Attribute to be inserted
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------
!***  Transfer INDX_CW
!----------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract INDX_CW from Parent ATM Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=EXP_STATE_ATM                        &  !<-- The Parent's ATM export state
                              ,name ='INDX_CW'                            &  !<-- Name of Attribute to extract
                              ,value=INDX_CW                              &  !<-- Put the extracted Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_CW into Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                 &  !<-- The Nesting Coupler's import state
                              ,name ='INDX_CW'                          &  !<-- The name of the Attribute to insert
                              ,value=INDX_CW                            &  !<-- The Attribute to be inserted
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF child_block
!
!-----------------------------------------------------------------------
!***  ALL PARENT TASKS NEED TO KNOW THE LOCAL SUBDOMAIN LIMITS OF EACH
!***  TASK ON THEIR CHILDREN.  
!-----------------------------------------------------------------------
!
      LIMITS(1)=ITS
      LIMITS(2)=ITE
      LIMITS(3)=JTS
      LIMITS(4)=JTE
!
!-----------------------------------------------------------------------
!***  CHILD TASKS SEND THEIR SUBDOMAIN LIMITS TO PARENT TASK 0.
!-----------------------------------------------------------------------
!
      IHANDLE_SEND=MPI_REQUEST_NULL
!
      IF(COMM_TO_MY_PARENT>0)THEN                                          !<-- Child tasks send to task 0 of parent
        CALL MPI_COMM_RANK(COMM_TO_MY_PARENT,MYPE,IERR)                    !<-- Obtain my rank
        CALL MPI_ISEND(LIMITS                                           &
                      ,4,MPI_INTEGER,0                                  &
                      ,MYPE,COMM_TO_MY_PARENT                           &
                      ,IHANDLE_SEND,IERR)
      ENDIF
!
!-----------------------------------------------------------------------
!***  RANK 0 PARENT TASKS RECV THEIR CHILDREN'S TASKS' SUBDOMAIN LIMITS.
!-----------------------------------------------------------------------
!
      IF(NUM_CHILDREN>0)THEN
!
        ALLOCATE(CTASK_LIMITS(1:NUM_CHILDREN))
!!!!    ALLOCATE(HANDLES     (1:NUM_CHILDREN))
!
        DO N=1,NUM_CHILDREN
          ID_DOM=CHILD_ID(N)
!
          ALLOCATE(CTASK_LIMITS(N)%LIMITS(1:4,1:FTASKS_DOMAIN(ID_DOM)),stat=IERR)  !<-- Pointer to hold each child task's subdomain limits
          IF(IERR/=0)WRITE(0,*)' Failed to allocate CTASK_LIMITS'
!
          CALL MPI_COMM_RANK(COMM_TO_MY_CHILDREN(N),MYPE,IERR)             !<-- Obtain the ranks of parent tasks 
!
          IF(MYPE==0)THEN
            DO NN=1,FTASKS_DOMAIN(ID_DOM)
              CALL MPI_RECV(CTASK_LIMITS(N)%LIMITS(1:4,NN),4,MPI_INTEGER &
                           ,NN-1,NN-1,COMM_TO_MY_CHILDREN(N),JSTAT,IERR)
            ENDDO
          ENDIF
!
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  DO NOT PROCEED UNTIL WE ARE CERTAIN THE PRECEDING COMMUNICATIONS
!***  COMPLETED SUCCESSFULLY.
!-----------------------------------------------------------------------
!
      IF(COMM_TO_MY_PARENT>0)THEN                                          !<-- Child tasks satisfy ISends
        CALL MPI_WAIT(IHANDLE_SEND,ISTAT,IERR)
      ENDIF
!
!-----------------------------------------------------------------------
!***  FOR EACH CHILD, PARENT TASK 0 SENDS ALL OF THE OTHER PARENT TASKS
!***  THE LOCAL SUBDOMAIN LIMITS OF EACH CHILD TASK SO THE PARENT TASKS
!***  WILL KNOW HOW TO DIVVY UP THE CHILD DOMAIN BOUNDARY DATA.
!-----------------------------------------------------------------------
!
      IF(NUM_CHILDREN>0)THEN
        DO N=1,NUM_CHILDREN
          ID_DOM=CHILD_ID(N)
!
          IF(MYPE==0)THEN
            DO NN=1,FTASKS_DOMAIN(ID_DOM)
              DO NX=2,FTASKS_DOMAIN(MY_DOMAIN_ID)
                CALL MPI_SEND(CTASK_LIMITS(N)%LIMITS(1:4,NN),4,MPI_INTEGER &
                             ,NX-1,NX-1,COMM_MY_DOMAIN,IERR)
              ENDDO
            ENDDO
!
          ELSE
            DO NN=1,FTASKS_DOMAIN(ID_DOM)
              CALL MPI_RECV(CTASK_LIMITS(N)%LIMITS(1:4,NN),4,MPI_INTEGER &
                           ,0,MYPE,COMM_MY_DOMAIN,JSTAT,IERR)
            ENDDO
          ENDIF
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_COUPLER_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_INTERP_SETUP(NUM_CHILDREN                 &
                                          ,MY_CHILDREN_ID               &
                                          ,IM_CHILD                     &
                                          ,JM_CHILD                     &
                                          ,FTASKS_DOMAIN                &
                                          ,CTASK_LIMITS                 &
                                          ,N_BLEND_H                    &
                                          ,N_BLEND_V                    &
                                          ,CF                           &
                                          ,ITS,ITE,JTS,JTE              &
                                          ,IDS,IDE,JDS,JDE )
!
!-----------------------------------------------------------------------
!
!***  ALLOCATE THREE PRIMARY INTERPOLATION QUANTITIES NEEDED BY 
!***  A PARENT DOMAIN TO GENERATE BOUNDARY DATA FOR ITS CHILDREN:
!
!    (1) Children's boundary index limits on each parent task;
!    (2) Parent I's and J's surrounding each child boundary point;
!    (3) Bilinear weights of each parent point surrounding each
!        child boundary point.   
!
!-----------------------------------------------------------------------
!
!***  ONLY PARENT TASKS EXECUTE THIS ROUTINE.
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: ITS,ITE,JTS,JTE                             &
                           ,IDS,IDE,JDS,JDE                             &
                           ,N_BLEND_H,N_BLEND_V,NUM_CHILDREN
!
      INTEGER,DIMENSION(1:NUM_CHILDREN),INTENT(IN) :: IM_CHILD          &
                                                     ,JM_CHILD
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: FTASKS_DOMAIN          &
                                                ,MY_CHILDREN_ID         
!
      TYPE(CHILD_TASKS),DIMENSION(:),INTENT(IN) :: CTASK_LIMITS
!
      TYPE(ESMF_Config),DIMENSION(1:NUM_CHILDREN),INTENT(INOUT) :: CF
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: N,N_CHILD_TASKS,NUM_CHILD_TASKS,THIS_CHILD_ID
!
      INTEGER :: EAST_LIMIT1 ,EAST_LIMIT2                               &
                ,WEST_LIMIT1 ,WEST_LIMIT2                               &
                ,NORTH_LIMIT1,NORTH_LIMIT2                              &
                ,SOUTH_LIMIT1,SOUTH_LIMIT2
!
      INTEGER :: RC,RC_SET 
!
      INTEGER,DIMENSION(1:NUM_CHILDREN) :: I_PARENT_START               &
                                          ,J_PARENT_START
!
      REAL,DIMENSION(1:NUM_CHILDREN) :: CHILD_PARENT_SPACE_RATIO
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_SET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS THAT HOLD THE FOUR H AND V PARENT POINTS
!***  THAT SURROUND EACH CHILD POINT IN THE CHILD'S BOUNDARY REGION.
!-----------------------------------------------------------------------
!
      ALLOCATE(PARENT_4_INDICES_H(1:NUM_CHILDREN))
      ALLOCATE(PARENT_4_INDICES_V(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS THAT HOLD THE WEIGHTS OF THE FOUR H AND V 
!***  PARENT POINTS THAT SURROUND EACH CHILD POINT IN THE CHILD'S 
!***  CHILD'S BOUNDARY REGION.
!-----------------------------------------------------------------------
!
      ALLOCATE(PARENT_4_WEIGHTS_H(1:NUM_CHILDREN))
      ALLOCATE(PARENT_4_WEIGHTS_V(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS THAT HOLD THE NUMBER OF CHILD TASKS 
!***  ON EACH SIDE OF THE CHILD BOUNDARIES THAT WILL BE SENT
!***  DATA FROM THE PARENT TASKS.
!-----------------------------------------------------------------------
!
      ALLOCATE(NUM_TASKS_SEND_H_S(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_H_N(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_H_W(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_H_E(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_V_S(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_V_N(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_V_W(1:NUM_CHILDREN)                       &
              ,NUM_TASKS_SEND_V_E(1:NUM_CHILDREN) )
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS THAT WILL HOLD THE RANKS OF ALL CHILD TASKS 
!***  ON EACH SIDE OF THE CHILD BOUNDARIES THAT WILL BE SENT DATA
!***  FROM THE PARENT TASKS.
!-----------------------------------------------------------------------
!
      ALLOCATE(CHILDTASK_BNDRY_H_RANKS(1:NUM_CHILDREN))
      ALLOCATE(CHILDTASK_BNDRY_V_RANKS(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS FOR STARTING/ENDING I's AND J's ON EACH
!***  PARENT TASK FOR EACH SIDE OF THE BOUNDARY.
!-----------------------------------------------------------------------
!
      ALLOCATE(CHILDTASK_H_SAVE(1:NUM_CHILDREN))
      ALLOCATE(CHILDTASK_V_SAVE(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  EXTRACT RELEVANT INFORMATION FROM THE CHILDREN'S CONFIGURE FILES.
!-----------------------------------------------------------------------
!
      child_loop_0: DO N=1,NUM_CHILDREN
!
        THIS_CHILD_ID=MY_CHILDREN_ID(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child_Interp_setup: Extract I of SW Mass Point on Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                    ,value =I_PARENT_START(N)           &  !<-- The variable filled (parent I of child's SW corner)
                                    ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)  
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Child_Init: Extract J of SW Mass Point on Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                    ,value =J_PARENT_START(N)           &  !<-- The variable filled (parent J of child's SW corner)
                                    ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INVERT THE PARENT-TO-CHILD SPACE RATIO FOR COMPUTATION.
!-----------------------------------------------------------------------
!
        CHILD_PARENT_SPACE_RATIO(N)=1./REAL(PARENT_CHILD_SPACE_RATIO(N))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE INDIVIDUAL POINTERS HOLDING THE FOUR H POINTS OF
!***  THE PARENT THAT SURROUND THIS CHILD'S BOUNDARY REGION H POINTS
!***  AND THE BILINEAR INTERPOLATION WEIGHTS OF THE FOUR PARENT POINTS
!***  SURROUNDING THOSE SAME CHILD POINTS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!   ***************************  NOTE  *****************************
!-----------------------------------------------------------------------
!     ALTHOUGH THE H POINTS IN THE NESTS' BOUNDARY REGION COVER ONLY
!     N_BLEND ROWS, WE ACTUALLY NEED TO DO HAVE THE NESTS' PD VALUES
!     ONE ROW FURTHER.  THAT IS BECAUSE WE ALSO NEED PD VALUES AT THE 
!     V POINTS IN THE NESTS' BOUNDARY REGION TO PERFORM THE PROPER
!     HYDROSTATIC UPDATING OF THE WINDS BY THE PARENTS THERE.  TO
!     DO THE 4-POINT AVERAGE NEEDED TO OBTAIN PD ON V POINTS, WE
!     OBVIOUSLY MUST HAVE THEM ON MASS POINTS ONE ROW BEYOND WHERE 
!     THEY ARE NEEDED FOR THE MASS POINTS ALONE.
!-----------------------------------------------------------------------
!   ***************************  NOTE  *****************************
!-----------------------------------------------------------------------
!
        SOUTH_LIMIT1=1
        SOUTH_LIMIT2=N_BLEND_H+1                                           !<-- Extend the region by 1 row for 4-point averaging of PD
!
        NORTH_LIMIT1=JM_CHILD(N)-N_BLEND_H                                 !<-- Extend the region by 1 row for 4-point averaging of PD
        NORTH_LIMIT2=JM_CHILD(N)
!
        WEST_LIMIT1=1
        WEST_LIMIT2=N_BLEND_H+1                                            !<-- Extend the region by 1 row for 4-point averaging of PD
!
        EAST_LIMIT1=IM_CHILD(N)-N_BLEND_H                                  !<-- Extend the region by 1 row for 4-point averaging of PD
        EAST_LIMIT2=IM_CHILD(N)
!
!--------------------------
!***  Parent point indices
!--------------------------
!
        ALLOCATE(PARENT_4_INDICES_H(N)%I_INDX_SBND(1:IM_CHILD(N)              &  !<-- Parent I's west/east of child south bndry H points
                                                  ,SOUTH_LIMIT1:SOUTH_LIMIT2  &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%I_INDX_NBND(1:IM_CHILD(N)              &  !<-- Parent I's west/east of child north bndry H points
                                                  ,NORTH_LIMIT1:NORTH_LIMIT2  &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%I_INDX_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Parent I's west/east of child west bndry H points
                                                  ,1:JM_CHILD(N)              &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%I_INDX_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Parent I's west/east of child east bndry H points
                                                  ,1:JM_CHILD(N)              &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%J_INDX_SBND(1:IM_CHILD(N)              &  !<-- Parent J's south/north of child south bndry H points
                                                  ,SOUTH_LIMIT1:SOUTH_LIMIT2  &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%J_INDX_NBND(1:IM_CHILD(N)              &  !<-- Parent J's south/north of child north bndry H points
                                                  ,NORTH_LIMIT1:NORTH_LIMIT2  &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%J_INDX_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Parent J's south/north of child west bndry H points
                                                  ,1:JM_CHILD(N)              &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_H(N)%J_INDX_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Parent J's south/north of child east bndry H points
                                                  ,1:JM_CHILD(N)              &
                                                  ,1:2))
!
!--------------------------
!***  Parent point weights
!--------------------------
!
        ALLOCATE(PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND(1:IM_CHILD(N)              &  !<-- Bilinear interpolation weights of parent points
                                                   ,SOUTH_LIMIT1:SOUTH_LIMIT2  &  !     surrounding child south bndry region H points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND(1:IM_CHILD(N)              &  !<-- Bilinear interpolation weights of parent points
                                                   ,NORTH_LIMIT1:NORTH_LIMIT2  &  !     surrounding child north bndry region H points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Bilinear interpolation weights of parent points
                                                   ,1:JM_CHILD(N)              &  !     surrounding child west bndry region H points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Bilinear interpolation weights of parent points
                                                   ,1:JM_CHILD(N)              &  !     surrounding child east bndry region H points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE INDIVIDUAL POINTERS HOLDING THE FOUR V POINTS OF
!***  THE PARENT THAT SURROUND THIS CHILD'S BOUNDARY REGION V POINTS.
!-----------------------------------------------------------------------
!
        SOUTH_LIMIT1=1
        SOUTH_LIMIT2=N_BLEND_V
!
        NORTH_LIMIT1=JM_CHILD(N)-1-N_BLEND_V+1
        NORTH_LIMIT2=JM_CHILD(N)-1
!
        WEST_LIMIT1=1
        WEST_LIMIT2=N_BLEND_V
!
        EAST_LIMIT1=IM_CHILD(N)-1-N_BLEND_V+1
        EAST_LIMIT2=IM_CHILD(N)-1
!
!--------------------------
!***  Parent point indices
!--------------------------
!
        ALLOCATE(PARENT_4_INDICES_V(N)%I_INDX_SBND(1:IM_CHILD(N)-1            &  !<-- Parent I's west/east of child south bndry V points.
                                                  ,SOUTH_LIMIT1:SOUTH_LIMIT2  &
                                                  ,1:2))                     
!
        ALLOCATE(PARENT_4_INDICES_V(N)%I_INDX_NBND(1:IM_CHILD(N)-1            &  !<-- Parent I's west/east of child north bndry V points.
                                                  ,NORTH_LIMIT1:NORTH_LIMIT2  &
                                                  ,1:2))                      
!
        ALLOCATE(PARENT_4_INDICES_V(N)%I_INDX_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Parent I's west/east of child west bndry V points.
                                                  ,1:JM_CHILD(N)-1            &
                                                  ,1:2))                       
!
        ALLOCATE(PARENT_4_INDICES_V(N)%I_INDX_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Parent I's west/east of child east bndry V points.
                                                  ,1:JM_CHILD(N)-1            &
                                                  ,1:2))                        
!
        ALLOCATE(PARENT_4_INDICES_V(N)%J_INDX_SBND(1:IM_CHILD(N)-1            &  !<-- Parent J's south/north of child south bndry V points.
                                                  ,SOUTH_LIMIT1:SOUTH_LIMIT2  &
                                                  ,1:2))                    
!
        ALLOCATE(PARENT_4_INDICES_V(N)%J_INDX_NBND(1:IM_CHILD(N)-1            &  !<-- Parent J's south/north of child north bndry V points.
                                                  ,NORTH_LIMIT1:NORTH_LIMIT2  &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_V(N)%J_INDX_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Parent J's south/north of child west bndry V points.
                                                  ,1:JM_CHILD(N)-1            &
                                                  ,1:2))
!
        ALLOCATE(PARENT_4_INDICES_V(N)%J_INDX_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Parent J's south/north of child east bndry V points.
                                                  ,1:JM_CHILD(N)-1            &
                                                  ,1:2))
!
!--------------------------
!***  Parent point weights
!--------------------------
!
        ALLOCATE(PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND(1:IM_CHILD(N)-1            &  !<-- Bilinear interpolation weights of parent points
                                                   ,SOUTH_LIMIT1:SOUTH_LIMIT2  &  !     surrounding child south bndry region V points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND(1:IM_CHILD(N)-1            &  !<-- Bilinear interpolation weights of parent points
                                                   ,NORTH_LIMIT1:NORTH_LIMIT2  &  !     surrounding child north bndry region V points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND(WEST_LIMIT1:WEST_LIMIT2    &  !<-- Bilinear interpolation weights of parent points
                                                   ,1:JM_CHILD(N)-1            &  !     surrounding child west bndry region V points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
        ALLOCATE(PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND(EAST_LIMIT1:EAST_LIMIT2    &  !<-- Bilinear interpolation weights of parent points
                                                   ,1:JM_CHILD(N)-1            &  !     surrounding child east bndry region V points.
                                                   ,1:4))                         !     1:4 indicates SW, SE, NW, NE of child point.
!
!-----------------------------------------------------------------------
!***  WHAT IS THE NUMBER OF FORECAST TASKS ON THE CHILD DOMAIN?
!-----------------------------------------------------------------------
!
        NUM_CHILD_TASKS=FTASKS_DOMAIN(THIS_CHILD_ID)
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS FOR STARTING/ENDING I's AND J's ON EACH
!***  PARENT TASK FOR EACH SIDE OF THE BOUNDARY.
!-----------------------------------------------------------------------
!
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_LO_SOUTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_HI_SOUTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_LO_NORTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_HI_NORTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_LO_WEST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_HI_WEST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_LO_EAST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_HI_EAST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER (1:NUM_CHILD_TASKS))
!
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_LO_SOUTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_HI_SOUTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_LO_NORTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_HI_NORTH   (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER(1:NUM_CHILD_TASKS)) 
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_LO_WEST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_HI_WEST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER(1:NUM_CHILD_TASKS)) 
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_LO_EAST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_HI_EAST    (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER(1:NUM_CHILD_TASKS))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS FOR THE CHILD TASK ID's THAT CONTAIN
!***  SEGMENTS OF THE CHILD BOUNDARY WITHIN A PARENT TASK FOR
!***  EACH SIDE OF THE BOUNDARY.
!-----------------------------------------------------------------------
!
        ALLOCATE(CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_H_RANKS(N)%NORTH(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_H_RANKS(N)%WEST (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_H_RANKS(N)%EAST (1:NUM_CHILD_TASKS))
!
        ALLOCATE(CHILDTASK_BNDRY_V_RANKS(N)%SOUTH(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_V_RANKS(N)%NORTH(1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_V_RANKS(N)%WEST (1:NUM_CHILD_TASKS))
        ALLOCATE(CHILDTASK_BNDRY_V_RANKS(N)%EAST (1:NUM_CHILD_TASKS))
!
!-----------------------------------------------------------------------
!
!***  NOW THE PARENT SETS UP QUANTITIES TO BE USED FOR GENERAL
!***  BILINEAR INTERPOLATION FROM THE PARENT TO ITS CHILDREN'S
!***  BOUNDARY REGIONS.  THESE QUANTITIES ARE:
!
!  (1a) The westernmost/eastermost I's of children's south/north
!       boundary region points on this parent task's subdomain.
!  (1b) The southernmost/northernmost J's of children's west/east
!       boundary region points on this parent task's subdomain.
!   (2) The I,J of the four parent points surrounding each
!       child's boundary region point.
!   (3) The bilinear interpolation weights for each of the four
!       parent points surrounding each child's boundary region 
!       point.
!-----------------------------------------------------------------------
!
!----------------------------------------------------
!***  Interpolation indices and weights for H points
!----------------------------------------------------
!
        CALL PARENT_TO_CHILD_INTERP_FACTORS('H_POINTS'                  &
                                           ,I_PARENT_START(N)           &
                                           ,J_PARENT_START(N)           &
                                           ,N_BLEND_H                   &
!
                                           ,IM_CHILD(N)                 &
                                           ,JM_CHILD(N)                 &
!
                                           ,NUM_CHILD_TASKS             &
                                           ,CTASK_LIMITS(N)%LIMITS      &
!                                                                                 ^
                                           ,CHILD_PARENT_SPACE_RATIO(N) & !       |
!                                                                         !       |
                                           ,ITS,ITE,JTS,JTE             & !       |  
                                           ,IDS,IDE,JDS,JDE             & !     INPUT
!                                                                          -----------------
                                         ,NUM_TASKS_SEND_H_S(N)         & !     OUTPUT
                                         ,NUM_TASKS_SEND_H_N(N)         & !       |
                                         ,NUM_TASKS_SEND_H_W(N)         & !       |
                                         ,NUM_TASKS_SEND_H_E(N)         & !       v
!                                                                         
                                   ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH      &
                                   ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH      &
                                   ,CHILDTASK_H_SAVE(N)%I_LO_NORTH      &
                                   ,CHILDTASK_H_SAVE(N)%I_HI_NORTH      &
                                   ,CHILDTASK_H_SAVE(N)%J_LO_WEST       &
                                   ,CHILDTASK_H_SAVE(N)%J_HI_WEST       &
                                   ,CHILDTASK_H_SAVE(N)%J_LO_EAST       & 
                                   ,CHILDTASK_H_SAVE(N)%J_HI_EAST       &
!
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH    &
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH    &
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%WEST     &
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%EAST     &
!
                                   ,PARENT_4_INDICES_H(N)%I_INDX_SBND   &
                                   ,PARENT_4_INDICES_H(N)%I_INDX_NBND   &
                                   ,PARENT_4_INDICES_H(N)%I_INDX_WBND   &
                                   ,PARENT_4_INDICES_H(N)%I_INDX_EBND   &
                                   ,PARENT_4_INDICES_H(N)%J_INDX_SBND   &
                                   ,PARENT_4_INDICES_H(N)%J_INDX_NBND   &
                                   ,PARENT_4_INDICES_H(N)%J_INDX_WBND   &
                                   ,PARENT_4_INDICES_H(N)%J_INDX_EBND   &
!
                                   ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND  &
                                   ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND  &
                                   ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND  &
                                   ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND)
!
!----------------------------------------------------
!***  Interpolation indices and weights for V points
!----------------------------------------------------
!
        CALL PARENT_TO_CHILD_INTERP_FACTORS('V_POINTS'                  &
                                           ,I_PARENT_START(N)           &
                                           ,J_PARENT_START(N)           &
                                           ,N_BLEND_V                   &
!
                                           ,IM_CHILD(N)                 &
                                           ,JM_CHILD(N)                 &
!
                                           ,NUM_CHILD_TASKS             &
                                           ,CTASK_LIMITS(N)%LIMITS      &
!
                                           ,CHILD_PARENT_SPACE_RATIO(N) & !       ^
!                                                                                 |
                                           ,ITS,ITE,JTS,JTE             & !       |  
                                           ,IDS,IDE,JDS,JDE             & !     INPUT
!                                                                          ----------------
                                         ,NUM_TASKS_SEND_V_S(N)         & !     OUTPUT
                                         ,NUM_TASKS_SEND_V_N(N)         & !       |
                                         ,NUM_TASKS_SEND_V_W(N)         & !       |
                                         ,NUM_TASKS_SEND_V_E(N)         & !       v
!                                                                         
                                   ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH      &
                                   ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH      &
                                   ,CHILDTASK_V_SAVE(N)%I_LO_NORTH      &
                                   ,CHILDTASK_V_SAVE(N)%I_HI_NORTH      &
                                   ,CHILDTASK_V_SAVE(N)%J_LO_WEST       &
                                   ,CHILDTASK_V_SAVE(N)%J_HI_WEST       &
                                   ,CHILDTASK_V_SAVE(N)%J_LO_EAST       & 
                                   ,CHILDTASK_V_SAVE(N)%J_HI_EAST       &
!
                                   ,CHILDTASK_BNDRY_V_RANKS(N)%SOUTH    &
                                   ,CHILDTASK_BNDRY_V_RANKS(N)%NORTH    &
                                   ,CHILDTASK_BNDRY_V_RANKS(N)%WEST     &
                                   ,CHILDTASK_BNDRY_V_RANKS(N)%EAST     &
!
                                   ,PARENT_4_INDICES_V(N)%I_INDX_SBND   &
                                   ,PARENT_4_INDICES_V(N)%I_INDX_NBND   &
                                   ,PARENT_4_INDICES_V(N)%I_INDX_WBND   &
                                   ,PARENT_4_INDICES_V(N)%I_INDX_EBND   &
                                   ,PARENT_4_INDICES_V(N)%J_INDX_SBND   &
                                   ,PARENT_4_INDICES_V(N)%J_INDX_NBND   &
                                   ,PARENT_4_INDICES_V(N)%J_INDX_WBND   &
                                   ,PARENT_4_INDICES_V(N)%J_INDX_EBND   &
!
                                   ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND  &
                                   ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND  &
                                   ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND  &
                                   ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND)
!
!-----------------------------------------------------------------------
!***  FOR V POINT VARIABLES, THE NUMBER OF POINTS TO BE TRANSFERRED
!***  FROM PARENTS TO THEIR CHILDREN'S BOUNDARIES IS THE SAME AS
!***  THE NUMBER OF COMPUTATION POINTS (NO EXTENSIONS AS IS NEEDED
!***  FOR PDB).
!-----------------------------------------------------------------------
!
        DO N_CHILD_TASKS=1,NUM_CHILD_TASKS
!
          CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER(N_CHILD_TASKS)=       &
            CHILDTASK_V_SAVE(N)%I_HI_SOUTH(N_CHILD_TASKS)
!
          CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER(N_CHILD_TASKS)=       &
            CHILDTASK_V_SAVE(N)%I_HI_NORTH(N_CHILD_TASKS)
!
          CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER(N_CHILD_TASKS)=        &
            CHILDTASK_V_SAVE(N)%J_HI_WEST(N_CHILD_TASKS)
!
          CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER(N_CHILD_TASKS)=        &
            CHILDTASK_V_SAVE(N)%J_HI_EAST(N_CHILD_TASKS)
!
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO child_loop_0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_INTERP_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_INTERP_FACTORS(FLAG_H_OR_V             &
                                               ,I_PARENT_START          &
                                               ,J_PARENT_START          &
                                               ,N_BLEND                 &
!
                                               ,IM_CHILD                &
                                               ,JM_CHILD                &        
!                                                                                 
                                               ,NUM_CHILD_TASKS         &         
                                               ,LIMITS                  &         
!                                                                          
                                               ,CHILD_PARENT_SPACE_RATIO & !     ^
!                                                                                |
                                               ,ITS,ITE,JTS,JTE         & !      |  
                                               ,IDS,IDE,JDS,JDE         & !    INPUT
!                                                                          --------------
                                               ,NUM_TASKS_SEND_S        & !    OUTPUT
                                               ,NUM_TASKS_SEND_N        & !      |
                                               ,NUM_TASKS_SEND_W        & !      |
                                               ,NUM_TASKS_SEND_E        & !      v
!                                                                         
                                               ,I_SAVE_LO_SOUTH         & 
                                               ,I_SAVE_HI_SOUTH         &
                                               ,I_SAVE_LO_NORTH         &
                                               ,I_SAVE_HI_NORTH         &
                                               ,J_SAVE_LO_WEST          &
                                               ,J_SAVE_HI_WEST          &
                                               ,J_SAVE_LO_EAST          &
                                               ,J_SAVE_HI_EAST          &
!
                                               ,LOCAL_TASK_RANK_S       &
                                               ,LOCAL_TASK_RANK_N       &
                                               ,LOCAL_TASK_RANK_W       &
                                               ,LOCAL_TASK_RANK_E       &
!
                                               ,I_INDX_SBND             &
                                               ,I_INDX_NBND             &
                                               ,I_INDX_WBND             &
                                               ,I_INDX_EBND             &
                                               ,J_INDX_SBND             &
                                               ,J_INDX_NBND             &
                                               ,J_INDX_WBND             &
                                               ,J_INDX_EBND             &
!
                                               ,WEIGHTS_SBND            &
                                               ,WEIGHTS_NBND            &
                                               ,WEIGHTS_WBND            &
                                               ,WEIGHTS_EBND )
!
!-----------------------------------------------------------------------
!***  PARENT COMPONENTS COMPUTE VARIOUS INDICES, WEIGHTS, ETC.
!***  NEEDED TO GENERATE BOUNDARY POINT DATA FOR THE GIVEN
!***  CHILD THROUGHOUT THE UPCOMING FORECAST.
!***  ONLY PARENT TASKS EXECUTE THIS ROUTINE.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: I_PARENT_START,J_PARENT_START               &  !<-- SW corner of nest lies on this I,J of parent
                           ,IM_CHILD,JM_CHILD                           &  !<-- Horizontal dimensions of nest domain
                           ,N_BLEND                                     &  !<-- Width (in rows) of boundary's blending region
                           ,NUM_CHILD_TASKS                                !<-- # of fcst tasks on the child's domain
!
      INTEGER,INTENT(IN) :: ITE,ITS,JTE,JTS                             &  !<-- Index limits on parent task subdomain
                           ,IDE,IDS,JDE,JDS                                !<-- Full dimensions of parent domain
!
      INTEGER,DIMENSION(1:4,1:NUM_CHILD_TASKS),INTENT(IN) :: LIMITS        !<-- ITS,ITE,JTS,JTE on each task of the child
!
      REAL,INTENT(IN) :: CHILD_PARENT_SPACE_RATIO                          !<-- Ratio of nest grid increment to parent's increment
!
      CHARACTER(*),INTENT(IN) :: FLAG_H_OR_V                               !<-- Are we dealing with H or V child boundary points?
!
      INTEGER,INTENT(OUT) :: NUM_TASKS_SEND_S                           &  !<-- # of child tasks with S bndry segments on this parent task
                            ,NUM_TASKS_SEND_N                           &  !<-- # of child tasks with N bndry segments on this parent task
                            ,NUM_TASKS_SEND_W                           &  !<-- # of child tasks with W bndry segments on this parent task
                            ,NUM_TASKS_SEND_E                              !<-- # of child tasks with E bndry segments on this parent task
!
      INTEGER,DIMENSION(NUM_CHILD_TASKS),INTENT(OUT) ::                 &
                                                     I_SAVE_LO_SOUTH    &  !<-- Starting I on child task containing south boundary segment
                                                    ,I_SAVE_HI_SOUTH    &  !<-- Ending I on child task containing south boundary segment
                                                    ,I_SAVE_LO_NORTH    &  !<-- Starting I on child task containing north boundary segment
                                                    ,I_SAVE_HI_NORTH    &  !<-- Ending I on child task containing north boundary segment
                                                    ,J_SAVE_LO_WEST     &  !<-- Starting J on child task containing west boundary segment
                                                    ,J_SAVE_HI_WEST     &  !<-- Ending J on child task containing west boundary segment
                                                    ,J_SAVE_LO_EAST     &  !<-- Starting J on child task containing east boundary segment
                                                    ,J_SAVE_HI_EAST     &  !<-- Ending J on child task containing east boundary segment
!
                                                   ,LOCAL_TASK_RANK_S   &  !<-- Child task ranks with S bndry on this parent task
                                                   ,LOCAL_TASK_RANK_N   &  !<-- Child task ranks with N bndry on this parent task
                                                   ,LOCAL_TASK_RANK_W   &  !<-- Child task ranks with W bndry on this parent task
                                                   ,LOCAL_TASK_RANK_E      !<-- Child task ranks with E bndry on this parent task
!
      INTEGER,DIMENSION(:,:,:),POINTER,INTENT(OUT) :: I_INDX_SBND       &  !<-- Parent I west/east of child south boundary point 
                                                     ,I_INDX_NBND       &  !<-- Parent I west/east of child north boundary point
                                                     ,I_INDX_WBND       &  !<-- Parent I west/east of child west boundary point
                                                     ,I_INDX_EBND       &  !<-- Parent I west/east of child east boundary point
                                                     ,J_INDX_SBND       &  !<-- Parent J south/north of child south boundary point
                                                     ,J_INDX_NBND       &  !<-- Parent J south/north of child north boundary point
                                                     ,J_INDX_WBND       &  !<-- Parent J south/north of child west boundary point
                                                     ,J_INDX_EBND          !<-- Parent J south/north of child east boundary point
!
      REAL,DIMENSION(:,:,:),POINTER,INTENT(OUT) :: WEIGHTS_SBND         &  !<-- Sbndry bilinear interp wghts for 4 surrounding parent points
                                                  ,WEIGHTS_NBND         &  !<-- Nbndry bilinear interp wghts for 4 surrounding parent points
                                                  ,WEIGHTS_WBND         &  !<-- Wbndry bilinear interp wghts for 4 surrounding parent points
                                                  ,WEIGHTS_EBND            !<-- Ebndry bilinear interp wghts for 4 surrounding parent points
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I_CHILD,J_CHILD,KOUNT_I,KOUNT_J                        &
                ,IM_END,JM_END,N,N_ADD,NC                               &
                ,NC_LAST_S,NC_LAST_N,NC_LAST_W,NC_LAST_E
!
      INTEGER,DIMENSION(1:NUM_CHILD_TASKS) :: I_LIMIT_LO                &
                                             ,I_LIMIT_HI                &
                                             ,ITS_CHILD                 &
                                             ,ITE_CHILD                 &
                                             ,J_LIMIT_LO                &
                                             ,J_LIMIT_HI                &
                                             ,JTS_CHILD                 &
                                             ,JTE_CHILD                 &
                                             ,NC_HOLD_S                 &
                                             ,NC_HOLD_N                 &
                                             ,NC_HOLD_W                 &
                                             ,NC_HOLD_E
!
      REAL :: EPS=1.E-4
!
      REAL :: ADD_INC,ARG1,ARG2                                         &
             ,R_ITS,R_ITE,R_IEND,R_JTS,R_JTE,R_JEND       
!
      REAL :: PARENT_I_CHILD_EBND,PARENT_I_CHILD_WBND                   &
             ,PARENT_J_CHILD_NBND,PARENT_J_CHILD_SBND                   &
             ,RATIO                                                     &
             ,REAL_I_PARENT,REAL_I_START                                &
             ,REAL_J_PARENT,REAL_J_START                                &
             ,RECIP_SUM
!
      REAL :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RATIO=CHILD_PARENT_SPACE_RATIO                                       !<-- Child-to-Parent gridspace ratio
!
!-----------------------------------------------------------------------
!***  CREATE THE REAL INDEX LIMITS ON THE PARENT GRID ACROSS WHICH
!***  THE CHILDREN'S BOUNDARY POINT VALUES WILL BE COMPUTED.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***                      !!!!! NOTE !!!!!
!
!***  FOR THE PURPOSE OF HANDLING CHILD BOUNDARIES, PARENT TASKS WILL
!***  "BEGIN" DIRECTLY ON THEIR SOUTHERNMOST/WESTERNMOST H OR V POINTS.
!***  EACH PARENT TASK COVERS THE GAP BETWEEN ITSELF AND THE NEXT TASK
!***  ON THE PARENT GRID IN EACH DIRECTION.
!***  THIS MEANS THAT IF A CHILD SOUTH BOUNDARY POINT LIES EXACTLY ON 
!***  A PARENT TASK POINT THAT ITSELF IS ON THE WESTERNMOST SIDE OF
!***  THAT PARENT TASK'S INTEGRATION SUBDOMAIN THEN THAT CHILD POINT
!***  WILL BE CONSIDERED TO LIE ON BOTH THAT PARENT TASK AND THE
!***  PARENT TASK TO THE WEST SIMPLY BECAUSE IT IS THE INTERSECTION
!***  OF THE REGIONS MANAGED BY BOTH OF THOSE PARENT TASKS.
!***  THIS SAME NOTION APPLIES FOR ALL OTHER DIRECTIONS AND SIDES.
!-----------------------------------------------------------------------
!
      R_ITS =REAL(ITS)-EPS                                                 !<-- REAL Istart of parent task's subdomain
      R_ITE =REAL(ITE)                                                     !<-- REAL Iend of parent task's subdomain
!
      R_JTS =REAL(JTS)-EPS                                                 !<-- REAL Jstart of parent task's subdomain
      R_JTE =REAL(JTE)                                                     !<-- REAL Jend of parent task's subdomain
!
!
!-----------------------------------------------------------------------
!***  BECAUSE EACH PARENT GRIDPOINT COVERS THE GAP TO THE NEXT PARENT
!***  GRIDPOINT AS EXPLAINED ABOVE, INCREASE THE SEARCH LIMIT FOR
!***  CHILD BOUNDARY POINTS.  THAT INCREASE WOULD BE 1 FOR BOTH H AND V
!***  BUT DUE TO THE NATURE OF THE B-GRID LAYOUT AND THE FACT THAT
!***  THE I INDEX OF CHILD V POINTS ON THE WEST BOUNDARY AND THE 
!***  J INDEX OF THE CHILD V POINTS ON THE SOUTH BOUNDARY HAVE SMALLER
!***  GRID INDEX VALUES IN TERMS OF THE PARENT INDICES, WE MUST SEARCH
!***  FOR CHILD H POINTS 1/2+0.5*(space_ratio) GRID INCREMENTS FURTHER
!***  THAN FOR CHILD V POINTS IN ORDER TO REACH THE SAME ACTUAL POSITION.
!-----------------------------------------------------------------------
!
!***  In this diagram the H's and V's are points on the parent task's
!***  subdomain while the h's and v's are points on a nest.  It shows
!***  how each parent point must look eastward.  The same goes for
!***  looking northwward.  A parent/nest ratio of 1/3 is used in this
!***  diagram.
!
!-----------------------------------------------------------------------
!
!     H           H           H
!
!
!
!           V           V
!                       v
!                 h   h   h   h
!
!     H           H           H
!
!     -----------> ----> ->
!           1       1/2 1/6
!     ^
!     |
! This parent
! gridpoint
! must reach
! cover area
! to next H.
! But V with
! same I as
! next H is
! 1.5 farther
! than this H.
! If nest v at
! 1.5 past this
! H is east
! bndry of nest
! then east h
! on bndry is
! 1+1/2+1/6.
! That is how
! far we must
! scan from
! this H.
!-----------------------------------------------------------------------
!
      IF(FLAG_H_OR_V=='H_POINTS')THEN
!       ADD_INC=1.5
        ADD_INC=1.5+0.5*RATIO+EPS

      ELSE
        ADD_INC=1.0
      ENDIF
!
!-----------------------------------------------------------------------
!
      DO N=1,NUM_CHILD_TASKS                                               
        I_SAVE_LO_SOUTH(N)=-1
        I_SAVE_LO_NORTH(N)=-1
        J_SAVE_LO_WEST (N)=-1
        J_SAVE_LO_EAST (N)=-1
!
        NC_HOLD_S(N)=-1
        NC_HOLD_N(N)=-1
        NC_HOLD_W(N)=-1
        NC_HOLD_E(N)=-1
      ENDDO
!
!-----------------------------------------------------------------------
!***  WHAT ARE THE CHILD I AND J INDEX LIMITS OF ANY SECTIONS OF ITS
!***  BOUNDARY THAT LIE WITHIN A PARENT TASK'S SUBDOMAIN?
!
!***  WHAT ARE THE INDICES OF THE FOUR PARENT GRIDPOINTS SURROUNDING
!***  EACH CHILD BOUNDARY POINT?
!
!***  WHAT ARE THE BILINEAR WEIGHTS ASSOCIATED WITH EACH OF THE FOUR
!***  SURROUNDING PARENT POINTS TO OBTAIN THE CHILD BOUNDARY POINT?
!
!***  THE PARENT WILL USE THESE PIECES OF INFORMATION TO INTERPOLATE
!***  FROM ITS GRID TO ITS CHILDREN'S BOUNDARY POINTS.
!----------------------------------------------------------------------- 
!
!-----------------------------------------------------------
!**********************  NOTE  *****************************
!-----------------------------------------------------------
!***  We assume that the WIDTH of the blending region of
!***  a child's boundary does NOT cross the border between
!***  two parent tasks' subdomains.
!-----------------------------------------------------------
!
      IF(FLAG_H_OR_V=='H_POINTS')THEN
        PARENT_J_CHILD_SBND=REAL(J_PARENT_START)                           !<-- J index of parent H for child's south H boundary
        PARENT_J_CHILD_NBND=PARENT_J_CHILD_SBND+(JM_CHILD-1)*RATIO         !<-- J index of parent H for child's north H boundary
        PARENT_I_CHILD_WBND=REAL(I_PARENT_START)                           !<-- I index of parent H for child's west H boundary
        PARENT_I_CHILD_EBND=PARENT_I_CHILD_WBND+(IM_CHILD-1)*RATIO         !<-- I index of parent H for child's east H boundary
        IM_END=IM_CHILD
        JM_END=JM_CHILD
        N_ADD=1                                                            !<-- Blending region along child's boundary
                                                                           !    increased by 1 row to allow 4-pt averaging of PD.
!
      ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
        PARENT_J_CHILD_SBND=REAL(J_PARENT_START)-0.5+RATIO*0.5             !<-- J index of parent V for child's south V boundary
        PARENT_J_CHILD_NBND=PARENT_J_CHILD_SBND+(JM_CHILD-2)*RATIO         !<-- J index of parent V for child's north V boundary
        PARENT_I_CHILD_WBND=REAL(I_PARENT_START)-0.5+RATIO*0.5             !<-- I index of parent V for child's west V boundary
        PARENT_I_CHILD_EBND=PARENT_I_CHILD_WBND+(IM_CHILD-2)*RATIO         !<-- I index of parent V for child's east V boundary
        IM_END=IM_CHILD-1
        JM_END=JM_CHILD-1
        N_ADD=0                                                            !<-- Blending region along child's boundary
                                                                           !    increased only for mass points (for PD averaging)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      DO N=1,NUM_CHILD_TASKS                                               !<-- Loop through forecast tasks on the child domain
        ITS_CHILD(N)=LIMITS(1,N)                                           !<-- ITS on this child task
        ITE_CHILD(N)=LIMITS(2,N)                                           !<-- ITE on this child task
        JTS_CHILD(N)=LIMITS(3,N)                                           !<-- JTS on this child task
        JTE_CHILD(N)=LIMITS(4,N)                                           !<-- JTE on this child task
        I_LIMIT_LO(N)=MAX(ITS_CHILD(N)-2,1)                                !<-- Starting I's for each child task on N/S bndries (2-pt halo)
!!!     I_LIMIT_HI(N)=MIN(ITE_CHILD(N)+2,IM_END)                           !<-- Ending I's for each child task on N/S bndries (2-pt halo)
        I_LIMIT_HI(N)=MIN(ITE_CHILD(N)+2+N_ADD,IM_END)                     !<-- Ending I's for each child task on N/S bndries (2-pt halo)
        J_LIMIT_LO(N)=MAX(JTS_CHILD(N)-2,1)                                !<-- Starting J's for each child task on W/E bndries (2-pt halo)
!!!     J_LIMIT_HI(N)=MIN(JTE_CHILD(N)+2,JM_END)                           !<-- Ending J's for each child task on W/E bndries (2-pt halo)
        J_LIMIT_HI(N)=MIN(JTE_CHILD(N)+2+N_ADD,JM_END)                     !<-- Ending J's for each child task on W/E bndries (2-pt halo)
      ENDDO
!
!-----------------------------------------------------
!-----------------------------------------------------
!***  Child's southern/northern boundaries
!
!***  Work eastward along these boundaries
!***  and save the basic indices and weights
!***  needed by the parent.
!-----------------------------------------------------
!-----------------------------------------------------
!
      NC_LAST_S      =-1
      NC_LAST_N      =-1
      NUM_TASKS_SEND_S=0
      NUM_TASKS_SEND_N=0
!
      ARG1=REAL(ITE)+ADD_INC
      ARG2=REAL(IDE)
      R_IEND=MIN(ARG1,ARG2)-EPS                                            !<-- REAL Iend of parent task's region for child N/S boundaries
      R_JEND=REAL(MIN(JTE+1,JDE))-EPS                                      !<-- Allow search for child points to go into
                                                                           !    the parent's halo
!-----------------------------------------------------
!
      REAL_I_START=PARENT_I_CHILD_WBND
!
!-----------------------------------------------------------------------
      i_loop: DO I_CHILD=1,IM_END                                          !<-- Loop over child I's across its South/North boundaries
!-----------------------------------------------------------------------
!
        REAL_I_PARENT=REAL_I_START+(I_CHILD-1)*RATIO                       !<-- Parent I index coinciding with child domain point
!     write(0,*)' real_i_start=',real_i_start,' real_i_parent=',real_i_parent,' r_its=',r_its,' r_ite=',r_ite
!
!       i_block: IF(REAL_I_PARENT>=R_ITS.AND.REAL_I_PARENT<=R_IEND)THEN    !<-- Column (I) of child's S/N bndry point lies on parent task?
        i_block: IF(REAL_I_PARENT>=R_ITS.AND.REAL_I_PARENT< R_IEND)THEN    !<-- Column (I) of child's S/N bndry point lies on parent task?
!     if(flag_h_or_v=='H_POINTS')then
!       write(0,*)' okay this H i_child=',i_child,' at parent i=',real_i_parent,' lies on parent task'
!     else
!       write(0,*)' okay this V i_child=',i_child,' at parent i=',real_i_parent,' lies on parent task'
!     endif
!
!-----------
!-----------
!***  South
!-----------
!-----------
!
          REAL_J_START=PARENT_J_CHILD_SBND
          KOUNT_J=0
!
!-----------------------------------------------------------------------
!     
          J_CHILD=1
!
!-----------------------------------------------------------------------
!***  Which child task contains this (I_CHILD,J_CHILD) point?
!***  Find out then save the I limits of that task's boundary
!***  segment on this parent task so that the parent task
!***  will know exactly which words to send to the child task.
!***  Also remember that the NMM boundary update routines go
!***  two points into the halo which means that two child
!***  tasks will share some boundary points on the parent.
!-----------------------------------------------------------------------
!
          child_ij_s: IF(REAL_J_START >=R_JTS.AND.REAL_J_START < R_JEND)THEN  !<-- Does parent task see this row of its child?
!
            DO NC=1,NUM_CHILD_TASKS                                        !<-- Loop through all tasks on child domain
!
              IF(I_CHILD>=I_LIMIT_LO(NC).AND.                           &  !<-- Does current child boundary point on this
                 I_CHILD<=I_LIMIT_HI(NC)                                &  !    parent task lie on child task "NC"?
                 .AND.                                                  &  !
                 J_CHILD>=JTS_CHILD(NC).AND.                            &  !
                 J_CHILD<=JTE_CHILD(NC))THEN
!
                IF(NC>NC_LAST_S)THEN                                       !<-- Encountered a new child task holding this S bndry point?
                  NUM_TASKS_SEND_S=NUM_TASKS_SEND_S+1                      !<-- Then increment the S bndry counter of the child tasks
                  LOCAL_TASK_RANK_S(NUM_TASKS_SEND_S)=NC-1                 !<-- Save this child task's local rank
                  NC_LAST_S=NC
                  NC_HOLD_S(NC)=NUM_TASKS_SEND_S
                ENDIF
!
                IF(I_SAVE_LO_SOUTH(NC_HOLD_S(NC))<0)THEN
                  I_SAVE_LO_SOUTH(NC_HOLD_S(NC))=I_CHILD                   !<-- Save the starting I index on this child task's segment
                ENDIF
                I_SAVE_HI_SOUTH(NC_HOLD_S(NC))=I_CHILD                     !<-- Save the ending I index on this child task's segment
!
              ENDIF
! 
            ENDDO
! 
!-----------------------------------------------------------------------
!
            j_south: DO J_CHILD=1,N_BLEND+N_ADD                            !<-- Blending region along child's southern boundary
!
              KOUNT_J=KOUNT_J+1
              REAL_J_PARENT=REAL_J_START+(KOUNT_J-1)*RATIO                 !<-- REAL parent J for this child's J
!
              I_INDX_SBND(I_CHILD,J_CHILD,1)=INT(REAL_I_PARENT+EPS)        !<-- Parent I west of child's south boundary point
              I_INDX_SBND(I_CHILD,J_CHILD,2)=INT(REAL_I_PARENT+EPS)+1      !<-- Parent I east of child's south boundary point
              J_INDX_SBND(I_CHILD,J_CHILD,1)=INT(REAL_J_PARENT+EPS)        !<-- Parent J south of child's south boundary point
              J_INDX_SBND(I_CHILD,J_CHILD,2)=INT(REAL_J_PARENT+EPS)+1      !<-- Parent J north of child's south boundary point
!
              WEIGHT_SW=(I_INDX_SBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (J_INDX_SBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_SE=(REAL_I_PARENT-I_INDX_SBND(I_CHILD,J_CHILD,1))* &
                        (J_INDX_SBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_NW=(I_INDX_SBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (REAL_J_PARENT-J_INDX_SBND(I_CHILD,J_CHILD,1))
              WEIGHT_NE=(REAL_I_PARENT-I_INDX_SBND(I_CHILD,J_CHILD,1))* &
                        (REAL_J_PARENT-J_INDX_SBND(I_CHILD,J_CHILD,1))
!
              RECIP_SUM=1./(WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE)
!
              WEIGHTS_SBND(I_CHILD,J_CHILD,INDX_SW)=WEIGHT_SW*RECIP_SUM    !<-- Bilin interp wght of parent point SW of child bndry point
              WEIGHTS_SBND(I_CHILD,J_CHILD,INDX_SE)=WEIGHT_SE*RECIP_SUM    !<-- Bilin interp wght of parent point SE of child bndry point
              WEIGHTS_SBND(I_CHILD,J_CHILD,INDX_NW)=WEIGHT_NW*RECIP_SUM    !<-- Bilin interp wght of parent point NW of child bndry point
              WEIGHTS_SBND(I_CHILD,J_CHILD,INDX_NE)=WEIGHT_NE*RECIP_SUM    !<-- Bilin interp wght of parent point NE of child bndry point
!
            ENDDO j_south
! 
!-----------------------------------------------------------------------
          ENDIF child_ij_s
!-----------------------------------------------------------------------
!
!
!-----------
!-----------
!***  North
!-----------
!-----------
!
          REAL_J_START=PARENT_J_CHILD_NBND-(N_BLEND-1+N_ADD)*RATIO         !<-- N_ADD accounts for the additional row for H points
          KOUNT_J=0
!
!-----------------------------------------------------------------------
!
          J_CHILD=JM_END-N_BLEND+1-N_ADD                                   !<-- Southernmost child J in north boundary
                                                                           !    blending region.
!-----------------------------------------------------------------------
!
          child_ij_n: IF(REAL_J_START >=R_JTS.AND.REAL_J_START < R_JEND)THEN  !<-- Does parent task see this row of its child?
!
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!***  Find the child tasks and their relevant limits
!***  along the child's northern boundary.
!-------------------------------------------------------------
!
            DO NC=1,NUM_CHILD_TASKS                                        !<-- Loop through all tasks on child domain
!
!     write(0,*)' real_j_start=',real_j_start,' nc=',nc
              IF(I_CHILD>=I_LIMIT_LO(NC).AND.                           &  !<-- Does current child boundary point on this
                 I_CHILD<=I_LIMIT_HI(NC)                                &  !    parent task lie on child task "NC"?
                 .AND.                                                  &  !
                 J_CHILD>=JTS_CHILD(NC).AND.                            &  !
                 J_CHILD<=JTE_CHILD(NC))THEN
!     write(0,*)' i_child=',i_child,' i_limit_lo=',i_limit_lo(nc),' i_limit_hi=',i_limit_hi(nc) &
!              ,' nc=',nc,' nc_last_n=',nc_last_n
!
                IF(NC>NC_LAST_N)THEN                                       !<-- Have we encountered a new child task holding this N bndry?
                  NUM_TASKS_SEND_N=NUM_TASKS_SEND_N+1                      !<-- Then increment the N bndry counter of the child tasks
                  LOCAL_TASK_RANK_N(NUM_TASKS_SEND_N)=NC-1                 !<-- Save this child task's local rank
                  NC_LAST_N=NC
                  NC_HOLD_N(NC)=NUM_TASKS_SEND_N
!     write(0,*)' inner block num_tasks_send_n=',num_tasks_send_n,' local_task_rank_n=',nc-1 &
!              ,' nc_last_n=',nc,' nc_hold_n=',num_tasks_send_n
                ENDIF
!
!     write(0,*)' i_save_lo_north=',i_save_lo_north(nc_hold_n(nc)),' nc_hold_n=',nc_hold_n(nc),' nc=',nc
                IF(I_SAVE_LO_NORTH(NC_HOLD_N(NC))<0)THEN
                  I_SAVE_LO_NORTH(NC_HOLD_N(NC))=I_CHILD                   !<-- Save the starting I index on this child task's segment
!     write(0,*)' save new i_save_lo_north=',i_save_lo_north(nc_hold_n(nc)),' nc_hold_n=',nc_hold_n(nc),' nc=',nc
                ENDIF
                I_SAVE_HI_NORTH(NC_HOLD_N(NC))=I_CHILD                     !<-- Save the ending I index on this child task's segment
! 
              ENDIF
! 
            ENDDO
! 
!-----------------------------------------------------------------------
!
            j_north: DO J_CHILD=JM_END-N_BLEND+1-N_ADD,JM_END              !<-- Blending region of child's northern boundary
!
!-----------------------------------------------------------------------
!
              KOUNT_J=KOUNT_J+1
              REAL_J_PARENT=REAL_J_START+(KOUNT_J-1)*RATIO                 !<-- REAL parent J for this child's J
!
              I_INDX_NBND(I_CHILD,J_CHILD,1)=INT(REAL_I_PARENT+EPS)        !<-- Parent I west of child's north boundary point
              I_INDX_NBND(I_CHILD,J_CHILD,2)=INT(REAL_I_PARENT+EPS)+1      !<-- Parent I east of child's north boundary point
              J_INDX_NBND(I_CHILD,J_CHILD,1)=INT(REAL_J_PARENT+EPS)        !<-- Parent J south of child's north boundary point
              J_INDX_NBND(I_CHILD,J_CHILD,2)=INT(REAL_J_PARENT+EPS)+1      !<-- Parent J north of child's north boundary point
!
              WEIGHT_SW=(I_INDX_NBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (J_INDX_NBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_SE=(REAL_I_PARENT-I_INDX_NBND(I_CHILD,J_CHILD,1))* &
                        (J_INDX_NBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_NW=(I_INDX_NBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (REAL_J_PARENT-J_INDX_NBND(I_CHILD,J_CHILD,1))
              WEIGHT_NE=(REAL_I_PARENT-I_INDX_NBND(I_CHILD,J_CHILD,1))* &
                        (REAL_J_PARENT-J_INDX_NBND(I_CHILD,J_CHILD,1))
!
              RECIP_SUM=1./(WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE)
!
              WEIGHTS_NBND(I_CHILD,J_CHILD,INDX_SW)=WEIGHT_SW*RECIP_SUM    !<-- Interp wght of parent point SW of child bndry point
              WEIGHTS_NBND(I_CHILD,J_CHILD,INDX_SE)=WEIGHT_SE*RECIP_SUM    !<-- Interp wght of parent point SE of child bndry point
              WEIGHTS_NBND(I_CHILD,J_CHILD,INDX_NW)=WEIGHT_NW*RECIP_SUM    !<-- Interp wght of parent point NW of child bndry point
              WEIGHTS_NBND(I_CHILD,J_CHILD,INDX_NE)=WEIGHT_NE*RECIP_SUM    !<-- Interp wght of parent point NE of child bndry point
!
            ENDDO j_north
!
!-----------------------------------------------------------------------
          ENDIF child_ij_n
!-----------------------------------------------------------------------
!
        ENDIF i_block
!
!-----------------------------------------------------------------------
!
      ENDDO i_loop
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------
!-----------------------------------------------------
!***  Child's western/eastern boundaries
!
!***  Work northward along these boundaries
!***  and save the basic indices and weights
!***  needed by the parent.
!-----------------------------------------------------
!-----------------------------------------------------
!
      NC_LAST_W      =-1
      NC_LAST_E      =-1
      NUM_TASKS_SEND_W=0                                                   !<-- Parent task sends to this many child tasks on W bndry
      NUM_TASKS_SEND_E=0                                                   !<-- Parent task sends to this many child tasks on E bndry
!
      ARG1=REAL(JTE)+ADD_INC
      ARG2=REAL(JDE)
      R_JEND=MIN(ARG1,ARG2)-EPS                                            !<-- REAL Jend of parent task's region for child W/E boundaries
!
      R_IEND=REAL(MIN(ITE+1,IDE))-EPS                                      !<-- Allow search for child points to go into 
                                                                           !    the parent's halo
!-----------------------------------------------------------------------
!
      REAL_J_START=PARENT_J_CHILD_SBND
!
!-----------------------------------------------------------------------
      j_loop: DO J_CHILD=1,JM_END                                          !<-- Loop through child J's across its W/E boundaries
!-----------------------------------------------------------------------

        REAL_J_PARENT=REAL_J_START+(J_CHILD-1)*RATIO                       !<-- Parent J index coinciding with child domain point
!
!       j_block: IF(REAL_J_PARENT>=R_JTS.AND.REAL_J_PARENT<=R_JEND)THEN    !<-- Row (J) of child's W/E bndry point lies on parent task?
        j_block: IF(REAL_J_PARENT>=R_JTS.AND.REAL_J_PARENT< R_JEND)THEN    !<-- Row (J) of child's W/E bndry point lies on parent task?
!
!----------
!----------
!***  West
!----------
!----------
!
          REAL_I_START=PARENT_I_CHILD_WBND
          KOUNT_I=0
!
!-----------------------------------------------------------------------
!
          I_CHILD=1
!
!-----------------------------------------------------------------------
!
          child_ij_w: IF(REAL_I_START >=R_ITS.AND.REAL_I_START < R_IEND)THEN  !<-- Does parent task see this column of its child?
!
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!***  Find the child tasks and their relevant limits
!***  along the child's western boundary.
!-------------------------------------------------------------
!
            DO NC=1,NUM_CHILD_TASKS                                        !<-- Loop through all tasks on child domain
!
              IF(J_CHILD>=J_LIMIT_LO(NC).AND.                           &  !<-- Does current child boundary point on this
                 J_CHILD<=J_LIMIT_HI(NC)                                &  !    parent task lie on child task "NC"?
                 .AND.                                                  &  !
                 I_CHILD>=ITS_CHILD(NC).AND.                            &  !
                 I_CHILD<=ITE_CHILD(NC))THEN
!
                IF(NC>NC_LAST_W)THEN                                       !<-- Have we encountered a new child task holding this W bndry?
                  NUM_TASKS_SEND_W=NUM_TASKS_SEND_W+1                      !<-- Then increment the W bndry counter of the child tasks
                  LOCAL_TASK_RANK_W(NUM_TASKS_SEND_W)=NC-1                 !<-- Save this child task's local rank
                  NC_LAST_W=NC
                  NC_HOLD_W(NC)=NUM_TASKS_SEND_W
                ENDIF
!
                IF(J_SAVE_LO_WEST(NC_HOLD_W(NC))<0)THEN
                  J_SAVE_LO_WEST(NC_HOLD_W(NC))=J_CHILD                    !<-- Save the starting J index on this child task's segment
                ENDIF
                J_SAVE_HI_WEST(NC_HOLD_W(NC))=J_CHILD                      !<-- Save the ending J index on this child task's segment
!
              ENDIF
!
            ENDDO
!
!-------------------------------------------------------------
            i_west: DO I_CHILD=1,N_BLEND+N_ADD                             !<-- Blending region of child's western boundary
!-------------------------------------------------------------
!
              KOUNT_I=KOUNT_I+1
              REAL_I_PARENT=REAL_I_START+(KOUNT_I-1)*RATIO                 !<-- REAL parent I for this child's I
!
              I_INDX_WBND(I_CHILD,J_CHILD,1)=INT(REAL_I_PARENT+EPS)        !<-- Parent I west of child's west boundary point
              I_INDX_WBND(I_CHILD,J_CHILD,2)=INT(REAL_I_PARENT+EPS)+1      !<-- Parent I east of child's west boundary point
              J_INDX_WBND(I_CHILD,J_CHILD,1)=INT(REAL_J_PARENT+EPS)        !<-- Parent J south of child's west boundary point
              J_INDX_WBND(I_CHILD,J_CHILD,2)=INT(REAL_J_PARENT+EPS)+1      !<-- Parent J north of child's west boundary point
!
              WEIGHT_SW=(I_INDX_WBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (J_INDX_WBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_SE=(REAL_I_PARENT-I_INDX_WBND(I_CHILD,J_CHILD,1))* &
                        (J_INDX_WBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_NW=(I_INDX_WBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (REAL_J_PARENT-J_INDX_WBND(I_CHILD,J_CHILD,1))
              WEIGHT_NE=(REAL_I_PARENT-I_INDX_WBND(I_CHILD,J_CHILD,1))* &
                        (REAL_J_PARENT-J_INDX_WBND(I_CHILD,J_CHILD,1))
!
              RECIP_SUM=1./(WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE)
!
              WEIGHTS_WBND(I_CHILD,J_CHILD,INDX_SW)=WEIGHT_SW*RECIP_SUM    !<-- Interp wght of parent point SW of child bndry point
              WEIGHTS_WBND(I_CHILD,J_CHILD,INDX_SE)=WEIGHT_SE*RECIP_SUM    !<-- Interp wght of parent point SE of child bndry point
              WEIGHTS_WBND(I_CHILD,J_CHILD,INDX_NW)=WEIGHT_NW*RECIP_SUM    !<-- Interp wght of parent point NW of child bndry point
              WEIGHTS_WBND(I_CHILD,J_CHILD,INDX_NE)=WEIGHT_NE*RECIP_SUM    !<-- Interp wght of parent point NE of child bndry point
!
            ENDDO i_west
!
!-----------------------------------------------------------------------
          ENDIF child_ij_w
!-----------------------------------------------------------------------
!
!----------
!----------
!***  East
!----------
!----------
!
          REAL_I_START=PARENT_I_CHILD_EBND-(N_BLEND-1+N_ADD)*RATIO
          KOUNT_I=0
!
!-----------------------------------------------------------------------
!***  Recall that we need an additional row of H points to allow 4-pt
!***  averaging of PD to V points.  We need only to search for the
!***  westernmost child J row of the east boundary blending region with
!***  the extra row because if that child I is on a parent task then
!***  all of the blending region must be on that task since we are 
!***  permitting the search to go into the parent tasks' haloes.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!      
          I_CHILD=IM_END-N_BLEND+1-N_ADD
!      
!-----------------------------------------------------------------------
!
          child_ij_e: IF(REAL_I_START >=R_ITS.AND.REAL_I_START < R_IEND)THEN  !<-- Does parent task see this column of its child?
!
!-----------------------------------------------------------------------
!-------------------------------------------------------------
!***  Find the child tasks and their relevant limits
!***  along the child's eastern boundary.
!-------------------------------------------------------------
!
            DO NC=1,NUM_CHILD_TASKS                                        !<-- Loop through all tasks on child domain
!
              IF(J_CHILD>=J_LIMIT_LO(NC).AND.                           &  !<-- Does current child boundary point on this
                 J_CHILD<=J_LIMIT_HI(NC)                                &  !    parent task lie on child task "NC"?
                 .AND.                                                  &  !
                 I_CHILD>=ITS_CHILD(NC).AND.                            &  !
                 I_CHILD<=ITE_CHILD(NC))THEN
!
                IF(NC>NC_LAST_E)THEN                                       !<-- Have we encountered a new child task holding this E bndry?
                  NUM_TASKS_SEND_E=NUM_TASKS_SEND_E+1                      !<-- Then increment the E bndry counter of the child tasks
                  LOCAL_TASK_RANK_E(NUM_TASKS_SEND_E)=NC-1                 !<-- Save this child task's local rank
                  NC_LAST_E=NC
                  NC_HOLD_E(NC)=NUM_TASKS_SEND_E
                ENDIF
!
                IF(J_SAVE_LO_EAST(NC_HOLD_E(NC))<0)THEN
                  J_SAVE_LO_EAST(NC_HOLD_E(NC))=J_CHILD                    !<-- Save the starting J index on this child task's segment
                ENDIF
                J_SAVE_HI_EAST(NC_HOLD_E(NC))=J_CHILD                      !<-- Save the ending J index on this child task's segment
!
              ENDIF
!
            ENDDO
!
!-----------------------------------------------------------------------
!
            i_east: DO I_CHILD=IM_END-N_BLEND+1-N_ADD,IM_END               !<-- Blending region of child's eastern boundary
!
              KOUNT_I=KOUNT_I+1
              REAL_I_PARENT=REAL_I_START+(KOUNT_I-1)*RATIO                 !<-- REAL parent I for this child's I
!
              I_INDX_EBND(I_CHILD,J_CHILD,1)=INT(REAL_I_PARENT+EPS)        !<-- Parent I west of child's east boundary point
              I_INDX_EBND(I_CHILD,J_CHILD,2)=INT(REAL_I_PARENT+EPS)+1      !<-- Parent I east of child's east boundary point
              J_INDX_EBND(I_CHILD,J_CHILD,1)=INT(REAL_J_PARENT+EPS)        !<-- Parent J south of child's east boundary point
              J_INDX_EBND(I_CHILD,J_CHILD,2)=INT(REAL_J_PARENT+EPS)+1      !<-- Parent J north of child's east boundary point
!
              WEIGHT_SW=(I_INDX_EBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (J_INDX_EBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_SE=(REAL_I_PARENT-I_INDX_EBND(I_CHILD,J_CHILD,1))* &
                        (J_INDX_EBND(I_CHILD,J_CHILD,2)-REAL_J_PARENT)
              WEIGHT_NW=(I_INDX_EBND(I_CHILD,J_CHILD,2)-REAL_I_PARENT)* &
                        (REAL_J_PARENT-J_INDX_EBND(I_CHILD,J_CHILD,1))
              WEIGHT_NE=(REAL_I_PARENT-I_INDX_EBND(I_CHILD,J_CHILD,1))* &
                        (REAL_J_PARENT-J_INDX_EBND(I_CHILD,J_CHILD,1))
!
              RECIP_SUM=1./(WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE)
!
              WEIGHTS_EBND(I_CHILD,J_CHILD,INDX_SW)=WEIGHT_SW*RECIP_SUM    !<-- Interp wght of parent point SW of child bndry point
              WEIGHTS_EBND(I_CHILD,J_CHILD,INDX_SE)=WEIGHT_SE*RECIP_SUM    !<-- Interp wght of parent point SE of child bndry point
              WEIGHTS_EBND(I_CHILD,J_CHILD,INDX_NW)=WEIGHT_NW*RECIP_SUM    !<-- Interp wght of parent point NW of child bndry point
              WEIGHTS_EBND(I_CHILD,J_CHILD,INDX_NE)=WEIGHT_NE*RECIP_SUM    !<-- Interp wght of parent point NE of child bndry point
!
            ENDDO i_east
!
!-----------------------------------------------------------------------
          ENDIF child_ij_e
!-----------------------------------------------------------------------
!
        ENDIF j_block
!
!-----------------------------------------------------------------------
!
      ENDDO j_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_INTERP_FACTORS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PRELIM_CHILD_INFO(IMP_STATE,EXP_STATE)
!
!-----------------------------------------------------------------------
!***  PARENTS SEND CHILDREN BASIC INFORMATION NEEDED FOR THE
!***  EXCHANGE OF BOUNDARY DATA DURING THE INTEGRATION.
!
!***  CHILDREN SEND THEIR SURFACE GEOPOTENTIAL TO PARENTS SO THAT
!***  THE PARENTS CAN PROPERLY BALANCE THEIR MASS DATA THAT IS
!***  INTERPOLATED TO CHILD GRIDPOINTS WHOSE TERRAIN IS DIFFERENT
!***  FROM THAT OF THE PARENT.
!-----------------------------------------------------------------------
!
!--------------
!*** Arguments
!--------------
!
      TYPE(ESMF_State),INTENT(IN)    :: IMP_STATE                           !<-- Parent-Child Coupler import state
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE                           !<-- Parent-Child Coupler export state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: ID_CHILDTASK,ID_DOM,IERR,ISTAT,KOUNT,LENGTH,MYPE       &
                ,N,N1,N2,NBASE,NT,NTX,NUM_WORDS,RC,RC_PRELIM
!
      INTEGER :: I_LO,I_HI,I_OFFSET,ILIM_HI,ILIM_LO                     &
                ,J_LO,J_HI,J_OFFSET,JLIM_HI,JLIM_LO
!
      INTEGER,DIMENSION(4) :: INFO=(/-1,0,0,0/)
      INTEGER,DIMENSION(4) :: TEMP
      INTEGER,DIMENSION(4,2) :: PARENT_INFO
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
      REAL,DIMENSION(:),ALLOCATABLE :: FIS_SEND
!
      TYPE(ESMF_Array) :: HOLD_ARRAY
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  PARENTS SEND THREE PIECES OF INFORMATION TO CHILD TASKS 
!***  ON CHILD DOMAIN BOUNDARIES SO THOSE CHILD TASKS WILL BE
!***  ABLE TO RECEIVE BOUNDARY DATA AND USE IT PROPERLY:
!
!     (1) The parent task's local rank
!     (2) The starting index on the child's boundary
!     (3) The ending index on the child's boundary
!
!***  THE CHILD TASK MUST BE ABLE TO KNOW IF THE DATA IT RECEIVES
!***  PERTAINS TO SOUTH BOUNDARY H OR V POINTS, NORTH BOUNDARY
!***  H OR V, POINTS, ETC.  THUS THE MPI TAG WILL INDICATE
!***  THE BOUNDARY'S SIDE AND VARIABLE TYPE.
!
!      101 --> South H
!      102 --> South V
!      103 --> North H
!      104 --> North V
!      105 --> West H  
!      106 --> West V
!      107 --> East H  
!      108 --> East V
!
!*** (THE CHILD TASKS KNOW WHICH SIDE OF THEIR DOMAIN'S BOUNDARY THEY
!***  ARE ON BUT SINCE A CORNER TASK IS ON MORE THAN ONE SIDE, THE
!***  CODE INDICATING THE SIDE IS USED FOR ALL CHILD TASKS.)
!-----------------------------------------------------------------------
!
      parent_sends: IF(NUM_CHILDREN>0)THEN                                 !<-- Consider the parent tasks
!
        child_loop1: DO N=1,NUM_CHILDREN                                   !<-- Parent loops through all its children
!
!-----------------------------------------------------------------------
!
          CALL MPI_COMM_RANK(COMM_TO_MY_CHILDREN(N),MYPE,IERR)             !<-- Obtain rank of parent task
          ID_DOM=MY_CHILDREN_ID(N)
!
!-------------
!***  South H
!-------------
!
          sh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_H_S(N)                               !<-- Look for a child task with south boundary H points
                ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ south boundary H points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NTX)              !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NTX)     !<-- Save the ending index of boundary segment on child
                  INFO(4)=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NTX)              !<-- Save the ending index of expanded boundary segment on child
!
                  CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Sbndry task
                               ,101                                     &  !<-- Tag for south boundary H points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE sh_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no south boundary H points
            CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,101                                           &  !<-- Tag for south boundary H points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO sh_loop
!
!-------------
!***  South V
!-------------
!
          sv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_V_S(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_V_S(N)                               !<-- Look for a child task with south boundary V points
                ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%SOUTH(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ south boundary V points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_V_SAVE(N)%I_LO_SOUTH(NTX)              !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_V_SAVE(N)%I_HI_SOUTH(NTX)              !<-- Save the ending index of boundary segment on child
!
                  CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Sbndry task
                               ,102                                     &  !<-- Tag for south boundary V points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE sv_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no south boundary H points
            CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,102                                           &  !<-- Tag for south boundary V points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO sv_loop
!
!-------------
!***  North H
!-------------
!
          nh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_H_N(N)                               !<-- Look for a child task with north boundary H points
                ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ north boundary H points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_H_SAVE(N)%I_LO_NORTH(NTX)              !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NTX)     !<-- Save the ending index of boundary segment on child
                  INFO(4)=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NTX)              !<-- Save the ending index of expanded boundary segment on child
!
                  CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Nbndry task
                               ,103                                     &  !<-- Tag for north boundary H points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE nh_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no north boundary H points
            CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,103                                           &  !<-- Tag for north boundary H points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO nh_loop
!
!-------------
!***  North V
!-------------
!
          nv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_V_N(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_V_N(N)                               !<-- Look for a child task with north boundary V points
                ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%NORTH(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ north boundary V points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_V_SAVE(N)%I_LO_NORTH(NTX)              !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_V_SAVE(N)%I_HI_NORTH(NTX)              !<-- Save the ending index of boundary segment on child
!
                  CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Nbndry task
                               ,104                                     &  !<-- Tag for north boundary V points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE nv_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no north boundary V points
            CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,104                                           &  !<-- Tag for north boundary V points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO nv_loop
!
!------------
!***  West H
!------------
!
          wh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_H_W(N)                               !<-- Look for a child task with west boundary H points
                ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ west boundary H points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_H_SAVE(N)%J_LO_WEST(NTX)               !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NTX)      !<-- Save the ending index of boundary segment on child
                  INFO(4)=CHILDTASK_H_SAVE(N)%J_HI_WEST(NTX)               !<-- Save the ending index of expanded boundary segment on child
!
                  CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Wbndry task
                               ,105                                     &  !<-- Tag for west boundary H points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE wh_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no west boundary H points
            CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,105                                           &  !<-- Tag for west boundary H points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO wh_loop
!
!------------
!***  West V
!------------
!
          wv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_V_W(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_V_W(N)                               !<-- Look for a child task with west boundary V points
                ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%WEST(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ west boundary V points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_V_SAVE(N)%J_LO_WEST(NTX)               !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_V_SAVE(N)%J_HI_WEST(NTX)               !<-- Save the ending index of boundary segment on child
!
                  CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Wbndry task
                               ,106                                     &  !<-- Tag for west boundary V points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE wv_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no west boundary V points
            CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,106                                           &  !<-- Tag for west boundary V points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO wv_loop
!
!------------
!***  East H
!------------
!
          eh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_H_E(N)                               !<-- Look for a child task with east boundary H points
                ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ east boundary H points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_H_SAVE(N)%J_LO_EAST(NTX)               !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NTX)      !<-- Save the ending index of boundary segment on child
                  INFO(4)=CHILDTASK_H_SAVE(N)%J_HI_EAST(NTX)               !<-- Save the ending index of expanded boundary segment on child
!
                  CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Ebndry task
                               ,107                                     &  !<-- Tag for east boundary H points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE eh_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no east boundary H points
            CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,107                                           &  !<-- Tag for east boundary H points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO eh_loop
!
!------------
!***  East V
!------------
!
          ev_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                         !<-- Each parent task loops through all tasks on each child
!
            INFO(1)=-1
!
            IF(NUM_TASKS_SEND_V_E(N)>0)THEN                                
              DO NTX=1,NUM_TASKS_SEND_V_E(N)                               !<-- Look for a child task with east boundary V points
                ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%EAST(NTX)
!              
                IF(NT==ID_CHILDTASK)THEN                                   !<-- If yes, we found a child task w/ east boundary V points
                  INFO(1)=MYPE                                             !<-- Save the parent task rank
                  INFO(2)=CHILDTASK_V_SAVE(N)%J_LO_EAST(NTX)               !<-- Save the starting index of boundary segment on child
                  INFO(3)=CHILDTASK_V_SAVE(N)%J_HI_EAST(NTX)               !<-- Save the ending index of boundary segment on child
!
                  CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK         &  !<-- Parent task sends the key data to the child Ebndry task
                               ,108                                     &  !<-- Tag for east boundary V points
                               ,COMM_TO_MY_CHILDREN(N),IERR)
!
                  CYCLE ev_loop                                            !<-- Move on to next child task
                ENDIF
              ENDDO
            ENDIF
!
            ID_CHILDTASK=NT                                                !<-- This child task has no east boundary V points
            CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK               &  !<-- So send this child task some dummy information
                         ,108                                           &  !<-- Tag for east boundary V points
                         ,COMM_TO_MY_CHILDREN(N),IERR)
!
          ENDDO ev_loop
!
!-----------------------------------------------------------------------
!
        ENDDO child_loop1
!
!-----------------------------------------------------------------------
!
      ENDIF parent_sends
!
!-----------------------------------------------------------------------
!***  THE CHILD TASKS RECEIVE THE KEY INFORMATION FROM THEIR
!***  PARENT TASKS.  AT THIS POINT THE CHILD TASKS DO NOT KNOW
!***  THE LOCAL RANKS OF THE PARENT TASKS THAT WILL BE SENDING
!***  INFORMATION TO THEM THUS MPI_ANY_SOURCE IS USED IN THE
!***  RECEIVE.  HOWEVER THIS MEANS THAT WHEN THERE ARE TWO
!***  PARENT TASKS SENDING TO A NEST TASK (RATHER THAN ONLY ONE)
!***  THEN THE TWO OVERLAP POINTS ON THE NEST BOUNDARY SEGMENT
!***  COMPUTED BY THE PARENT TASKS WILL ULTIMATELY HAVE VALUES
!***  DEPENDING ON WHICH PARENT TASK'S PRELIMINARY INFORMATION
!***  IS RECEIVED LAST IN THE MPI_ANY_SOURCE RECV BELOW.  SINCE
!***  THE VALUES IN THOSE OVERLAP POINTS ARE NOT BIT IDENTICAL
!***  THEN ANY SUCCESSIVE RUNS CAN HAVE SLIGHTLY DIFFERENT ANSWERS.
!***  TO AVOID THAT HAPPENING WHEN TWO PARENT TASKS ARE SENDING
!***  THE CHILD TASK WILL RECEIVE THEIR RANKS AND THEN PUT THEM
!***  IN ASCENDING ORDER SO THAT ALL SUBSEQENT UPDATES OF THE
!***  NEST BOUNDARY OVERLAP POINTS ARE ALWAYS DONE IN THE SAME
!***  WAY REGARDLESS OF THE ORDER THE PRELIMINARY INFORMATION
!***  IS RECEIVED WITH MPI_ANY_SOURCE.
!-----------------------------------------------------------------------
!
      child_recvs: IF(COMM_TO_MY_PARENT>0)THEN                             !<-- Only children will receive
!
        ID_DOM=ID_PARENTS(MY_DOMAIN_ID)                                    !<-- Domain ID of this child's parent
        ALLOCATE(PARENT_TASK(1:FTASKS_DOMAIN(ID_DOM)))
!
!-------------
!***  South H
!-------------
!
        KOUNT=0
        INDX_MIN_H%SOUTH= 1000000
        INDX_MAX_H%SOUTH=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,4                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,101                                             &  !<-- Tag used for south boundary H points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task sent key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,4
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,4                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%SOUTH_H%ID_SOURCE   =PARENT_INFO(1,N)           !<-- Rank of parent task that will send Sboundary H data segment
            PARENT_TASK(N)%SOUTH_H%INDX_START  =PARENT_INFO(2,N)           !<-- Istart on child grid of the boundary data segment
            PARENT_TASK(N)%SOUTH_H%INDX_END    =PARENT_INFO(3,N)           !<-- Iend on child grid of the boundary data segment
            PARENT_TASK(N)%SOUTH_H%INDX_END_EXP=PARENT_INFO(4,N)           !<-- Iend on child grid of the expanded boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
            LENGTH=(3*LM+1)*NBASE
!
            PARENT_TASK(N)%SOUTH_H%LENGTH=LENGTH                           !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%SOUTH_H%STRING(1:LENGTH))              !<-- Sboundary H datastring to be received from parent task
!
            INDX_MIN_H%SOUTH=MIN(INDX_MIN_H%SOUTH,PARENT_INFO(2,N))        !<-- Starting child I for union of parent task segments sent  
            INDX_MAX_H%SOUTH=MAX(INDX_MAX_H%SOUTH,PARENT_INFO(3,N))        !<-- Ending child I for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_H%SOUTH=KOUNT
!!!     LENGTH_BND_SEG_H%SOUTH=INDX_MAX_H%SOUTH-INDX_MIN_H%SOUTH+1
!
        south_h: IF(NUM_PARENT_TASKS_SENDING_H%SOUTH>0)THEN                    !<-- Does this child task recv Sboundary H data from parent?
!
          ALLOCATE(PDB_S(INDX_MIN_H%SOUTH:INDX_MAX_H%SOUTH,1:N_BLEND_H))       !<-- Full PDB south H boundary segment on this child task
          ALLOCATE( TB_S(INDX_MIN_H%SOUTH:INDX_MAX_H%SOUTH,1:N_BLEND_H,1:LM))  !<-- Full TB south H boundary segment on this child task
          ALLOCATE( QB_S(INDX_MIN_H%SOUTH:INDX_MAX_H%SOUTH,1:N_BLEND_H,1:LM))  !<-- Full QB south H boundary segment on this child task
          ALLOCATE(CWB_S(INDX_MIN_H%SOUTH:INDX_MAX_H%SOUTH,1:N_BLEND_H,1:LM))  !<-- Full CWB south H boundary segment on this child task
!
          ILIM_LO=INDX_MIN_H%SOUTH
          ILIM_HI=INDX_MAX_H%SOUTH
          JLIM_LO=1
          JLIM_HI=N_BLEND_H
!
          LENGTH=(3*LM+1)*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_SOUTH_H(1:LENGTH),stat=ISTAT)                  !<-- 1-D combined H-point data on child task's Sbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_SOUTH_H'               &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_SOUTH_H'               &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_SOUTH_H'               &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_SOUTH_H'               &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF south_h
!
!-------------
!***  South V
!-------------
!
        KOUNT=0
        INDX_MIN_V%SOUTH= 1000000
        INDX_MAX_V%SOUTH=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,3                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,102                                             &  !<-- Tag used for south boundary V points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,3
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,3                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%SOUTH_V%ID_SOURCE =PARENT_INFO(1,N)             !<-- Rank of parent task that will send Sboundary V data segment
            PARENT_TASK(N)%SOUTH_V%INDX_START=PARENT_INFO(2,N)             !<-- Istart on child grid of the boundary data segment
            PARENT_TASK(N)%SOUTH_V%INDX_END  =PARENT_INFO(3,N)             !<-- Iend on child grid of the boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
            LENGTH=2*LM*NBASE
!
            PARENT_TASK(N)%SOUTH_V%LENGTH=LENGTH                           !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%SOUTH_V%STRING(1:LENGTH))              !<-- Sboundary V datastring to be received from parent task
!
            INDX_MIN_V%SOUTH=MIN(INDX_MIN_V%SOUTH,PARENT_INFO(2,N))        !<-- Starting child I for union of parent task segments sent
            INDX_MAX_V%SOUTH=MAX(INDX_MAX_V%SOUTH,PARENT_INFO(3,N))        !<-- Ending child I for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_V%SOUTH=KOUNT
!!!     LENGTH_BND_SEG_V%SOUTH=INDX_MAX_V%SOUTH-INDX_MIN_V%SOUTH+1
!
        south_v: IF(NUM_PARENT_TASKS_SENDING_V%SOUTH>0)THEN                !<-- Does this child task recv any Sboundary V data from parent?
!
          ALLOCATE(UB_S(INDX_MIN_V%SOUTH:INDX_MAX_V%SOUTH,1:N_BLEND_H,1:LM))
          ALLOCATE(VB_S(INDX_MIN_V%SOUTH:INDX_MAX_V%SOUTH,1:N_BLEND_H,1:LM))
!
          ILIM_LO=INDX_MIN_V%SOUTH
          ILIM_HI=INDX_MAX_V%SOUTH
          JLIM_LO=1
          JLIM_HI=N_BLEND_V
!
          LENGTH=2*LM*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_SOUTH_V(1:LENGTH),stat=ISTAT)                  !<-- 1-D combined V-point data on child task's Sbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_SOUTH_V'               &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_SOUTH_V'               &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_SOUTH_V'               &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_SOUTH_V'               &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF south_v
!
!-------------
!***  North H
!-------------
!
        KOUNT=0
        INDX_MIN_H%NORTH= 1000000
        INDX_MAX_H%NORTH=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,4                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,103                                             &  !<-- Tag used for north boundary H points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,4
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,4                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%NORTH_H%ID_SOURCE   =PARENT_INFO(1,N)           !<-- Rank of parent task that will send Nboundary H data segment
            PARENT_TASK(N)%NORTH_H%INDX_START  =PARENT_INFO(2,N)           !<-- Istart on child grid of the boundary data segment
            PARENT_TASK(N)%NORTH_H%INDX_END    =PARENT_INFO(3,N)           !<-- Iend on child grid of the boundary data segment
            PARENT_TASK(N)%NORTH_H%INDX_END_EXP=PARENT_INFO(4,N)           !<-- Iend on child grid of the expanded boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
            LENGTH=(3*LM+1)*NBASE
!
            PARENT_TASK(N)%NORTH_H%LENGTH=LENGTH                           !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%NORTH_H%STRING(1:LENGTH))              !<-- Nboundary H datastring to be received from parent task
!
            INDX_MIN_H%NORTH=MIN(INDX_MIN_H%NORTH,PARENT_INFO(2,N))        !<-- Starting child I for union of parent task segments sent
            INDX_MAX_H%NORTH=MAX(INDX_MAX_H%NORTH,PARENT_INFO(3,N))        !<-- Ending child I for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_H%NORTH=KOUNT
!!!     LENGTH_BND_SEG_H%NORTH=INDX_MAX_H%NORTH-INDX_MIN_H%NORTH+1
!
        north_h: IF(NUM_PARENT_TASKS_SENDING_H%NORTH>0)THEN                    !<-- Does this child task recv Nboundary H data from parent?
!
          ALLOCATE(PDB_N(INDX_MIN_H%NORTH:INDX_MAX_H%NORTH,1:N_BLEND_H))
          ALLOCATE( TB_N(INDX_MIN_H%NORTH:INDX_MAX_H%NORTH,1:N_BLEND_H,1:LM))
          ALLOCATE( QB_N(INDX_MIN_H%NORTH:INDX_MAX_H%NORTH,1:N_BLEND_H,1:LM))
          ALLOCATE(CWB_N(INDX_MIN_H%NORTH:INDX_MAX_H%NORTH,1:N_BLEND_H,1:LM))
!
!
          ILIM_LO=INDX_MIN_H%NORTH
          ILIM_HI=INDX_MAX_H%NORTH
          JLIM_LO=1
          JLIM_HI=N_BLEND_H
!
          LENGTH=(3*LM+1)*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_NORTH_H(1:LENGTH),stat=ISTAT)                  !<-- 1-D combined H-point data on child task's Nbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_NORTH_H'               &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_NORTH_H'               &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_NORTH_H'               &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_NORTH_H'               &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF north_h
!
!-------------
!***  North V
!-------------
!
        KOUNT=0
        INDX_MIN_V%NORTH= 1000000
        INDX_MAX_V%NORTH=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,3                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,104                                             &  !<-- Tag used for north boundary V points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,3
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,3                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%NORTH_V%ID_SOURCE =PARENT_INFO(1,N)             !<-- Rank of parent task that will send Nboundary V data segment
            PARENT_TASK(N)%NORTH_V%INDX_START=PARENT_INFO(2,N)             !<-- Istart on child grid of the boundary data segment
            PARENT_TASK(N)%NORTH_V%INDX_END  =PARENT_INFO(3,N)             !<-- Iend on child grid of the boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
            LENGTH=2*LM*NBASE
!
            PARENT_TASK(N)%NORTH_V%LENGTH=LENGTH                           !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%NORTH_V%STRING(1:LENGTH))              !<-- Nboundary V datastring to be received from parent task
!
            INDX_MIN_V%NORTH=MIN(INDX_MIN_V%NORTH,PARENT_INFO(2,N))        !<-- Starting child I for union of parent task segments sent
            INDX_MAX_V%NORTH=MAX(INDX_MAX_V%NORTH,PARENT_INFO(3,N))        !<-- Ending child I for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_V%NORTH=KOUNT
!!!     LENGTH_BND_SEG_V%NORTH=INDX_MAX_V%NORTH-INDX_MIN_V%NORTH+1
!
        north_v: IF(NUM_PARENT_TASKS_SENDING_V%NORTH>0)THEN                !<-- Does this child task recv any Nboundary V data from parent?
!
          ALLOCATE(UB_N(INDX_MIN_V%NORTH:INDX_MAX_V%NORTH,1:N_BLEND_V,1:LM))
          ALLOCATE(VB_N(INDX_MIN_V%NORTH:INDX_MAX_V%NORTH,1:N_BLEND_V,1:LM))
!
          ILIM_LO=INDX_MIN_V%NORTH
          ILIM_HI=INDX_MAX_V%NORTH
          JLIM_LO=1
          JLIM_HI=N_BLEND_V
!
          LENGTH=2*LM*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_NORTH_V(1:LENGTH),stat=ISTAT)                  !<-- 1-D combined V-point data on child task's Nbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_NORTH_V'               &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_NORTH_V'               &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_NORTH_V'               &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_NORTH_V'               &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF north_v
!
!------------
!***  West H
!------------
!
        KOUNT=0
        INDX_MIN_H%WEST= 1000000
        INDX_MAX_H%WEST=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,4                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,105                                             &  !<-- Tag used for west boundary H points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,4
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,4                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%WEST_H%ID_SOURCE   =PARENT_INFO(1,N)            !<-- Rank of parent task that will send Wboundary H data segment
            PARENT_TASK(N)%WEST_H%INDX_START  =PARENT_INFO(2,N)            !<-- Jstart on child grid of the boundary data segment
            PARENT_TASK(N)%WEST_H%INDX_END    =PARENT_INFO(3,N)            !<-- Jend on child grid of the boundary data segment
            PARENT_TASK(N)%WEST_H%INDX_END_EXP=PARENT_INFO(4,N)            !<-- Jend on child grid of the expanded boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
            LENGTH=(3*LM+1)*NBASE
!
            PARENT_TASK(N)%WEST_H%LENGTH=LENGTH                            !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%WEST_H%STRING(1:LENGTH))               !<-- Wboundary H datastring to be received from parent task
!
            INDX_MIN_H%WEST=MIN(INDX_MIN_H%WEST,PARENT_INFO(2,N))          !<-- Starting child J for union of parent task segments sent
            INDX_MAX_H%WEST=MAX(INDX_MAX_H%WEST,PARENT_INFO(3,N))          !<-- Ending child J for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_H%WEST=KOUNT
!!!     LENGTH_BND_SEG_H%WEST=INDX_MAX_H%WEST-INDX_MIN_H%WEST+1
!
        west_h: IF(NUM_PARENT_TASKS_SENDING_H%WEST>0)THEN                  !<-- Does this child task recv Wboundary H data from parent?
!
          ALLOCATE(PDB_W(1:N_BLEND_H,INDX_MIN_H%WEST:INDX_MAX_H%WEST))
          ALLOCATE( TB_W(1:N_BLEND_H,INDX_MIN_H%WEST:INDX_MAX_H%WEST,1:LM))
          ALLOCATE( QB_W(1:N_BLEND_H,INDX_MIN_H%WEST:INDX_MAX_H%WEST,1:LM))
          ALLOCATE(CWB_W(1:N_BLEND_H,INDX_MIN_H%WEST:INDX_MAX_H%WEST,1:LM))
!
          ILIM_LO=1
          ILIM_HI=N_BLEND_H
          JLIM_LO=INDX_MIN_H%WEST
          JLIM_HI=INDX_MAX_H%WEST
!
          LENGTH=(3*LM+1)*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_WEST_H(1:LENGTH),stat=ISTAT)                   !<-- 1-D combined H-point data on child task's Wbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_WEST_H'                &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_WEST_H'                &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_WEST_H'                &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_WEST_H'                &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF west_h
!
!------------
!***  West V
!------------
!
        KOUNT=0
        INDX_MIN_V%WEST= 1000000
        INDX_MAX_V%WEST=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,3                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,106                                             &  !<-- Tag used for west boundary V points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,3
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,3                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%WEST_V%ID_SOURCE =PARENT_INFO(1,N)              !<-- Rank of parent task that will send Wboundary V data segment
            PARENT_TASK(N)%WEST_V%INDX_START=PARENT_INFO(2,N)              !<-- Jstart on child grid of the boundary data segment
            PARENT_TASK(N)%WEST_V%INDX_END  =PARENT_INFO(3,N)              !<-- Jend on child grid of the boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
            LENGTH=2*LM*NBASE
!
            PARENT_TASK(N)%WEST_V%LENGTH=LENGTH                            !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%WEST_V%STRING(1:LENGTH))               !<-- Wboundary V datastring to be received from parent task
!
            INDX_MIN_V%WEST=MIN(INDX_MIN_V%WEST,PARENT_INFO(2,N))          !<-- Starting child J for union of parent task segments sent
            INDX_MAX_V%WEST=MAX(INDX_MAX_V%WEST,PARENT_INFO(3,N))          !<-- Ending child J for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_V%WEST=KOUNT
!!!     LENGTH_BND_SEG_V%WEST=INDX_MAX_V%WEST-INDX_MIN_V%WEST+1
!
        west_v: IF(NUM_PARENT_TASKS_SENDING_V%WEST>0)THEN                  !<-- Does this child task recv Wboundary V data from parent?
!
          ALLOCATE(UB_W(1:N_BLEND_H,INDX_MIN_V%WEST:INDX_MAX_V%WEST,1:LM))
          ALLOCATE(VB_W(1:N_BLEND_H,INDX_MIN_V%WEST:INDX_MAX_V%WEST,1:LM))
!
          ILIM_LO=1
          ILIM_HI=N_BLEND_H
          JLIM_LO=INDX_MIN_V%WEST
          JLIM_HI=INDX_MAX_V%WEST
!
          LENGTH=2*LM*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_WEST_V(1:LENGTH),stat=ISTAT)                   !<-- 1-D combined V-point data on child task's Wbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_WEST_V'                &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_WEST_V'                &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_WEST_V'                &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_WEST_V'                &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF west_v
!
!------------
!***  East H
!------------
!
        KOUNT=0
        INDX_MIN_H%EAST= 1000000
        INDX_MAX_H%EAST=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,4                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,107                                             &  !<-- Tag used for east boundary H points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,4
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,4                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
            PARENT_TASK(N)%EAST_H%ID_SOURCE   =PARENT_INFO(1,N)            !<-- Rank of parent task that will send Eboundary H data segment
            PARENT_TASK(N)%EAST_H%INDX_START  =PARENT_INFO(2,N)            !<-- Jstart on child grid of the boundary data segment
            PARENT_TASK(N)%EAST_H%INDX_END    =PARENT_INFO(3,N)            !<-- Jend on child grid of the boundary data segment
            PARENT_TASK(N)%EAST_H%INDX_END_EXP=PARENT_INFO(4,N)            !<-- Jend on child grid of the expanded boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
            LENGTH=(3*LM+1)*NBASE
!
            PARENT_TASK(N)%EAST_H%LENGTH=LENGTH                            !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%EAST_H%STRING(1:LENGTH))               !<-- Eboundary H datastring to be received from parent task
!
            INDX_MIN_H%EAST=MIN(INDX_MIN_H%EAST,PARENT_INFO(2,N))          !<-- Starting child J for union of parent task segments sent
            INDX_MAX_H%EAST=MAX(INDX_MAX_H%EAST,PARENT_INFO(3,N))          !<-- Ending child J for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_H%EAST=KOUNT
!!!     LENGTH_BND_SEG_H%EAST=INDX_MAX_H%EAST-INDX_MIN_H%EAST+1
!
        east_h: IF(NUM_PARENT_TASKS_SENDING_H%EAST>0)THEN                  !<-- Does this child task recv Eboundary H data from parent?
!
          ALLOCATE(PDB_E(1:N_BLEND_H,INDX_MIN_H%EAST:INDX_MAX_H%EAST))
          ALLOCATE( TB_E(1:N_BLEND_H,INDX_MIN_H%EAST:INDX_MAX_H%EAST,1:LM))
          ALLOCATE( QB_E(1:N_BLEND_H,INDX_MIN_H%EAST:INDX_MAX_H%EAST,1:LM))
          ALLOCATE(CWB_E(1:N_BLEND_H,INDX_MIN_H%EAST:INDX_MAX_H%EAST,1:LM))
!
          ILIM_LO=1
          ILIM_HI=N_BLEND_H
          JLIM_LO=INDX_MIN_H%EAST
          JLIM_HI=INDX_MAX_H%EAST
!
          LENGTH=(3*LM+1)*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_EAST_H(1:LENGTH),stat=ISTAT)                   !<-- 1-D combined H-point data on child task's Ebndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_EAST_H'                &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_EAST_H'                &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_EAST_H'                &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_EAST_H'                &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF east_h
!
!------------
!***  East V
!------------
!
        KOUNT=0
        INDX_MIN_V%EAST= 1000000
        INDX_MAX_V%EAST=-1000000
!
        DO N=1,FTASKS_DOMAIN(ID_DOM)                                       !<-- Child task loops through its parent's tasks
!
          CALL MPI_RECV(INFO                                            &  !<-- Receive data packet from each parent task
                       ,3                                               &  !<-- # of words in data packet
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Accept data from any parent task that is sending
                       ,108                                             &  !<-- Tag used for east boundary V points
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between child and parent
                       ,STATUS                                          &  !<-- Status of Recv
                       ,IERR)
!
          IF(INFO(1)>=0)THEN                                               !<-- If yes, this parent task has key preliminary bndry info
            KOUNT=KOUNT+1
            DO N1=1,3
              PARENT_INFO(N1,KOUNT)=INFO(N1)                               !<-- Save the data from that parent task
            ENDDO
          ENDIF
!
        ENDDO
!
        IF(KOUNT==2)THEN                                                   !<-- Nest task recvs data from two parent tasks
          IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                        !<-- Data recvd from 'out of order' parent tasks
!
            DO N1=1,3                                                      !<-- Save parent data in order of ascending task IDs
              TEMP(N1)         =PARENT_INFO(N1,1)                          !
              PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                          !
              PARENT_INFO(N1,2)=TEMP(N1)                                   !<--
            ENDDO
!
          ENDIF
        ENDIF
!
        IF(KOUNT>0)THEN
          DO N=1,KOUNT
!
            PARENT_TASK(N)%EAST_V%ID_SOURCE =PARENT_INFO(1,N)              !<-- Rank of parent task that will send Eboundary V data segment
            PARENT_TASK(N)%EAST_V%INDX_START=PARENT_INFO(2,N)              !<-- Jstart on child grid of the boundary data segment
            PARENT_TASK(N)%EAST_V%INDX_END  =PARENT_INFO(3,N)              !<-- Jend on child grid of the boundary data segment
!
            NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
            LENGTH=2*LM*NBASE
!
            PARENT_TASK(N)%EAST_V%LENGTH=LENGTH                            !<-- # of words in this parent task's datastring of child bndry
!
            ALLOCATE(PARENT_TASK(N)%EAST_V%STRING(1:LENGTH))               !<-- Eboundary V datastring to be received from parent task
!
            INDX_MIN_V%EAST=MIN(INDX_MIN_V%EAST,PARENT_INFO(2,N))          !<-- Starting child J for union of parent task segments sent
            INDX_MAX_V%EAST=MAX(INDX_MAX_V%EAST,PARENT_INFO(3,N))          !<-- Ending child J for union of parent task segments sent
          ENDDO
        ENDIF
!
        NUM_PARENT_TASKS_SENDING_V%EAST=KOUNT
!!!     LENGTH_BND_SEG_V%EAST=INDX_MAX_V%EAST-INDX_MIN_V%EAST+1
!
        east_v: IF(NUM_PARENT_TASKS_SENDING_V%EAST>0)THEN                  !<-- Does this child task recv Eboundary V data from parent?
!
          ALLOCATE(UB_E(1:N_BLEND_H,INDX_MIN_V%EAST:INDX_MAX_V%EAST,1:LM))
          ALLOCATE(VB_E(1:N_BLEND_H,INDX_MIN_V%EAST:INDX_MAX_V%EAST,1:LM))
!
          ILIM_LO=1
          ILIM_HI=N_BLEND_H
          JLIM_LO=INDX_MIN_V%EAST
          JLIM_HI=INDX_MAX_V%EAST
!
          LENGTH=2*LM*(ILIM_HI-ILIM_LO+1)*(JLIM_HI-JLIM_LO+1)
          ALLOCATE(BOUND_1D_EAST_V(1:LENGTH),stat=ISTAT)                   !<-- 1-D combined V-point data on child task's Ebndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_LO_EAST_V'                &  !<-- Name of the boundary array's lower I limit
                                ,value=ILIM_LO                         &  !<-- The boundary array's lower I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='ILIM_HI_EAST_V'                &  !<-- Name of the boundary array's upper I limit
                                ,value=ILIM_HI                         &  !<-- The boundary array's upper I limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_LO_EAST_V'                &  !<-- Name of the boundary array's lower J limit
                                ,value=JLIM_LO                         &  !<-- The boundary array's lower J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                       &  !<-- The Parent_Child Coupler export state
                                ,name ='JLIM_HI_EAST_V'                &  !<-- Name of the boundary array's upper J limit
                                ,value=JLIM_HI                         &  !<-- The boundary array's upper J limit
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF east_v
!
!-----------------------------------------------------------------------
!
      ENDIF child_recvs
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  NOW THE CHILDREN SEND THEIR SFC GEOPOTENTIAL TO THEIR PARENTS
!***  SO THE PARENTS CAN PROPERLY BALANCE THEIR OWN DATA THAT THEY
!***  INTERPOLATE TO CHILD GRIDPOINTS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      child_sends: IF(COMM_TO_MY_PARENT>0)THEN                             !<-- Select the nests' tasks
!
!-----------------------------------------------------------------------
!***  EXTRACT THE SFC GEOPOTENTIAL FROM THE COUPLER'S IMPORT STATE.
!***  IF THIS CHILD DOMAIN IS ALSO A PARENT THEN IT ALREADY EXTRACTED
!***  ITS FIS IN child_block of PARENT_CHILD_CPL_INITIALIZE BUT WE
!***  NOW EXTRACT FIS AGAIN IN CASE THIS CHILD DOMAIN IS NOT A PARENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract FIS from Parent-Child Coupler Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state   =IMP_STATE                           &  !<-- The parent-child coupler import state
                          ,itemName='FIS'                               &  !<-- Extract FIS
                          ,array   =HOLD_ARRAY                          &  !<-- Put the extracted Array here
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract FIS from ESMF Array in Parent-Child Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ArrayGet(array    =HOLD_ARRAY                         &  !<-- Array that holds the data pointer
                          ,localDe  =0                                  &
                          ,farrayPtr=FIS                                &  !<-- Put the pointer here
                          ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_PRELIM)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        CALL MPI_COMM_RANK(COMM_TO_MY_PARENT,MYPE,IERR)                    !<-- Obtain rank of this child task
!
!-----------------------------------------------------------------------
!
!!!     I_OFFSET=IMS+NHALO                                                 !<-- Offset of I in unloaded FIS vs. original FIS
!!!     J_OFFSET=JMS+NHALO                                                 !<-- Offset of J in unloaded FIS vs. original FIS
        I_OFFSET=IMS-1+NHALO                                               !<-- Offset of I in unloaded FIS vs. original FIS
        J_OFFSET=JMS-1+NHALO                                               !<-- Offset of J in unloaded FIS vs. original FIS
!
!------------------------------
!***  Child South Boundary FIS
!------------------------------
!
        IF(NUM_PARENT_TASKS_SENDING_H%SOUTH>0)THEN                         !<-- Child tasks know which parent tasks compute their BCs
!
          DO N=1,NUM_PARENT_TASKS_SENDING_H%SOUTH                          !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
            I_LO=PARENT_TASK(N)%SOUTH_H%INDX_START                         !<-- Starting I of child covered by parent task N
            I_HI=PARENT_TASK(N)%SOUTH_H%INDX_END_EXP                       !<-- Ending I of child for expanded area covered by parent task N
!!!         NUM_WORDS=(I_HI-I_LO+1)*N_BLEND_H                              !<-- # of child points covered by parent task N
            NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                          !<-- # of child points covered by parent task N
            ALLOCATE(FIS_SEND(1:NUM_WORDS))                                !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
            KOUNT=0
!
!!!         DO N2=1,N_BLEND_H
            DO N2=1,N_BLEND_H+1                                            !<-- Extra row for 4-pt interp of PD to V pts
            DO N1=I_LO,I_HI
              KOUNT=KOUNT+1
              FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                 !<-- 1-D FIS of child points covered by parent task N
            ENDDO
            ENDDO
!
            CALL MPI_SEND(FIS_SEND                                      &  !<-- Send FIS data to parent task N
                         ,NUM_WORDS                                     &   
                         ,MPI_REAL                                      &
                         ,PARENT_TASK(N)%SOUTH_H%ID_SOURCE              &
                         ,MYPE                                          &  
                         ,COMM_TO_MY_PARENT                             &
                         ,IERR)
!
            DEALLOCATE(FIS_SEND)
!
          ENDDO
!
        ENDIF
!
!------------------------------
!***  Child North Boundary FIS
!------------------------------
!
        IF(NUM_PARENT_TASKS_SENDING_H%NORTH>0)THEN                         !<-- Child tasks know which parent tasks compute their BCs
!                                                                               
          DO N=1,NUM_PARENT_TASKS_SENDING_H%NORTH                          !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
            I_LO=PARENT_TASK(N)%NORTH_H%INDX_START                         !<-- Starting I of child covered by parent task N
            I_HI=PARENT_TASK(N)%NORTH_H%INDX_END_EXP                       !<-- Ending I of child for expanded area covered by parent task N
!!!         NUM_WORDS=(I_HI-I_LO+1)*N_BLEND_H                              !<-- # of child points covered by parent task N
            NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                          !<-- # of child points covered by parent task N
            ALLOCATE(FIS_SEND(1:NUM_WORDS))                                !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
            KOUNT=0
!
!!!         DO N2=JTE-N_BLEND_H+1,JTE
            DO N2=JTE-N_BLEND_H  ,JTE                                      !<-- Extra row for 4-pt interp of PD to V pts
            DO N1=I_LO,I_HI
              KOUNT=KOUNT+1
              FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                 !<-- 1-D FIS of child points covered by parent task N
            ENDDO
            ENDDO
!
            CALL MPI_SEND(FIS_SEND                                      &  !<-- Send FIS data to parent task N
                         ,NUM_WORDS                                     &
                         ,MPI_REAL                                      &
                         ,PARENT_TASK(N)%NORTH_H%ID_SOURCE              &
                         ,MYPE                                          &
                         ,COMM_TO_MY_PARENT                             &
                         ,IERR)
!
            DEALLOCATE(FIS_SEND)
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        I_OFFSET=IMS-1+NHALO                                               !<-- Offset of I in unloaded FIS vs. original FIS
        J_OFFSET=JMS-1+NHALO                                               !<-- Offset of J in unloaded FIS vs. original FIS
!
!-----------------------------
!***  Child West Boundary FIS
!-----------------------------
!
        IF(NUM_PARENT_TASKS_SENDING_H%WEST>0)THEN                          !<-- Child tasks know which parent tasks compute their BCs
!                                                                          
          DO N=1,NUM_PARENT_TASKS_SENDING_H%WEST                           !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
            J_LO=PARENT_TASK(N)%WEST_H%INDX_START                          !<-- Starting J of child covered by parent task N
            J_HI=PARENT_TASK(N)%WEST_H%INDX_END_EXP                        !<-- Ending J of child for expanded area covered by parent task N
!!!         NUM_WORDS=(J_HI-J_LO+1)*N_BLEND_H                              !<-- # of child points covered by parent task N
            NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                          !<-- # of child points covered by parent task N
            ALLOCATE(FIS_SEND(1:NUM_WORDS))                                !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
            KOUNT=0
!
            DO N2=J_LO,J_HI
!!!         DO N1=1,N_BLEND_H  
            DO N1=1,N_BLEND_H+1                                            !<-- Extra row for 4-pt interp of PD to V pts
              KOUNT=KOUNT+1
              FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                 !<-- 1-D FIS of child points covered by parent task N
            ENDDO
            ENDDO
!
            CALL MPI_SEND(FIS_SEND                                      &  !<-- Send FIS data to parent task N
                         ,NUM_WORDS                                     &
                         ,MPI_REAL                                      &
                         ,PARENT_TASK(N)%WEST_H%ID_SOURCE               &
                         ,MYPE                                          &
                         ,COMM_TO_MY_PARENT                             &
                         ,IERR)
!
            DEALLOCATE(FIS_SEND)
          ENDDO
!
        ENDIF
!
!-----------------------------
!***  Child East Boundary FIS
!-----------------------------
!
        IF(NUM_PARENT_TASKS_SENDING_H%EAST>0)THEN                          !<-- Child tasks know which parent tasks compute their BCs
!
          DO N=1,NUM_PARENT_TASKS_SENDING_H%EAST                           !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
            J_LO=PARENT_TASK(N)%EAST_H%INDX_START                          !<-- Starting J of child covered by parent task N
            J_HI=PARENT_TASK(N)%EAST_H%INDX_END_EXP                        !<-- Ending J of child for expanded area covered by parent task N
!!!         NUM_WORDS=(J_HI-J_LO+1)*N_BLEND_H                              !<-- # of child points covered by parent task N
            NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                          !<-- # of child points covered by parent task N
            ALLOCATE(FIS_SEND(1:NUM_WORDS))                                !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
            KOUNT=0
!
            DO N2=J_LO,J_HI
!!!         DO N1=ITE-N_BLEND_H+1,ITE
            DO N1=ITE-N_BLEND_H  ,ITE                                      !<-- Extra row for 4-pt interp of PD to V pts
              KOUNT=KOUNT+1
              FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                 !<-- 1-D FIS of child points covered by parent task N
            ENDDO
            ENDDO
!
            CALL MPI_SEND(FIS_SEND                                      &  !<-- Send FIS data to parent task N
                         ,NUM_WORDS                                     &
                         ,MPI_REAL                                      &
                         ,PARENT_TASK(N)%EAST_H%ID_SOURCE               &
                         ,MYPE                                          &
                         ,COMM_TO_MY_PARENT                             &
                         ,IERR)
!
            DEALLOCATE(FIS_SEND)
          ENDDO
!
        ENDIF
!
      ENDIF child_sends
!
!-----------------------------------------------------------------------
!***  FINALLY THE PARENT TASKS RECEIVE SFC GEOPOTENTIALS FROM THE
!***  CHILD BOUNDARY TASKS.
!-----------------------------------------------------------------------
!
      parent_recvs: IF(NUM_CHILDREN>0)THEN                                 !<-- Select the parents' tasks
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE POINTERS THAT WILL HOLD THE SURFACE GEOPOTENTIAL
!***  OF CHILD TASKS ON EACH SIDE OF THE CHILD BOUNDARIES THAT WILL
!***  BE SENT TO THE APPROPRIATE PARENT TASKS.
!-----------------------------------------------------------------------
!
        ALLOCATE(FIS_CHILD_SOUTH(1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_NORTH(1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_WEST (1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_EAST (1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!
        child_loop2: DO N=1,NUM_CHILDREN
!
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child South Boundary FIS
!------------------------------
!
          IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                  !<-- This parent task covers some child Sboundary H points
!
            ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1:NUM_TASKS_SEND_H_S(N)))    !<-- FIS data slot for each Sbndry child task on parent task
!
            DO NTX=1,NUM_TASKS_SEND_H_S(N)                                 !<-- Loop through those particular child tasks 
!                                                                          
              I_LO=CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NTX)                     !<-- Starting I of child point bndry segment on this parent task
              I_HI=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NTX)                     !<-- Ending I of child point bndry segment on this parent task
!!!           NUM_WORDS=(I_HI-I_LO+1)*N_BLEND_H                            !<-- # of child points in its bndry segment on parent task
              NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                        !<-- # of child points in its bndry segment on parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
              ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))    !<-- FIS on child S bndry covered by this parent task
!
              CALL MPI_RECV(FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA          &  !<-- Recv FIS values from Sbndry child task NTX
                           ,NUM_WORDS                                   &  !<-- # of FIS values
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)       &  !<-- The child task sending
                           ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)       &  !<-- Tag
                           ,COMM_TO_MY_CHILDREN(N)                      &  !<-- MPI communicator
                           ,STATUS                                      &
                           ,IERR)
!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child North Boundary FIS
!------------------------------
!
          IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                  !<-- This parent task covers some child Nboundary H points
!
            ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1:NUM_TASKS_SEND_H_N(N)))    !<-- FIS data slot for each Nbndry child task on parent task
!
            DO NTX=1,NUM_TASKS_SEND_H_N(N)                                 !<-- Loop through those particular child tasks 
!                                                                          
              I_LO=CHILDTASK_H_SAVE(N)%I_LO_NORTH(NTX)                     !<-- Starting I of child point bndry segment on this parent task
              I_HI=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NTX)                     !<-- Ending I of child point bndry segment on this parent task
!!!           NUM_WORDS=(I_HI-I_LO+1)*N_BLEND_H                            !<-- # of child points in its bndry segment on parent task
              NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                        !<-- # of child points in its bndry segment on parent task
              ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))    !<-- FIS on child N bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
              CALL MPI_RECV(FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA          &  !<-- Recv FIS values from Nbndry child task NTX
                           ,NUM_WORDS                                   &  !<-- # of FIS values
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)       &  !<-- The child task sending
                           ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)       &  !<-- Tag
                           ,COMM_TO_MY_CHILDREN(N)                      &  !<-- MPI communicator
                           ,STATUS                                      &
                           ,IERR)


!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child West Boundary FIS
!-----------------------------
!
          IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                  !<-- This parent task covers some child Wboundary H points
!
            ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1:NUM_TASKS_SEND_H_W(N)))     !<-- FIS data slot for each Wbndry child task on parent task
!
            DO NTX=1,NUM_TASKS_SEND_H_W(N)                                 !<-- Loop through those particular child tasks 
!                                                                          
              J_LO=CHILDTASK_H_SAVE(N)%J_LO_WEST(NTX)                      !<-- Starting J of child point bndry segment on this parent task
              J_HI=CHILDTASK_H_SAVE(N)%J_HI_WEST(NTX)                      !<-- Ending J of child point bndry segment on this parent task
!!!           NUM_WORDS=(J_HI-J_LO+1)*N_BLEND_H                            !<-- # of child points in its bndry segment on parent task
              NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                        !<-- # of child points in its bndry segment on parent task
              ALLOCATE(FIS_CHILD_WEST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))     !<-- FIS on child W bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
              CALL MPI_RECV(FIS_CHILD_WEST(N)%TASKS(NTX)%DATA           &  !<-- Recv FIS values from Wbndry child task NTX
                           ,NUM_WORDS                                   &  !<-- # of FIS values
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)        &  !<-- The child task sending
                           ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)        &  !<-- Tag
                           ,COMM_TO_MY_CHILDREN(N)                      &  !<-- MPI communicator
                           ,STATUS                                      &
                           ,IERR)
!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child East Boundary FIS
!-----------------------------
!
          IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                  !<-- This parent task covers some child Eboundary H points
!
            ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1:NUM_TASKS_SEND_H_E(N)))     !<-- FIS data slot for each Ebndry child task on parent task
!
            DO NTX=1,NUM_TASKS_SEND_H_E(N)                                 !<-- Loop through those particular child tasks 
!                                                                          
              J_LO=CHILDTASK_H_SAVE(N)%J_LO_EAST(NTX)                      !<-- Starting J of child point bndry segment on this parent task
              J_HI=CHILDTASK_H_SAVE(N)%J_HI_EAST(NTX)                      !<-- Ending J of child point bndry segment on this parent task
!!!           NUM_WORDS=(J_HI-J_LO+1)*N_BLEND_H                            !<-- # of child points in its bndry segment on parent task
              NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                        !<-- # of child points in its bndry segment on parent task
              ALLOCATE(FIS_CHILD_EAST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))     !<-- FIS on child E bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
              CALL MPI_RECV(FIS_CHILD_EAST(N)%TASKS(NTX)%DATA           &  !<-- Recv FIS values from Ebndry child task NTX
                           ,NUM_WORDS                                   &  !<-- # of FIS values
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)        &  !<-- The child task sending
                           ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)        &  !<-- Tag
                           ,COMM_TO_MY_CHILDREN(N)                      &  !<-- MPI communicator
                           ,STATUS                                      &
                           ,IERR)
!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO child_loop2
!
      ENDIF parent_recvs
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PRELIM_CHILD_INFO
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_UPDATE_CHILD_PSFC(FIS,PD,T,Q                    &
                                         ,PT,PDTOP,SG1,SG2              &
                                         ,IMS,IME,JMS,JME               &
                                         ,NLEV                          &
!
                                         ,NUM_CHILD_TASKS               &
                                         ,LIMITS                        &
!
                                         ,FIS_CHILD_SBND                &
                                         ,FIS_CHILD_NBND                &
                                         ,FIS_CHILD_WBND                &
                                         ,FIS_CHILD_EBND                &
!
                                         ,LOCAL_TASK_RANKS_S            &
                                         ,LOCAL_TASK_RANKS_N            &
                                         ,LOCAL_TASK_RANKS_W            &
                                         ,LOCAL_TASK_RANKS_E            &
!
                                         ,NUM_TASKS_SEND_SBND           &
                                         ,NUM_TASKS_SEND_NBND           &
                                         ,NUM_TASKS_SEND_WBND           &
                                         ,NUM_TASKS_SEND_EBND           &
!
                                         ,I_INDX_PARENT_SBND            &
                                         ,I_INDX_PARENT_NBND            &
                                         ,I_INDX_PARENT_WBND            &
                                         ,I_INDX_PARENT_EBND            &
                                         ,J_INDX_PARENT_SBND            &
                                         ,J_INDX_PARENT_NBND            &
                                         ,J_INDX_PARENT_WBND            &
                                         ,J_INDX_PARENT_EBND            &
!                                                                        
                                         ,I_LO_SOUTH                    &
                                         ,I_HI_SOUTH                    &
                                         ,I_HI_SOUTH_TRANSFER           &
                                         ,I_LO_NORTH                    &
                                         ,I_HI_NORTH                    &
                                         ,I_HI_NORTH_TRANSFER           &
                                         ,J_LO_WEST                     &
                                         ,J_HI_WEST                     &
                                         ,J_HI_WEST_TRANSFER            &
                                         ,J_LO_EAST                     &
                                         ,J_HI_EAST                     &
                                         ,J_HI_EAST_TRANSFER            &
!                                                                        
                                         ,WEIGHT_SBND                   & 
                                         ,WEIGHT_NBND                   & 
                                         ,WEIGHT_WBND                   & 
                                         ,WEIGHT_EBND                   & 
!                                                                                ^
                                         ,N_BLEND                       & !      |
                                         ,IM_CHILD_X                    & !      |
                                         ,JM_CHILD_X                    & !    INPUT
!                                                                           -----------
                                         ,CHILD_H_SBND                  & !    OUTPUT
                                         ,CHILD_H_NBND                  & !      |
                                         ,CHILD_H_WBND                  & !      |
                                         ,CHILD_H_EBND                  & !      v
!
                                         ,PDB_SBND                      & !
                                         ,PDB_NBND                      & !
                                         ,PDB_WBND                      & !
                                         ,PDB_EBND )                      !
!
!-----------------------------------------------------------------------
!***  GIVEN A CHILD'S ACTUAL SURFACE GEOPOTENTIAL, GENERATE A NEW 
!***  VALUE OF PD FOR THE CHILD BOUNDARY POINTS BASED ON THE 
!***  SURROUNDING PARENT POINTS.
!-----------------------------------------------------------------------
!
!-----------
!***  Input
!-----------
!
      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME                             &  !<-- Parent task's memory limits
                           ,IM_CHILD_X,JM_CHILD_X                       &  !<-- Index limits of the nest domain
                           ,N_BLEND                                     &  !<-- # of domain boundary blending rows 
                           ,NLEV                                        &  !<-- # of vertical levels in parent array
                           ,NUM_CHILD_TASKS                             &  !<-- # of fcst tasks on this child
                           ,NUM_TASKS_SEND_SBND                         &  !<-- # of child tasks with Sboundary regions on parent task
                           ,NUM_TASKS_SEND_NBND                         &  !<-- # of child tasks with Nboundary regions on parent task
                           ,NUM_TASKS_SEND_WBND                         &  !<-- # of child tasks with Wboundary regions on parent task
                           ,NUM_TASKS_SEND_EBND                            !<-- # of child tasks with Eboundary regions on parent task
!
!!!   INTEGER,DIMENSION(1:NUM_CHILD_TASKS),INTENT(IN) ::                &
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: LOCAL_TASK_RANKS_S     &  !<-- Ranks of this child's Sbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_N     &  !<-- Ranks of this child's Nbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_W     &  !<-- Ranks of this child's Wbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_E        !<-- Ranks of this child's Ebndry fcst tasks
!
      INTEGER,DIMENSION(1:4,1:NUM_CHILD_TASKS),INTENT(IN) :: LIMITS        !<-- ITS,ITE,JTS,JTE on each task of the child
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: I_LO_SOUTH             &  !<-- Starting I of Sbndry region on child tasks
                                                ,I_HI_SOUTH             &  !<-- Ending I of Sbndry region on child tasks
                                                ,I_HI_SOUTH_TRANSFER    &  !<-- Ending I of Sbndry region for transfer to child
                                                ,I_LO_NORTH             &  !<-- Starting I of Nbndry region on child tasks
                                                ,I_HI_NORTH             &  !<-- Ending I of Nbndry region on child tasks
                                                ,I_HI_NORTH_TRANSFER    &  !<-- Ending I of Nbndry region for transfer to child
                                                ,J_LO_WEST              &  !<-- Starting J of Wbndry region on child tasks
                                                ,J_HI_WEST              &  !<-- Ending J of Wbndry region on child tasks
                                                ,J_HI_WEST_TRANSFER     &  !<-- Ending J of Wbndry region for transfer to child
                                                ,J_LO_EAST              &  !<-- Starting J of Ebndry region on child tasks
                                                ,J_HI_EAST              &  !<-- Ending J of Ebndry region on child tasks
                                                ,J_HI_EAST_TRANSFER        !<-- Ending J of Ebndry region for transfer to child
!
      INTEGER,DIMENSION(:,:,:),POINTER,INTENT(IN) :: I_INDX_PARENT_SBND &  !<-- I indices of parent points west/east of child Sbndry point
                                                    ,I_INDX_PARENT_NBND &  !<-- I indices of parent points west/east of child Nbndry point
                                                    ,I_INDX_PARENT_WBND &  !<-- I indices of parent points west/east of child Wbndry point
                                                    ,I_INDX_PARENT_EBND &  !<-- I indices of parent points west/east of child Ebndry point
                                                    ,J_INDX_PARENT_SBND &  !<-- J indices of parent points south/north of child Sbndry point
                                                    ,J_INDX_PARENT_NBND &  !<-- J indices of parent points south/north of child Nbndry point
                                                    ,J_INDX_PARENT_WBND &  !<-- J indices of parent points south/north of child Wbndry point
                                                    ,J_INDX_PARENT_EBND    !<-- J indices of parent points south/north of child Ebndry point
!
      REAL,INTENT(IN) :: PT                                             &  !<-- Top pressure of model domain (Pa)
                        ,PDTOP                                             !<-- Pressure at top of sigma domain (Pa)
!
      REAL,DIMENSION(:),POINTER,INTENT(IN) :: SG1                       &  !<-- Interface sigmas, pressure domain
                                             ,SG2                          !<-- Interface sigmas, sigma domain
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS                 &  !<-- Parent FIS
                                                   ,PD                     !<-- Parent PD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) :: T            &  !<-- Parent sensible temperature (K)
                                                          ,Q               !<-- Parent specific humidity (kg/kg)
!
      REAL,DIMENSION(:,:,:),POINTER,INTENT(IN) :: WEIGHT_SBND           &  !<-- Bilinear interp weights of the 4 parent points around 
                                                 ,WEIGHT_NBND           &  !    each point on child's boundary sides (S,N,W,E).
                                                 ,WEIGHT_WBND           &  !    
                                                 ,WEIGHT_EBND              !    
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: FIS_CHILD_SBND &  !<-- Sfc geopot on Sbndry points of each child task
                                                        ,FIS_CHILD_NBND &  !<-- Sfc geopot on Nbndry points of each child task
                                                        ,FIS_CHILD_WBND &  !<-- Sfc geopot on Wbndry points of each child task
                                                        ,FIS_CHILD_EBND    !<-- Sfc geopot on Ebndry points of each child task
!
!------------
!***  Output
!------------
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(OUT) :: CHILD_H_SBND  &  !<-- All H point data for child Sbndry to be sent by parent
                                                         ,CHILD_H_NBND  &  !<-- All H point data for child Nbndry to be sent by parent
                                                         ,CHILD_H_WBND  &  !<-- All H point data for child Wbndry to be sent by parent
                                                         ,CHILD_H_EBND  &  !<-- All H point data for child Ebndry to be sent by parent
!
                                                         ,PDB_SBND      &  !<-- Child boundary PD (Pa) on child domain Sbndry
                                                         ,PDB_NBND      &  !<-- Child boundary PD (Pa) on child domain Nbndry
                                                         ,PDB_WBND      &  !<-- Child boundary PD (Pa) on child domain Wbndry
                                                         ,PDB_EBND         !<-- Child boundary PD (Pa) on child domain Ebndry
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,L                                                  &
                ,I_EAST,I_WEST,J_SOUTH,J_NORTH                          &
                ,I_START,I_END,J_START,J_END                            &
                ,I_START_TRANSFER,I_END_TRANSFER                        &
                ,J_START_TRANSFER,J_END_TRANSFER                        &
                ,KOUNT_PTS                                              &
                ,KOUNT_TRANSFER                                         &
                ,N_SIDE,NUM_TASKS_SEND,NTX                              &
                ,RC
!
      INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_PARENT_BND             &
                                         ,J_INDX_PARENT_BND
!
      REAL :: COEFF_1,COEFF_2,D_LNP_DFI,FIS_CHILD                       &
             ,LOG_P1_PARENT,PDTOP_PT,PHI_DIFF,PSFC_CHILD                &
             ,Q_INTERP,T_INTERP
!
      REAL :: WGHT_NE,WGHT_NW,WGHT_SE,WGHT_SW
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: PX
!
      REAL,DIMENSION(:,:),ALLOCATABLE :: LOG_PBOT                       &
                                        ,LOG_PTOP
!
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: PHI_INTERP                   &
                                          ,PINT_INTERP
!
      REAL,DIMENSION(:,:,:),POINTER :: WEIGHT_BND
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER :: CHILD_BOUND_H             &
                                             ,FIS_CHILD_BND             &
                                             ,PDB
!
      integer,dimension(8) :: values
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  LOOP THROUGH THE FOUR SIDES OF THE NEST DOMAIN BOUNDARY (S,N,W,E).
!***  WE USE SOME DUMMY VARIABLES/POINTERS GENERICALLY FOR ALL FOUR
!***  OF THE SIDES.
!-----------------------------------------------------------------------
!
      loop_sides: DO N_SIDE=1,4                                              !<-- Loop through the 4 lateral boundaries (S,N,W,E)
!
!-----------------------------------------------------------------------
!
        IF(N_SIDE==1)THEN
          IF(NUM_TASKS_SEND_SBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Sbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_SBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_SBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_SBND
          WEIGHT_BND=>WEIGHT_SBND
          FIS_CHILD_BND=>FIS_CHILD_SBND
          PDB=>PDB_SBND
          CHILD_BOUND_H=>CHILD_H_SBND
!
        ELSEIF(N_SIDE==2)THEN
          IF(NUM_TASKS_SEND_NBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Nbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_NBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_NBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_NBND
          WEIGHT_BND=>WEIGHT_NBND
          FIS_CHILD_BND=>FIS_CHILD_NBND
          PDB=>PDB_NBND
          CHILD_BOUND_H=>CHILD_H_NBND
!
        ELSEIF(N_SIDE==3)THEN
          IF(NUM_TASKS_SEND_WBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Wbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_WBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_WBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_WBND
          WEIGHT_BND=>WEIGHT_WBND
          FIS_CHILD_BND=>FIS_CHILD_WBND
          PDB=>PDB_WBND
          CHILD_BOUND_H=>CHILD_H_WBND
!
        ELSEIF(N_SIDE==4)THEN
          IF(NUM_TASKS_SEND_EBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Ebndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_EBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_EBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_EBND
          WEIGHT_BND=>WEIGHT_EBND
          FIS_CHILD_BND=>FIS_CHILD_EBND
          PDB=>PDB_EBND
          CHILD_BOUND_H=>CHILD_H_EBND
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        child_task_loop: DO NTX=1,NUM_TASKS_SEND                            !<-- Fill bndry data for each child task on the child bndry
                                                                            !    segment seen by this parent task.
!-----------------------------------------------------------------------
!
!----------------------------------------------
!***  South boundary limits on this child task
!----------------------------------------------
!
          IF(N_SIDE==1)THEN
            I_START=I_LO_SOUTH(NTX)
            I_END  =I_HI_SOUTH(NTX)
            J_START=1
!           J_END  =N_BLEND
            J_END  =N_BLEND+1                                               !<-- Extend by one row to allow 4-pt averaging of PD to V pts
!
            I_START_TRANSFER=I_START
            I_END_TRANSFER  =I_HI_SOUTH_TRANSFER(NTX)                       !<-- Extra PD points for 4-pt average are not transferred
!
            J_START_TRANSFER=J_START
            J_END_TRANSFER  =J_END-1                                        !<-- The extra row of PD is not transferred to the nests
!
!----------------------------------------------
!***  North boundary limits on this child task
!----------------------------------------------
!
          ELSEIF(N_SIDE==2)THEN
            I_START=I_LO_NORTH(NTX)
            I_END  =I_HI_NORTH(NTX)
!           J_START=JM_CHILD_X-N_BLEND+1
            J_START=JM_CHILD_X-N_BLEND                                      !<-- Extend by one row to allow 4-pt averaging of PD to V pts
            J_END  =JM_CHILD_X
!
            I_START_TRANSFER=I_START
            I_END_TRANSFER  =I_HI_NORTH_TRANSFER(NTX)                       !<-- Extra PD points for 4-pt average are not transferred
!
            J_START_TRANSFER=J_START+1                                      !<-- The extra row of PD is not transferred to the nests
            J_END_TRANSFER  =J_END
!
!---------------------------------------------
!***  West boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==3)THEN
            I_START=1
!           I_END  =N_BLEND
            I_END  =N_BLEND+1                                               !<-- Extend by one row to allow 4-pt averaging of PD to V pts
            J_START=J_LO_WEST(NTX)
            J_END  =J_HI_WEST(NTX)
!
            I_START_TRANSFER=I_START
            I_END_TRANSFER  =I_END-1                                        !<-- The extra row of PD is not transferred to the nests
!
            J_START_TRANSFER=J_START
            J_END_TRANSFER  =J_HI_WEST_TRANSFER(NTX)                        !<-- Extra PD points for 4-pt average are not transferred
!
!---------------------------------------------
!***  East boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==4)THEN
!           I_START=IM_CHILD_X-N_BLEND+1
            I_START=IM_CHILD_X-N_BLEND                                      !<-- Extend by one row to allow 4-pt averaging of PD to V pts
            I_END  =IM_CHILD_X
            J_START=J_LO_EAST(NTX)
            J_END  =J_HI_EAST(NTX)
!
            I_START_TRANSFER=I_START+1                                      !<-- The extra row of PD is not transferred to the nests
            I_END_TRANSFER  =I_END  
!
            J_START_TRANSFER=J_START
            J_END_TRANSFER  =J_HI_EAST_TRANSFER(NTX)                        !<-- Extra PD points for 4-pt average are not transferred
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE WORKING NEST ARRAYS VALID FOR THE CURRENT CHILD TASK
!***  ON THE CHILD'S GRID.
!-----------------------------------------------------------------------
!
          ALLOCATE(PINT_INTERP(I_START:I_END,J_START:J_END,1:LM+1))
          ALLOCATE( PHI_INTERP(I_START:I_END,J_START:J_END,1:LM+1))
          ALLOCATE(   LOG_PTOP(I_START:I_END,J_START:J_END))
          ALLOCATE(   LOG_PBOT(I_START:I_END,J_START:J_END))
!
!-----------------------------------------------------------------------
!***  COMPUTE PARENT HEIGHTS OF LAYER INTERFACES AT THE FOUR POINTS
!***  SURROUNDING EACH CHILD BOUNDARY POINT.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST THE BOTTOM LAYER (L=NLEV)
!-----------------------------------------------------------------------
!
          DO J=J_START,J_END                                               !<-- J limits of child task bndry region on parent task
          DO I=I_START,I_END                                               !<-- I limits of child task bndry region on parent task
!
            I_WEST =I_INDX_PARENT_BND(I,J,1)                               !<-- Parent I index on or west of child's boundary point
            I_EAST =I_INDX_PARENT_BND(I,J,2)                               !<-- Parent I index east of child's boundary point
            J_SOUTH=J_INDX_PARENT_BND(I,J,1)                               !<-- Parent J index on or south of child's boundary point
            J_NORTH=J_INDX_PARENT_BND(I,J,2)                               !<-- Parent J index north of child's boundary point
!
            WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                                !<-- Bilinear weight for parent's point SW of nest's point
            WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                                !<-- Bilinear weight for parent's point SE of nest's point
            WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                                !<-- Bilinear weight for parent's point NW of nest's point
            WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                                !<-- Bilinear weight for parent's point NE of nest's point
!
            PX(I_WEST,J_SOUTH)=PD(I_WEST,J_SOUTH)+PT                 !<-- Sfc pressure on parent point SW of nest point
            PX(I_EAST,J_SOUTH)=PD(I_EAST,J_SOUTH)+PT                 !<-- Sfc pressure on parent point SE of nest point
            PX(I_WEST,J_NORTH)=PD(I_WEST,J_NORTH)+PT                 !<-- Sfc pressure on parent point NW of nest point
            PX(I_EAST,J_NORTH)=PD(I_EAST,J_NORTH)+PT                 !<-- Sfc pressure on parent point NE of nest point
!
            PINT_INTERP(I,J,LM+1)=WGHT_SW*PX(I_WEST,J_SOUTH)            &  !<-- Parent's surface pressure interp'd to this child's 
                                 +WGHT_SE*PX(I_EAST,J_SOUTH)            &  !    gridpoint (I,J) along child's boundary 
                                 +WGHT_NW*PX(I_WEST,J_NORTH)            &  !    on child task NTX.
                                 +WGHT_NE*PX(I_EAST,J_NORTH)
!
            LOG_PBOT(I,J)=LOG(PINT_INTERP(I,J,LM+1))                       !<-- Log of parent's horizontally interpolated sfc pressure 
                                                                           !    at child boundary point (I,J)
!
            PHI_INTERP(I,J,LM+1)=WGHT_SW*FIS(I_WEST,J_SOUTH)            &  !<-- Parent sfc geoptential interp'd to nest bndry point (I,J)
                                +WGHT_SE*FIS(I_EAST,J_SOUTH)            &
                                +WGHT_NW*FIS(I_WEST,J_NORTH)            &
                                +WGHT_NE*FIS(I_EAST,J_NORTH)
!
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  NOW THAT WE HAVE THE PARENT'S SFC PRESSURE AND SFC GEOPOTENTIAL
!***  AT THE CHILD BOUNDARY POINTS (I,J), COMPUTE THE INTERFACE HEIGHTS
!***  BASED ON THE HORIZONTALLY INTERPOLATED INTERFACE PRESSURE AND 
!***  THE T AND Q.
!-----------------------------------------------------------------------
!
          DO L=NLEV,1,-1                                                   !<-- Work upward to obtain interface geopotentials
!
            PDTOP_PT=PDTOP*SG1(L)+PT
!
            DO J=J_START,J_END                                             !<-- J limits of child task bndry region on parent task
            DO I=I_START,I_END                                             !<-- I limits of child task bndry region on parent task
!
              I_WEST =I_INDX_PARENT_BND(I,J,1)                             !<-- Parent I index on or west of child's boundary point
              I_EAST =I_INDX_PARENT_BND(I,J,2)                             !<-- Parent I index east of child's boundary point
              J_SOUTH=J_INDX_PARENT_BND(I,J,1)                             !<-- Parent J index on or south of child's boundary point
              J_NORTH=J_INDX_PARENT_BND(I,J,2)                             !<-- Parent J index north of child's boundary point
!
              PX(I_WEST,J_SOUTH)=SG2(L)*PD(I_WEST,J_SOUTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point SW of nest point
              PX(I_EAST,J_SOUTH)=SG2(L)*PD(I_EAST,J_SOUTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point SE of nest point
              PX(I_WEST,J_NORTH)=SG2(L)*PD(I_WEST,J_NORTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point NW of nest point
              PX(I_EAST,J_NORTH)=SG2(L)*PD(I_EAST,J_NORTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point NE of nest point
!
              WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                              !<-- Bilinear weight for parent's point SW of child's point
              WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                              !<-- Bilinear weight for parent's point SE of child's point
              WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                              !<-- Bilinear weight for parent's point NW of child's point
              WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                              !<-- Bilinear weight for parent's point NE of child's point
!
              PINT_INTERP(I,J,L)=WGHT_SW*PX(I_WEST,J_SOUTH)             &  !<-- Top interface pressure interp'd to child gridpoint
                                +WGHT_SE*PX(I_EAST,J_SOUTH)             &  !    along child's boundary for child task NTX
                                +WGHT_NW*PX(I_WEST,J_NORTH)             &
                                +WGHT_NE*PX(I_EAST,J_NORTH)
!
              LOG_PTOP(I,J)=LOG(PINT_INTERP(I,J,L))                        !<-- Log of parent (top) interface pressure at child bndry point
!
              T_INTERP=WGHT_SW*T(I_WEST,J_SOUTH,L)                      &  !<-- T interp'd to child gridpoint along child's
                      +WGHT_SE*T(I_EAST,J_SOUTH,L)                      &  !    boundary for child task NTX
                      +WGHT_NW*T(I_WEST,J_NORTH,L)                      &
                      +WGHT_NE*T(I_EAST,J_NORTH,L)
!
              Q_INTERP=WGHT_SW*Q(I_WEST,J_SOUTH,L)                      &  !<-- Q interp'd to child gridpoint along child's
                      +WGHT_SE*Q(I_EAST,J_SOUTH,L)                      &  !    boundary for child task NTX
                      +WGHT_NW*Q(I_WEST,J_NORTH,L)                      &
                      +WGHT_NE*Q(I_EAST,J_NORTH,L)
!
              PHI_INTERP(I,J,L)=PHI_INTERP(I,J,L+1)                     &  !<-- Top interface geopotl of parent at child gridpoint (I,J)
                               +R_D*T_INTERP*(1.+P608*Q_INTERP)         &
                               *(LOG_PBOT(I,J)-LOG_PTOP(I,J))
!
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  !<--- Move Log(Ptop) to bottom of next model layer up
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  ! 
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  ! 
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  !<--
!
            ENDDO
            ENDDO
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  USE THE CHILD'S ACTUAL SFC GEOPOTENTIAL TO DERIVE THE VALUE OF
!***  PD AT THE CHILD BOUNDARY POINTS BASED ON THE PARENT'S HEIGHTS 
!***  AND PRESSURES AT ITS (THE PARENT'S) LAYER INTERFACES OVER THE
!***  CHILD'S BOUNDARY POINTS.
!
!***  IF THE CHILD'S TERRAIN IS LOWER THAN THE VALUE OF THE PARENT'S
!***  TERRAIN INTERPOLATED TO THE CHILD POINT THEN EXTRAPOLATE THE
!***  PARENT'S INTERPOLATED SFC PRESSURE DOWN TO THE CHILD'S TERRAIN
!***  QUADRATICALLY.
!-----------------------------------------------------------------------
!
          KOUNT_PTS=0
          KOUNT_TRANSFER=0
!
!-----------------------------------------------------------------------
          core_loop: DO J=J_START,J_END                                    !<-- J limits of child task bndry region on parent task
          DO I=I_START,I_END                                               !<-- I limits of child task bndry region on parent task
!-----------------------------------------------------------------------
!
            KOUNT_PTS=KOUNT_PTS+1
            FIS_CHILD=FIS_CHILD_BND(NTX)%DATA(KOUNT_PTS)
!
            IF(FIS_CHILD<PHI_INTERP(I,J,LM+1))THEN                         !<-- Child's terrain lies below parent's (interp'd to child)
!
              COEFF_1=(PINT_INTERP(I,J,LM+1)-PINT_INTERP(I,J,LM))       &  !<-- Coefficient for linear term
                     /(PHI_INTERP(I,J,LM+1)-PHI_INTERP(I,J,LM))
!
              COEFF_2=(COEFF_1                                          &  !<-- Coefficient for quadratic term
                      -(PINT_INTERP(I,J,LM)-PINT_INTERP(I,J,LM-1))      &  !
                      /(PHI_INTERP(I,J,LM)-PHI_INTERP(I,J,LM-1)))       &  !
                      /(PHI_INTERP(I,J,LM+1)-PHI_INTERP(I,J,LM-1))         !
!
              PHI_DIFF=FIS_CHILD-PHI_INTERP(I,J,LM+1)                      !<-- Diff between nest sfc geopotential and parent's interp'd
!
              PSFC_CHILD=COEFF_2*PHI_DIFF*PHI_DIFF                      &  !<-- Parent pressure at child's boundary terrain surface
                        +COEFF_1*PHI_DIFF                               &
                        +PINT_INTERP(I,J,LM+1)
!
            ELSE                                                           !<-- Child's terrain is at or above parent's
              DO L=LM,1,-1
                IF(FIS_CHILD<PHI_INTERP(I,J,L))THEN
                  LOG_P1_PARENT=LOG(PINT_INTERP(I,J,L+1))                  !<-- Log pressure on bottom of parent model layer
                  D_LNP_DFI=(LOG_P1_PARENT-LOG(PINT_INTERP(I,J,L)))     &  !<-- d[ln(p)]/d[fi] in parent layer L
                           /(PHI_INTERP(I,J,L+1)-PHI_INTERP(I,J,L))
                  PSFC_CHILD=EXP(LOG_P1_PARENT                          &  !<-- Parent pressure at child's boundary terrain surface
                                +D_LNP_DFI*(FIS_CHILD-PHI_INTERP(I,J,L+1)))
                  EXIT
                ENDIF
              ENDDO
            ENDIF
!
            PDB(NTX)%DATA(KOUNT_PTS)=PSFC_CHILD-PT                         !<-- Parent's approximation of child's PD on child's boundary
!
!-----------------------------------------------------------------------
!***  SAVE ONLY THE PD's THAT ARE ON THE NESTS' BOUNDARY REGION'S 
!***  MASS POINTS INTO THE DATASTRING THAT WILL BE TRANSFERRED FROM
!***  THE PARENT TO THE NESTS.  THE EXTRA POINTS IN THE PDB POINTER
!***  ARE FOR 4-PT AVERAGING ONTO V POINTS.
!***  RECALL THAT THE CHILD_BOUND_* POINTER IS THE 1-D STRING OF ALL
!***  DATA SENT FROM PARENT TASKS TO CHILD BOUNDARY TASKS AND PDB IS
!***  THE FIRST DATA IN THAT STRING.
!-----------------------------------------------------------------------
!
            IF(I>=I_START_TRANSFER.AND.I<=I_END_TRANSFER                &
                  .AND.                                                 &
               J>=J_START_TRANSFER.AND.J<=J_END_TRANSFER)THEN
!
              KOUNT_TRANSFER=KOUNT_TRANSFER+1
              CHILD_BOUND_H(NTX)%DATA(KOUNT_TRANSFER)=                  &
                        PDB(NTX)%DATA(KOUNT_PTS)
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO
          ENDDO core_loop
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(PINT_INTERP)
          DEALLOCATE(PHI_INTERP)
          DEALLOCATE(LOG_PTOP)  
          DEALLOCATE(LOG_PBOT)  
!
!-----------------------------------------------------------------------
!
        ENDDO child_task_loop
!
!-----------------------------------------------------------------------
!
      ENDDO loop_sides
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_UPDATE_CHILD_PSFC
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_UPDATE_CHILD_BNDRY(VBL_PARENT                   &
                                          ,PD,PT,PDTOP                  &
                                          ,PSGML1,SGML2,SG1,SG2         &
!
                                          ,PD_SBND                      &
                                          ,PD_NBND                      &
                                          ,PD_WBND                      &
                                          ,PD_EBND                      &
!
                                          ,IMS,IME,JMS,JME              &
                                          ,NLEV                         &
                                          ,N_REMOVE                     &
!
                                          ,NUM_TASKS_SEND_SBND          &
                                          ,NUM_TASKS_SEND_NBND          &
                                          ,NUM_TASKS_SEND_WBND          &
                                          ,NUM_TASKS_SEND_EBND          &
!
                                          ,I_INDX_PARENT_SBND           &
                                          ,I_INDX_PARENT_NBND           &
                                          ,I_INDX_PARENT_WBND           &
                                          ,I_INDX_PARENT_EBND           &
                                          ,J_INDX_PARENT_SBND           &
                                          ,J_INDX_PARENT_NBND           &
                                          ,J_INDX_PARENT_WBND           &
                                          ,J_INDX_PARENT_EBND           &
!                                                                       
                                          ,I_LO_SOUTH                   &
                                          ,I_HI_SOUTH                   &
                                          ,I_HI_SOUTH_TRANSFER          &
                                          ,I_LO_NORTH                   &
                                          ,I_HI_NORTH                   &
                                          ,I_HI_NORTH_TRANSFER          &
                                          ,J_LO_WEST                    &
                                          ,J_HI_WEST                    &
                                          ,J_HI_WEST_TRANSFER           &
                                          ,J_LO_EAST                    &
                                          ,J_HI_EAST                    &
                                          ,J_HI_EAST_TRANSFER           &
!
                                          ,WEIGHT_SBND                  &
                                          ,WEIGHT_NBND                  &
                                          ,WEIGHT_WBND                  &
                                          ,WEIGHT_EBND                  &
!                                                                               ^
                                          ,N_BLEND                      & !     |
                                          ,IM_CHILD_X                   & !     |
                                          ,JM_CHILD_X                   & !   INPUT
!                                                                           ----------
                                          ,VBL_CHILD_SBND               & !   OUTPUT
                                          ,VBL_CHILD_NBND               & !     |
                                          ,VBL_CHILD_WBND               & !     |
                                          ,VBL_CHILD_EBND )               !     v
!
!-----------------------------------------------------------------------
!***  PARENT TASKS INTERPOLATE THEIR VALUES OF VARIABLES
!***  TO CHILD GRID POINTS. 
!-----------------------------------------------------------------------
!
!-----------
!***  Input
!-----------
!
      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME                             &  !<-- Parent task's memory limits
                           ,IM_CHILD_X,JM_CHILD_X                       &  !<-- Index limits of the nest domain
                           ,N_BLEND                                     &  !<-- # of domain boundary blending rows
                           ,NLEV                                        &  !<-- # of vertical levels in parent array
                           ,N_REMOVE                                    &  !<-- # of rows to ignore on north/east sides (H=>0;V=>1)
!
                           ,NUM_TASKS_SEND_SBND                         &  !<-- # of child tasks with Sbndry points on parent task
                           ,NUM_TASKS_SEND_NBND                         &  !<-- # of child tasks with Nbndry points on parent task
                           ,NUM_TASKS_SEND_WBND                         &  !<-- # of child tasks with Wbndry points on parent task
                           ,NUM_TASKS_SEND_EBND                            !<-- # of child tasks with Ebndry points on parent task
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: I_LO_SOUTH             &  !<-- Starting I of Sbndry region on child tasks
                                                ,I_HI_SOUTH             &  !<-- Ending I of Sbndry region on child tasks
                                                ,I_HI_SOUTH_TRANSFER    &  !<-- Ending I of Sbndry for transfer to child
                                                ,I_LO_NORTH             &  !<-- Starting I of Nbndry region on child tasks
                                                ,I_HI_NORTH             &  !<-- Ending I of Nbndry region on child tasks
                                                ,I_HI_NORTH_TRANSFER    &  !<-- Ending I of Nbndry for transfer to child
                                                ,J_LO_WEST              &  !<-- Starting J of Wbndry region on child tasks
                                                ,J_HI_WEST              &  !<-- Ending J of Wbndry region on child tasks
                                                ,J_HI_WEST_TRANSFER     &  !<-- Ending J of Wbndry for transfer to child
                                                ,J_LO_EAST              &  !<-- Starting J of Ebndry region on child tasks
                                                ,J_HI_EAST              &  !<-- Ending J of Ebndry region on child tasks
                                                ,J_HI_EAST_TRANSFER        !<-- Ending J of Ebndry for transfer to child
!
      INTEGER,DIMENSION(:,:,:),POINTER,INTENT(IN) :: I_INDX_PARENT_SBND &  !<-- I indices of parent points west/east of child Sbndry point
                                                    ,I_INDX_PARENT_NBND &  !<-- I indices of parent points west/east of child Nbndry point
                                                    ,I_INDX_PARENT_WBND &  !<-- I indices of parent points west/east of child Wbndry point
                                                    ,I_INDX_PARENT_EBND &  !<-- I indices of parent points west/east of child Ebndry point
                                                    ,J_INDX_PARENT_SBND &  !<-- J indices of parent points south/north of child Sbndry point
                                                    ,J_INDX_PARENT_NBND &  !<-- J indices of parent points south/north of child Nbndry point
                                                    ,J_INDX_PARENT_WBND &  !<-- J indices of parent points south/north of child Wbndry point
                                                    ,J_INDX_PARENT_EBND    !<-- J indices of parent points south/north of child Ebndry point
!
      REAL,INTENT(IN) :: PT                                             &  !<-- Top pressure of model domain (Pa)
                        ,PDTOP                                             !<-- Pressure at top of sigma domain (Pa)
!
      REAL,DIMENSION(:),POINTER,INTENT(IN) :: PSGML1                    &  !<-- Midlayer pressures, pressure domain
                                             ,SGML2                     &  !<-- Midlayer sigmas, sigma domain
                                             ,SG1                       &  !<-- Interface sigmas, pressure domain
                                             ,SG2                          !<-- Interface sigmas, sigma domain
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PD                     !<-- Parent PD

!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) :: VBL_PARENT      !<-- Parent mass point variable
!
      REAL,DIMENSION(:,:,:),POINTER,INTENT(IN) :: WEIGHT_SBND           &  !<-- Bilinear interp weights of the 4 parent points around 
                                                 ,WEIGHT_NBND           &  !    each point on child's boundaries.
                                                 ,WEIGHT_WBND           &  !   
                                                 ,WEIGHT_EBND              !  
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: PD_SBND        &  !<-- Boundary region PD (Pa) (column mass in sigma domain)
                                                        ,PD_NBND        &  !    on the four sides of the child boundary.
                                                        ,PD_WBND        &  !
                                                        ,PD_EBND           !
!
!------------
!***  Output
!------------
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(OUT) :: VBL_CHILD_SBND &  !<-- Mass variable in child bndry region as computed
                                                         ,VBL_CHILD_NBND &  !<-- by parent.
                                                         ,VBL_CHILD_WBND &  !
                                                         ,VBL_CHILD_EBND    !
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,L                                                  &
                ,I_EAST,I_WEST,J_SOUTH,J_NORTH                          &
                ,I_START,I_START_EXPAND                                 &
                ,I_END,I_END_EXPAND                                     &
                ,J_START,J_START_EXPAND                                 &
                ,J_END,J_END_EXPAND                                     &
                ,KNT_PTS,KNT_PTS_X                                      &
                ,LOC_1,LOC_2                                            &
                ,N_ADD,N_EXP,N_SIDE,N_STRIDE,NTX                        &
                ,NUM_LEVS_SPLINE,NUM_TASKS_SEND                         &
                ,RC
!
      INTEGER,DIMENSION(:,:,:),POINTER :: I_INDX_PARENT_BND             &
                                         ,J_INDX_PARENT_BND
!
      REAL :: COEFF_1,DELP_EXTRAP                                       &
             ,PDTOP_PT,R_DELP
!
      REAL :: WGHT_NE,WGHT_NW,WGHT_SE,WGHT_SW
!
      REAL,DIMENSION(1:LM) :: PMID_CHILD
!
      REAL,DIMENSION(1:LM+1) :: P_INPUT                                 &
                               ,SEC_DERIV                               &
                               ,VBL_INPUT
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: PX
!
      REAL,DIMENSION(:,:),ALLOCATABLE :: PINT_INTERP_HI                 &
                                        ,PINT_INTERP_LO                 &
                                        ,PMID_INTERP                    &
                                        , VBL_INTERP
!
      REAL,DIMENSION(:),POINTER :: VBL_COL_CHILD
!
      REAL,DIMENSION(:,:,:),POINTER :: WEIGHT_BND
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER :: PDB                       &
                                             ,VBL_CHILD_BND
!
      integer :: nnn,lll
      integer,dimension(8) :: values
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N_EXP=1-N_REMOVE                                                     !<-- Handles expansion of PDB range (H->1; V->0)
!
!-----------------------------------------------------------------------
!***  LOOP THROUGH THE FOUR SIDES OF THE NEST DOMAIN BOUNDARY (S,N,W,E).
!***  WE USE SOME DUMMY VARIABLES/POINTERS GENERICALLY FOR ALL FOUR
!***  OF THE SIDES.
!-----------------------------------------------------------------------
!
      loop_sides: DO N_SIDE=1,4                                            
!
!-----------------------------------------------------------------------
!
        IF(N_SIDE==1)THEN
          IF(NUM_TASKS_SEND_SBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Sbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_SBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_SBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_SBND
          WEIGHT_BND=>WEIGHT_SBND
          PDB=>PD_SBND
          VBL_CHILD_BND=>VBL_CHILD_SBND
!
        ELSEIF(N_SIDE==2)THEN
          IF(NUM_TASKS_SEND_NBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Nbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_NBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_NBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_NBND
          WEIGHT_BND=>WEIGHT_NBND
          PDB=>PD_NBND
          VBL_CHILD_BND=>VBL_CHILD_NBND
!
        ELSEIF(N_SIDE==3)THEN
          IF(NUM_TASKS_SEND_WBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Wbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_WBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_WBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_WBND
          WEIGHT_BND=>WEIGHT_WBND
          PDB=>PD_WBND
          VBL_CHILD_BND=>VBL_CHILD_WBND
!
        ELSEIF(N_SIDE==4)THEN
          IF(NUM_TASKS_SEND_EBND==0)CYCLE                                    !<-- Move on if parent task sees no part of nest's Ebndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_EBND
          I_INDX_PARENT_BND=>I_INDX_PARENT_EBND
          J_INDX_PARENT_BND=>J_INDX_PARENT_EBND
          WEIGHT_BND=>WEIGHT_EBND
          PDB=>PD_EBND
          VBL_CHILD_BND=>VBL_CHILD_EBND
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        child_task_loop: DO NTX=1,NUM_TASKS_SEND                           !<-- Fill bndry data for each child task on the child bndry
                                                                           !    segment seen by this parent task.
!-----------------------------------------------------------------------
!
!----------------------------------------------
!***  South boundary limits on this child task
!----------------------------------------------
!
          IF(N_SIDE==1)THEN
            I_START =I_LO_SOUTH(NTX)
            I_END   =I_HI_SOUTH_TRANSFER(NTX)
            J_START =1
            J_END   =N_BLEND
!
            I_START_EXPAND=I_START                                         !<-- Expanded limits for extra row of PD for 4-pt averaging
            I_END_EXPAND=I_HI_SOUTH(NTX)                                   !
            J_START_EXPAND=J_START                                         !
            J_END_EXPAND=J_END+N_EXP                                       !<--
!
!----------------------------------------------
!***  North boundary limits on this child task
!----------------------------------------------
!
          ELSEIF(N_SIDE==2)THEN
            I_START =I_LO_NORTH(NTX)
            I_END   =I_HI_NORTH_TRANSFER(NTX)
            J_START =JM_CHILD_X-N_BLEND+1-N_REMOVE
            J_END   =JM_CHILD_X-N_REMOVE
!
            I_START_EXPAND=I_START                                         !<-- Expanded limits for extra row of PD for 4-pt averaging
            I_END_EXPAND=I_HI_NORTH(NTX)                                   !
            J_START_EXPAND=J_START-N_EXP                                   !
            J_END_EXPAND=J_END                                             !<--
!
!---------------------------------------------
!***  West boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==3)THEN
            I_START =1
            I_END   =N_BLEND
            J_START =J_LO_WEST(NTX)
            J_END   =J_HI_WEST_TRANSFER(NTX)
!
            I_START_EXPAND=I_START                                         !<-- Expanded limits for extra row of PD for 4-pt averaging
            I_END_EXPAND=I_END+N_EXP                                       !
            J_START_EXPAND=J_START                                         !
            J_END_EXPAND=J_HI_WEST(NTX)                                    !<--
!
!---------------------------------------------
!***  East boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==4)THEN
            I_START =IM_CHILD_X-N_BLEND+1-N_REMOVE
            I_END   =IM_CHILD_X-N_REMOVE
            J_START =J_LO_EAST(NTX)
            J_END   =J_HI_EAST_TRANSFER(NTX)
!
            I_START_EXPAND=I_START-N_EXP                                   !<-- Expanded limits for extra row of PD for 4-pt averaging
            I_END_EXPAND=I_END                                             !
            J_START_EXPAND=J_START                                         !
            J_END_EXPAND=J_HI_EAST(NTX)                                    !<--
!
          ENDIF
!
!-----------------------------------------------------------------------
!
          N_STRIDE=(I_END-I_START+1)*(J_END-J_START+1)                     !<-- # of pts, this side, this child task's bndry region
!
          ALLOCATE(PMID_INTERP(1:N_STRIDE,1:LM))
          ALLOCATE(VBL_INTERP (1:N_STRIDE,1:LM))
!
          ALLOCATE(PINT_INTERP_HI(I_START:I_END,J_START:J_END))
          ALLOCATE(PINT_INTERP_LO(I_START:I_END,J_START:J_END))
!
!-----------------------------------------------------------------------
!***  WE NEED THE MID-LAYER PRESSURE VALUES IN THE PARENT LAYERS 
!***  OVER THE CHILD BOUNDARY POINT LOCATIONS SINCE THOSE ARE
!***  REQUIRED FOR THE VERTICAL INTERPOLATION OF VARIABLES
!***  TO THE MID-LAYERS IN THE CHILD.
!***  COMPUTE THE INTERFACE PRESSURES OF THE PARENT LAYERS 
!***  THEN TAKE THE MEANS TO GET THE MID-LAYER VALUES.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  START WITH THE BOTTOM LAYER (L=NLEV).
!-----------------------------------------------------------------------
!
          DO J=J_START,J_END                                               !<-- J limits of child task bndry region on parent task
          DO I=I_START,I_END                                               !<-- I limits of child task bndry region on parent task
!
            I_WEST =I_INDX_PARENT_BND(I,J,1)                               !<-- Parent I index on or west of child's boundary point
            I_EAST =I_INDX_PARENT_BND(I,J,2)                               !<-- Parent I index east of child's boundary point
            J_SOUTH=J_INDX_PARENT_BND(I,J,1)                               !<-- Parent J index on or south of child's boundary point
            J_NORTH=J_INDX_PARENT_BND(I,J,2)                               !<-- Parent J index north of child's boundary point
!
            WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                                !<-- Bilinear weight for parent's point SW of nest's point
            WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                                !<-- Bilinear weight for parent's point SE of nest's point
            WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                                !<-- Bilinear weight for parent's point NW of nest's point
            WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                                !<-- Bilinear weight for parent's point NE of nest's point
!
            PX(I_WEST,J_SOUTH)=PD(I_WEST,J_SOUTH)+PT                       !<-- Sfc pressure on parent point SW of nest point
            PX(I_EAST,J_SOUTH)=PD(I_EAST,J_SOUTH)+PT                       !<-- Sfc pressure on parent point SE of nest point
            PX(I_WEST,J_NORTH)=PD(I_WEST,J_NORTH)+PT                       !<-- Sfc pressure on parent point NW of nest point
            PX(I_EAST,J_NORTH)=PD(I_EAST,J_NORTH)+PT                       !<-- Sfc pressure on parent point NE of nest point
!
            PINT_INTERP_LO(I,J)=WGHT_SW*PX(I_WEST,J_SOUTH)              &  !<-- Parent's surface pressure interp'd to this child's 
                               +WGHT_SE*PX(I_EAST,J_SOUTH)              &  !    gridpoint (I,J) along child's boundary for
                               +WGHT_NW*PX(I_WEST,J_NORTH)              &  !    child task NTX.
                               +WGHT_NE*PX(I_EAST,J_NORTH)
!
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  NOW COMPUTE THOSE MID-LAYER PRESSURES IN THE PARENT LAYERS
!***  AS WELL AS THE VALUES OF THE PARENT'S VARIABLE AT THOSE
!***  PRESSURE LEVELS.
!-----------------------------------------------------------------------
!
          DO L=NLEV,1,-1                                                   !<-- Work upward to get geopotentials on child layer interfaces
!
            KNT_PTS=0
            PDTOP_PT=SG1(L)*PDTOP+PT
!
            DO J=J_START,J_END                                             !<-- J limits of child task bndry region on parent task
            DO I=I_START,I_END                                             !<-- I limits of child task bndry region on parent task
!
              I_WEST =I_INDX_PARENT_BND(I,J,1)                             !<-- Parent I index on or west of child's boundary point
              I_EAST =I_INDX_PARENT_BND(I,J,2)                             !<-- Parent I index east of child's boundary point
              J_SOUTH=J_INDX_PARENT_BND(I,J,1)                             !<-- Parent J index on or south of child's boundary point
              J_NORTH=J_INDX_PARENT_BND(I,J,2)                             !<-- Parent J index north of child's boundary point
!
              PX(I_WEST,J_SOUTH)=SG2(L)*PD(I_WEST,J_SOUTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point SW of nest point
              PX(I_EAST,J_SOUTH)=SG2(L)*PD(I_EAST,J_SOUTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point SE of nest point
              PX(I_WEST,J_NORTH)=SG2(L)*PD(I_WEST,J_NORTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point NW of nest point
              PX(I_EAST,J_NORTH)=SG2(L)*PD(I_EAST,J_NORTH)+PDTOP_PT        !<-- Pressure, top of layer L, parent point NE of nest point
!
              WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                              !<-- Bilinear weight for parent's point SW of child's point
              WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                              !<-- Bilinear weight for parent's point SE of child's point
              WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                              !<-- Bilinear weight for parent's point NW of child's point
              WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                              !<-- Bilinear weight for parent's point NE of child's point
!
              PINT_INTERP_HI(I,J)=WGHT_SW*PX(I_WEST,J_SOUTH)            &  !<-- Top interface pressure interp'd to child gridpoint
                                 +WGHT_SE*PX(I_EAST,J_SOUTH)            &  !    in child's boundary region on child task NTX
                                 +WGHT_NW*PX(I_WEST,J_NORTH)            &
                                 +WGHT_NE*PX(I_EAST,J_NORTH)
!
              KNT_PTS=KNT_PTS+1
!
              PMID_INTERP(KNT_PTS,L)=0.5*(PINT_INTERP_HI(I,J)           &  !<-- Parent midlayer pressure interp'd to child gridpoint
                                         +PINT_INTERP_LO(I,J))             !    in child's boundary region of child task NTX
!
              VBL_INTERP(KNT_PTS,L)=WGHT_SW                             &  !<-- Parent variable interp'd to child gridpoint 
                                    *VBL_PARENT(I_WEST,J_SOUTH,L)       &  !    in child's boundary region of child task NTX.
                                   +WGHT_SE                             &  !  
                                    *VBL_PARENT(I_EAST,J_SOUTH,L)       &  !
                                   +WGHT_NW                             &  !  
                                    *VBL_PARENT(I_WEST,J_NORTH,L)       &  !
                                   +WGHT_NE                             &  !
                                    *VBL_PARENT(I_EAST,J_NORTH,L)          !<--
!
              PINT_INTERP_LO(I,J)=PINT_INTERP_HI(I,J)
!
            ENDDO
            ENDDO
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE VALUES OF THE VARIABLE AT MID-LAYERS OF THE CHILD
!***  OVER NEST BOUNDARY REGION POINTS BASED ON THE PARENT'S
!***  INTERPOLATED VALUES.
!-----------------------------------------------------------------------
!
          KNT_PTS  =0
          KNT_PTS_X=0
          LOC_1=0
          N_ADD=(LM-1)*N_STRIDE
!
!-----------------------------------------------------------------------
!***  RECALL THAT THE PDB_H POINTER CONTAINS DATA FOR ONE EXTRA ROW
!***  BEYOND THE CHILD BOUNDARY ON ALL SIDES SINCE WE NEED THE
!***  ABILITY TO DO 4-PT AVERAGES OF PD TO V POINTS.   THE PDB_V 
!***  POINTER DID NOT NEED THE EXTRA ROW AND SO DOES NOT HAVE IT.
!***  HERE WE MUST USE PDB BUT THE VALUES FOR THE PROGNOSTIC VARIABLES
!***  ARE ONLY GENERATED ON THE TRUE BOUNDARY POINTS.  TO ADDRESS 
!***  PDB CORRECTLY WE LOOP OVER THE POTENTIALLY EXPANDED BOUNDARY
!***  REGION BUT CYCLE WHEN WE ARE AT POINTS NOT WITHIN THE TRUE
!***  BOUNDARY REGION OF THE CHILD.
!-----------------------------------------------------------------------
!
          main_loop: DO J=J_START_EXPAND,J_END_EXPAND                      !<-- J limits of expanded child task bndry region on parent task
          DO I=I_START_EXPAND,I_END_EXPAND                                 !<-- Child's expanded south boundary I values on this parent task
!
!-----------------------------------------------------------------------
!
            KNT_PTS_X=KNT_PTS_X+1                                          !<-- Location in 1-D datastring of child PDB
!
!-----------------------------------------------------------------------
!
            IF(I<I_START.OR.I>I_END                                     &
               .OR.                                                     &
               J<J_START.OR.J>J_END)CYCLE
!
            KNT_PTS=KNT_PTS+1                                              !<-- Location in 1-D datastring of child bndry variable
!
!-----------------------------------------------------------------------
!***  MIDLAYER PRESSURES FOR THE NEST POINTS.
!-----------------------------------------------------------------------
!
            DO L=1,LM                                 
              PMID_CHILD(L)=PSGML1(L)+SGML2(L)*PDB(NTX)%DATA(KNT_PTS_X)    ! <-- Nest midlayer pressure
            ENDDO
! 
!-----------------------------------------------------------------------
!***  USE SPLINE INTERPOLATION TO MOVE VARIABLES FROM THEIR 
!***  VERTICAL LOCATIONS IN THE COLUMN AFTER HORIZONTAL INTERPOLATION
!***  FROM THE SURROUNDING PARENT POINTS TO CHILD BOUNDARY POINT LEVELS.
!***  THE TARGET LOCATIONS ARE THE NEW MIDLAYER PRESSURES IN THE
!***  NEST BOUNDARY POINT COLUMNS BASED ON THE NEW SURFACE PRESSURE 
!***  FOR THE NEST'S TERRAIN.
!
!***  IF THE TARGET LOCATION LIES BELOW THE MIDDLE OF THE LOWEST PARENT
!***  LAYER IN THE NEWLY CREATED CHILD COLUMN THEN EXTRAPOLATE LINEARLY
!***  IN PRESSURE TO OBTAIN A VALUE AT THE LOWEST CHILD MID-LAYER AND 
!***  FILL IN THE REMAINING 'UNDERGROUND' LEVELS USING THE CALL TO 
!***  'SPLINE' JUST AS WITH ALL THE OTHER HIGHER LEVELS.
!-----------------------------------------------------------------------
!
            DO L=1,LM                                                      !<-- Extract mid-layer values of parent in nest column
              P_INPUT  (L)=PMID_INTERP(KNT_PTS,L)
              VBL_INPUT(L)= VBL_INTERP(KNT_PTS,L)
            ENDDO
!
            LOC_1=LOC_1+1
            LOC_2=LOC_1+N_ADD
            VBL_COL_CHILD=>VBL_CHILD_BND(NTX)%DATA(LOC_1:LOC_2:N_STRIDE)   !<-- Point working column pointer into 1-D horizontal
!                                                                               output pointer for this variable.
            NUM_LEVS_SPLINE=LM
!
            IF(PMID_CHILD(LM)>P_INPUT(LM))THEN                             !<-- Nest's lowest mid-layer is lower than parent's
              NUM_LEVS_SPLINE=LM+1                                         !    so add another input level at nest's lowest
              P_INPUT(LM+1)=PMID_CHILD(LM)                                 !    mid-layer.
              R_DELP=1./(P_INPUT(LM)-P_INPUT(LM-1))
              DELP_EXTRAP=PMID_CHILD(LM)-P_INPUT(LM)
!
              COEFF_1=(VBL_INPUT(LM)-VBL_INPUT(LM-1))*R_DELP
              VBL_INPUT(LM+1)=VBL_INPUT(LM)+COEFF_1*DELP_EXTRAP            !<-- Create extrapolated value at nest's lowest mid-layer
!                                                                               in input pointer.
            ENDIF
!
            DO L=1,LM+1
              SEC_DERIV(L)=0.                                              !<-- Initialize 2nd derivatives of the spline to zero.
            ENDDO
!
            CALL SPLINE(NUM_LEVS_SPLINE                                 &  !<-- # of input levels
                       ,P_INPUT                                         &  !<-- Input mid-layer pressures
                       ,VBL_INPUT                                       &  !<-- Input mid-layer mass variable value
                       ,SEC_DERIV                                       &  !<-- Specified 2nd derivatives (=0) at parent points
                       ,LM                                              &  !<-- # of child mid-layers to interpolate to
                       ,PMID_CHILD                                      &  !<-- Child mid-layer pressures to interpolate to
                       ,VBL_COL_CHILD)                                     !<-- Child mid-layer variable value returned
!
!-----------------------------------------------------------------------
!
          ENDDO
          ENDDO main_loop
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(PMID_INTERP)
          DEALLOCATE( VBL_INTERP)
!
          DEALLOCATE(PINT_INTERP_HI)
          DEALLOCATE(PINT_INTERP_LO)
!
!-----------------------------------------------------------------------
!
        ENDDO child_task_loop
!
!-----------------------------------------------------------------------
!
      ENDDO loop_sides
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_UPDATE_CHILD_BNDRY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PRESSURE_ON_NEST_BNDRY_V(PD_H                          &
                                         ,IMS,IME,JMS,JME               &
!
                                         ,PDB_S_H                       &
                                         ,PDB_N_H                       &
                                         ,PDB_W_H                       &
                                         ,PDB_E_H                       &
!
                                         ,NUM_TASKS_SEND_SBND           &
                                         ,NUM_TASKS_SEND_NBND           &
                                         ,NUM_TASKS_SEND_WBND           &
                                         ,NUM_TASKS_SEND_EBND           &
!
                                         ,I_LO_SOUTH_V                  &
                                         ,I_HI_SOUTH_V                  &
                                         ,I_LO_NORTH_V                  &
                                         ,I_HI_NORTH_V                  &
                                         ,J_LO_WEST_V                   &
                                         ,J_HI_WEST_V                   &
                                         ,J_LO_EAST_V                   &
                                         ,J_HI_EAST_V                   & 
!
                                         ,I_LO_SOUTH_H                  &
                                         ,I_HI_SOUTH_H                  &
                                         ,I_LO_NORTH_H                  &
                                         ,I_HI_NORTH_H                  &
                                         ,J_LO_WEST_H                   &
                                         ,J_HI_WEST_H                   &
                                         ,J_LO_EAST_H                   &
                                         ,J_HI_EAST_H                   & 
!                                                                               ^
                                         ,N_BLEND_V                     & !     |
                                         ,IM_CHILD                      & !     |
                                         ,JM_CHILD                      & !     |
                                                                          !     |
                                         ,INC_FIX_N                     & !   INPUT
!                                                                           ---------
                                         ,PD_V                          & !   OUTPUT
                                         ,PDB_S_V                       & !     |
                                         ,PDB_N_V                       & !     |
                                         ,PDB_W_V                       & !     v
                                         ,PDB_E_V )                       !  
!
!-----------------------------------------------------------------------
!***  USE 4-PT HORIZONTAL INTERPOLATION TO COMPUTE PD ON V POINTS
!***  OF THE PARENT DOMAIN AND OF THE NEST BOUNDARY GIVEN THOSE 
!***  VALUES ON H POINTS.
!-----------------------------------------------------------------------
!
!-----------
!***  Input
!-----------
!
      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME                             &  !<-- Parent task's memory limits
                           ,IM_CHILD,JM_CHILD                           &  !<-- Nest domain limits
                           ,INC_FIX_N                                   &  !<-- Increment for selecting nest tasks for averaging H to V
                           ,N_BLEND_V                                   &  !<-- V rows in nest's boundary region
                           ,NUM_TASKS_SEND_SBND                         &  !<-- # of child tasks with Sbndry V points on parent task
                           ,NUM_TASKS_SEND_NBND                         &  !<-- # of child tasks with Nbndry V points on parent task
                           ,NUM_TASKS_SEND_WBND                         &  !<-- # of child tasks with Wbndry V points on parent task
                           ,NUM_TASKS_SEND_EBND                            !<-- # of child tasks with Ebndry V points on parent task
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: I_LO_SOUTH_V           &  !<-- Starting I of Sbndry region V points on child tasks
                                                ,I_HI_SOUTH_V           &  !<-- Ending I of Sbndry region V points on child tasks
                                                ,I_LO_NORTH_V           &  !<-- Starting I of Nbndry region V points on child tasks
                                                ,I_HI_NORTH_V           &  !<-- Ending I of Nbndry region V points on child tasks
                                                ,J_LO_WEST_V            &  !<-- Starting J of Wbndry region V points on child tasks
                                                ,J_HI_WEST_V            &  !<-- Ending J of Wbndry region V points on child tasks
                                                ,J_LO_EAST_V            &  !<-- Starting J of Ebndry region V points on child tasks
                                                ,J_HI_EAST_V               !<-- Ending J of Ebndry region V points on child tasks
!
      INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: I_LO_SOUTH_H           &  !<-- Starting I of Sbndry region H points on child tasks
                                                ,I_HI_SOUTH_H           &  !<-- Ending I of Sbndry region H points on child tasks
                                                ,I_LO_NORTH_H           &  !<-- Starting I of Nbndry region H points on child tasks
                                                ,I_HI_NORTH_H           &  !<-- Ending I of Nbndry region H points on child tasks
                                                ,J_LO_WEST_H            &  !<-- Starting J of Wbndry region H points on child tasks
                                                ,J_HI_WEST_H            &  !<-- Ending J of Wbndry region H points on child tasks
                                                ,J_LO_EAST_H            &  !<-- Starting J of Ebndry region H points on child tasks
                                                ,J_HI_EAST_H               !<-- Ending J of Ebndry region H points on child tasks
!
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PD_H                   !<-- Parent PD (Pa) (column mass in sigma domain) on H points
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: PDB_S_H        &  !<-- Boundary region PD (Pa) (column mass in sigma domain)
                                                        ,PDB_N_H        &  !    on mass points on the four sides of the child boundary.
                                                        ,PDB_W_H        &  !
                                                        ,PDB_E_H           !<--
!
!------------
!***  Output
!------------
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: PD_V                  !<-- Parent PD (Pa) (column mass in sigma domain) on V points
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(OUT) :: PDB_S_V       &  !<-- Child boundary PD (Pa) on child domain Sbndry V points
                                                         ,PDB_N_V       &  !<-- Child boundary PD (Pa) on child domain Nbndry V points
                                                         ,PDB_W_V       &  !<-- Child boundary PD (Pa) on child domain Wbndry V points
                                                         ,PDB_E_V          !<-- Child boundary PD (Pa) on child domain Ebndry V points
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: DIFF_START,DIFF_START_PTS                              &
                ,I,J                                                    &
                ,I_ADD_H,I_INC_H                                        &
                ,I_START_H,I_START_V                                    &
                ,I_END_H,I_END_V                                        &
                ,I_HI_H,I_LO_H                                          &
                ,J_START_H,J_START_V                                    &
                ,J_END_H,J_END_V                                        &
                ,KNT_H,KNT_V                                            &
                ,N_ADD_H,N_OFFSET_H,N_SIDE                              &
                ,NTX,NTX_H,NUM_TASKS_SEND
!
      INTEGER,DIMENSION(4) :: NTX_ADD
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER :: PDB_H                     &
                                             ,PDB_V
!
      integer,dimension(8) :: values
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      DO N_SIDE=1,4
        NTX_ADD(N_SIDE)=0
      ENDDO
!
!-----------------------------------------------------------------------
!***  FIRST OBTAIN PD ON THE PARENT'S V POINTS.
!-----------------------------------------------------------------------
!
      DO J=JMS,JME-1
      DO I=IMS,IME-1
        PD_V(I,J)=0.25*(PD_H(I,J)+PD_H(I+1,J)+PD_H(I,J+1)+PD_H(I+1,J+1))
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOW AVERAGE VALUES OF PD ON AND ADJACENT TO THE NESTS' BOUNDARY
!***  TO OBTAIN PD ON THE NESTS' BOUNDARY V POINTS.
!-----------------------------------------------------------------------
!
      loop_sides: DO N_SIDE=1,4                                            !<-- Loop through the 4 lateral boundaries (S,N,W,E)
!
!-----------------------------------------------------------------------
!
        IF(N_SIDE==1)THEN
          IF(NUM_TASKS_SEND_SBND==0)CYCLE                                  !<-- Move on if parent task sees no part of nest's Sbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_SBND                               !<-- # of Sbndry nest tasks to which parent sends PDB on V
          PDB_H=>PDB_S_H                                                   !<-- String of Sbndry PDB on H for this nest task
          PDB_V=>PDB_S_V                                                   !<-- String of Sbndry PDB on V for this nest task
!
!
        ELSEIF(N_SIDE==2)THEN
          IF(NUM_TASKS_SEND_NBND==0)CYCLE                                  !<-- Move on if parent task sees no part of nest's Nbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_NBND                               !<-- # of Nbndry nest tasks to which parent sends PDB on V
          PDB_H=>PDB_N_H                                                   !<-- String of Nbndry PDB on H for this nest task
          PDB_V=>PDB_N_V                                                   !<-- String of Nbndry PDB on V for this nest task
!
        ELSEIF(N_SIDE==3)THEN
          IF(NUM_TASKS_SEND_WBND==0)CYCLE                                  !<-- Move on if parent task sees no part of nest's Wbndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_WBND                               !<-- # of Wbndry nest tasks to which parent sends PDB on V
          PDB_H=>PDB_W_H                                                   !<-- String of Wbndry PDB on H for this nest task
          PDB_V=>PDB_W_V                                                   !<-- String of Wbndry PDB on V for this nest task
!
        ELSEIF(N_SIDE==4)THEN
          IF(NUM_TASKS_SEND_EBND==0)CYCLE                                  !<-- Move on if parent task sees no part of nest's Ebndry
!
          NUM_TASKS_SEND=NUM_TASKS_SEND_EBND                               !<-- # of Ebndry nest tasks to which parent sends PDB on V
          PDB_H=>PDB_E_H                                                   !<-- String of Ebndry PDB on H for this nest task
          PDB_V=>PDB_E_V                                                   !<-- String of Ebndry PDB on V for this nest task
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        child_task_loop: DO NTX=1,NUM_TASKS_SEND                           !<-- Compute PD on V points for each child task with
                                                                           !    bndry segments seen by this parent task.
!-----------------------------------------------------------------------
!
!----------------------------------------------
!***  South boundary limits on this child task
!----------------------------------------------
!
          IF(N_SIDE==1)THEN
            I_START_V=I_LO_SOUTH_V(NTX)                                    !<-- I index of first Sbndry V point on this child task
            I_END_V  =I_HI_SOUTH_V(NTX)                                    !<-- I index of last Sbndry V point on this child task
            J_START_V=1                                                    !<-- J index of first Sbndry V point on this child task
            J_END_V  =N_BLEND_V                                            !<-- J index of last Sbndry V point on this child task
!
            I_START_H=I_LO_SOUTH_H(NTX)                                    !<-- I index of 1st Sbndry H point west of 1st Sbndry V point
            I_END_H  =I_HI_SOUTH_H(NTX)                                    !<-- I index of last Sbndry H point east of last Sbndry V point
            J_START_H=J_START_V                                            !<-- J index of 1st Sbndry H point west of 1st Sbndry V point
            J_END_H  =J_END_V                                              !<-- J index of last Sbndry H point east of last Sbndry V point
!
            IF(NTX==1.AND.NUM_TASKS_SEND>1                              &
                     .AND.I_START_V==I_LO_SOUTH_H(NTX+1)+INC_FIX_N      &
                     .AND.I_END_V>=I_END_H)THEN
              NTX_ADD(N_SIDE)=1
            ENDIF
!
            NTX_H=NTX+NTX_ADD(N_SIDE)                                      !<-- Reference parent task for PDB on nest H points
            DIFF_START=I_START_V-I_LO_SOUTH_H(NTX_H)
            I_HI_H=I_HI_SOUTH_H(NTX_H)
            I_LO_H=I_LO_SOUTH_H(NTX_H)
            N_ADD_H=I_HI_H-I_LO_H+1
!
!----------------------------------------------
!***  North boundary limits on this child task
!----------------------------------------------
!
          ELSEIF(N_SIDE==2)THEN
            I_START_V=I_LO_NORTH_V(NTX)                                    !<-- I index of first Nbndry V point on this child task
            I_END_V  =I_HI_NORTH_V(NTX)                                    !<-- I index of last Nbndry V point on this child task
            J_START_V=JM_CHILD-N_BLEND_V                                   !<-- J index of first Nbndry V point on this child task
            J_END_V  =JM_CHILD-1                                           !<-- J index of last Nbndry V point on this child task
!
            I_START_H=I_LO_NORTH_H(NTX)                                    !<-- I index of 1st Nbndry H point west of 1st Nbndry V point
            I_END_H  =I_HI_NORTH_H(NTX)                                    !<-- I index of last Nbndry H point east of last Nbndry V point
            J_START_H=J_START_V                                            !<-- J index of 1st Nbndry H point west of 1st Nbndry V point
            J_END_H  =J_END_V                                              !<-- J index of last Nbndry H point east of last Nbndry V point
!
            IF(NTX==1.AND.NUM_TASKS_SEND>1                              &
                     .AND.I_START_V==I_LO_NORTH_H(NTX+1)+INC_FIX_N      &
                     .AND.I_END_V>=I_END_H)THEN
              NTX_ADD(N_SIDE)=1
            ENDIF
            NTX_H=NTX+NTX_ADD(N_SIDE)                                      !<-- Reference parent task for PDB on nest H points
            DIFF_START=I_START_V-I_LO_NORTH_H(NTX_H)
            I_HI_H=I_HI_NORTH_H(NTX_H)
            I_LO_H=I_LO_NORTH_H(NTX_H)
            N_ADD_H=I_HI_H-I_LO_H+1
!
!----------------------------------------------
!***  West boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==3)THEN
            I_START_V=1                                                    !<-- I index of first Wbndry V point on this child task
            I_END_V  =N_BLEND_V                                            !<-- I index of last Wbndry V point on this child task
            J_START_V=J_LO_WEST_V(NTX)                                     !<-- J index of first Wbndry V point on this child task
            J_END_V  =J_HI_WEST_V(NTX)                                     !<-- J index of last Wbndry V point on this child task
!
            I_START_H=I_START_V                                            !<-- I index of 1st Wbndry H point west of 1st Wbndry V point
            I_END_H  =I_END_V+1                                            !<-- I index of last Wbndry H point east of last Wbndry V point
            J_START_H=J_LO_WEST_H(NTX)                                     !<-- J index of 1st Wbndry H point west of 1st Wbndry V point
            J_END_H  =J_HI_WEST_H(NTX)                                     !<-- J index of last Wbndry H point east of last Wbndry V point
!
            IF(NTX==1.AND.NUM_TASKS_SEND>1                              &
                     .AND.J_START_V==J_LO_WEST_H(NTX+1)+INC_FIX_N       &
                     .AND.J_END_V>=J_END_H)THEN
              NTX_ADD(N_SIDE)=1
            ENDIF
!
            NTX_H=NTX+NTX_ADD(N_SIDE)                                      !<-- Reference parent task for PDB on nest H points
            DIFF_START=J_START_V-J_LO_WEST_H(NTX_H)
            DIFF_START_PTS=DIFF_START*(N_BLEND_H+1)
            N_ADD_H=N_BLEND_H+1
!
!---------------------------------------------
!***  East boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==4)THEN
            I_START_V=IM_CHILD-N_BLEND_V                                   !<-- I index of first Ebndry V point on this child task
            I_END_V  =IM_CHILD-1                                           !<-- I index of last Ebndry V point on this child task
            J_START_V=J_LO_EAST_V(NTX)                                     !<-- J index of first Ebndry V point on this child task
            J_END_V  =J_HI_EAST_V(NTX)                                     !<-- J index of last Ebndry V point on this child task
!
            I_START_H=I_START_V                                            !<-- I index of 1st Ebndry H point west of 1st Ebndry V point
            I_END_H  =I_END_V+1                                            !<-- I index of last Ebndry H point east of last Ebndry V point
            J_START_H=J_LO_EAST_H(NTX)                                     !<-- J index of 1st Ebndry H point west of 1st Ebndry V point
            J_END_H  =J_HI_EAST_H(NTX)                                     !<-- J index of last Ebndry H point east of last Ebndry V point
!
            IF(NTX==1.AND.NUM_TASKS_SEND>1                              &
                     .AND.J_START_V==J_LO_EAST_H(NTX+1)+INC_FIX_N       &
                     .AND.J_END_V>=J_END_H)THEN
              NTX_ADD(N_SIDE)=1
            ENDIF
!
            NTX_H=NTX+NTX_ADD(N_SIDE)                                      !<-- Reference parent task for PDB on nest H points
            DIFF_START=J_START_V-J_LO_EAST_H(NTX_H)
            DIFF_START_PTS=DIFF_START*(N_BLEND_H+1)
            N_ADD_H=N_BLEND_H+1
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  READY TO AVERAGE NEST BOUNDARY PD ON H TO PD ON V.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RECALL THAT THE PDB_H POINTER CONTAINS DATA FOR ONE EXTRA ROW
!***  BEYOND THE CHILD BOUNDARY ON ALL SIDES SINCE WE NEED THE
!***  ABILITY TO DO 4-PT AVERAGES OF PD TO V POINTS.
!***  HERE WE MUST USE PDB_H BUT THE VALUES FOR PDB_V (AND THUS U,V)
!***  ARE ONLY GENERATED ON THE TRUE BOUNDARY POINTS.  TO ADDRESS PDB_H
!***  CORRECTLY WE MUST TAKE INTO ACCOUNT THOSE EXTRA POINTS IN PDB_H
!***  AS WE MARCH THROUGH THE V POINT LOCATIONS.
!
!***  ALSO WE MUST BE AWARE THAT THE NEST BOUNDARY V SEGMENTS THAT ARE
!***  BEING CONSIDERED FOR EACH NEST TASK MUST CORRESPOND TO THE SAME
!***  NEST TASK H POINT VALUES OF PDB.  BECAUSE OF OVERLAP THAT CAN
!***  OCCUR DUE TO THE SEGMENTS OF PDB ON H BEING LARGER THAN PBD ON
!***  V WE MUST EXPLICITLY CHECK TO BE SURE NEST PDB ON V IS BEING
!***  FILLED BY NEST PDB ON H FOR THE SAME NEST TASK.
!-----------------------------------------------------------------------
!
          KNT_V=0
!
          DO J=J_START_V,J_END_V
            IF(N_SIDE<=2)THEN
              N_OFFSET_H=(DIFF_START+1)*(J-J_START_V+1)-1
            ELSEIF(N_SIDE>=3)THEN
              N_OFFSET_H=(J-J_START_V)*N_BLEND_H+DIFF_START_PTS
            ENDIF
!
            DO I=I_START_V,I_END_V
              KNT_V=KNT_V+1
              KNT_H=KNT_V+N_OFFSET_H
              PDB_V(NTX)%DATA(KNT_V)=(PDB_H(NTX_H)%DATA(KNT_H)            &
                                     +PDB_H(NTX_H)%DATA(KNT_H+1)          &
                                     +PDB_H(NTX_H)%DATA(KNT_H+N_ADD_H)    &
                                     +PDB_H(NTX_H)%DATA(KNT_H+N_ADD_H+1)) &
                                     *0.25
            ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDDO child_task_loop
!
!-----------------------------------------------------------------------
!
      ENDDO loop_sides
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PRESSURE_ON_NEST_BNDRY_V
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE CHILD_DATA_FROM_STRING(LENGTH_DATA                     &
                                       ,DATASTRING                      &
                                       ,ILIM_LO,ILIM_HI                 &
                                       ,JLIM_LO,JLIM_HI                 &
                                       ,I_START,I_END                   &
                                       ,J_START,J_END                   &
                                       ,PDB,TB,QB,CWB                   &
                                       ,UB,VB )
!
!-----------------------------------------------------------------------
!***  EXTRACT VARIABLES FOR A NEST'S BOUNDARY FROM THE DATASTRING
!***  RECEIVED BY THE CHILD FROM ITS PARENT.  A CHILD TASK MIGHT 
!***  RECEIVE SEGMENTS OF BOUNDARY DATA FROM TWO PARENT TASKS IN 
!***  WHICH CASE THE PIECES ARE BE COMBINED.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,INTENT(IN) :: LENGTH_DATA                                 &  !<-- # of words in datastring
                           ,I_START                                     &  !<-- Starting I of data segment in string on child's grid
                           ,I_END                                       &  !<-- Ending I of data segment in string on child's grid
                           ,J_START                                     &  !<-- Starting J of data segment in string on child's grid
                           ,J_END                                       &  !<-- Ending J of data segment in string on child's grid
                           ,ILIM_LO                                     &  !<-- Lower I limit of full boundary segment 
                           ,ILIM_HI                                     &  !<-- Upper I limit of full boundary segment 
                           ,JLIM_LO                                     &  !<-- Lower J limit of full boundary segment 
                           ,JLIM_HI                                        !<-- Upper J limit of full boundary segment 
!
      REAL,DIMENSION(:),POINTER,INTENT(IN) :: DATASTRING                   !<-- The string of boundary data from the parent
!      
      REAL,DIMENSION(ILIM_LO:ILIM_HI                                    &
                    ,JLIM_LO:JLIM_HI),INTENT(OUT),OPTIONAL :: PDB          !<-- PD for segment of the child boundary
!      
      REAL,DIMENSION(ILIM_LO:ILIM_HI                                    &
                    ,JLIM_LO:JLIM_HI                                    &
                    ,      1:LM     ),INTENT(OUT),OPTIONAL :: TB        &  !<-- Temperature for segment of the child boundary
                                                             ,QB        &  !<-- Specific humidity for segment of the child boundary
                                                             ,CWB       &  !<-- Cloud condensate for segment of the child boundary
                                                             ,UB        &  !<-- U wind for segment of the child boundary
                                                             ,VB           !<-- V wind for segment of the child boundary
!      
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,K,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  EXTRACT THE APPROPRIATE BOUNDARY VARIABLES FROM THE DATASTRING
!***  PROVIDED BY THE PARENT.
!***  HERE WE TRANSFER FROM THE 1-D STRING TO THE 2-D AND 3-D ARRAYS.
!***  THIS EASILY ALLOWS FOR PROPER COMBINING OF SEPARATE STRINGS
!***  WHOSE ENDS MAY OVERLAP THAT ARE ARRIVING FROM DIFFERENT PARENT
!***  TASKS.  TO EXPORT THESE BOUNDARY ARRAYS OUT OF THE COUPLER
!***  THOUGH THE DATA MUST BE PUT BACK INTO 1-D.
!-----------------------------------------------------------------------
!
      N=0
!
      IF(PRESENT(TB))THEN                                                  !<-- Datastring for mass variables has been sent in
!
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          PDB(I,J)=DATASTRING(N)
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          TB(I,J,K)=DATASTRING(N)
        ENDDO
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          QB(I,J,K)=DATASTRING(N)
        ENDDO
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          CWB(I,J,K)=DATASTRING(N)
        ENDDO
        ENDDO
        ENDDO
!
      ELSEIF(PRESENT(UB))THEN                                              !<-- Datastring for wind variables has been sent in
!
        DO K=1,LM
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          UB(I,J,K)=DATASTRING(N)
        ENDDO
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=J_START,J_END
        DO I=I_START,I_END
          N=N+1
          VB(I,J,K)=DATASTRING(N)
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CHILD_DATA_FROM_STRING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE EXPORT_CHILD_BOUNDARY(PDB                              &
                                      ,TB                               &
                                      ,QB                               &
                                      ,CWB                              &
!
                                      ,UB                               &
                                      ,VB                               &
!
                                      ,ILIM_LO,ILIM_HI                  &
                                      ,JLIM_LO,JLIM_HI                  &
!
                                      ,DATA_NAME                        &
!
                                      ,DATA_EXP                         &
!
                                      ,EXPORT_STATE )
!
!-----------------------------------------------------------------------
!***  LOAD THE CHILD BOUNDARY VALUES RECEIVED FROM THE PARENT 
!***  INTO THE PARENT-CHILD COUPLER EXPORT STATE.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,INTENT(IN) :: ILIM_LO                                     &  !<-- Lower I limit of full boundary segment
                           ,ILIM_HI                                     &  !<-- Upper I limit of full boundary segment
                           ,JLIM_LO                                     &  !<-- Lower J limit of full boundary segment
                           ,JLIM_HI                                        !<-- Upper J limit of full boundary segment
!      
      REAL,DIMENSION(ILIM_LO:ILIM_HI                                    &
                    ,JLIM_LO:JLIM_HI),INTENT(IN),OPTIONAL :: PDB           !<-- PD for segment of the child boundary
!      
      REAL,DIMENSION(ILIM_LO:ILIM_HI                                    &
                    ,JLIM_LO:JLIM_HI                                    &
                    ,      1:LM     ),INTENT(IN),OPTIONAL :: TB         &  !<-- Temperature for segment of the child boundary
                                                            ,QB         &  !<-- Specific humidity for segment of the child boundary
                                                            ,CWB        &  !<-- Cloud condensate for segment of the child boundary
                                                            ,UB         &  !<-- U wind for segment of the child boundary
                                                            ,VB            !<-- V wind for segment of the child boundary
!
      CHARACTER(*),INTENT(IN) :: DATA_NAME                                 !<-- Name used for each child task's boundary segment
!
      REAL,DIMENSION(:),POINTER,INTENT(INOUT) :: DATA_EXP                  !<-- Combined boundary segment data on child task
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXPORT_STATE                       !<-- The Parent-Child Coupler export state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,K
      INTEGER :: ISTAT,RC,RC_EXP_BNDRY
      INTEGER,SAVE :: NN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_EXP_BNDRY=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  THE CHILDREN'S FINAL BOUNDARY DATA WILL BE LOADED INTO THE
!***  PARENT-CHILD COUPLER'S EXPORT STATE AS 1-D ARRAYS (Attributes)
!***  SINCE THEY ARE NOT SPANNING THE CHILDRENS' ESMF Grids (as Fields).
!***  BUT BECAUSE THEY ARE ATTRIBUTES AND NOT FIELDS, THEY WILL NEED
!***  TO BE RESET EVERY TIME.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      main: IF(PRESENT(PDB))THEN                                           !<-- If true then H point variables were sent in
!-----------------------------------------------------------------------
!
        NN=0                                                               !<-- The data always begins with H points
!
        DO J=JLIM_LO,JLIM_HI
        DO I=ILIM_LO,ILIM_HI
          NN=NN+1
          DATA_EXP(NN)=PDB(I,J)                                            !<-- Insert complete PDB into the 1-D H boundary data
        ENDDO
        ENDDO
!
        DO K=1,LM
        DO J=JLIM_LO,JLIM_HI
        DO I=ILIM_LO,ILIM_HI
          DATA_EXP(NN+1)= TB(I,J,K)                                        !<-- Insert complete TB into the 1-D H boundary data
          DATA_EXP(NN+2)= QB(I,J,K)                                        !<-- Insert complete QB into the 1-D H boundary data
          DATA_EXP(NN+3)=CWB(I,J,K)                                        !<-- Insert complete CWB into the 1-D H boundary data
          NN=NN+3
        ENDDO
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry Data "//DATA_NAME//" Into Coupler Export"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                  &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                     &  !<-- Name of the children's new boundary H data
                              ,count    =NN                            &  !<-- # of words in the data 
                              ,valueList=DATA_EXP                      &  !<-- The children's new boundary H data
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_EXP_BNDRY)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ELSEIF(PRESENT(UB))THEN  
!
!-----------------------------------------------------------------------
!
        NN=0
!
        DO K=1,LM
        DO J=JLIM_LO,JLIM_HI
        DO I=ILIM_LO,ILIM_HI
          DATA_EXP(NN+1)=UB(I,J,K)
          DATA_EXP(NN+2)=VB(I,J,K)
          NN=NN+2
        ENDDO
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry Data "//DATA_NAME//" Into Coupler Export"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                    &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                       &  !<-- Name of the children's new boundary V data
                              ,count    =NN                              &  !<-- # of words in the data 
                              ,valueList=DATA_EXP                        &  !<-- The children's new boundary V data
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_EXP_BNDRY)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF main
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EXPORT_CHILD_BOUNDARY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE SPLINE(NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW)
!-----------------------------------------------------------------------
!
!     ******************************************************************
!     *                                                                *
!     *  This is a one-dimensional cubic spline fitting routine        *
!     *  programmed for a small scalar machine.                        *
!     *                                                                *
!     *  Programmer: Z. Janjic, Yugoslav Fed. Hydromet. Inst., Beograd *
!     *                                                                *
!     *  NOLD - Number of given values of the function.  Must be >= 3. *
!     *  XOLD - Locations of the points at which the values of the     *
!     *         function are given.  Must be in ascending order.       *
!     *  YOLD - The given values of the function at the points XOLD.   *
!     *  Y2   - The second derivatives at the points XOLD.  If natural *
!     *         spline is fitted Y2(1)=0 and Y2(nold)=0. Must be       *
!     *         specified.                                             *
!     *  NNEW - Number of values of the function to be calculated.     *
!     *  XNEW - Locations of the points at which the values of the     *
!     *         function are calculated.  XNEW(K) must be >= XOLD(1)   *
!     *         and <= XOLD(NOLD).                                     *
!     *  YNEW - The values of the function to be calculated.           *
!     *                                                                *
!     ******************************************************************
!
!-----------------------------------------------------------------------
!***  Arguments
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NNEW,NOLD
!
      REAL,DIMENSION(1:NOLD),INTENT(IN) :: XOLD,YOLD
      REAL,DIMENSION(1:NNEW),INTENT(IN) :: XNEW
!
      REAL,DIMENSION(1:NOLD),INTENT(INOUT) :: Y2
!
      REAL,DIMENSION(1:NNEW),INTENT(OUT) :: YNEW
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER :: K,K1,K2,KOLD,NOLDM1
!
      REAL :: AK,BK,CK,DEN,DX,DXC,DXL,DXR,DYDXL,DYDXR,RDX,RTDXC         &
             ,X,XK,XSQ,Y2K,Y2KP1
!
      REAL,DIMENSION(1:NOLD-2) :: P,Q
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NOLDM1=NOLD-1
!
      DXL=XOLD(2)-XOLD(1)
      DXR=XOLD(3)-XOLD(2)
      DYDXL=(YOLD(2)-YOLD(1))/DXL
      DYDXR=(YOLD(3)-YOLD(2))/DXR
      RTDXC=0.5/(DXL+DXR)
!
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
      Q(1)=-RTDXC*DXR
!
      IF(NOLD==3) GO TO 700
!
!-----------------------------------------------------------------------
      K=3
!
  100 CONTINUE
      DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
!
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
!
      K=K+1
      IF(K<NOLD) GO TO 100
!
!-----------------------------------------------------------------------
!
  700 CONTINUE
      K=NOLDM1
!
  200 CONTINUE
      Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
!
      K=K-1
      IF(K>1) GO TO 200
!
!-----------------------------------------------------------------------
!
      K1=1
!
  300 CONTINUE
      XK=XNEW(K1)
!
      DO 400 K2=2,NOLD
        IF(XOLD(K2)<=XK) GO TO 400
        KOLD=K2-1
        GO TO 450
  400 CONTINUE
!
      YNEW(K1)=YOLD(NOLD)
      GO TO 600
!
  450 CONTINUE
      IF(K1==1)   GO TO 500
      IF(K==KOLD) GO TO 550
!
  500 CONTINUE
      K=KOLD
!
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
!
      AK=0.1666667*RDX*(Y2KP1-Y2K)
      BK=0.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-0.1666667*DX*(Y2KP1+Y2K+Y2K)
!
  550 CONTINUE
      X=XK-XOLD(K)
      XSQ=X*X
!
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
!
  600 CONTINUE
      K1=K1+1
!
      IF(K1<=NNEW) GO TO 300
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SPLINE
!
!-----------------------------------------------------------------------
!      
      END MODULE MODULE_PARENT_CHILD_CPL_COMP
!
!-----------------------------------------------------------------------
